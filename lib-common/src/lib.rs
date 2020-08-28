/// lib-common -- shared functionality
use std::path::Path;

use regex::Regex;

use bio_types::genome::Interval;
use rust_htslib::{bam, bam::Read, bcf};

pub mod stats;

/// Global error type.
#[derive(thiserror::Error, Debug)]
pub enum Error {
    /// A command line option is missing.
    #[error("missing command line argument")]
    OptionMissing(),
    /// Inconsistent input files.
    #[error("inconsistent input files")]
    InconsistentInput(),
    /// The output file already exists.
    #[error("output file already exists")]
    OutputFileExists(),
    /// Incorrect cluster setting name.
    #[error("unknown cluster setting name")]
    UnknownClusterSettingName(),
    /// Problem with file I/O.
    #[error("problem with I/O")]
    Io {
        #[from]
        source: std::io::Error,
        // TODO: add experimental backtrace feature?
    },
    /// Problem with htslib
    #[error("problem with BCF file access")]
    HtslibBcfError {
        #[from]
        source: bcf::errors::Error, // TODO: add experimental backtrace feature?
    },
    /// Problem with htslib
    #[error("problem with BAM file access")]
    HtslibBamError {
        #[from]
        source: bam::errors::Error, // TODO: add experimental backtrace feature?
    },
    /// Problem with string conversion
    #[error("problem with string conversion")]
    StrUtf8Error {
        #[from]
        source: std::str::Utf8Error, // TODO: add experimental backtrace feature?
    },
    /// Problem with string conversion
    #[error("problem with string conversion")]
    StringUtf8Error {
        #[from]
        source: std::string::FromUtf8Error, // TODO: add experimental backtrace feature?
    },
}

/// Return `bam::Format` for the given filename.
pub fn guess_bam_format(filename: &str) -> bam::Format {
    if filename.ends_with(".bam") {
        bam::Format::BAM
    } else {
        bam::Format::SAM
    }
}

#[derive(Debug)]
pub struct BcfFormatInfo {
    /// Guessed format.
    pub format: bcf::Format,
    /// Guessed compression status.
    pub uncompressed: bool,
}

/// Return `bcf::Format` for the given filename.
pub fn guess_bcf_format(filename: &str) -> BcfFormatInfo {
    if filename.ends_with(".bcf") {
        return BcfFormatInfo {
            format: bcf::Format::BCF,
            uncompressed: false,
        };
    } else if filename.ends_with(".vcf.gz") {
        return BcfFormatInfo {
            format: bcf::Format::VCF,
            uncompressed: false,
        };
    } else {
        return BcfFormatInfo {
            format: bcf::Format::VCF,
            uncompressed: true,
        };
    }
}

/// Build new `bcf::Header`.
pub fn build_vcf_header(template: &bcf::header::HeaderView) -> Result<bcf::Header, bcf::Error> {
    let mut header = bcf::Header::new();

    // Copy over sequence records.
    for record in template.header_records() {
        match record {
            bcf::header::HeaderRecord::Contig { key: _, values } => {
                header.push_record(
                    format!(
                        "##contig=<ID={},length={}>",
                        values
                            .get("ID")
                            .expect("source contig header does not have ID"),
                        values
                            .get("length")
                            .expect("source contig header does not have length")
                    )
                    .as_bytes(),
                );
            }
            _ => (),
        }
    }

    // Fields: ALT, INFO, FORMAT.
    let alts = vec![
        ("DEL", "Deletion"),
        ("DUP", "Duplication"),
        ("INV", "Inversion"),
    ];
    for (id, desc) in alts {
        header.push_record(format!("##ALT=<ID={},length={}>", &id, &desc).as_bytes());
    }
    let infos = vec![
        ("SVTYPE", "1", "String", "Type of structural variant"),
        ("CHR2", "1", "String", "Chromosome of end coordinate"),
        ("END", "1", "Integer", "End position of linear SV"),
        ("END2", "1", "Integer", "End position of BND"),
        ("STRANDS", "1", "String", "Breakpoint strandedness"),
        ("SVLEN", "1", "Integer", "SV length"),
        ("ALGORITHMS", ".", "String", "Source algorithms"),
    ];
    for (id, number, type_, desc) in infos {
        header.push_record(
            format!(
                "##INFO=<ID={},Number={},Type={},Description={}>",
                &id, &number, &type_, &desc
            )
            .as_bytes(),
        );
    }
    let formats = vec![
        ("GT", "1", "String", "Genotype"),
        ("delly", "1", "Integer", "Called by Delly"),
    ];
    for (id, number, type_, desc) in formats {
        header.push_record(
            format!(
                "##FORMAT=<ID={},Number={},Type={},Description={}>",
                &id, &number, &type_, &desc
            )
            .as_bytes(),
        );
    }

    // Add samples.
    for name in template.samples() {
        header.push_sample(name);
    }

    Ok(header)
}

/// Generate list of all contigs from BAM header.
pub fn build_chroms_bam(
    header: &bam::HeaderView,
    regex: Option<String>,
) -> Result<Vec<Interval>, Error> {
    let mut result: Vec<Interval> = Vec::new();

    let re = if let Some(regex) = regex {
        Some(Regex::new(&regex).unwrap())
    } else {
        None
    };

    for (no, name) in header.target_names().iter().enumerate() {
        let name = std::str::from_utf8(name)?.to_string();
        let len = header.target_len(no as u32).unwrap() as u64;
        let do_push = match &re {
            Some(re) => re.is_match(&name),
            None => true,
        };
        if do_push {
            result.push(Interval::new(name, 0..len));
        }
    }

    Ok(result)
}

/// Parse @RG lane into triple (id, sm).
fn parse_line_rg(line: String) -> Option<(String, String)> {
    let line_split = line.split("\t");
    let mut id: Option<String> = None;
    let mut sm: Option<String> = None;
    for s in line_split {
        let token: Vec<&str> = s.split(":").collect();
        if token.len() >= 2 {
            match token[0] {
                "ID" => {
                    id = Some(token[1].to_string());
                }
                "SM" => {
                    sm = Some(token[1].to_string());
                }
                _ => (),
            }
        }
    }

    match (id, sm) {
        (Some(id), Some(sm)) => Some((id, sm)),
        _ => None,
    }
}

/// Construct sample list from BAM header.
fn samples_from_header(header: &Vec<u8>) -> Vec<String> {
    let text = String::from_utf8(header.to_vec()).unwrap();
    let mut samples = Vec::new();

    for line in text.lines() {
        if line.starts_with("@RG") {
            match parse_line_rg(line.to_string()) {
                Some((_id, sm)) => {
                    if !samples.contains(&sm) {
                        samples.push(sm.clone());
                    }
                }
                None => (),
            }
        }
    }

    samples
}

/// Construct sample list from BAM file.
pub fn samples_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<String>, Error> {
    let path = path.as_ref().to_str().unwrap();
    let reader = bam::Reader::from_path(path)?;

    Ok(samples_from_header(&Vec::from(reader.header().as_bytes())))
}

/// Enumeration for calling algorithms.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Algorithm {
    /// Delly2
    Delly,
}

#[cfg(test)]
mod tests {
    use super::*;
    use matches::assert_matches;

    #[test]
    fn test_guess_bam_format() {
        assert_matches!(guess_bam_format("ex.sam"), bam::Format::SAM);
        assert_matches!(guess_bam_format("ex.bam"), bam::Format::BAM);
        assert_matches!(guess_bam_format("ex.xxx"), bam::Format::SAM);
    }

    #[test]
    fn test_guess_bcf_format() {
        assert_matches!(guess_bcf_format("ex.vcf").format, bcf::Format::VCF);
        assert_eq!(guess_bcf_format("ex.vcf").uncompressed, true);
        assert_matches!(guess_bcf_format("ex.vcf.gz").format, bcf::Format::VCF);
        assert_eq!(guess_bcf_format("ex.vcf.gz").uncompressed, false);
        assert_matches!(guess_bcf_format("ex.bcf").format, bcf::Format::BCF);
        assert_eq!(guess_bcf_format("ex.bcf").uncompressed, false);
        assert_matches!(guess_bcf_format("ex.xxx").format, bcf::Format::VCF);
        assert_eq!(guess_bcf_format("ex.xxx").uncompressed, true);
    }
}
