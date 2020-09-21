use rust_htslib::{bcf, bcf::Read};

use super::error::Error;

#[derive(Debug)]
pub struct FormatInfo {
    /// Guessed format.
    pub format: bcf::Format,
    /// Guessed compression status.
    pub uncompressed: bool,
}

/// Return `FormatInfo` for the given filename.
pub fn guess_bcf_format(filename: &str) -> FormatInfo {
    if filename.ends_with(".bcf") {
        FormatInfo {
            format: bcf::Format::BCF,
            uncompressed: false,
        }
    } else if filename.ends_with(".vcf.gz") {
        FormatInfo {
            format: bcf::Format::VCF,
            uncompressed: false,
        }
    } else {
        FormatInfo {
            format: bcf::Format::VCF,
            uncompressed: true,
        }
    }
}

/// Build new `bcf::Header`.
pub fn build_vcf_header(template: &bcf::header::HeaderView) -> Result<bcf::Header, Error> {
    let mut header = bcf::Header::new();

    // Copy over sequence records.
    for record in template.header_records() {
        if let bcf::header::HeaderRecord::Contig { key: _, values } = record {
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
    }

    // Fields: ALT, INFO, FORMAT.
    let alts = vec![
        ("DEL", "Deletion"),
        ("DUP", "Duplication"),
        ("INV", "Inversion"),
        ("CNV", "Copy number variant"),
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
        ("cnmops", "1", "Integer", "Called by cnMOPS"),
        ("delly", "1", "Integer", "Called by Delly"),
        ("manta", "1", "Integer", "Called by Manta"),
        ("PR", "1", "Float", "Paired read evidence"),
        ("SR", "1", "Float", "Split read evidence"),
        ("RD", "1", "Float", "Read depth evidence"),
        ("BF", "1", "Float", "B allele frequency evidence"),
        ("VL", "1", "Integer", "SNV count left of CNV region"),
        ("VM", "1", "Integer", "SNV count within CNV region"),
        ("VR", "1", "Integer", "SNV count right of CNV region"),
        ("ROH", "1", "Integer", "Run of homozygosity"),
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

pub fn collect_contigs(reader: &bcf::IndexedReader) -> Result<Vec<String>, Error> {
    let mut result: Vec<String> = Vec::new();
    let header = reader.header();
    for rid in 0..header.contig_count() {
        result.push(std::str::from_utf8(header.rid2name(rid)?)?.to_string());
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use matches::assert_matches;

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
