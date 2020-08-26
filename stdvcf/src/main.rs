/// stdvcf -- Extract and standardize tool output VCF to standardized VCF.
use std::fs;
use std::path::Path;

use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, LevelFilter};
use rust_htslib::{bcf, bcf::Read};

use lib_common::guess_bcf_format;
use lib_config::Config;

/// Global error type.
#[derive(thiserror::Error, Debug)]
enum Error {
    /// A command line option is missing.
    #[error("missing command line argument")]
    OptionMissing(),
    /// The output file already exists.
    #[error("output file already exists")]
    OutputFileExists(),
    /// Problem with file I/O.
    #[error("problem with I/O")]
    Io {
        #[from]
        source: std::io::Error,
        // TODO: add experimental backtrace feature?
    },
    /// Problem with htslib
    #[error("problem with BCF file access")]
    Htslib {
        #[from]
        source: bcf::errors::Error, // TODO: add experimental backtrace feature?
    },
    /// Problem with string conversion
    #[error("problem with string conversion")]
    Utf8Error {
        #[from]
        source: std::str::Utf8Error, // TODO: add experimental backtrace feature?
    },
}

/// Command line options
#[derive(Debug)]
struct Options {
    /// Verbosity level
    verbosity: u64,
    /// Path to configuration file to use,
    path_config: Option<String>,
    /// Path to input file.
    path_input: String,
    /// Path to output file.
    path_output: String,
    /// Overwrite output file.
    overwrite: bool,
}

impl Options {
    pub fn from_arg_matches<'a>(matches: &ArgMatches<'a>) -> Result<Self, Error> {
        Ok(Options {
            verbosity: matches.occurrences_of("v"),
            path_config: matches.value_of("config").map(|s| s.to_string()),
            path_input: match matches.value_of("input") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            path_output: match matches.value_of("output") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            overwrite: matches.occurrences_of("overwrite") > 0,
        })
    }
}

/// Build new `bcf::Header`.
fn build_header(template: &bcf::header::HeaderView) -> Result<bcf::Header, Error> {
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
        ("SVLEN", "1", "String", "SV length"),
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

/// Enumeration for calling algorithms.
enum Algorithm {
    /// Delly2
    Delly,
}

/// Summarize a single record.
fn summarize_record(
    src: &mut bcf::Record,
    src_header: &bcf::header::HeaderView,
    dst: &mut bcf::Record,
    dst_header: &bcf::header::HeaderView,
    config: &Config,
) -> Result<bool, Error> {
    let caller = Algorithm::Delly;

    if config.stdvcf_apply_filters
        && src
            .filters()
            .map(|id| String::from_utf8(src.header().id_to_name(id)).unwrap())
            .any(|f| f != "PASS")
    {
        return Ok(false);
    }

    dst.set_rid(Some(
        dst_header.name2rid(src_header.rid2name(src.rid().unwrap())?)?,
    ));
    dst.set_pos(src.pos());
    dst.set_alleles(&src.alleles())?;

    // INFO/END2
    let end2 = match src.info(b"END2").integer() {
        Ok(Some(ends)) => ends[0],
        _ => match src.info(b"END").integer() {
            Ok(Some(ends)) => ends[0],
            _ => panic!("Could not read END or END2"),
        },
    };
    dst.push_info_integer(b"END2", &[end2])?;
    // INFO/CHR2
    dst.push_info_string(b"CHR2", &src.info(b"CHR2").string()?.unwrap())?;
    // INFO/SVTYPE
    let sv_type = String::from(std::str::from_utf8(
        src.info(b"SVTYPE").string()?.unwrap()[0],
    )?);
    if sv_type == "INS" {
        return Ok(false); // don't care about INS
    }
    dst.push_info_string(b"SVTYPE", &[sv_type.as_bytes()])?;
    // INFO/STRANDS
    let raw_strands = std::str::from_utf8(src.info(b"CT").string()?.unwrap()[0])?;
    let strands = match raw_strands {
        "5to3" => "-+",
        "3to5" => "+-",
        "5to5" => "--",
        "3to3" => "++",
        "NtoN" => {
            if sv_type == "INS" {
                "+-"
            } else {
                panic!("improper strands")
            }
        }
        _ => panic!("improper strands"),
    };
    dst.push_info_string(b"STRANDS", &[strands.as_bytes()])?;
    // INFO/SVLEN
    let svtype = std::str::from_utf8(src.info(b"SVTYPE").string()?.unwrap()[0])?;
    if svtype == "BND" {
        dst.push_info_integer(b"SVLEN", &[-1])?;
    } else {
        dst.push_info_integer(b"SVLEN", &[end2 - (src.pos() as i32)])?;
    }
    // INFO/ALGORITHMS
    dst.push_info_string(b"ALGORITHMS", &["delly".as_bytes()])?;

    // FORMAT/GT
    let sample_count = dst_header.sample_count() as usize;
    let mut gts: Vec<i32> = Vec::new();
    Vec::from(src.format(b"GT").integer()?)[..sample_count]
        .iter()
        .for_each(|xs| gts.extend_from_slice(xs));
    dst.push_format_integer(b"GT", &gts)?;
    // FORMAT/delly
    match caller {
        Algorithm::Delly => {
            // For each sample, collect wehther there is any called allele.
            let called: Vec<i32> = (0..sample_count)
                .map(|i| {
                    let ith_genotype: bcf::record::Genotype =
                        src.genotypes().unwrap().get(i as usize);
                    let any_var: bool = ith_genotype.iter().any(|gt_allele| match gt_allele {
                        bcf::record::GenotypeAllele::Unphased(i)
                        | bcf::record::GenotypeAllele::Phased(i) => (*i != 0),
                        _ => false,
                    });
                    if any_var {
                        1
                    } else {
                        0
                    }
                })
                .collect();
            dst.push_format_integer(b"delly", &called)?;
        }
    }

    Ok(true)
}

/// Main entry point after parsing command line and loading options.
fn perform_extraction(options: &Options, config: &Config) -> Result<(), Error> {
    info!("Starting to extract BCF file...");
    let mut reader = bcf::Reader::from_path(&options.path_input)?;
    let header = build_header(reader.header())?;
    let guessed = guess_bcf_format(&options.path_output);
    let mut writer = bcf::Writer::from_path(
        &options.path_output,
        &header,
        guessed.uncompressed,
        guessed.format,
    )?;

    let spinner_style = ProgressStyle::default_spinner()
        .tick_chars(".oO@* ")
        .template("{prefix:.bold.dim} {spinner} {msg} {elapsed}");
    let spinner = if options.verbosity == 0 {
        let spinner = ProgressBar::new_spinner();
        spinner.set_style(spinner_style.clone());
        spinner.enable_steady_tick(100);
        Some(spinner)
    } else {
        None
    };

    let mut counter = 0;
    let mut buffer_read = reader.empty_record();
    loop {
        if !reader.read(&mut buffer_read)? {
            break; // done
        }
        let mut buffer_write = writer.empty_record();
        if summarize_record(
            &mut buffer_read,
            reader.header(),
            &mut buffer_write,
            writer.header(),
            &config,
        )? {
            writer.write(&buffer_write)?;
        }

        counter += 1;
        if counter % 1_000 == 0 {
            if let Some(spinner) = &spinner {
                spinner.set_message(&format!(
                    "currently at {}:{}",
                    &std::str::from_utf8(reader.header().rid2name(buffer_read.rid().unwrap())?)?,
                    buffer_read.pos() + 1
                ));
            }
        }
    }
    if let Some(spinner) = &spinner {
        spinner.finish();
    }

    info!("Done extracting BCF file...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("snappysv-stdvcf")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Extract and standardize records from tool VCF files")
        .args(&[
            Arg::from_usage("-v... 'Increase verbosity'"),
            Arg::from_usage("--overwrite 'Allow overwriting of output file'"),
            Arg::from_usage("-c, --config=[FILE] 'Sets a custom config file'"),
            Arg::from_usage("<input> 'input file to read from'"),
            Arg::from_usage("<output> 'output file to write to'"),
        ])
        .get_matches();
    let options = Options::from_arg_matches(&matches)?;

    // Output file must not exist yet.
    if Path::new(&options.path_output).exists() && !options.overwrite {
        return Err(Error::OutputFileExists());
    }

    // Setup logging verbosity.
    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{} [{}] {}",
                chrono::Local::now().format("[%Y-%m-%d %H:%M:%S]"),
                record.level(),
                message
            ))
        })
        .level(if matches.is_present("verbose") {
            LevelFilter::Debug
        } else {
            LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();
    info!("Starting snappysv-stdvcf");
    info!("options: {:?}", &options);

    // Parse further settings from configuration file.
    let config: Config = match &options.path_config {
        None => toml::from_str("").unwrap(),
        Some(path_config) => {
            debug!("Loading config file: {}", &path_config);
            let contents = fs::read_to_string(&path_config)?;
            toml::from_str(&contents).unwrap()
        }
    };
    info!("options: {:?}", &config);

    // Perform the record extraction.
    perform_extraction(&options, &config)?;

    info!("All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `perform_extraction()` and compares the result.
    fn _perform_extraction_and_test(
        tmp_dir: &TempDir,
        path_input: &str,
        path_expected: &str,
        config_str: &str,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());
        let options = super::Options {
            verbosity: 1, // disable progress bar
            path_config: None,
            path_input: String::from(path_input),
            path_output: path_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str(config_str).unwrap();

        super::perform_extraction(&options, &config)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }

    #[test]
    fn test_stdvcf_delly2_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-delly-filter.vcf",
            "./src/tests/data/ex-delly-filter.expected.vcf",
            "stdvcf_apply_filters = true",
        )?;
        Ok(())
    }

    #[test]
    fn test_stdvcf_delly2_nofilter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-delly-nofilter.vcf",
            "./src/tests/data/ex-delly-nofilter.expected.vcf",
            "stdvcf_apply_filters = false",
        )?;
        Ok(())
    }
}
