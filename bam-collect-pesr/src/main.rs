/// bam-collect-pesr -- Collect paired end and split read evidence from BAM.
use std::fs;
use std::path::Path;

use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, LevelFilter};
use rust_htslib::{bam, bam::Read};

use lib_common::bam::library::{
    estimate_library_insert_size, is_discordant_pair, is_split_read_left, is_split_read_right,
    LibraryProperties,
};
use lib_common::error::Error;
use lib_common::read_evidence;
use lib_config::Config;

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

fn extract_evidence(
    record: &bam::Record,
    reader: &bam::IndexedReader,
    config: &Config,
    lib_properties: &LibraryProperties,
) -> Result<Vec<read_evidence::Record>, Error> {
    let cigar = record.cigar();
    let mut result = Vec::new();

    let cl = is_split_read_left(record, &cigar, config);
    let cr = is_split_read_right(record, &cigar, config);

    if cl || cr {
        result.push(read_evidence::Record::SplitRead {
            read_name: std::str::from_utf8(record.qname())?.to_string(),
            is_first: record.is_first_in_template(),
            contig: std::str::from_utf8(reader.header().tid2name(record.tid() as u32))?.to_string(),
            start: record.pos(),
            end: record.cigar().end_pos(),
            clipped_sides: match (cl, cr) {
                (true, true) => read_evidence::Sides::Both,
                (true, false) => read_evidence::Sides::Left,
                (false, true) => read_evidence::Sides::Right,
                _ => panic!("clipped record not clipped?"),
            },
        })
    }

    if is_discordant_pair(record, lib_properties) {
        result.push(read_evidence::Record::PairedRead {
            read_name: std::str::from_utf8(record.qname())?.to_string(),
            is_first1: record.is_first_in_template(),
            contig1: std::str::from_utf8(reader.header().tid2name(record.tid() as u32))?
                .to_string(),
            start1: record.pos(),
            end1: record.cigar().end_pos(),
            strand1: if record.is_reverse() {
                read_evidence::Strand::Forward
            } else {
                read_evidence::Strand::Reverse
            },
            contig2: if record.is_paired() && record.mtid() >= 0 {
                Some(
                    std::str::from_utf8(reader.header().tid2name(record.mtid() as u32))?
                        .to_string(),
                )
            } else {
                None
            },
            start2: if record.is_paired() {
                Some(record.mpos())
            } else {
                None
            },
            strand2: if record.is_paired() {
                Some(if record.is_mate_reverse() {
                    read_evidence::Strand::Forward
                } else {
                    read_evidence::Strand::Reverse
                })
            } else {
                None
            },
            tlen: if record.is_paired() && record.insert_size() >= 0 {
                Some(record.insert_size())
            } else {
                None
            },
        })
    }

    Ok(result)
}

/// Perform extraction of paired read/split read signal.
fn perform_collection(
    options: &Options,
    config: &Config,
    lib_properties: &LibraryProperties,
) -> Result<(), Error> {
    // We can handle the paired read/split read signal collection with a single scan.
    info!("Starting to scan BAM file...");

    // We will scan the BAM file contig wise.
    let mut reader = bam::IndexedReader::from_path(&options.path_input)?;
    let target_count = reader.header().target_count() as usize;
    if config.htslib_io_threads > 0 {
        reader.set_threads(config.htslib_io_threads)?;
    }

    // Evidence is written to an extended BED3 file.
    let mut writer = read_evidence::Writer::from_path(&options.path_output)?;

    for target_id in 0..target_count {
        let target_name = std::str::from_utf8(reader.header().target_names()[target_id])
            .unwrap()
            .to_string();
        if options.verbosity > 0 {
            info!(
                "Starting to scan target #{}/{}: {}",
                target_id + 1,
                target_count,
                &target_name
            );
        }
        let target_len = reader.header().target_len(target_id as u32).unwrap() as i64;

        reader.fetch(target_id as u32, 0, target_len as u64)?;

        let progress_bar =
            if options.verbosity == 0 {
                let prog_bar = ProgressBar::new((target_len / 1_000) as u64);
                prog_bar.set_style(ProgressStyle::default_bar()
                .template(
                    "scanning {msg:.green.bold} [{elapsed_precise}] [{wide_bar:.cyan/blue}] \
                    {pos:>7}/{len:7} Kbp {elapsed}/{eta}")
                .progress_chars("=>-"));
                prog_bar.set_message(&target_name);
                Some(prog_bar)
            } else {
                None
            };

        let mut buffer = bam::Record::new();
        let mut counter: usize = 0;
        loop {
            if !reader.read(&mut buffer)? {
                break;
            } else {
                for evidence in extract_evidence(&buffer, &reader, &config, &lib_properties)? {
                    writer.write(&evidence)?;
                }

                counter += 1;
                if counter % 10_000 == 0 {
                    if let Some(prog_bar) = &progress_bar {
                        prog_bar.set_position(buffer.pos() as u64 / 10000);
                    }
                }
            }
        }

        if let Some(prog_bar) = &progress_bar {
            prog_bar.finish();
        }
    }
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("snappysv-bam-collect-pesr")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Collect paired end and split read evidence from BAM")
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
    info!("Starting snappysv-bam-collect-pesr");
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

    // Estimate the library insert size.
    let lib_properties = estimate_library_insert_size(&options.path_input, &config)?;
    // Run the collection algorithm.
    perform_collection(&options, &config, &lib_properties)?;

    info!("All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `perform_collection()` and compares the result.
    fn _perform_collection_and_test(
        tmp_dir: &TempDir,
        path_input: &str,
        path_expected: &str,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.tsv").to_str().unwrap());
        let options = super::Options {
            verbosity: 1, // disable progress bar
            path_config: None,
            path_input: String::from(path_input),
            path_output: path_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str("").unwrap();
        let library_properties = super::LibraryProperties {
            max_rlen: 100,
            median_isize: 300.0,
            std_dev_isize: 10.0,
            max_normal_isize: 330,
        };

        super::perform_collection(&options, &config, &library_properties)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }

    #[test]
    fn identify_pairs_by_hard_clipping_leading() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-hard-leading.sorted.bam",
            "./src/tests/data/ex-pe-hard-leading.expected.tsv",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_hard_clipping_trailing() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-hard-trailing.sorted.bam",
            "./src/tests/data/ex-pe-hard-trailing.expected.tsv",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_soft_clipping_negative() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-soft-neg.sorted.bam",
            "./src/tests/data/ex-pe-soft-neg.expected.tsv",
        )?;
        Ok(())
    }
    #[test]
    fn identify_pairs_by_soft_clipping_positive() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-soft-pos.sorted.bam",
            "./src/tests/data/ex-pe-soft-pos.expected.tsv",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_different_tid() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-tid.sorted.bam",
            "./src/tests/data/ex-pe-tid.expected.tsv",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_large_tlen() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-tlen.sorted.bam",
            "./src/tests/data/ex-pe-tlen.expected.tsv",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_orientation() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-orient.sorted.bam",
            "./src/tests/data/ex-pe-orient.expected.tsv",
        )?;
        Ok(())
    }
}
