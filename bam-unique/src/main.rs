/// bam-unique -- Remove duplicate rows from BAM files (as created by scanbam)
use std::collections::HashSet;
use std::fs;
use std::path::Path;

use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use log::{debug, info, LevelFilter};
use rust_htslib::{bam, bam::Read};

use lib_common::bam::guess_bam_format;
use lib_common::error::Error;
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

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
struct RecordIdentifier {
    first: bool,
    tid: i32,
    pos: i64,
}

impl RecordIdentifier {
    fn from_record(record: &bam::Record) -> Self {
        Self {
            first: record.is_first_in_template(),
            tid: record.tid(),
            pos: record.pos(),
        }
    }
}

/// Write unique reads from chunk into the writer.
fn write_unique(chunk: &[bam::Record], writer: &mut bam::Writer) -> Result<(), bam::errors::Error> {
    let mut seen = HashSet::new();

    for record in chunk.iter() {
        let record_id = RecordIdentifier::from_record(record);
        if !seen.contains(&record_id) {
            writer.write(&record)?;
            seen.insert(record_id);
        }
    }

    Ok(())
}

/// Main entry point after parsing command line and loading options.
fn perform_filtration(options: &Options, config: &Config) -> Result<(), bam::errors::Error> {
    info!("Starting to scan BAM file...");

    let mut reader = bam::Reader::from_path(&options.path_input)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(
        &options.path_output,
        &header,
        guess_bam_format(&options.path_output),
    )?;
    if config.htslib_io_threads > 0 {
        reader.set_threads(config.htslib_io_threads)?;
        writer.set_threads(config.htslib_io_threads)?;
    }

    let mut chunk: Vec<bam::Record> = Vec::new();
    let mut buffer = bam::Record::new();
    loop {
        if !reader.read(&mut buffer)? {
            break;
        } else if chunk.is_empty() {
            chunk.push(buffer.clone());
        } else {
            if buffer.qname() != chunk[0].qname() {
                write_unique(&chunk, &mut writer)?;
                chunk.clear();
            }
            chunk.push(buffer.clone());
        }
    }
    if !chunk.is_empty() {
        write_unique(&chunk, &mut writer)?;
    }

    info!("Done scanning BAM file...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("snappysv-uniqbam")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Remove duplicate rows from BAM files (as created by scanbam)")
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
    info!("Starting snappysv-uniqbam");
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

    // Perform the block-wise filtration.
    perform_filtration(&options, &config)?;

    info!("All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `perform_filtration()` and compares the result.
    fn _perform_filtration_and_test(
        tmp_dir: &TempDir,
        path_input: &str,
        path_expected: &str,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.sam").to_str().unwrap());
        let options = super::Options {
            verbosity: 1, // disable progress bar
            path_config: None,
            path_input: String::from(path_input),
            path_output: path_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str("").unwrap();

        super::perform_filtration(&options, &config)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }

    #[test]
    fn perform_filtration_one_pair() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_filtration_and_test(
            &tmp_dir,
            "./src/tests/data/ex-duplicates.bam",
            "./src/tests/data/ex-duplicates.expected.sam",
        )?;
        Ok(())
    }

    #[test]
    fn perform_filtration_two_pairs() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_filtration_and_test(
            &tmp_dir,
            "./src/tests/data/ex-duplicates2.bam",
            "./src/tests/data/ex-duplicates2.expected.sam",
        )?;
        Ok(())
    }
}
