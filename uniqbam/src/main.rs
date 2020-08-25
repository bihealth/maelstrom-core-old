/// uniqbam -- Remove duplicate rows from BAM files (as created by scanbam)
use std::fs;
use std::path::Path;

use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use log::{debug, info, LevelFilter};
use rust_htslib::bam;

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
    #[error("problem with BAM file access")]
    Htslib {
        #[from]
        source: bam::errors::Error, // TODO: add experimental backtrace feature?
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

/// Main entry point after parsing command line and loading options.
fn perform_filtration(_options: &Options, _config: &Config) -> Result<(), bam::errors::Error> {
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
    info!("Starting snappysv-scanbam");
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
mod tests {}
