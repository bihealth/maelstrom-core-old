use std::fs;

use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use log::{debug, info, warn, LevelFilter};
use rust_htslib::{bam, bam::Read};
use serde::Deserialize;
use thiserror::Error;

fn default_lib_estimation_sample_size() -> usize {
    100_000
}

fn default_library_cutoff_deviation() -> f64 {
    0.1
}

fn default_library_cutoff_sd_mult() -> f64 {
    7.0
}

fn default_lib_estimation_sd_mult() -> f64 {
    3.0
}
/// Program configuration, from config file.
#[derive(Deserialize, Debug)]
struct Config {
    /// Number of records to read for estimating the library length.
    #[serde(default = "default_lib_estimation_sample_size")]
    lib_estimation_sample_size: usize,

    /// Library deviation to use for cut-off.
    #[serde(default = "default_library_cutoff_deviation")]
    library_cutoff_deviation: f64,

    /// Standard deviation cutoff.
    #[serde(default = "default_library_cutoff_sd_mult")]
    library_cutoff_sd_mult: f64,

    /// Standard deviation multiplicator
    #[serde(default = "default_lib_estimation_sd_mult")]
    lib_estimation_sd_mult: f64,
}

/// Global error type.
#[derive(Error, Debug)]
enum ScanbamError {
    /// A command line option is missing.
    #[error("missing command line argument")]
    OptionMissing(),
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
}

impl Options {
    pub fn from_arg_matches<'a>(matches: &ArgMatches<'a>) -> Result<Self, ScanbamError> {
        Ok(Options {
            verbosity: matches.occurrences_of("v"),
            path_config: matches.value_of("config").map(|s| s.to_string()),
            path_input: match matches.value_of("input") {
                Some(x) => String::from(x),
                None => return Err(ScanbamError::OptionMissing()),
            },
            path_output: match matches.value_of("output") {
                Some(x) => String::from(x),
                None => return Err(ScanbamError::OptionMissing()),
            },
        })
    }
}

/// Library properties.
#[derive(Debug)]
struct LibraryProperties {
    /// Median insert size
    median_isize: f64,
    /// Inset size standard deviation
    std_dev_isize: f64,
    /// Maximal normal insert size
    max_normal_isize: f64,
}

/// Estimate the library insert size.
///
/// Current main limitation: PE read, no artifact filter.
fn estimate_library_insert_size(
    options: &Options,
    config: &Config,
) -> Result<LibraryProperties, bam::errors::Error> {
    info!(
        "reading {} records to estimate insert size...",
        config.lib_estimation_sample_size
    );

    let mut reader = bam::Reader::from_path(&options.path_input)?;
    let mut insert_sizes = Vec::new();
    for r in reader.records() {
        let record = r?;
        if record.is_paired()
            && record.is_proper_pair()
            && record.insert_size() > 0
            && record.tid() == record.mtid()
            && !record.is_unmapped()
            && !record.is_mate_unmapped()
            && !record.is_secondary()
            && !record.is_supplementary()
            && !record.is_quality_check_failed()
            && !record.is_duplicate()
        {
            insert_sizes.push(record.insert_size());
            if insert_sizes.len() > config.lib_estimation_sample_size {
                break;
            }
        }
    }

    if insert_sizes.len() == 0 {
        panic!("Found no reads in input file!");
    } else if insert_sizes.len() < config.lib_estimation_sample_size {
        warn!(
            "Only found {} records instead of {}",
            insert_sizes.len(),
            config.lib_estimation_sample_size
        );
    }

    insert_sizes.sort();
    let median: f64 = insert_sizes[insert_sizes.len() / 2] as f64;
    let delta: f64 = config.library_cutoff_deviation * config.library_cutoff_sd_mult * median;
    let cutoff_max: f64 = median + delta;
    let cutoff_min: f64 = median - delta;
    let cutoff_min: f64 = if (cutoff_min < 0.0) || (cutoff_max < cutoff_min) {
        0.0
    } else {
        cutoff_min
    };

    let mut count = 0;
    let mut variance: f64 = 0.0;
    for i in insert_sizes {
        let i = i as f64;
        if i >= cutoff_min && i <= cutoff_max {
            variance += (i - median) * (i - median);
            count += 1;
        }
    }
    let std_dev = (variance / (count as f64)).sqrt();
    let max_normal = median + config.lib_estimation_sd_mult * std_dev;

    let result = LibraryProperties {
        median_isize: median,
        std_dev_isize: std_dev,
        max_normal_isize: max_normal,
    };

    info!("library properties: {:?}", &result);

    Ok(result)
}

fn main() -> Result<(), ScanbamError> {
    // Setup command line parser and parse options.
    let matches = App::new("snappysv-scanbam")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Scan BAM file for discordant and clipped reads")
        .args(&[
            Arg::from_usage("-v... 'Increase verbosity'"),
            Arg::from_usage("-c, --config=[FILE] 'Sets a custom config file'"),
            Arg::from_usage("<input> 'input file to read from'"),
            Arg::from_usage("<output> 'output file to write to'"),
        ])
        .get_matches();
    let options = Options::from_arg_matches(&matches)?;

    // Setup logging verbosity.
    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
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

    // Run the actual code.
    let lib_properties = estimate_library_insert_size(&options, &config)?;

    Ok(())
}
