/// bam-scan -- Scan BAM file for discordant and clipped reads
use std::fs;
use std::path::Path;
use std::str;

use bloom::{BloomFilter, ASMS};
use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, LevelFilter};
use rust_htslib::{bam, bam::Read};

use lib_common::bam::guess_bam_format;
use lib_common::bam::library::{estimate_library_insert_size, is_interesting, LibraryProperties};
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

struct ExtractionState {
    buffer: bam::RecordBuffer,
    writer: bam::Writer,
    bloom: BloomFilter,
}

/// Extract reads from the buffer in the current window.
fn extract_reads_from_current_window(
    progress_bar: &Option<ProgressBar>,
    state: &mut ExtractionState,
    window_end: i64,
    config: &Config,
    lib_properties: &LibraryProperties,
    pass: i32,
) -> Result<(), Error> {
    let mut pos = 0;
    for record in state.buffer.iter() {
        pos = record.pos();
        if pos > window_end {
            break;
        }
        if state.bloom.contains(&record.qname())
            || is_interesting(&record, &lib_properties, &config)
        {
            debug!(
                "Writing {}/{}",
                str::from_utf8(&record.qname()).unwrap(),
                record.is_first_in_template()
            );
            state.writer.write(&record)?;
            state.bloom.insert(&record.qname());
        }
    }

    if let Some(prog_bar) = progress_bar {
        if pass == 2 && pos >= 0 {
            prog_bar.set_position((pos / 1_000) as u64);
        }
    }

    Ok(())
}

/// Extract reads from the path to the BAM file in the options (path_input) to the output BAM file.
fn extract_reads(
    options: &Options,
    config: &Config,
    lib_properties: &LibraryProperties,
) -> Result<(), Error> {
    // We will fetch records from window of size buffer_length into a bam::RecordBuffer and shift
    // the window with a span of of window_overlap.  This will ensure that all records with a
    // normal insert size are within the same window.
    //
    // We will then iterate over this window twice.  First, we capture the read names of all
    // reads showing a PE/SR signal.  We then iterate a second time.  We collect the reads
    // indicating a PE/SR signal or that were already present in the output file.  The second
    // iteration is needed to capture pairs where the second read shows the split read signal.
    info!("Starting to scan BAM file...");
    let window_size = config.sliding_window_size;
    let window_overlap =
        config.sliding_window_margin + lib_properties.max_normal_isize + lib_properties.max_rlen;

    let mut reader = bam::IndexedReader::from_path(&options.path_input)?;
    let mut header = bam::Header::new();
    for header_line in std::str::from_utf8(reader.header().as_bytes())
        .unwrap()
        .split('\n')
    {
        if !header_line.is_empty() {
            let header_line = String::from(header_line);
            if header_line.starts_with("@HD") {
                header.push_record(&bam::header::HeaderRecord::new(b"HD\tVN:1.6\tSO:unordered"));
            } else if header_line.starts_with("@CO") {
                header.push_comment(header_line[4..].as_bytes());
            } else {
                header.push_record(&bam::header::HeaderRecord::new(header_line[1..].as_bytes()));
            }
        }
    }

    let header_view = bam::HeaderView::from_header(&header);
    let target_count = header_view.target_count() as usize;
    let mut writer = bam::Writer::from_path(
        &options.path_output,
        &header,
        guess_bam_format(&options.path_output),
    )?;
    if config.htslib_io_threads > 0 {
        reader.set_threads(config.htslib_io_threads)?;
        writer.set_threads(config.htslib_io_threads)?;
    }

    let mut state = ExtractionState {
        writer,
        buffer: bam::RecordBuffer::new(reader, true),
        bloom: BloomFilter::with_rate(
            config.bloom_false_positive_rate,
            config.bloom_expected_read_count,
        ),
    };

    for target_id in 0..target_count {
        let target_name = str::from_utf8(header_view.target_names()[target_id]).unwrap();
        if options.verbosity > 0 {
            info!(
                "Starting to scan target #{}/{}: {}",
                target_id + 1,
                target_count,
                &target_name
            );
        }
        let target_len = header_view.target_len(target_id as u32).unwrap() as i64;
        let window_count = (target_len + window_size - 1) / window_size;
        let mut window_no = 0;

        let progress_bar = if options.verbosity == 0 {
            let prog_bar = ProgressBar::new((target_len / 1_000) as u64);
            prog_bar.set_style(ProgressStyle::default_bar()
                .template("scanning {msg:.green.bold} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos:>7}/{len:7} Kbp {elapsed}/{eta}")
                .progress_chars("=>-"));
            prog_bar.set_message(&target_name);
            Some(prog_bar)
        } else {
            None
        };

        while window_no < window_count {
            let window_start = if window_no == 0 {
                0
            } else {
                window_no * window_size - window_overlap
            };
            let window_end = (window_no + 1) * window_size;
            let target_name = &header_view.target_names()[target_id];
            if options.verbosity > 0 {
                info!(
                    "Scanning target's window #{}/{}: {}:{}-{}",
                    window_no + 1,
                    window_count,
                    str::from_utf8(target_name).unwrap(),
                    window_start,
                    window_end
                );
            }
            state
                .buffer
                .fetch(target_name, window_start as u64, window_end as u64)?;
            for pass in 1..=2 {
                extract_reads_from_current_window(
                    &progress_bar,
                    &mut state,
                    window_end,
                    &config,
                    &lib_properties,
                    pass,
                )?;
            }
            window_no += 1;
        }

        if let Some(prog_bar) = progress_bar {
            prog_bar.finish();
        }
    }

    info!("Done scanning BAM file...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("snappysv-scanbam")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Scan BAM file for discordant and clipped reads")
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

    // Estimate the library insert size.
    let lib_properties = estimate_library_insert_size(&options.path_input, &config)?;
    // Run the extraction algorithm.
    extract_reads(&options, &config, &lib_properties)?;

    info!("All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `extract_reads()` and compares the result.
    fn _extract_reads_and_test(
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
        let library_properties = super::LibraryProperties {
            max_rlen: 100,
            median_isize: 300.0,
            std_dev_isize: 10.0,
            max_normal_isize: 330,
        };

        super::extract_reads(&options, &config, &library_properties)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }

    #[test]
    fn identify_pairs_by_hard_clipping_leading() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-hard-leading.sorted.bam",
            "./src/tests/data/ex-pe-hard-leading.expected.sam",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_hard_clipping_trailing() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-hard-trailing.sorted.bam",
            "./src/tests/data/ex-pe-hard-trailing.expected.sam",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_soft_clipping_negative() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-soft-neg.sorted.bam",
            "./src/tests/data/ex-pe-soft-neg.expected.sam",
        )?;
        Ok(())
    }
    #[test]
    fn identify_pairs_by_soft_clipping_positive() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-soft-pos.sorted.bam",
            "./src/tests/data/ex-pe-soft-pos.expected.sam",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_different_tid() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-tid.sorted.bam",
            "./src/tests/data/ex-pe-tid.expected.sam",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_large_tlen() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-tlen.sorted.bam",
            "./src/tests/data/ex-pe-tlen.expected.sam",
        )?;
        Ok(())
    }

    #[test]
    fn identify_pairs_by_orientation() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _extract_reads_and_test(
            &tmp_dir,
            "./src/tests/data/ex-pe-orient.sorted.bam",
            "./src/tests/data/ex-pe-orient.expected.sam",
        )?;
        Ok(())
    }
}
