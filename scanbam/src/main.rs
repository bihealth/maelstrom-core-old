/// scanbam -- Scan BAM file for discordant and clipped reads
use std::fs;
use std::path::Path;
use std::str;

use bloom::{BloomFilter, ASMS};
use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn, LevelFilter};
use rust_htslib::{bam, bam::Read};

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

/// Library properties.
#[derive(Debug)]
struct LibraryProperties {
    /// Maximal read length.
    max_rlen: i64,
    /// Median insert size
    median_isize: f64,
    /// Inset size standard deviation
    std_dev_isize: f64,
    /// Maximal normal insert size
    max_normal_isize: i64,
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
        max_rlen: 150, // TODO: fixme
        median_isize: median,
        std_dev_isize: std_dev,
        max_normal_isize: max_normal.ceil() as i64,
    };

    info!("library properties: {:?}", &result);

    Ok(result)
}

/// Determine whether the record shows PE or SR signal.
fn is_interesting(
    record: &bam::Record,
    lib_properties: &LibraryProperties,
    config: &Config,
) -> bool {
    // We need to extract the CIGAR information for split read analysis.
    let cigar = &record.cigar_cached().unwrap();

    // Check for "split read alignment" that is marked as supplementary by BWA-MEM (if
    // configured, allow having supplementary alignments masked as secondary).
    if (record.is_supplementary()
        || (config.supplementary_masked_as_secondary && record.is_secondary()))
        && (cigar.leading_hardclips() > 0 || cigar.trailing_hardclips() > 0)
    {
        return true;
    }

    // Early exit if not interesting for paired end or split read analysis.
    if !record.is_paired()
        || record.is_quality_check_failed()
        || record.is_duplicate()
        || record.is_unmapped()
        || record.is_mate_unmapped()
        || (!config.supplementary_masked_as_secondary && record.is_secondary())
    {
        // Singletons, QC-fail, and duplicate records cannot be properly used, same as pairs with
        // unaligned reads.  In the future, one-end anchored reads might be considered, though.
        return false;
    }

    // Paired read signal is quite easy to identify.
    if (record.tid() >= 0 && record.mtid() >= 0 && record.tid() != record.mtid())
        || record.insert_size() > lib_properties.max_normal_isize
        || record.is_reverse() == record.is_mate_reverse()
    {
        // Discordantly aligning pairs are immediately interesting.
        return true;
    }

    // Check for split-read anchors (alignment is primary and has sufficient number of
    // soft-clipped bases).
    if !record.is_secondary() && !record.is_supplementary() {
        if cigar.leading_softclips() >= config.min_clipped_bases
            || cigar.leading_hardclips() >= config.min_clipped_bases
            || cigar.trailing_softclips() >= config.min_clipped_bases
            || cigar.trailing_hardclips() >= config.min_clipped_bases
        {
            return true;
        }
    }

    // Fall-through => not interesting.
    false
}

/// Extract reads from the buffer in the current window.
fn extract_reads_from_current_window(
    progress_bar: &Option<ProgressBar>,
    buffer: &mut bam::RecordBuffer,
    window_end: i64,
    writer: &mut bam::Writer,
    bloom: &mut BloomFilter,
    config: &Config,
    lib_properties: &LibraryProperties,
    pass: i32,
) -> Result<(), Error> {
    let mut pos = 0;
    for record in buffer.iter() {
        pos = record.pos();
        if pos > window_end {
            break;
        }
        if bloom.contains(&record.qname()) || is_interesting(&record, &lib_properties, &config) {
            println!(
                "Writing {}/{}",
                str::from_utf8(&record.qname()).unwrap(),
                record.is_first_in_template()
            );
            writer.write(&record)?;
            bloom.insert(&record.qname());
        }
    }

    if let Some(bar) = progress_bar {
        if pass == 2 && pos >= 0 {
            bar.set_position((pos / 1_000) as u64);
        }
    }

    Ok(())
}

/// Return `bam::Format` for the given filename.
fn guess_bam_format(filename: &str) -> bam::Format {
    if filename.ends_with(".bam") {
        bam::Format::BAM
    } else {
        bam::Format::SAM
    }
}

/// Extract reads from the path to the BAM file in the options (path_input) to the output BAM file.
fn extract_reads(
    options: &Options,
    config: &Config,
    lib_properties: &LibraryProperties,
) -> Result<(), Error> {
    // Construct bloom filter; used for collecting reads to read.
    let mut bloom = BloomFilter::with_rate(
        config.bloom_false_positive_rate,
        config.bloom_expected_read_count,
    );

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
        .split("\n")
    {
        if header_line.len() > 0 {
            let header_line = String::from(header_line);
            if header_line.starts_with("@HD") {
                header.push_record(&bam::header::HeaderRecord::new(
                    "HD\tVN:1.6\tSO:unordered".as_bytes(),
                ));
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
    let mut buffer = bam::RecordBuffer::new(reader, true);

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
            let bar = ProgressBar::new((target_len / 1_000) as u64);
            bar.set_style(ProgressStyle::default_bar()
                .template("scanning {msg:.green.bold} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos:>7}/{len:7} KB {elapsed}/{eta}")
                .progress_chars("=>-"));
            bar.set_message(&target_name);
            Some(bar)
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
            buffer.fetch(target_name, window_start as u64, window_end as u64)?;
            for pass in 1..=2 {
                extract_reads_from_current_window(
                    &progress_bar,
                    &mut buffer,
                    window_end,
                    &mut writer,
                    &mut bloom,
                    &config,
                    &lib_properties,
                    pass,
                )?;
            }
            window_no += 1;
        }

        if let Some(bar) = progress_bar {
            bar.finish();
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
    let lib_properties = estimate_library_insert_size(&options, &config)?;
    // Run the exraction algorithm.
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
