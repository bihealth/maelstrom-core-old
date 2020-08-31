use super::super::error::Error;

use log::{info, warn};
use rust_htslib::{bam, bam::Read};

use lib_config::Config;

/// Library properties.
#[derive(Debug)]
pub struct LibraryProperties {
    /// Maximal read length.
    pub max_rlen: i64,
    /// Median insert size
    pub median_isize: f64,
    /// Inset size standard deviation
    pub std_dev_isize: f64,
    /// Maximal normal insert size
    pub max_normal_isize: i64,
}

/// Estimate the library insert size.
///
/// Current main limitation: PE read, no artifact filter.
pub fn estimate_library_insert_size(
    path_input: &str,
    config: &Config,
) -> Result<LibraryProperties, Error> {
    info!(
        "reading {} records to estimate insert size...",
        config.lib_estimation_sample_size
    );

    let mut reader = bam::Reader::from_path(&path_input)?;
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

    if insert_sizes.is_empty() {
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

/// Return if split read clipped on left side.
pub fn is_split_read_left(
    record: &bam::Record,
    cigar: &bam::record::CigarStringView,
    config: &Config,
) -> bool {
    ((record.is_supplementary()
        || (config.supplementary_masked_as_secondary && record.is_secondary()))
        && (cigar.leading_hardclips() >= config.min_clipped_bases))
        || (!record.is_supplementary()
            && !record.is_secondary()
            && (cigar.leading_hardclips() >= config.min_clipped_bases
                || cigar.leading_softclips() >= config.min_clipped_bases))
}

/// Return if split read clipped on right side.
pub fn is_split_read_right(
    record: &bam::Record,
    cigar: &bam::record::CigarStringView,
    config: &Config,
) -> bool {
    ((record.is_supplementary()
        || (config.supplementary_masked_as_secondary && record.is_secondary()))
        && (cigar.trailing_hardclips() >= config.min_clipped_bases))
        || (!record.is_supplementary()
            && !record.is_secondary()
            && (cigar.trailing_hardclips() >= config.min_clipped_bases
                || cigar.trailing_softclips() >= config.min_clipped_bases))
}

/// Return if pair is discordant.
pub fn is_discordant_pair(record: &bam::Record, lib_properties: &LibraryProperties) -> bool {
    record.is_paired() && (record.tid() >= 0 && record.mtid() >= 0 && record.tid() != record.mtid())
        || record.insert_size().abs() > lib_properties.max_normal_isize
        || record.is_reverse() == record.is_mate_reverse()
}

/// Determine whether the record shows PE or SR signal.
pub fn is_interesting(
    record: &bam::Record,
    lib_properties: &LibraryProperties,
    config: &Config,
) -> bool {
    // We need to extract the CIGAR information for split read analysis.
    let cigar = &record.cigar_cached().unwrap();

    // Early exit if not interesting at all.
    if record.is_quality_check_failed()
        || record.is_duplicate()
        || record.is_unmapped()
        || record.is_mate_unmapped()
        || (!config.supplementary_masked_as_secondary && record.is_secondary())
    {
        return false;
    }

    is_split_read_left(&record, &cigar, &config)
        || is_split_read_right(&record, &cigar, &config)
        || is_discordant_pair(record, lib_properties)
}
