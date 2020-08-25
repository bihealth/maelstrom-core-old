/// lib-config -- shared configuration.
use serde::Deserialize;

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

fn default_sliding_window_margin() -> i64 {
    200
}

fn default_sliding_window_size() -> i64 {
    10_000
}

fn default_bloom_false_positive_rate() -> f32 {
    return 0.0001;
}

fn default_bloom_expected_read_count() -> u32 {
    return 100_000_000;
}

fn default_min_clipped_bases() -> i64 {
    return 20;
}

fn default_supplementary_masked_as_secondary() -> bool {
    return true;
}

fn default_htslib_io_threads() -> usize {
    return 0;
}

/// Program configuration, from config file.
#[derive(Deserialize, Debug)]
pub struct Config {
    /// Number of records to read for estimating the library length.
    #[serde(default = "default_lib_estimation_sample_size")]
    pub lib_estimation_sample_size: usize,

    /// Library deviation to use for cut-off.
    #[serde(default = "default_library_cutoff_deviation")]
    pub library_cutoff_deviation: f64,

    /// Standard deviation cutoff.
    #[serde(default = "default_library_cutoff_sd_mult")]
    pub library_cutoff_sd_mult: f64,

    /// Standard deviation multiplicator.
    #[serde(default = "default_lib_estimation_sd_mult")]
    pub lib_estimation_sd_mult: f64,

    /// Window margin for fetching.
    #[serde(default = "default_sliding_window_margin")]
    pub sliding_window_margin: i64,

    /// Window size for fetching.
    #[serde(default = "default_sliding_window_size")]
    pub sliding_window_size: i64,

    /// Bloom filter FPR.
    #[serde(default = "default_bloom_false_positive_rate")]
    pub bloom_false_positive_rate: f32,

    /// Expected number of reads; for bloom filter.
    #[serde(default = "default_bloom_expected_read_count")]
    pub bloom_expected_read_count: u32,

    /// Number of clipped bases (and maybe split aligned bases) that are interesting.
    #[serde(default = "default_min_clipped_bases")]
    pub min_clipped_bases: i64,

    /// Whether to expect supplementary alignments to be masked as secondary.  This is/was common
    /// when older Picard versions were used with BWA-MEM that did not support the supplementary
    /// alignment flag yet.
    #[serde(default = "default_supplementary_masked_as_secondary")]
    pub supplementary_masked_as_secondary: bool,

    /// Number of I/O threads to use.
    #[serde(default = "default_htslib_io_threads")]
    pub htslib_io_threads: usize,
}
