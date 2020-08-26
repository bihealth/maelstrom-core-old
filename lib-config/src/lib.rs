/// lib-config -- shared configuration.
use serde::Deserialize;

/// Clustering configuration for clusvcf.
#[derive(Deserialize, Debug, Clone)]
pub struct ClusterSettings {
    /// Reciprocal overlap by size.
    pub reciprocal_overlap: f32,
    /// Maximal breakpoint distance.
    pub max_bp_distance: Option<i64>,
    /// Require strands to match.
    pub match_strands: bool,
    /// Require sv_type to match.
    pub match_sv_type: bool,
    /// Sample overlap to require, if any.
    pub sample_overlap: Option<f32>,
}

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
    0.0001
}

fn default_bloom_expected_read_count() -> u32 {
    100_000_000
}

fn default_min_clipped_bases() -> i64 {
    20
}

fn default_supplementary_masked_as_secondary() -> bool {
    true
}

fn default_htslib_io_threads() -> usize {
    0
}

fn default_stdvcf_apply_filters() -> bool {
    true
}

fn default_clusvcf_presets_per_tool_pe_sr() -> ClusterSettings {
    ClusterSettings {
        reciprocal_overlap: 0.1,
        max_bp_distance: Some(300),
        match_strands: true,
        match_sv_type: true,
        sample_overlap: None,
    }
}

fn default_clusvcf_presets_per_tool_doc() -> ClusterSettings {
    ClusterSettings {
        reciprocal_overlap: 0.8,
        max_bp_distance: None,
        match_strands: true,
        match_sv_type: true,
        sample_overlap: None,
    }
}

/// Program configuration, from config file.
#[derive(Deserialize, Debug, Clone)]
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

    /// Whether or not to interpret FILTER values in stdvcf.
    #[serde(default = "default_stdvcf_apply_filters")]
    pub stdvcf_apply_filters: bool,

    /// Cluster preset for per-tool (PE/SR) aggregation.
    #[serde(default = "default_clusvcf_presets_per_tool_pe_sr")]
    pub clusvcf_presets_per_tool_pe_sr: ClusterSettings,

    /// Cluster preset for per-tool (DoC) aggregation.
    #[serde(default = "default_clusvcf_presets_per_tool_doc")]
    pub clusvcf_presets_per_tool_doc: ClusterSettings,
}
