/// lib-config -- shared configuration.
use serde::Deserialize;

fn default_min_mapq() -> Option<u8> {
    Some(0)
}
fn default_min_unclipped() -> Option<f32> {
    Some(0.6)
}
fn default_count_kind() -> String {
    String::from("coverage")
}
fn default_window_length() -> usize {
    100
}

/// Configuration for bam-collect-doc.
#[derive(Deserialize, Debug, Clone)]
pub struct DepthOfCoverageConfig {
    /// Minimal MAPQ value for an alignment to be considered.
    #[serde(default = "default_min_mapq")]
    pub min_mapq: Option<u8>,
    /// Minimal proportion of unclipped bases for an alignment to be considered.
    #[serde(default = "default_min_unclipped")]
    pub min_unclipped: Option<f32>,
    /// The counts to generate, one of "coverage", and "fragments".
    #[serde(default = "default_count_kind")]
    pub count_kind: String,
    /// The window length,
    #[serde(default = "default_window_length")]
    pub window_length: usize,
}

fn default_reciprocal_overlap() -> f32 {
    0.8
}
fn default_max_bp_distance() -> Option<i64> {
    Some(300)
}
fn default_match_strands() -> bool {
    true
}
fn default_match_sv_type() -> bool {
    true
}

fn default_sv_type_out() -> Option<String> {
    None
}

fn default_sample_overlap() -> Option<f32> {
    Some(0.5)
}

/// Clustering configuration for vcf-cluster.
#[derive(Deserialize, Debug, Clone)]
pub struct ClusterSettings {
    /// Reciprocal overlap by size.
    #[serde(default = "default_reciprocal_overlap")]
    pub reciprocal_overlap: f32,
    /// Maximal breakpoint distance.
    #[serde(default = "default_max_bp_distance")]
    pub max_bp_distance: Option<i64>,
    /// Require strands to match.
    #[serde(default = "default_match_strands")]
    pub match_strands: bool,
    /// Require sv_type to match.
    #[serde(default = "default_match_sv_type")]
    pub match_sv_type: bool,
    /// The SV type to write out when ignoring it.
    #[serde(default = "default_sv_type_out")]
    pub sv_type_out: Option<String>,
    /// Sample overlap to require, if any.
    #[serde(default = "default_sample_overlap")]
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

fn default_vcf_cluster_presets_per_tool_pesr() -> ClusterSettings {
    ClusterSettings {
        reciprocal_overlap: 0.1,
        max_bp_distance: Some(300),
        match_strands: true,
        match_sv_type: true,
        sv_type_out: None,
        sample_overlap: None,
    }
}

fn default_vcf_cluster_presets_per_tool_doc() -> ClusterSettings {
    ClusterSettings {
        reciprocal_overlap: 0.8,
        max_bp_distance: None,
        match_strands: false,
        match_sv_type: false,
        sv_type_out: Some("CNV".to_string()),
        sample_overlap: None,
    }
}

fn default_collect_doc_config() -> DepthOfCoverageConfig {
    DepthOfCoverageConfig {
        min_mapq: default_min_mapq(),
        min_unclipped: default_min_unclipped(),
        window_length: default_window_length(),
        count_kind: default_count_kind(),
    }
}

fn default_annotate_read_evidence_max_dist() -> i64 {
    1_000
}

fn default_blocked_regions_bed() -> Option<String> {
    None
}

fn default_doc_annotation_min_bins() -> usize {
    10
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
    #[serde(default = "default_vcf_cluster_presets_per_tool_pesr")]
    pub vcf_cluster_presets_per_tool_pesr: ClusterSettings,

    /// Cluster preset for per-tool (DoC) aggregation.
    #[serde(default = "default_vcf_cluster_presets_per_tool_doc")]
    pub vcf_cluster_presets_per_tool_doc: ClusterSettings,

    /// Preset for depth of coverage extraction.
    #[serde(default = "default_collect_doc_config")]
    pub collect_doc_config: DepthOfCoverageConfig,

    /// Number of bases to look for read evidence to annotate.
    #[serde(default = "default_annotate_read_evidence_max_dist")]
    pub annotate_read_evidence_max_dist: i64,

    /// Optionally, a BED file with blocked regions.
    #[serde(default = "default_blocked_regions_bed")]
    pub blocked_regions_bed: Option<String>,

    /// Minimal number of bins for DoC annotation.
    #[serde(default = "default_doc_annotation_min_bins")]
    pub doc_annotation_min_bins: usize,
}
