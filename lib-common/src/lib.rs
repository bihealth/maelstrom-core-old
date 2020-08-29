/// lib-common -- shared functionality
pub mod bam;
pub mod bcf;
pub mod error;
pub mod read_evidence;
pub mod stats;
pub mod sv;

/// Enumeration for calling algorithms.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Algorithm {
    /// Delly2
    Delly,
}
