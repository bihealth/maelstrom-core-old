/// lib-common -- shared functionality
pub mod bam;
pub mod bcf;
pub mod error;
pub mod read_evidence;
pub mod stats;
pub mod sv;

use bio_types::genome::Interval;
use regex::Regex;

use error::Error;

/// Enumeration for calling algorithms.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Algorithm {
    /// Delly2
    Delly,
    /// Manta
    Manta,
}

/// Parse string into region.
pub fn parse_region(s: &str) -> Result<Interval, Error> {
    let re = Regex::new(r"([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)").unwrap();
    if let Some(c) = re.captures(s) {
        let chrom = c.get(1).unwrap().as_str().to_string();
        let a = c.get(2).unwrap().as_str().parse().unwrap();
        let b = c.get(3).unwrap().as_str().parse().unwrap();
        Ok(Interval::new(chrom, a..b))
    } else {
        Err(Error::InvalidRegion())
    }
}
