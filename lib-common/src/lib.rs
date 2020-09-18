/// lib-common -- shared functionality
pub mod bam;
pub mod bcf;
pub mod error;
pub mod read_evidence;
pub mod stats;
pub mod sv;
use log::info;

use core::hash::Hash;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::genome::Interval;
use bio_types::strand::NoStrand;
use regex::Regex;

use error::Error;

/// Enumeration for calling algorithms.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Algorithm {
    /// Delly2
    Delly,
    /// Manta
    Manta,
    // cnMOPS
    CNMOPS,
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

/// Load BED file file and retur as AnnotMap
pub fn bed_to_annot_map<R>(
    path: &str,
    contig_map: &HashMap<String, R>,
) -> Result<AnnotMap<R, ()>, Error>
where
    R: Hash + Eq + Clone + std::borrow::ToOwned<Owned = R>,
{
    info!("Loading BED file {}", &path);
    let mut result = AnnotMap::new();

    let mut count = 0;
    for line in io::BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let arr: Vec<&str> = line.split('\t').collect();
        if arr.len() < 3 {
            return Err(error::Error::InvalidBEDFile(format!(
                "Unexpected number of fields in {} (must have >= 3)",
                &line
            )));
        }

        let start = arr[1].parse::<isize>()?;
        let end = arr[2].parse::<isize>()?;

        #[allow(clippy::or_fun_call)]
        let rid = contig_map
            .get(arr[0])
            .ok_or(Error::InvalidBEDFile(format!("Unknown contig: {}", arr[0])))?;

        let loc = Contig::new(
            rid.to_owned(),
            start,
            (end - start) as usize,
            NoStrand::Unknown,
        );
        result.insert_at((), &loc);
        count += 1;
    }
    info!("=> {} entries", count);

    Ok(result)
}
