/// Analysis of the reference sequence for GC content and gaps.
use std::path::Path;

use bio::io::fasta;
use bio::utils::Text;

use lib_common::error::Error;
use log::{debug, info};

/// Some statistics on the reference.
pub struct ReferenceStats {
    /// GC content of windows.
    pub gc_content: Vec<f32>,
    /// Whether or not the window contains a gap (`N`).
    pub has_gap: Vec<bool>,
}

impl ReferenceStats {
    /// Create a new statistics from a given path.
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        chrom: &str,
        window_length: usize,
    ) -> Result<Self, Error> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => Ok(Self::build(p, chrom, window_length)?),
            _ => Err(Error::InvalidPath()),
        }
    }

    /// Internal builder function.
    fn build(path: &str, chrom: &str, window_length: usize) -> Result<Self, Error> {
        info!("Loading GC content and gap (is-N) status...");

        let seq = Self::load_seq(path, chrom)?;

        let (gc_content, has_gap) = Self::look_at_chars(&seq, window_length)?;

        Ok(ReferenceStats {
            gc_content,
            has_gap,
        })
    }

    /// Load sequence from `path`.
    fn load_seq(path: &str, chrom: &str) -> Result<Text, Error> {
        debug!("Loading reference sequence {}: {}...", path, chrom);

        debug!("Opening indexed FASTA reader");
        let mut ref_reader = fasta::IndexedReader::from_file(&path)?;

        debug!("Reading reference {}", chrom);
        let mut seq = Text::new();
        ref_reader.fetch_all(chrom)?;
        ref_reader.read(&mut seq)?;
        let end = seq.len();
        debug!("Reference seq {} has {} characters", chrom, end);

        Ok(seq)
    }

    /// Perform analysis for GC ratio and "has N" flags.
    fn look_at_chars(seq: &[u8], window_length: usize) -> Result<(Vec<f32>, Vec<bool>), Error> {
        debug!("Analyzing sequence composition...");

        let num_buckets = ((seq.len() + window_length - 1) / window_length) as usize;

        // Count, GC characters, non-N chracters, and record "is N" flags.
        let mut gc_count = vec![0 as i32; num_buckets];
        let mut non_n_count = vec![0 as i32; num_buckets];
        let mut has_gap = vec![false; num_buckets];

        // Count GC chars and establish gap status.
        debug!("Counting GC and N characters...");
        for (i, c) in seq.iter().enumerate() {
            let bucket = i / window_length;

            match *c as char {
                'g' | 'c' | 'G' | 'C' => {
                    gc_count[bucket] += 1;
                    non_n_count[bucket] += 1;
                }
                'a' | 't' | 'A' | 'T' => {
                    non_n_count[bucket] += 1;
                }
                'n' | 'N' => {
                    has_gap[bucket] = true;
                }
                _ => (),
            }
        }

        // Compute GC content from counts.
        debug!("Computing GC content from GC counts...");
        let mut gc_ratio = vec![0 as f32; num_buckets];
        for (i, count) in gc_count.iter().enumerate() {
            gc_ratio[i] = if non_n_count[i] == 0 {
                std::f32::NAN
            } else {
                (*count as f32) / (non_n_count[i] as f32)
            };
        }

        Ok((gc_ratio, has_gap))
    }
}
