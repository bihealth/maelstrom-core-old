/// Code for aggregating BAM records in windows and target regions.
use std::cmp::min;

use lib_config::DepthOfCoverageConfig;

use rust_htslib::{bam, bam::Read};

use bio_types::genome::{AbstractInterval, Interval};
use lib_common::error::Error;
use lib_common::stats::Stats;

/// Struct for the result in one region (window or target).
#[derive(Debug)]
pub struct AggregationStats {
    /// Coverage in the window to use after the "coverage" step.
    pub cov: f32,
    /// Optional standard deviation of coverage in the window.
    pub cov_sd: Option<f32>,

    /// Start of the region on the contig.
    pub start: usize,
    /// End of the region on the contig.
    pub end: usize,

    /// Mean MAPQ in region.
    pub mean_mapq: f32,
}

/// Struct with common information for aggregator.
#[derive(Debug)]
pub struct BaseAggregator {
    /// Configuration.
    config: DepthOfCoverageConfig,
    /// Information about the contig to process.
    contig: Interval,

    /// Number of processed records.
    num_processed: u32,
    /// Number of skipped records.
    num_skipped: u32,
}

impl BaseAggregator {
    pub fn contig_length(&self) -> usize {
        (self.contig.range().end - self.contig.range().start) as usize
    }
}

/// Trait for alignment aggregation from BAM files.
pub trait BamRecordAggregator {
    /// Put all `fetch()`ed records from `reader` into the aggregator.
    fn put_fetched_records(
        &mut self,
        reader: &mut bam::IndexedReader,
        prog: &dyn std::ops::Fn(i64) -> (),
    ) -> Result<(), Error>;

    /// Return statistics for the given window.
    fn get_stats(&self, window_id: usize) -> AggregationStats;

    /// Return number of windows/targets.
    fn num_regions(&self) -> usize;

    /// Number of processed records.
    fn num_processed(&self) -> u32;

    /// Number of skipped records.
    fn num_skipped(&self) -> u32;
}

/// Struct for aggregating fragment counts in a genome-wide fashion.
#[derive(Debug)]
pub struct FragmentsAggregator {
    /// Common information from all aggregators.
    base: BaseAggregator,

    /// Per-window fragment counts (each bin is a `u32`).
    counters: Vec<u32>,
    /// Sum of MAPQ values.
    mapq_sums: Vec<u64>,
}

impl FragmentsAggregator {
    /// Construct new aggregator with the given BAM `header`.  This information is
    /// necessary to appropriately allocate buffers for all samples in the header.
    pub fn new(config: DepthOfCoverageConfig, contig: Interval) -> FragmentsAggregator {
        let base = BaseAggregator {
            config,
            contig,
            num_processed: 0,
            num_skipped: 0,
        };
        let contig_length = base.contig_length();
        let num_bins = (contig_length + base.config.window_length - 1) / base.config.window_length;

        FragmentsAggregator {
            base,
            counters: vec![0; num_bins],
            mapq_sums: vec![0; num_bins],
        }
    }

    fn put_bam_record(&mut self, record: &bam::Record) {
        if !self.skip_mapq(record)
            && !self.skip_flags(record)
            && !self.skip_discordant(record)
            && !self.skip_clipping(record)
            && !self.skip_paired_and_all_but_leftmost(record)
        {
            self.base.num_processed += 1;

            let fragment_center = if record.is_paired() {
                record.pos() + record.insert_size() / 2
            } else {
                record.pos() + record.cigar().end_pos()
            } as u32;

            let window_length = self.base.config.window_length;
            let bin = fragment_center as usize / window_length;
            self.counters[bin] += 1;
            self.mapq_sums[bin] += record.mapq() as u64;
        }
    }

    // Skip `record` based on `MAPQ`?
    fn skip_mapq(&self, record: &bam::Record) -> bool {
        if let Some(min_mapq) = self.base.config.min_mapq {
            record.mapq() < min_mapq
        } else {
            false
        }
    }

    // Skip `record` because of flags.
    fn skip_flags(&self, record: &bam::Record) -> bool {
        record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
            || record.is_quality_check_failed()
    }

    /// Skip `record` because of discordant alignment?
    fn skip_discordant(&self, record: &bam::Record) -> bool {
        if !record.is_paired() {
            false // unpaired cannot be discordant
        } else {
            !record.is_proper_pair()
        }
    }

    // Skip `record` because of too much clipping?
    fn skip_clipping(&self, record: &bam::Record) -> bool {
        if let Some(min_unclipped) = self.base.config.min_unclipped {
            let cigar = record.cigar();
            ((record.seq().len() as i64
                - cigar.leading_hardclips()
                - cigar.leading_softclips()
                - cigar.trailing_hardclips()
                - cigar.trailing_softclips()) as f32)
                / (record.seq().len() as f32)
                < min_unclipped
        } else {
            false
        }
    }

    // Skip `record` because it is paired and not the leftmost fragment (or discordant).
    fn skip_paired_and_all_but_leftmost(&self, record: &bam::Record) -> bool {
        record.is_paired() && (record.insert_size() <= 0)
    }
}

impl BamRecordAggregator for FragmentsAggregator {
    fn put_fetched_records(
        &mut self,
        reader: &mut bam::IndexedReader,
        prog: &dyn std::ops::Fn(i64) -> (),
    ) -> Result<(), Error> {
        let mut counter = 0;
        let mut record = bam::Record::new();
        loop {
            if !reader.read(&mut record)? {
                break;
            } else {
                self.put_bam_record(&record);
            }

            if counter % 10_000 == 0 {
                prog(record.pos());
            }
            counter += 1;
        }

        Ok(())
    }

    fn get_stats(&self, window_id: usize) -> AggregationStats {
        let window_length = self.base.config.window_length;

        AggregationStats {
            cov: self.counters[window_id] as f32,
            cov_sd: None,

            start: window_id * window_length,
            end: min(self.base.contig_length(), (window_id + 1) * window_length),
            mean_mapq: if self.counters[window_id] == 0 {
                0.0
            } else {
                (self.mapq_sums[window_id] as f64 / self.counters[window_id] as f64) as f32
            },
        }
    }

    fn num_regions(&self) -> usize {
        (self.base.contig_length() + self.base.config.window_length - 1)
            / self.base.config.window_length
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }
}

// Bin for coverage aggregation.
#[derive(Debug, Clone)]
pub struct CoverageBin {
    /// Mean coverage.
    pub cov_mean: f32,
    /// Coverage standard deviation.
    pub cov_stddev: f32,
    /// Mean MAPQ, weighted by aligned bases.
    pub mapq_mean: f32,
}

impl CoverageBin {
    fn new() -> Self {
        CoverageBin {
            cov_mean: 0_f32,
            cov_stddev: 0_f32,
            mapq_mean: 0_f32,
        }
    }
}

// Struct for aggregating as coverage.
#[derive(Debug)]
pub struct CoverageAggregator {
    /// Common information from all aggregators.
    pub base: BaseAggregator,
    /// Number of bases for each base in the current bin for the one sample in the BAM file.
    pub coverage: Vec<CoverageBin>,

    /// Base-wise coverage information for the current window.
    pub depths: Vec<usize>,
    /// Sum of MAPQ values.
    pub mapqs: Vec<usize>,
}

impl CoverageAggregator {
    pub fn new(config: DepthOfCoverageConfig, contig: Interval) -> CoverageAggregator {
        let base = BaseAggregator {
            config,
            contig,
            num_processed: 0,
            num_skipped: 0,
        };
        let contig_length = base.contig_length();
        let window_length = base.config.window_length;
        let num_bins = (contig_length + window_length - 1) / window_length;

        CoverageAggregator {
            base,
            coverage: vec![CoverageBin::new(); num_bins],
            depths: vec![0; window_length],
            mapqs: vec![0; window_length],
        }
    }

    fn push_window(&mut self, window_id: usize) {
        let depth_sum: usize = self.depths.iter().sum();
        let mapq_sum: usize = self.mapqs.iter().sum();
        let depths: Vec<f64> = self.depths.iter().map(|x| *x as f64).collect();
        self.coverage[window_id] = CoverageBin {
            cov_mean: (&depths).mean() as f32,
            cov_stddev: (&depths).std_dev() as f32,
            mapq_mean: if depth_sum == 0 {
                0.0
            } else {
                ((mapq_sum as f64) / (depth_sum as f64)) as f32
            },
        };
        let window_length = self.base.config.window_length as usize;
        self.depths = vec![0; window_length];
        self.mapqs = vec![0; window_length];
    }
}

impl BamRecordAggregator for CoverageAggregator {
    /// Put all `fetch()`ed records from `reader` into the aggregator.
    fn put_fetched_records(
        &mut self,
        reader: &mut bam::IndexedReader,
        prog: &dyn std::ops::Fn(i64) -> (),
    ) -> Result<(), Error> {
        let window_length = self.base.config.window_length as usize;

        // Iterate over all pileups
        let mut window_id = None;
        for (counter, pileup) in reader.pileup().enumerate() {
            let pileup = pileup.unwrap();
            let pos = pileup.pos() as usize;

            // On window change, push window to result.
            let next_window_id = pos / window_length;
            if let Some(window_id) = window_id {
                if window_id != next_window_id {
                    self.push_window(window_id);
                }
            }
            window_id = Some(next_window_id);

            // Compute depth of "valid" reads (note that only single/first-read coverage)
            // is computed.
            let mapqs = pileup
                .alignments()
                .filter(|alignment| {
                    let record = alignment.record();
                    !record.is_secondary()
                        && !record.is_duplicate()
                        && !record.is_supplementary()
                        && !record.is_duplicate()
                        && !record.is_quality_check_failed()
                        && self
                            .base
                            .config
                            .min_mapq
                            .map(|min_mapq| record.mapq() >= min_mapq)
                            .unwrap_or(true)
                        && (!record.is_paired() || record.is_proper_pair())
                })
                .map(|alignment| alignment.record().mapq())
                .collect::<Vec<u8>>();
            self.depths[pos % window_length] = mapqs.len();
            self.mapqs[pos % window_length] += mapqs.iter().map(|x| *x as usize).sum::<usize>();

            if counter % 10_000 == 0 {
                prog(pos as i64);
            }
        }

        if let Some(window_id) = window_id {
            self.push_window(window_id);
        }

        Ok(())
    }

    fn get_stats(&self, window_id: usize) -> AggregationStats {
        let window_length = self.base.config.window_length;
        AggregationStats {
            cov: self.coverage[window_id as usize].cov_mean,
            cov_sd: Some(self.coverage[window_id as usize].cov_stddev),
            mean_mapq: self.coverage[window_id as usize].mapq_mean,

            start: window_id * window_length,
            end: min(self.base.contig_length(), (window_id + 1) * window_length),
        }
    }

    fn num_regions(&self) -> usize {
        (self.base.contig_length() + self.base.config.window_length - 1)
            / self.base.config.window_length
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }
}
