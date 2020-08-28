/// vcf-cluster -- Cluster VCF files with structural variants.
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::iter::FromIterator;
use std::path::Path;

use bio::data_structures::interval_tree::IntervalTree;
use clap::{App, Arg, ArgMatches};
use disjoint_sets::UnionFind;
use git_version::git_version;
// use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, error, info, LevelFilter};
use rust_htslib::{bcf, bcf::Read};

use lib_common::{build_vcf_header, guess_bcf_format, Error};
use lib_config::{ClusterSettings, Config};

/// Command line options
#[derive(Debug)]
struct Options {
    /// Verbosity level
    verbosity: u64,
    /// Path to configuration file to use,
    path_config: Option<String>,
    /// Path to input files.
    paths_input: Vec<String>,
    /// Path to output file.
    path_output: String,
    /// Overwrite output file.
    overwrite: bool,
    /// Cluster setting name.
    setting: String,
}

impl Options {
    pub fn from_arg_matches<'a>(matches: &ArgMatches<'a>) -> Result<Self, Error> {
        Ok(Options {
            verbosity: matches.occurrences_of("v"),
            path_config: matches.value_of("config").map(|s| s.to_string()),
            paths_input: match matches.values_of("input") {
                Some(xs) => xs.map(String::from).collect(),
                None => return Err(Error::OptionMissing()),
            },
            path_output: match matches.value_of("output") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            overwrite: matches.occurrences_of("overwrite") > 0,
            setting: match matches.value_of("setting") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
        })
    }
}

/// Representation of a structural variant record suitable for clustering.
#[derive(Debug, Clone)]
struct StandardizedRecord {
    chrom: String,
    pos: i64,
    reference: String,
    alt: String,
    filters: Vec<String>,
    chrom2: String,
    end2: i64,
    sv_type: String,
    strands: String,
    sv_len: i64,
    algorithms: Vec<String>,
    samples: Vec<String>,
    gts: Vec<String>,
    called_by: Vec<Vec<String>>,
}

impl PartialEq for StandardizedRecord {
    fn eq(&self, other: &Self) -> bool {
        self.chrom == other.chrom
            && self.pos == other.pos
            && self.chrom2 == other.chrom2
            && self.end2 == other.end2
            && self.reference == other.reference
            && self.alt == other.alt
    }
}

impl Eq for StandardizedRecord {}

impl Ord for StandardizedRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.pos, self.end2, &self.sv_type, &other.strands).cmp(&(
            other.pos,
            other.end2,
            &other.sv_type,
            &other.strands,
        ))
    }
}

impl PartialOrd for StandardizedRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl StandardizedRecord {
    fn new() -> Self {
        Self {
            chrom: String::new(),
            pos: 0,
            reference: String::new(),
            alt: String::new(),
            filters: Vec::new(), // TODO: remove?
            chrom2: String::new(),
            end2: 0,
            sv_type: String::new(),
            strands: String::new(),
            sv_len: 0,
            algorithms: Vec::new(),
            samples: Vec::new(),
            gts: Vec::new(),
            called_by: Vec::new(),
        }
    }

    fn from_bcf_record(record: &mut bcf::Record) -> Result<Self, Error> {
        let sample_count = record.header().sample_count() as usize;

        let mut called_by: Vec<Vec<String>> = vec![vec![]; sample_count];
        for algorithm in vec!["delly".to_string()] {
            if let Ok(arr) = record.format(b"delly").integer() {
                for i in 0..sample_count {
                    if arr[i][0] != 0 {
                        called_by[i].push(algorithm.clone());
                    }
                }
            }
        }

        Ok(Self {
            chrom: String::from_utf8(record.header().rid2name(record.rid().unwrap())?.to_vec())?,
            pos: record.pos(),
            reference: String::from_utf8(record.alleles()[0].to_vec())?,
            alt: String::from_utf8(record.alleles()[1].to_vec())?,
            filters: record
                .filters()
                .map(|id| String::from_utf8(record.header().id_to_name(id)).unwrap())
                .collect(),
            chrom2: String::from_utf8(record.info(b"CHR2").string()?.unwrap()[0].to_vec())?,
            end2: record.info(b"END2").integer()?.unwrap()[0] as i64,
            sv_type: String::from_utf8(record.info(b"SVTYPE").string()?.unwrap()[0].to_vec())?,
            strands: String::from_utf8(record.info(b"STRANDS").string()?.unwrap()[0].to_vec())?,
            sv_len: record.info(b"SVLEN").integer()?.unwrap()[0] as i64,
            algorithms: record
                .info(b"ALGORITHMS")
                .string()?
                .unwrap()
                .iter()
                .map(|s| String::from_utf8(s.to_vec()).unwrap())
                .collect(),
            samples: record
                .header()
                .samples()
                .iter()
                .map(|s| String::from_utf8(s.to_vec()).unwrap())
                .collect(),
            gts: (0..sample_count)
                .map(|i| {
                    record
                        .genotypes()
                        .unwrap()
                        .get(i as usize)
                        .iter()
                        .map(|gt_allele| match gt_allele {
                            bcf::record::GenotypeAllele::Unphased(i)
                            | bcf::record::GenotypeAllele::Phased(i) => format!("{}", i),
                            _ => ".".to_string(),
                        })
                        .collect::<Vec<_>>()
                        .join("-")
                })
                .collect(),
            called_by,
        })
    }

    fn interval(&self) -> std::ops::Range<i64> {
        self.extended_interval(0)
    }

    fn extended_interval(&self, delta: i64) -> std::ops::Range<i64> {
        let mut tmp = if self.sv_type == "BND" {
            self.pos..(self.pos + 1)
        } else {
            self.pos..self.end2
        };
        if tmp.start > tmp.end {
            std::mem::swap(&mut tmp.start, &mut tmp.end);
        }
        tmp.start -= delta;
        tmp.end += delta;
        tmp
    }
}

fn overlap(lhs: &std::ops::Range<i64>, rhs: &std::ops::Range<i64>) -> std::ops::Range<i64> {
    if lhs.end <= rhs.start || rhs.end <= lhs.start {
        #[allow(clippy::reversed_empty_ranges)]
        std::ops::Range { start: 0, end: 0 }
    } else {
        std::ops::Range {
            start: std::cmp::max(lhs.start, rhs.start),
            end: std::cmp::min(lhs.end, rhs.end),
        }
    }
}

fn range_len(r: &std::ops::Range<i64>) -> i64 {
    r.end - r.start
}

fn collect_contigs(reader: &bcf::IndexedReader) -> Result<Vec<String>, Error> {
    let mut result: Vec<String> = Vec::new();
    let header = reader.header();
    for rid in 0..header.contig_count() {
        result.push(std::str::from_utf8(header.rid2name(rid)?)?.to_string());
    }
    Ok(result)
}

/// Check contigs between all input files.
fn check_contigs(paths_input: &[String]) -> Result<(), Error> {
    info!("Checking input files for having contigs...");
    if paths_input.is_empty() {
        error!("No input files given (you should never reach here)!");
        return Err(Error::InconsistentInput());
    }

    let reader = bcf::IndexedReader::from_path(&paths_input[0])?;
    let first_contigs = collect_contigs(&reader)?;
    let from_first: HashSet<String> = HashSet::from_iter(first_contigs.iter().cloned());

    for path in paths_input {
        let reader = bcf::IndexedReader::from_path(&path)?;
        let from_current = HashSet::from_iter(collect_contigs(&reader)?.iter().cloned());
        if from_first != from_current {
            return Err(Error::InconsistentInput());
        }
    }

    info!("All input files have consistent contigs.");
    Ok(())
}

/// Extract StandardizedRecords from VCF file at the given path.
fn load_std_records(path: &str, contig_name: &str) -> Result<Vec<StandardizedRecord>, Error> {
    let mut result: Vec<StandardizedRecord> = Vec::new();
    let mut reader = bcf::IndexedReader::from_path(path)?;
    let mut record = reader.empty_record();

    let rid = reader.header().name2rid(contig_name.as_bytes())?;
    if reader.fetch(rid, 0, 1_000_000_000).is_ok() {
        loop {
            if !reader.read(&mut record)? {
                break; // done
            }
            result.push(StandardizedRecord::from_bcf_record(&mut record)?);
        }
    };

    Ok(result)
}

/// Check whether the two records are to be merged.
fn is_overlap_okay(
    lhs: &StandardizedRecord,
    rhs: &StandardizedRecord,
    config: &ClusterSettings,
) -> bool {
    debug!("Comparing \n  {:?}\nand\n  {:?}\n--", lhs, rhs);
    if lhs.sv_type != rhs.sv_type {
        panic!("Comparing SVs of different type!");
    }
    let ovl = overlap(&lhs.interval(), &rhs.interval());
    let ovl_lhs = (range_len(&ovl) as f32) / (range_len(&lhs.interval()) as f32);
    let ovl_rhs = (range_len(&ovl) as f32) / (range_len(&rhs.interval()) as f32);
    let rec_ovl = ovl_lhs.min(ovl_rhs);
    debug!(
        "len ovl: {} || {} || {} || {}",
        ovl_lhs, ovl_rhs, rec_ovl, config.reciprocal_overlap,
    );
    if rec_ovl < config.reciprocal_overlap {
        return false;
    }
    if let Some(max_bp_distance) = config.max_bp_distance {
        debug!(
            "bp distance: {} || {} || {}",
            (lhs.interval().start - rhs.interval().start).abs(),
            (lhs.interval().end - rhs.interval().end).abs(),
            max_bp_distance
        );
        if (lhs.interval().start - rhs.interval().start).abs() > max_bp_distance
            || (lhs.interval().end - rhs.interval().end).abs() > max_bp_distance
        {
            return false;
        }
    }
    debug!("strands: {} || {}", lhs.strands, rhs.strands);
    if config.match_strands && lhs.strands != rhs.strands {
        return false;
    }
    if let Some(sample_overlap) = config.sample_overlap {
        let lhs_samples: HashSet<String> = HashSet::from_iter(lhs.samples.iter().cloned());
        let rhs_samples: HashSet<String> = HashSet::from_iter(rhs.samples.iter().cloned());
        let int_samples: HashSet<String> = lhs_samples
            .intersection(&rhs_samples)
            .to_owned()
            .cloned()
            .collect();
        let ovl_lhs_samples = (int_samples.len() as f32) / (lhs_samples.len() as f32);
        let ovl_rhs_samples = (int_samples.len() as f32) / (rhs_samples.len() as f32);
        debug!(
            "sample ovl: {} || {} || {}",
            ovl_lhs_samples, ovl_rhs_samples, sample_overlap
        );
        if ovl_lhs_samples.min(ovl_rhs_samples) < sample_overlap {
            return false;
        }
    }
    debug!("Overlap!");

    true
}

/// Cluster the SVs (all are assumed to be on the same contig).
fn cluster_records(
    options: &Options,
    config: &Config,
    records: &[StandardizedRecord],
) -> Result<Vec<StandardizedRecord>, Error> {
    let mut result = Vec::new();

    let cluster_settings = match options.setting.as_str() {
        "per_tool_pesr" => config.clusvcf_presets_per_tool_pesr.clone(),
        "per_tool_doc" => config.clusvcf_presets_per_tool_doc.clone(),
        _ => return Err(Error::UnknownClusterSettingName()),
    };

    for sv_type in &["DEL", "DUP", "INV", "BND"] {
        // Build interval tree for efficient overlap detection.
        let mut ids = Vec::new();
        let mut tree = IntervalTree::new();
        for (i, record) in records.iter().enumerate() {
            if record.sv_type == *sv_type {
                tree.insert(
                    record.extended_interval(cluster_settings.max_bp_distance.unwrap_or(0)),
                    ids.len(),
                );
                ids.push(i);
            }
        }

        // Use union-find data structure for joining overlapping records into one cluster.
        let mut uf = UnionFind::new(ids.len());
        for record_idx in 0..ids.len() {
            let record = &records[ids[record_idx]];
            debug!("find {:?}", record.interval());
            for ovl_it in tree.find(record.interval()) {
                let other_idx = *ovl_it.data();
                let other = &records[ids[other_idx]];

                if is_overlap_okay(&record, &other, &cluster_settings) {
                    debug!("union: {} {}", record_idx, other_idx);
                    uf.union(record_idx, other_idx);
                }
            }
        }

        // Gather clusters of record ids.
        let mut cluster: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (record_idx, record_id) in ids.iter().enumerate() {
            let set_id = uf.find(record_idx);
            cluster
                .entry(set_id)
                .or_insert_with(Vec::new)
                .push(*record_id);
        }

        // Collapse each cluster.
        for record_ids in cluster.values() {
            debug!("collapsing: {:?}", &record_ids);
            let mut record = StandardizedRecord::new();
            record.chrom = records[record_ids[0]].chrom.clone();
            record.reference = String::from("N");
            record.alt = records[record_ids[0]].alt.clone();
            record.chrom2 = records[record_ids[0]].chrom2.clone();
            record.sv_type = records[record_ids[0]].sv_type.clone();
            record.strands = records[record_ids[0]].strands.clone();

            for record_id in record_ids {
                let mut algorithms = records[*record_id].algorithms.clone();
                record.algorithms.append(&mut algorithms);
            }
            record.algorithms.sort();
            record.algorithms.dedup();

            let mut sample_idx: HashMap<String, usize> = HashMap::new();
            for record_id in record_ids {
                let l_record = &records[*record_id];
                for (l_idx, sample) in l_record.samples.iter().enumerate() {
                    if let Some(&idx) = sample_idx.get(sample) {
                        record.gts[idx] = match std::cmp::max(
                            record.gts[idx].matches('1').count(),
                            l_record.gts[l_idx].matches('1').count(),
                        ) {
                            0 => "0/0".to_string(),
                            1 => "0/1".to_string(),
                            _ => "1/1".to_string(),
                        };
                        let mut called_by = l_record.called_by[l_idx].clone();
                        record.called_by[idx].append(&mut called_by);
                        record.called_by[idx].sort();
                        record.called_by[idx].dedup();
                    } else {
                        sample_idx.insert(sample.clone(), record.samples.len());
                        record.samples.push(sample.clone());
                        record.gts.push(l_record.gts[l_idx].clone());
                        record.called_by.push(l_record.called_by[l_idx].clone());
                    }
                }
            }

            let mut pos: Vec<_> = record_ids.iter().map(|id| records[*id].pos).collect();
            pos.sort();
            debug!("pos = {:?}", &pos);
            record.pos = if pos.len() % 2 == 1 {
                pos[pos.len() / 2]
            } else {
                (pos[pos.len() / 2] + pos[pos.len() / 2 - 1]) / 2
            };
            if *sv_type == "BND" {
                record.end2 = record.pos + 1;
                record.sv_len = -1;
            } else {
                let mut end2: Vec<_> = record_ids.iter().map(|id| records[*id].end2).collect();
                end2.sort();
                debug!("end2 = {:?}", &end2);
                record.end2 = if end2.len() % 2 == 1 {
                    end2[end2.len() / 2]
                } else {
                    (end2[end2.len() / 2] + end2[end2.len() / 2 - 1]) / 2
                };
                record.sv_len = record.end2 - record.pos;
            }

            result.push(record);
        }
    }

    result.sort();

    Ok(result)
}

/// Write out clustered records.
fn write_records(records: &[StandardizedRecord], writer: &mut bcf::Writer) -> Result<(), Error> {
    let sample_count = writer.header().sample_count();

    for r in records {
        let mut record = writer.empty_record();
        record.set_rid(Some(writer.header().name2rid(r.chrom.as_bytes())?));
        record.set_pos(r.pos);
        let alleles: Vec<&[u8]> = vec![&r.reference.as_bytes(), &r.alt.as_bytes()];
        record.set_alleles(&alleles)?;
        for filter in &r.filters {
            record.push_filter(writer.header().name_to_id(&filter.as_bytes())?);
        }
        record.push_info_integer(b"END2", vec![r.end2 as i32].as_slice())?;
        record.push_info_string(b"CHR2", vec![r.chrom2.as_bytes()].as_slice())?;
        record.push_info_string(b"SVTYPE", vec![r.sv_type.as_bytes()].as_slice())?;
        record.push_info_string(b"STRANDS", vec![r.strands.as_bytes()].as_slice())?;
        record.push_info_integer(b"SVLEN", vec![r.sv_len as i32].as_slice())?;
        let algorithms: Vec<&[u8]> = r.algorithms.iter().map(|a| a.as_bytes()).collect();
        record.push_info_string(b"ALGORITHMS", algorithms.as_slice())?;

        // Prepare buffers for FORMAT fields.
        let mut called_by: HashMap<String, Vec<i32>> = HashMap::new();
        for algorithm in &r.algorithms {
            called_by.insert(algorithm.clone(), vec![0; sample_count as usize]);
        }

        // Write information from "sparse" StandardizedRecord into BCF record.
        let mut genotypes = vec![0; 2 * sample_count as usize];
        for (i, sample) in r.samples.iter().enumerate() {
            let sample_idx = *writer.header().sample_to_id(sample.as_bytes())? as usize;

            genotypes[2 * sample_idx] = (r.gts[i].starts_with('1') as i32 + 1) << 1;
            genotypes[2 * sample_idx + 1] = (r.gts[i].ends_with('1') as i32 + 1) << 1;
            for algorithm in &r.called_by[i] {
                if let Some(arr) = called_by.get_mut(algorithm) {
                    arr[sample_idx] = 1;
                }
            }
        }

        // Write FORMAT fields into output record.
        record
            .push_format_integer(b"GT", &genotypes)
            .expect("could not write genotype");
        for (algorithm, is_called) in &called_by {
            record.push_format_integer(algorithm.as_bytes(), &is_called)?;
        }

        // Actually write out output record.
        writer.write(&record)?;
    }

    Ok(())
}

/// Main entry point after parsing command line and loading options.
fn perform_clustering(options: &Options, config: &Config) -> Result<(), Error> {
    info!("Starting to cluster BCF file...");

    // Check that the contigs are the same for all input files.
    check_contigs(&options.paths_input)?;

    // Open the output files (piggy-back reading of contig names from first input file).
    let (contig_names, mut writer) = {
        let reader = bcf::IndexedReader::from_path(&options.paths_input[0])?;
        let header = {
            let mut header = build_vcf_header(reader.header())?;
            let mut seen: HashSet<String> = HashSet::from_iter(
                reader
                    .header()
                    .samples()
                    .iter()
                    .map(|&x| String::from_utf8(x.to_vec()).unwrap()),
            );
            for path in &options.paths_input[1..] {
                let reader = bcf::IndexedReader::from_path(path)?;
                let unseen: Vec<String> = reader
                    .header()
                    .samples()
                    .iter()
                    .map(|&x| String::from_utf8(x.to_vec()).unwrap())
                    .filter(|s| !seen.contains(s))
                    .collect();
                for name in &unseen {
                    seen.insert(name.clone());
                    header.push_sample(name.as_bytes());
                }
            }

            header
        };
        let guessed = guess_bcf_format(&options.path_output);

        (
            collect_contigs(&reader)?,
            bcf::Writer::from_path(
                &options.path_output,
                &header,
                guessed.uncompressed,
                guessed.format,
            )?,
        )
    };

    for contig_name in contig_names {
        // Extract StandardizedRecord data for all contigs from all files.
        let records = {
            let mut records = Vec::new();
            for path in &options.paths_input {
                let mut result = load_std_records(&path, &contig_name)?;
                records.append(&mut result);
            }
            records
        };

        // Cluster collected StandardizedRecords on contig.
        let records = cluster_records(&options, &config, &records)?;

        // Write out clusters for current contig.
        write_records(&records, &mut writer)?;
    }

    info!("Done clustering BCF file...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("snappysv-stdvcf")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Extract and standardize records from tool VCF files")
        .args(&[
            Arg::from_usage("-v... 'Increase verbosity'"),
            Arg::from_usage("--overwrite 'Allow overwriting of output file'"),
            Arg::from_usage("-c, --config=[FILE] 'Sets a custom config file'"),
            Arg::from_usage(
                "-s, --setting=<SETTING> 'Use cluster settings name, one of \
                {per_tool_pesr,per_tool_doc}, default per_tool_pesr'",
            ),
            Arg::from_usage("<input>... 'input file to read from'"),
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
    info!("Starting snappysv-clusvcf");
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

    // Perform the record extraction.
    perform_clustering(&options, &config)?;

    info!("All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `perform_clustering()` and compares the result.
    fn _perform_clustering_and_test(
        tmp_dir: &TempDir,
        paths_input: &Vec<String>,
        path_expected: &str,
        config_str: &str,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());
        let options = super::Options {
            setting: "per_tool_pesr".to_string(),
            verbosity: 1, // disable progress bar
            path_config: None,
            paths_input: paths_input.clone(),
            path_output: path_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str(config_str).unwrap();

        super::perform_clustering(&options, &config)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }

    #[test]
    fn test_clusvcf_delly2_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_clustering_and_test(
            &tmp_dir,
            &vec![
                "./src/tests/data/ex-delly-1.vcf.gz".to_string(),
                "./src/tests/data/ex-delly-2.vcf.gz".to_string(),
            ],
            "./src/tests/data/ex-delly.expected.vcf",
            "stdvcf_apply_filters = true",
        )?;
        Ok(())
    }
}
