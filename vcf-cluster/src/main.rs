use std::cmp::Ordering;
/// vcf-cluster -- Cluster VCF files with structural variants.
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::iter::FromIterator;
use std::path::Path;

use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Interval};
use clap::{App, Arg, ArgMatches};
use disjoint_sets::UnionFind;
use git_version::git_version;
// use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, error, info, LevelFilter};
use rust_htslib::{bcf, bcf::Read};

use lib_common::bcf::{build_vcf_header, collect_contigs, guess_bcf_format};
use lib_common::error::Error;
use lib_common::parse_region;
use lib_common::sv::*;
use lib_config::{ClusterSettings, Config};

/// Command line options
#[derive(Debug)]
struct Options {
    /// Verbosity level
    verbosity: u64,
    /// List of regions to call.
    regions: Option<Vec<Interval>>,
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
            regions: matches
                .value_of("regions")
                .map(|s| {
                    let x: Result<Vec<Interval>, Error> =
                        s.split(',').map(|t| parse_region(&t)).collect();
                    x
                })
                .transpose()?,
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
            setting: matches
                .value_of("setting")
                .unwrap_or("per_tool_pesr")
                .to_string(),
        })
    }
}

fn range_len(r: &std::ops::Range<i64>) -> i64 {
    r.end - r.start
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
fn load_std_records(path: &str, region: &Interval) -> Result<Vec<StandardizedRecord>, Error> {
    let mut result: Vec<StandardizedRecord> = Vec::new();
    let mut reader = bcf::IndexedReader::from_path(path)?;
    let mut record = reader.empty_record();

    let rid = reader.header().name2rid(region.contig().as_bytes())?;
    if reader
        .fetch(rid, region.range().start, region.range().end)
        .is_ok()
    {
        while reader.read(&mut record)? {
            result.push(StandardizedRecord::from_bcf_record(&mut record)?);
            debug!("=> {:?}", &result.last());
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
    if config.match_sv_type && lhs.sv_type != rhs.sv_type {
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

/// Cluster the SVs for the given SV type.
fn cluster_records_sv_type(
    records: &[StandardizedRecord],
    cluster_settings: &ClusterSettings,
    sv_type: &str,
    result: &mut Vec<StandardizedRecord>,
) -> Result<(), Error> {
    // Build interval tree for efficient overlap detection.
    let mut ids = Vec::new();
    let mut tree = IntervalTree::new();
    for (i, record) in records.iter().enumerate() {
        if !cluster_settings.match_sv_type || record.sv_type == sv_type {
            debug!("inserting...");
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
        record.alt = if let Some(sv_type) = &cluster_settings.sv_type_out {
            format!("<{}>", &sv_type)
        } else {
            records[record_ids[0]].alt.clone()
        };
        record.chrom2 = records[record_ids[0]].chrom2.clone();
        record.sv_type = if let Some(sv_type) = &cluster_settings.sv_type_out {
            sv_type.clone()
        } else {
            records[record_ids[0]].sv_type.clone()
        };
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
        if sv_type == "BND" {
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

    Ok(())
}

/// Cluster the SVs (all are assumed to be on the same contig).
fn cluster_records(
    options: &Options,
    config: &Config,
    records: &[StandardizedRecord],
) -> Result<Vec<StandardizedRecord>, Error> {
    let mut result = Vec::new();

    let cluster_settings = match options.setting.as_str() {
        "per_tool_pesr" => config.vcf_cluster_presets_per_tool_pesr.clone(),
        "per_tool_doc" => config.vcf_cluster_presets_per_tool_doc.clone(),
        _ => return Err(Error::UnknownClusterSettingName()),
    };

    if cluster_settings.match_sv_type {
        for sv_type in &["DEL", "DUP", "INV", "BND"] {
            cluster_records_sv_type(records, &cluster_settings, &sv_type, &mut result)?;
        }
    } else {
        cluster_records_sv_type(records, &cluster_settings, &"<IGNORED>", &mut result)?;
    }

    result.sort();

    Ok(result)
}

/// Write out clustered records.
fn write_records(
    records: &[StandardizedRecord],
    writer: &mut bcf::Writer,
    offset: usize,
) -> Result<usize, Error> {
    for (i, r) in records.iter().enumerate() {
        let mut record = writer.empty_record();
        r.update_bcf_record(&mut record)?;
        record.set_id(format!("SV{:08}", offset + i + 1).as_bytes())?;

        // Actually write out output record.
        writer.write(&record)?;
    }

    Ok(offset + records.len())
}

/// Main entry point after parsing command line and loading options.
fn perform_clustering(options: &Options, config: &Config) -> Result<(), Error> {
    info!("Starting to cluster BCF file...");

    // Check that the contigs are the same for all input files.
    check_contigs(&options.paths_input)?;

    // Open the output files (piggy-back reading of contig names from first input file).
    let (regions, mut writer) = {
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

        let regions = if let Some(regions) = &options.regions {
            regions.clone()
        } else {
            collect_contigs(&reader)?
                .iter()
                .map(|name| Interval::new(name.clone(), 0..10_000_000_000))
                .collect()
        };

        (
            regions,
            bcf::Writer::from_path(
                &options.path_output,
                &header,
                guessed.uncompressed,
                guessed.format,
            )?,
        )
    };

    let mut offset = 0;
    for region in &regions {
        // Extract StandardizedRecord data for all contigs from all files.
        let records = {
            let mut records = Vec::new();
            for path in &options.paths_input {
                let mut result = load_std_records(&path, region)?;
                records.append(&mut result);
            }
            records
        };

        // Cluster collected StandardizedRecords on contig.
        let records = cluster_records(&options, &config, &records)?;

        // Write out clusters for current contig.
        offset = write_records(&records, &mut writer, offset)?;
    }

    info!("Done clustering BCF file...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("maelstrom-vcf-cluster")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Extract and standardize records from tool VCF files")
        .args(&[
            Arg::from_usage("-v... 'Increase verbosity'"),
            Arg::from_usage("--overwrite 'Allow overwriting of output file'"),
            Arg::from_usage("-c, --config=[FILE] 'Sets a custom config file'"),
            Arg::from_usage("-r, --regions=[REGIONS] 'comma-separated list of regions'"),
            Arg::from_usage(
                "-s, --setting=[SETTING] 'Use cluster settings name, one of \
                {per_tool_pesr,per_tool_doc}, default per_tool_pesr'",
            ),
            Arg::from_usage("<input>... 'input file to read from'"),
            Arg::from_usage("<output> 'output file to write to'"),
        ])
        .get_matches();
    let options = Options::from_arg_matches(&matches)?;

    // Output file must not exist yet.
    if options.path_output != "-"
        && options.path_output != "/dev/stdout"
        && Path::new(&options.path_output).exists()
        && !options.overwrite
    {
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
        .level(match matches.occurrences_of("verbose").cmp(&1) {
            Ordering::Less => LevelFilter::Info,
            Ordering::Equal => LevelFilter::Debug,
            Ordering::Greater => LevelFilter::Trace,
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();
    info!("Starting maelstrom-vcf-cluster");
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
    use super::Interval;
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `perform_clustering()` and compares the result.
    fn _perform_clustering_and_test(
        tmp_dir: &TempDir,
        setting: &str,
        paths_input: &Vec<String>,
        path_expected: &str,
        config_str: &str,
        regions: &Option<Vec<Interval>>,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());
        let options = super::Options {
            setting: setting.to_string(),
            verbosity: 1, // disable progress bar
            regions: regions.clone(),
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
    fn test_vcf_cluster_delly2_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_clustering_and_test(
            &tmp_dir,
            "per_tool_pesr",
            &vec![
                "./src/tests/data/ex-delly-1.vcf.gz".to_string(),
                "./src/tests/data/ex-delly-2.vcf.gz".to_string(),
            ],
            "./src/tests/data/ex-delly.expected.vcf",
            "stdvcf_apply_filters = true",
            &None,
        )?;
        Ok(())
    }

    #[test]
    fn test_vcf_cluster_cnmops() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_clustering_and_test(
            &tmp_dir,
            "per_tool_doc",
            &vec![
                "./src/tests/data/ex-cnmops-1.vcf.gz".to_string(),
                "./src/tests/data/ex-cnmops-2.vcf.gz".to_string(),
            ],
            "./src/tests/data/ex-cnmops.expected.vcf",
            "",
            &None,
        )?;
        Ok(())
    }
}
