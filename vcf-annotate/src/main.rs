/// vcf-annotate -- Create annotations for VCF file with SVs.
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::annot::contig::SeqContigStranded;
use bio_types::annot::loc::Loc;
use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::NoStrand;
use bio_types::strand::ReqStrand;
use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::sorted;
use log::{debug, info, LevelFilter};
use rust_htslib::{bcf, bcf::Read};

use lib_common::bcf::collect_contigs;
use lib_common::bed_to_annot_map;
use lib_common::error::Error;
use lib_common::parse_region;
use lib_common::read_evidence;
use lib_common::read_evidence::Sides;
use lib_common::read_evidence::Strand;
use lib_common::stats::Stats;
use lib_common::sv;
use lib_config::Config;

/// Command line options
#[derive(Debug)]
struct Options {
    /// Verbosity level
    verbosity: u64,
    /// List of regions to call.
    regions: Option<Vec<Interval>>,
    /// Path to configuration file to use,
    path_config: Option<String>,
    /// Path to input PE/SR evidence file.
    path_pesr_evidence: Option<String>,
    /// Path to input DoC evidence file.
    path_doc_evidence: Option<String>,
    /// Path to input small variant VCF.
    path_snv_vcf: Option<String>,
    /// Name of sample to annotate for.
    sample: String,
    /// Path to input VCF file.
    path_input: String,
    /// Path to output file.
    prefix_output: String,
    /// Overwrite output file(s).
    overwrite: bool,
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
            path_pesr_evidence: matches
                .value_of("path-pesr-evidence")
                .map(|s| s.to_string()),
            path_doc_evidence: matches.value_of("path-doc-evidence").map(|s| s.to_string()),
            path_snv_vcf: matches.value_of("path-snv-vcf").map(|s| s.to_string()),
            sample: match matches.value_of("sample") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            path_input: match matches.value_of("input") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            prefix_output: match matches.value_of("output") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            overwrite: matches.occurrences_of("overwrite") > 0,
        })
    }
}

/// Code to implement matrix I/O.
mod mtx_io {
    use std::io::Write;

    pub struct Writer {
        handle: std::fs::File,
    }

    impl Writer {
        pub fn from(path: &str, sample: &str) -> Result<Self, std::io::Error> {
            let mut handle = std::fs::File::create(path)?;
            handle.write_all(format!("#sv_id\t{}\n", &sample).as_bytes())?;
            Ok(Self { handle })
        }

        pub fn write(&mut self, sv_id: &str, value: f64) -> Result<(), std::io::Error> {
            self.handle
                .write_all(format!("{}\t{}\n", &sv_id, value).as_bytes())?;
            Ok(())
        }
    }
}

/// Load all evidence records for the given region.
fn load_read_evidence(
    annot_map: &mut AnnotMap<String, read_evidence::Record>,
    interval: &Interval,
    path_pesr_evidence: &str,
    options: &Options,
    blocked: &Option<AnnotMap<String, ()>>,
) -> Result<usize, Error> {
    let mut skipped = 0;

    info!("loading evidence for {:?}...", interval);
    let mut reader = read_evidence::IndexedReader::from_path(&path_pesr_evidence)?;
    if !reader.fetch(
        &interval.contig(),
        interval.range().start,
        interval.range().end,
    )? {
        return Ok(skipped);
    }

    let spinner_style = ProgressStyle::default_spinner()
        .tick_chars(".oO@* ")
        .template("{prefix:.bold.dim} {spinner} {msg} {elapsed}");
    let spinner = if options.verbosity == 0 {
        let spinner = ProgressBar::new_spinner();
        spinner.set_style(spinner_style);
        spinner.enable_steady_tick(100);
        Some(spinner)
    } else {
        None
    };

    let mut counter: usize = 0;

    while let Some(record) = reader.read_record()? {
        debug!("record = {:?}", &record);
        let location = match &record {
            read_evidence::Record::PairedRead {
                contig1,
                start1,
                end1,
                ..
            } => Contig::new(
                contig1.clone(),
                *start1 as isize,
                (end1 - start1) as usize,
                NoStrand::Unknown,
            ),
            read_evidence::Record::SplitRead {
                contig, start, end, ..
            } => Contig::new(
                contig.clone(),
                *start as isize,
                (end - start) as usize,
                NoStrand::Unknown,
            ),
        };
        let skip = if let Some(blocked) = blocked {
            blocked.find(&location).next().is_some()
        } else {
            false
        };
        if !skip {
            annot_map.insert_at(record, &location);
        } else {
            skipped += 1;
            debug!("skipping!");
        }

        counter += 1;
        if counter % 1_000 == 0 {
            if let Some(spinner) = &spinner {
                spinner.set_message(&format!("currently at {}", &location.contig(),));
            }
        }
    }
    if let Some(spinner) = &spinner {
        spinner.finish_and_clear();
    }
    info!("... done loading evidence for {:?}", &interval);

    Ok(skipped)
}

#[derive(Debug, PartialEq, PartialOrd, Eq, Ord)]
struct ReadEvidenceCount {
    sv_id: String,
    pe_count: usize,
    sr_count: usize,
}

/// Return PR/SR read names.
fn fetch_read_evidence(
    query: &SeqContigStranded,
    read_evidence: &AnnotMap<String, read_evidence::Record>,
    blocked: &Option<AnnotMap<String, ()>>,
) -> (HashSet<String>, HashSet<(String, bool)>) {
    debug!("  fetch_read_evidence");
    let mut prs = HashSet::new();
    let mut srs = HashSet::new();

    debug!("    finding {:?}", &query);
    debug!("    #found = {}", read_evidence.find(query).count());

    for record in read_evidence.find(query) {
        match record.data() {
            read_evidence::Record::PairedRead {
                read_name,
                strand1,
                contig1,
                start1,
                contig2,
                start2,
                ..
            } => {
                debug!(
                    "    [[PR]] read_name = {:?}, strand1 = {:?}",
                    &read_name, strand1
                );
                match (query.strand(), strand1) {
                    (ReqStrand::Forward, Strand::Forward)
                    | (ReqStrand::Reverse, Strand::Reverse) => {
                        let skip = if let Some(blocked) = blocked {
                            let location1 = Contig::new(
                                contig1.clone(),
                                *start1 as isize,
                                1,
                                NoStrand::Unknown,
                            );
                            let skip1 = blocked.find(&location1).next().is_some();
                            let skip2 = if let Some(start2) = start2 {
                                let location2 = Contig::new(
                                    contig2.clone().unwrap(),
                                    *start2 as isize,
                                    1,
                                    NoStrand::Unknown,
                                );
                                blocked.find(&location2).next().is_some()
                            } else {
                                false
                            };
                            skip1 || skip2
                        } else {
                            false
                        };
                        if !skip {
                            prs.insert(read_name.clone());
                        }
                    }
                    _ => (), // ignored; no strand match
                }
            }
            read_evidence::Record::SplitRead {
                read_name,
                clipped_sides,
                is_first,
                contig,
                start,
                ..
            } => {
                debug!(
                    "    [[SR]] read_name = {:?}, clipped_sides = {:?}",
                    &read_name, clipped_sides
                );
                match (query.strand(), clipped_sides) {
                    (ReqStrand::Forward, Sides::Both)
                    | (ReqStrand::Reverse, Sides::Both)
                    | (ReqStrand::Forward, Sides::Right)
                    | (ReqStrand::Reverse, Sides::Left) => {
                        let skip = if let Some(blocked) = blocked {
                            let location =
                                Contig::new(contig.clone(), *start as isize, 1, NoStrand::Unknown);
                            blocked.find(&location).next().is_some()
                        } else {
                            false
                        };
                        if !skip {
                            srs.insert((read_name.clone(), *is_first));
                        }
                    }
                    _ => (), // ignored; no side match
                }
            }
        }
    }

    (prs, srs)
}

/// Count PE and SR evidence.
fn count_evidence(
    left: &SeqContigStranded,
    right: &SeqContigStranded,
    read_evidence: &AnnotMap<String, read_evidence::Record>,
    blocked: &Option<AnnotMap<String, ()>>,
) -> (usize, usize) {
    let (left_prs, left_srs) = fetch_read_evidence(left, read_evidence, blocked);
    let (right_prs, right_srs) = fetch_read_evidence(right, read_evidence, blocked);
    debug!("count_evidence");
    debug!("  left_prs = {:?}", &left_prs);
    debug!("  left_srs = {:?}", &left_srs);
    debug!("  right_prs = {:?}", &right_prs);
    debug!("  right_srs = {:?}", &right_srs);
    debug!("------");
    debug!(
        "  both prs = {:?}",
        &left_prs
            .intersection(&right_prs)
            .cloned()
            .collect::<String>()
    );
    debug!(
        "  both prs = {:?}",
        &left_prs.intersection(&right_prs).count()
    );

    (
        left_prs.intersection(&right_prs).count(),
        left_srs.intersection(&right_srs).count(),
    )
}

/// Perform PE/SR annotation of SV.
fn annotate_pesr(
    options: &Options,
    config: &Config,
    read_evidence: &AnnotMap<String, read_evidence::Record>,
    region: &Interval,
    blocked: &Option<AnnotMap<String, ()>>,
) -> Result<Vec<ReadEvidenceCount>, Error> {
    fn f(a: isize, b: isize) -> isize {
        if b > a {
            0
        } else {
            a - b
        }
    }

    let mut reader = bcf::IndexedReader::from_path(&options.path_input)?;
    let res = reader.fetch(
        reader.header().name2rid(region.contig().as_bytes())?,
        region.range().start,
        region.range().end,
    );
    match res {
        Err(rust_htslib::bcf::errors::Error::Seek { .. }) => return Ok(vec![]),
        Err(_) => res?,
        Ok(_) => (),
    }

    let mut record = reader.empty_record();
    let mut result = Vec::new();
    while reader.read(&mut record)? {
        let sv_id = std::str::from_utf8(&record.id())?.to_string();
        debug!(">>> sv_id = {}", &sv_id);
        let record = sv::StandardizedRecord::from_bcf_record(&mut record)?;
        // Define a couple of shortcuts to make the matching code below more compact.
        use ReqStrand::{Forward, Reverse};
        let delta = config.annotate_read_evidence_max_dist as usize;
        let deltai = delta as isize;
        let chrom = &record.chrom;
        let chrom2 = &record.chrom2;
        let pos = record.pos as isize;
        let end2 = record.end2 as isize;
        let mut pe_count = 0;
        let mut sr_count = 0;

        let search_wheres = match (&record.sv_type[..], &record.strands[..]) {
            ("DEL", _) => vec![(
                SeqContigStranded::new(chrom.clone(), f(pos, deltai), delta, Forward),
                SeqContigStranded::new(chrom2.clone(), end2, delta, Reverse),
            )],
            ("DUP", _) => vec![(
                SeqContigStranded::new(chrom.clone(), pos, delta, Reverse),
                SeqContigStranded::new(chrom2.clone(), end2 - deltai, delta, Forward),
            )],
            ("INV", _) => vec![
                (
                    SeqContigStranded::new(chrom.clone(), f(pos, deltai), delta, Forward),
                    SeqContigStranded::new(chrom2.clone(), f(end2, deltai), delta, Forward),
                ),
                (
                    SeqContigStranded::new(chrom.clone(), pos, delta, Reverse),
                    SeqContigStranded::new(chrom2.clone(), f(end2, deltai), delta, Reverse),
                ),
            ],
            ("BND", "--") | ("BND", "++") => vec![(
                SeqContigStranded::new(chrom.clone(), pos, delta, Reverse),
                SeqContigStranded::new(chrom2.clone(), end2, delta, Reverse),
            )],
            ("BND", "-+") => vec![(
                SeqContigStranded::new(chrom.clone(), pos, delta, Reverse),
                SeqContigStranded::new(chrom2.clone(), f(end2, deltai), delta, Reverse),
            )],
            ("BND", "+-") => vec![(
                SeqContigStranded::new(chrom.clone(), f(pos, deltai), delta, Forward),
                SeqContigStranded::new(chrom2.clone(), end2, delta, Reverse),
            )],
            _ => panic!(format!(
                "Unknown SV/strands combination: {}/{}",
                &record.sv_type, &record.strands
            )),
        };

        debug!("Search where: {:?}", &search_wheres);

        for (left, right) in &search_wheres {
            debug!(">>>>> searching: {}/{}", &left, &right);
            let (pe, sr) = count_evidence(left, right, read_evidence, blocked);
            debug!(">>>>> pe = {}, sr = {}", pe, sr);
            pe_count += pe;
            sr_count += sr;
        }

        result.push(ReadEvidenceCount {
            sv_id,
            pe_count,
            sr_count,
        });
    }

    Ok(result)
}

/// Coverage evidence information for a structural variant.
#[derive(Debug, PartialEq, PartialOrd)]
struct CoverageEvidence {
    /// ID of the called SV
    sv_id: String,
    /// Normalized coverage
    norm_cov: f64,
}

struct DocWindow {
    cov: f64,
    mapq: f64,
}

/// Perform DoC annotation of SV.
fn annotate_doc(
    options: &Options,
    config: &Config,
    region: &Interval,
    median_doc: f64,
) -> Result<Vec<CoverageEvidence>, Error> {
    let mut reader = bcf::IndexedReader::from_path(&options.path_input)?;
    let res = reader.fetch(
        reader.header().name2rid(region.contig().as_bytes())?,
        region.range().start,
        region.range().end,
    );
    match res {
        Err(rust_htslib::bcf::errors::Error::Seek { .. }) => return Ok(vec![]),
        Err(_) => res?,
        Ok(_) => (),
    }

    let mut doc_reader = bcf::IndexedReader::from_path(options.path_doc_evidence.clone().unwrap())?;
    let mut doc_record = doc_reader.empty_record();

    let mut record = reader.empty_record();
    let mut result = Vec::new();
    while reader.read(&mut record)? {
        let sv_id = std::str::from_utf8(&record.id())?.to_string();
        debug!(">>> sv_id = {}", &sv_id);
        let record = sv::StandardizedRecord::from_bcf_record(&mut record)?;

        let norm_cov =
            if record.sv_type == "DEL" || record.sv_type == "DUP" || record.sv_type == "CNV" {
                let rid: u32 = doc_reader.header().name2rid(record.chrom.as_bytes())?;
                if doc_reader
                    .fetch(rid, record.pos as u64, record.end2 as u64)
                    .is_ok()
                {
                    let mut doc_windows = Vec::new();
                    while doc_reader.read(&mut doc_record)? {
                        doc_windows.push(DocWindow {
                            cov: doc_record.format(b"RCV").float()?[0][0].into(),
                            mapq: doc_record.format(b"MQ").float()?[0][0].into(),
                        });
                    }

                    let q0_windows = doc_windows.iter().filter(|w| w.mapq >= 1.0).count();
                    let do_filter = q0_windows >= config.doc_annotation_min_bins;
                    let covs: Vec<f64> = doc_windows
                        .iter()
                        .filter(|w| !do_filter || w.mapq < 1.0)
                        .map(|w| w.cov)
                        .collect();

                    if covs.is_empty() {
                        0.0_f64
                    } else {
                        covs.median() / median_doc
                    }
                } else {
                    std::f64::NAN
                }
            } else {
                std::f64::NAN
            };

        result.push(CoverageEvidence { sv_id, norm_cov });
    }

    Ok(result)
}

/// Perform BAF annotation of SV.
fn annotate_baf(
    _options: &Options,
    _config: &Config,
    _region: &Interval,
) -> Result<Vec<(String, f64)>, Error> {
    panic!("BAF analysis not implemented yet!");
    // Ok(vec![])
}

/// Load DoC from file and compute median.
fn load_doc_median(path: &str) -> Result<f64, Error> {
    let mut reader = bcf::Reader::from_path(path)?;
    let mut covs: Vec<f64> = Vec::new();
    let mut record = reader.empty_record();
    while reader.read(&mut record)? {
        covs.push(record.format(b"RCV").float()?[0][0].into());
    }
    Ok((&covs).median())
}

/// Main entry point after parsing command line and loading options.
fn perform_annotation(options: &Options, config: &Config) -> Result<(), Error> {
    info!("Starting to annotate variants for sample...");

    let mut pe_writer = options
        .path_pesr_evidence
        .clone()
        .map(|_p| {
            mtx_io::Writer::from(
                &format!("{}.pe.tsv", &options.prefix_output),
                &options.sample,
            )
        })
        .transpose()?;
    let mut sr_writer = options
        .path_pesr_evidence
        .clone()
        .map(|_p| {
            mtx_io::Writer::from(
                &format!("{}.sr.tsv", &options.prefix_output),
                &options.sample,
            )
        })
        .transpose()?;
    let median_doc = if let Some(path_doc_evidence) = &options.path_doc_evidence {
        info!("Computing median depth of coverage (DoC)...");
        let median_doc = load_doc_median(path_doc_evidence)?;
        info!("... median DoC is {}", median_doc);
        Some(median_doc)
    } else {
        None
    };
    let mut doc_writer = options
        .path_doc_evidence
        .clone()
        .map(|_p| {
            mtx_io::Writer::from(
                &format!("{}.doc.tsv", &options.prefix_output),
                &options.sample,
            )
        })
        .transpose()?;
    let mut baf_writer = options
        .path_snv_vcf
        .clone()
        .map(|_p| {
            mtx_io::Writer::from(
                &format!("{}.baf.tsv", &options.prefix_output),
                &options.sample,
            )
        })
        .transpose()?;

    let reader = bcf::IndexedReader::from_path(&options.path_input)?;

    let contigs = collect_contigs(&reader)?;

    let blocked = if let Some(blocked_regions_bed) = &config.blocked_regions_bed {
        let contigs: HashMap<String, String> =
            contigs.iter().map(|s| (s.clone(), s.clone())).collect();
        Some(bed_to_annot_map(blocked_regions_bed, &contigs)?)
    } else {
        None
    };

    let regions = if let Some(regions) = &options.regions {
        regions.clone()
    } else {
        contigs
            .iter()
            .map(|name| Interval::new(name.clone(), 0..10_000_000_000))
            .collect()
    };

    let mut skipped = 0;
    let read_evidence = if let Some(path_pesr_evidence) = &options.path_pesr_evidence {
        info!("Loading read-based evidence...");
        let mut read_evidence: AnnotMap<String, read_evidence::Record> = AnnotMap::new();
        for region in &regions {
            debug!("region = {:?}", &region);
            skipped += load_read_evidence(
                &mut read_evidence,
                region,
                path_pesr_evidence,
                &options,
                &blocked,
            )?;
        }
        debug!("evidence: {:?}", &read_evidence);
        Some(read_evidence)
    } else {
        None
    };
    info!("Skipped {} blocked evidence records", skipped);

    info!("Processing regions/contigs...");
    for region in &regions {
        info!("Processing contig {:?}", region);

        if let (Some(pe_writer), Some(sr_writer)) = (&mut pe_writer, &mut sr_writer) {
            for ReadEvidenceCount {
                sv_id,
                pe_count,
                sr_count,
            } in sorted(
                annotate_pesr(
                    &options,
                    &config,
                    &read_evidence.as_ref().unwrap(),
                    &region,
                    &blocked,
                )?
                .iter(),
            ) {
                pe_writer.write(&sv_id, *pe_count as f64)?;
                sr_writer.write(&sv_id, *sr_count as f64)?;
            }
        }
        if let Some(doc_writer) = &mut doc_writer {
            for CoverageEvidence { sv_id, norm_cov } in
                annotate_doc(&options, &config, &region, median_doc.unwrap())?.iter()
            {
                doc_writer.write(&sv_id, *norm_cov)?;
            }
        }
        if let Some(baf_writer) = &mut baf_writer {
            for (sv_id, score) in annotate_baf(&options, &config, &region)?.iter() {
                baf_writer.write(&sv_id, *score)?;
            }
        }
    }

    info!("Done annotating variants...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("maelstrom-vcf-annotate")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Create annotations for VCF file with SVs.")
        .args(&[
            Arg::from_usage("-v... 'Increase verbosity'"),
            Arg::from_usage("--overwrite 'Allow overwriting of output file'"),
            Arg::from_usage("-c, --config=[FILE] 'Sets a custom config file'"),
            Arg::from_usage("-r, --regions=[REGIONS] 'comma-separated list of regions'"),
            Arg::from_usage("--path-pesr-evidence=[FILE] 'Path to PE/SR evidence file'"),
            Arg::from_usage("--path-doc-evidence=[FILE] 'Path to DoC evidence file'"),
            Arg::from_usage("--path-snv-vcf=[FILE] 'Path to BAF evidence file'"),
            Arg::from_usage("-s, --sample=<SAMPLE> 'Set sample to analyze'"),
            Arg::from_usage("<input> 'input VCF file to read from'"),
            Arg::from_usage("<output> 'output prefix; suffixes: .{pe,sr,baf,doc}.tsv"),
        ])
        .get_matches();
    let options = Options::from_arg_matches(&matches)?;

    // Output files must not exist yet.
    if options.path_pesr_evidence.is_some()
        && Path::new(&format!("{}.pesr.tsv", options.prefix_output)).exists()
        && !options.overwrite
    {
        return Err(Error::OutputFileExists());
    }
    if options.path_doc_evidence.is_some()
        && Path::new(&format!("{}.doc.tsv", options.prefix_output)).exists()
        && !options.overwrite
    {
        return Err(Error::OutputFileExists());
    }
    if options.path_snv_vcf.is_some()
        && Path::new(&format!("{}.baf.tsv", options.prefix_output)).exists()
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
        .level(if matches.is_present("v") {
            LevelFilter::Debug
        } else {
            LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();
    info!("Starting maelstrom-vcf-annotate");
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

    // Perform the record annotation.
    perform_annotation(&options, &config)?;

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
    fn _perform_annotation_and_test(
        tmp_dir: &TempDir,
        sample: &str,
        path_pesr_evidence: Option<String>,
        path_doc_evidence: Option<String>,
        path_snv_vcf: Option<String>,
        path_input: &str,
        path_expected_prefix: &str,
        regions: &Option<Vec<Interval>>,
        config_text: &str,
    ) -> Result<(), super::Error> {
        let prefix_output = String::from(tmp_dir.path().join("out").to_str().unwrap());
        let options = super::Options {
            verbosity: 1, // disable progress bar
            regions: regions.clone(),
            path_config: None,
            path_pesr_evidence: path_pesr_evidence,
            path_doc_evidence: path_doc_evidence,
            path_snv_vcf: path_snv_vcf,
            sample: sample.to_string(),
            path_input: path_input.to_string(),
            prefix_output: prefix_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str(config_text).unwrap();

        super::perform_annotation(&options, &config)?;

        if options.path_pesr_evidence.is_some() {
            assert_eq!(
                fs::read_to_string(format!("{}.pe.tsv", path_expected_prefix)).unwrap(),
                fs::read_to_string(format!("{}.pe.tsv", prefix_output)).unwrap(),
                "paired-end",
            );
            assert_eq!(
                fs::read_to_string(format!("{}.sr.tsv", path_expected_prefix)).unwrap(),
                fs::read_to_string(format!("{}.sr.tsv", prefix_output)).unwrap(),
                "split-read",
            );
        }
        if options.path_doc_evidence.is_some() {
            assert_eq!(
                fs::read_to_string(format!("{}.doc.tsv", path_expected_prefix)).unwrap(),
                fs::read_to_string(format!("{}.doc.tsv", prefix_output)).unwrap(),
                "depth-of-coverage",
            );
        }
        if options.path_snv_vcf.is_some() {
            assert_eq!(
                fs::read_to_string(format!("{}.baf.tsv", path_expected_prefix)).unwrap(),
                fs::read_to_string(format!("{}.baf.tsv", prefix_output)).unwrap(),
                "b-allele-fraction",
            );
        }

        Ok(())
    }

    #[test]
    fn test_vcf_cluster_delly2_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_annotation_and_test(
            &tmp_dir,
            "sample-1",
            Some(String::from("./src/tests/data/ex-delly-pesr.tsv.gz")),
            Some(String::from("./src/tests/data/ex-delly-doc.vcf.gz")),
            None, // Some(String::from("./src/tests/data/ex-delly-snv.vcf.gz")),
            "./src/tests/data/ex-delly-svs.vcf.gz",
            "./src/tests/data/ex-delly.expected",
            &None,
            "",
        )?;
        Ok(())
    }

    #[test]
    fn test_vcf_cluster_delly2_filter_with_blocked() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_annotation_and_test(
            &tmp_dir,
            "sample-1",
            Some(String::from("./src/tests/data/ex-delly-pesr.tsv.gz")),
            Some(String::from("./src/tests/data/ex-delly-doc.vcf.gz")),
            None, // Some(String::from("./src/tests/data/ex-delly-snv.vcf.gz")),
            "./src/tests/data/ex-delly-svs.vcf.gz",
            "./src/tests/data/ex-delly.expected-blocked",
            &None,
            "blocked_regions_bed = \"./src/tests/data/ex-delly-blocked.bed\"",
        )?;
        Ok(())
    }
}
