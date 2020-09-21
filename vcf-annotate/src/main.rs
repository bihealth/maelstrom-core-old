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
use log::{debug, info, LevelFilter};
use rust_htslib::{bcf, bcf::Read};

use lib_common::bcf::{build_vcf_header, collect_contigs, guess_bcf_format};
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
    path_output: String,
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
            path_output: match matches.value_of("output") {
                Some(x) => String::from(x),
                None => return Err(Error::OptionMissing()),
            },
            overwrite: matches.occurrences_of("overwrite") > 0,
        })
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
        let slack = config.annotate_read_evidence_slack as usize;
        let slacki = config.annotate_read_evidence_slack as isize;
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
                SeqContigStranded::new(chrom.clone(), f(pos, deltai), delta + slack, Forward),
                SeqContigStranded::new(chrom2.clone(), end2 - slacki, delta + slack, Reverse),
            )],
            ("DUP", _) => vec![(
                SeqContigStranded::new(chrom.clone(), pos - slacki, delta + slack, Reverse),
                SeqContigStranded::new(chrom2.clone(), end2 - deltai, delta + slack, Forward),
            )],
            ("INV", _) => vec![
                (
                    SeqContigStranded::new(chrom.clone(), f(pos, deltai), delta + slack, Forward),
                    SeqContigStranded::new(chrom2.clone(), f(end2, deltai), delta + slack, Forward),
                ),
                (
                    SeqContigStranded::new(chrom.clone(), pos - slacki, delta + slack, Reverse),
                    SeqContigStranded::new(chrom2.clone(), f(end2, deltai), delta + slack, Reverse),
                ),
            ],
            ("BND", "--") | ("BND", "++") => vec![(
                SeqContigStranded::new(chrom.clone(), pos - slacki, delta + slack, Reverse),
                SeqContigStranded::new(chrom2.clone(), end2 - slacki, delta + slack, Reverse),
            )],
            ("BND", "-+") => vec![(
                SeqContigStranded::new(chrom.clone(), pos - slacki, delta + slack, Reverse),
                SeqContigStranded::new(chrom2.clone(), f(end2, deltai), delta + slack, Reverse),
            )],
            ("BND", "+-") => vec![(
                SeqContigStranded::new(chrom.clone(), f(pos, deltai), delta + slack, Forward),
                SeqContigStranded::new(chrom2.clone(), end2 - slacki, delta + slack, Reverse),
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
                let annotation_doc_baf_limit = config.annotation_doc_baf_limit as i64;
                let (start, end) = if record.end2 - record.pos > annotation_doc_baf_limit as i64 {
                    let shift =
                        (record.end2 - record.pos - annotation_doc_baf_limit as i64) as u64 / 2;
                    (record.pos as u64 + shift, record.end2 as u64 - shift)
                } else {
                    (record.pos as u64, record.end2 as u64)
                };

                let rid: u32 = doc_reader.header().name2rid(record.chrom.as_bytes())?;
                if doc_reader.fetch(rid, start, end).is_ok() {
                    let mut doc_windows = Vec::new();
                    while doc_reader.read(&mut doc_record)? {
                        doc_windows.push(DocWindow {
                            cov: doc_record.format(b"RCV").float()?[0][0].into(),
                            mapq: doc_record.format(b"MQ").float()?[0][0].into(),
                        });
                    }

                    let min_mapq = 55.0;
                    let q0_windows = doc_windows.iter().filter(|w| w.mapq >= min_mapq).count();
                    let do_filter = q0_windows >= config.doc_annotation_min_bins;
                    let covs: Vec<f64> = doc_windows
                        .iter()
                        .filter(|w| !do_filter || w.mapq >= min_mapq)
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

/// Het SNV frequency / count evidence
#[derive(Debug, PartialEq, PartialOrd)]
struct SNVEvidence {
    /// Number of SNVs in left flanking region.
    snvs_left: i32,
    /// Number of SNVs within region.
    snvs_within: i32,
    /// Number of SNVs in right flanking region.
    snvs_right: i32,
    /// BAF value
    baf_mean: f32,
}

impl SNVEvidence {
    fn new() -> Self {
        Self {
            snvs_left: 0,
            snvs_within: 0,
            snvs_right: 0,
            baf_mean: std::f32::NAN,
        }
    }
}

/// Perform BAF annotation of SV.
fn annotate_snv(
    options: &Options,
    config: &Config,
    region: &Interval,
) -> Result<Vec<SNVEvidence>, Error> {
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

    let mut baf_reader = bcf::IndexedReader::from_path(options.path_snv_vcf.clone().unwrap())?;
    let mut baf_record = baf_reader.empty_record();

    let mut record = reader.empty_record();
    let mut result = Vec::new();
    while reader.read(&mut record)? {
        let record = sv::StandardizedRecord::from_bcf_record(&mut record)?;

        let snv_evidence = if record.sv_type == "DEL"
            || record.sv_type == "DUP"
            || record.sv_type == "CNV"
        {
            let rid: u32 = baf_reader.header().name2rid(record.chrom.as_bytes())?;
            let annotation_doc_baf_limit = config.annotation_doc_baf_limit as i64;

            let (start, end) = if record.end2 - record.pos > annotation_doc_baf_limit as i64 {
                let shift = (record.end2 - record.pos - annotation_doc_baf_limit as i64) as u64 / 2;
                (record.pos as u64 + shift, record.end2 as u64 - shift)
            } else {
                (record.pos as u64, record.end2 as u64)
            };
            let length = end - start;

            // Count SNVs left of CNV.
            let mut snvs_left: i32 = 0;
            if baf_reader.fetch(rid, start - length, start).is_ok() {
                while baf_reader.read(&mut baf_record)? {
                    let genotype: bcf::record::Genotype = baf_record.genotypes().unwrap().get(0);
                    let gt0 = genotype.get(0).unwrap().index().unwrap_or(0) as usize;
                    let gt1 = genotype.get(1).unwrap().index().unwrap_or(0) as usize;
                    if baf_record.alleles().len() == 2  // only biallelic SNVs
                            && baf_record.alleles()[0].len() == 1
                            && baf_record.alleles()[1].len() == 1
                            && gt0 != gt1
                    {
                        snvs_left += 1;
                    }
                }
            }

            // Count SNVs right of CNV.
            let mut snvs_right: i32 = 0;
            if baf_reader.fetch(rid, start - length, start).is_ok() {
                while baf_reader.read(&mut baf_record)? {
                    let genotype: bcf::record::Genotype = baf_record.genotypes().unwrap().get(0);
                    let gt0 = genotype.get(0).unwrap().index().unwrap_or(0) as usize;
                    let gt1 = genotype.get(1).unwrap().index().unwrap_or(0) as usize;
                    if baf_record.alleles().len() == 2  // only biallelic SNVs
                            && baf_record.alleles()[0].len() == 1
                            && baf_record.alleles()[1].len() == 1
                            && gt0 != gt1
                    {
                        snvs_right += 1;
                    }
                }
            }

            if baf_reader.fetch(rid, start, end).is_ok() {
                let mut bafs = Vec::new();
                while baf_reader.read(&mut baf_record)? {
                    let genotype: bcf::record::Genotype = baf_record.genotypes().unwrap().get(0);
                    let gt0 = genotype.get(0).unwrap().index().unwrap_or(0) as usize;
                    let gt1 = genotype.get(1).unwrap().index().unwrap_or(0) as usize;
                    if baf_record.alleles().len() == 2  // only biallelic SNVs
                            && baf_record.alleles()[0].len() == 1
                            && baf_record.alleles()[1].len() == 1
                            && gt0 != gt1
                    {
                        let a0 = baf_record.format(b"AD").integer()?[0][gt0] as f64;
                        let a1 = baf_record.format(b"AD").integer()?[0][gt1] as f64;
                        debug!("{} -- {} // {}", baf_record.pos(), a0, a1);
                        bafs.push(a1 / (a0 + a1));
                    }
                }

                SNVEvidence {
                    snvs_left,
                    snvs_right,
                    snvs_within: bafs.len() as i32,
                    baf_mean: if bafs.is_empty() {
                        std::f32::NAN
                    } else {
                        bafs.mean() as f32
                    },
                }
            } else {
                SNVEvidence::new()
            }
        } else {
            SNVEvidence::new()
        };

        result.push(snv_evidence);
    }

    Ok(result)
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

/// Write annotated variants.
fn write_annotated(
    options: &Options,
    _config: &Config,
    region: &Interval,
    read_evidence: &Option<Vec<ReadEvidenceCount>>,
    doc_evidence: &Option<Vec<CoverageEvidence>>,
    baf_evidence: &Option<Vec<SNVEvidence>>,
    writer: &mut bcf::Writer,
) -> Result<(), Error> {
    let mut reader = bcf::IndexedReader::from_path(&options.path_input)?;
    let res = reader.fetch(
        reader.header().name2rid(region.contig().as_bytes())?,
        region.range().start,
        region.range().end,
    );
    match res {
        Err(rust_htslib::bcf::errors::Error::Seek { .. }) => return Ok(()),
        Err(_) => res?,
        Ok(_) => (),
    }

    let mut idx = 0;
    let mut record = reader.empty_record();
    while reader.read(&mut record)? {
        writer.translate(&mut record);

        if let Some(evidence) = read_evidence {
            let elem = evidence.get(idx).unwrap();
            record.push_format_float(b"PR", &[elem.pe_count as f32])?;
            record.push_format_float(b"SR", &[elem.sr_count as f32])?;
        }
        if let Some(evidence) = doc_evidence {
            let elem = evidence.get(idx).unwrap();
            record.push_format_float(b"RD", &[elem.norm_cov as f32])?;
        }
        if let Some(evidence) = baf_evidence {
            let elem = evidence.get(idx).unwrap();
            record.push_format_integer(b"VL", &[elem.snvs_left as i32])?;
            record.push_format_integer(b"VM", &[elem.snvs_within as i32])?;
            record.push_format_integer(b"VR", &[elem.snvs_right as i32])?;
            if !elem.baf_mean.is_nan() {
                record.push_format_float(b"BF", &[elem.baf_mean as f32])?;
            }

            let length = record.info(b"SVLEN").integer()?.unwrap()[0];
            let min_snvs = if length >= 100_000 {
                50
            } else {
                5 * (length / 10_000)
            };
            if (elem.snvs_right < min_snvs || elem.snvs_left < min_snvs)
                && elem.snvs_within < min_snvs
            {
                record.push_format_integer(b"ROH", &[1])?;
            }
        }

        writer.write(&record)?;

        idx += 1;
    }

    Ok(())
}

/// Main entry point after parsing command line and loading options.
fn perform_annotation(options: &Options, config: &Config) -> Result<(), Error> {
    info!("Starting to annotate variants for sample...");

    let median_doc = if let Some(path_doc_evidence) = &options.path_doc_evidence {
        info!("Computing median depth of coverage (DoC)...");
        let median_doc = load_doc_median(path_doc_evidence)?;
        info!("... median DoC is {}", median_doc);
        Some(median_doc)
    } else {
        None
    };

    let reader = bcf::IndexedReader::from_path(&options.path_input)?;
    let header = build_vcf_header(reader.header())?;
    let guessed = guess_bcf_format(&options.path_output);
    let mut writer = bcf::Writer::from_path(
        &options.path_output,
        &header,
        guessed.uncompressed,
        guessed.format,
    )?;

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

        let read_evidence = read_evidence
            .as_ref()
            .map(|re| annotate_pesr(&options, &config, &re, &region, &blocked))
            .transpose()?;
        let doc_evidence = options
            .path_doc_evidence
            .as_ref()
            .map(|_| annotate_doc(&options, &config, &region, median_doc.unwrap()))
            .transpose()?;
        let baf_evidence = options
            .path_snv_vcf
            .as_ref()
            .map(|_| annotate_snv(&options, &config, &region))
            .transpose()?;

        write_annotated(
            &options,
            &config,
            &region,
            &read_evidence,
            &doc_evidence,
            &baf_evidence,
            &mut writer,
        )?;
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
            Arg::from_usage("<output> 'output VCF file"),
        ])
        .get_matches();
    let options = Options::from_arg_matches(&matches)?;

    // Output files must not exist yet.
    if options.path_pesr_evidence.is_some()
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
        path_expected: &str,
        regions: &Option<Vec<Interval>>,
        config_text: &str,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());
        let options = super::Options {
            verbosity: 1, // disable progress bar
            regions: regions.clone(),
            path_config: None,
            path_pesr_evidence: path_pesr_evidence,
            path_doc_evidence: path_doc_evidence,
            path_snv_vcf: path_snv_vcf,
            sample: sample.to_string(),
            path_input: path_input.to_string(),
            path_output: path_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str(config_text).unwrap();

        super::perform_annotation(&options, &config)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(path_output).unwrap(),
        );

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
            Some(String::from("./src/tests/data/ex-delly-snvs.vcf.gz")),
            "./src/tests/data/ex-delly-svs.vcf.gz",
            "./src/tests/data/ex-delly.expected.vcf",
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
            Some(String::from("./src/tests/data/ex-delly-snvs.vcf.gz")),
            "./src/tests/data/ex-delly-svs.vcf.gz",
            "./src/tests/data/ex-delly.expected-blocked.vcf",
            &None,
            "blocked_regions_bed = \"./src/tests/data/ex-delly-blocked.bed\"",
        )?;
        Ok(())
    }
}
