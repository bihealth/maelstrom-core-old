/// vcf-standardize -- Extract and standardize tool output VCF to standardized VCF.
use std::fs;
use std::path::Path;

use bio_types::genome::{AbstractInterval, Interval};
use clap::{App, Arg, ArgMatches};
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn, LevelFilter};
use regex::Regex;
use rust_htslib::{bcf, bcf::Read};

use lib_common::bcf::{build_vcf_header, collect_contigs, guess_bcf_format};
use lib_common::error::Error;
use lib_common::{parse_region, Algorithm};
use lib_config::Config;

/// Command line options
#[derive(Debug)]
struct Options {
    /// Verbosity level
    verbosity: u64,
    /// The SV calling tool used for the input file.
    tool: Algorithm,
    /// List of regions to call.
    regions: Option<Vec<Interval>>,
    /// Path to configuration file to use,
    path_config: Option<String>,
    /// Path to input file.
    path_input: String,
    /// Path to output file.
    path_output: String,
    /// Overwrite output file.
    overwrite: bool,
}

impl Options {
    pub fn from_arg_matches<'a>(matches: &ArgMatches<'a>) -> Result<Self, Error> {
        Ok(Options {
            verbosity: matches.occurrences_of("v"),
            tool: match matches.value_of("tool").ok_or(Error::OptionMissing())? {
                "cnmops" => Algorithm::CNMOPS,
                "delly" => Algorithm::Delly,
                "manta" => Algorithm::Manta,
                _ => return Err(Error::OptionMissing()),
            },
            regions: matches
                .value_of("regions")
                .map(|s| {
                    let x: Result<Vec<Interval>, Error> =
                        s.split(',').map(|t| parse_region(&t)).collect();
                    x
                })
                .transpose()?,
            path_config: matches.value_of("config").map(|s| s.to_string()),
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

/// cnMOPS-specific part of summarizing a single BCF record.
fn summarize_record_cnmops(
    src: &mut bcf::Record,
    _src_header: &bcf::header::HeaderView,
    dst: &mut bcf::Record,
    dst_header: &bcf::header::HeaderView,
) -> Result<bool, Error> {
    dst.set_alleles(&src.alleles())?;

    // INFO/END2
    let end2 = src.info(b"END").integer().expect("No END?").unwrap()[0];
    dst.push_info_integer(b"END2", &[end2])?;
    // INFO/CHR2
    dst.push_info_string(
        b"CHR2",
        &[&dst_header.rid2name(dst.rid().expect("No REF?"))?],
    )?;
    // INFO/SVTYPE
    let sv_type = String::from(std::str::from_utf8(
        src.info(b"SVTYPE").string()?.unwrap()[0],
    )?);
    dst.push_info_string(b"SVTYPE", &[sv_type.as_bytes()])?;

    // NB: no INFO/STRANDS
    // INFO/SVLEN
    dst.push_info_integer(b"SVLEN", &[end2 - (src.pos() as i32)])?;
    // INFO/ALGORITHMS
    dst.push_info_string(b"ALGORITHMS", &[b"cnmops"])?;

    // FORMAT/GT
    let sample_count = dst_header.sample_count() as usize;
    let mut gts: Vec<i32> = Vec::new();
    src.format(b"GT").integer()?[..sample_count]
        .iter()
        .for_each(|xs| gts.extend_from_slice(xs));
    dst.push_format_integer(b"GT", &gts)?;

    // FORMAT/cnmops
    // For each sample, collect wehther there is any called allele.
    let called: Vec<i32> = (0..sample_count)
        .map(|i| {
            let ith_genotype: bcf::record::Genotype = src.genotypes().unwrap().get(i as usize);
            let any_var: bool = ith_genotype.iter().any(|gt_allele| match gt_allele {
                bcf::record::GenotypeAllele::Unphased(i)
                | bcf::record::GenotypeAllele::Phased(i) => (*i != 0),
                _ => false,
            });
            if any_var {
                1
            } else {
                0
            }
        })
        .collect();
    dst.push_format_integer(b"cnmops", &called)?;

    Ok(true)
}

/// Delly-specific part of summarizing a single BCF record.
fn summarize_record_delly(
    src: &mut bcf::Record,
    _src_header: &bcf::header::HeaderView,
    dst: &mut bcf::Record,
    dst_header: &bcf::header::HeaderView,
) -> Result<bool, Error> {
    dst.set_alleles(&src.alleles())?;

    // INFO/END2
    let end2 = match src.info(b"END2").integer() {
        Ok(Some(ends)) => ends[0],
        _ => match src.info(b"END").integer() {
            Ok(Some(ends)) => ends[0],
            _ => panic!("Could not read END or END2"),
        },
    };
    dst.push_info_integer(b"END2", &[end2])?;
    // INFO/CHR2
    dst.push_info_string(b"CHR2", &src.info(b"CHR2").string()?.unwrap())?;
    // INFO/SVTYPE
    let sv_type = String::from(std::str::from_utf8(
        src.info(b"SVTYPE").string()?.unwrap()[0],
    )?);
    if sv_type == "INS" {
        return Ok(false); // don't care about INS
    }
    dst.push_info_string(b"SVTYPE", &[sv_type.as_bytes()])?;

    // INFO/STRANDS
    let raw_strands = std::str::from_utf8(src.info(b"CT").string()?.unwrap()[0])?;
    let strands = match raw_strands {
        "5to3" => "-+",
        "3to5" => "+-",
        "5to5" => "--",
        "3to3" => "++",
        "NtoN" => {
            if sv_type == "INS" {
                "+-"
            } else {
                panic!("improper strands")
            }
        }
        _ => panic!("improper strands"),
    };
    dst.push_info_string(b"STRANDS", &[strands.as_bytes()])?;
    // INFO/SVLEN
    let svtype = std::str::from_utf8(src.info(b"SVTYPE").string()?.unwrap()[0])?;
    if svtype == "BND" {
        dst.push_info_integer(b"SVLEN", &[-1])?;
    } else {
        dst.push_info_integer(b"SVLEN", &[end2 - (src.pos() as i32)])?;
    }
    // INFO/ALGORITHMS
    dst.push_info_string(b"ALGORITHMS", &[b"delly"])?;

    // FORMAT/GT
    let sample_count = dst_header.sample_count() as usize;
    let mut gts: Vec<i32> = Vec::new();
    src.format(b"GT").integer()?[..sample_count]
        .iter()
        .for_each(|xs| gts.extend_from_slice(xs));
    dst.push_format_integer(b"GT", &gts)?;

    // FORMAT/delly
    // For each sample, collect wehther there is any called allele.
    let called: Vec<i32> = (0..sample_count)
        .map(|i| {
            let ith_genotype: bcf::record::Genotype = src.genotypes().unwrap().get(i as usize);
            let any_var: bool = ith_genotype.iter().any(|gt_allele| match gt_allele {
                bcf::record::GenotypeAllele::Unphased(i)
                | bcf::record::GenotypeAllele::Phased(i) => (*i != 0),
                _ => false,
            });
            if any_var {
                1
            } else {
                0
            }
        })
        .collect();
    dst.push_format_integer(b"delly", &called)?;

    Ok(true)
}

/// Parse Breakend BND position string.
fn parse_bnd_pos(alt: &str) -> Result<(&[u8], i32), Error> {
    let re = Regex::new(r".*[\[\]](.+):(.+)[\[\]].*").unwrap();
    let cap = re
        .captures(alt)
        .unwrap_or_else(|| panic!("Could not match breakend: {}", &alt));

    Ok((
        cap.get(1).unwrap().as_str().as_bytes(),
        cap.get(2).unwrap().as_str().parse()?,
    ))
}

/// Parse Breakend BND position string to strands.
fn parse_bnd_strands(alt: &str) -> String {
    if alt.ends_with('[') {
        "+-".to_string()
    } else if alt.ends_with(']') {
        "++".to_string()
    } else if alt.starts_with(']') {
        "-+".to_string()
    } else if alt.starts_with('[') {
        "--".to_string()
    } else {
        panic!("Unexpected alt string: {}", &alt);
    }
}

/// Delly-specific part of summarizing a single BCF record.
fn summarize_record_manta(
    src: &mut bcf::Record,
    src_header: &bcf::header::HeaderView,
    dst: &mut bcf::Record,
    dst_header: &bcf::header::HeaderView,
) -> Result<bool, Error> {
    // TODO: skip "mated" records

    // Obtain the SVTYPE, normalize, skip in case of insertions.
    let sv_type = std::str::from_utf8(src.info(b"SVTYPE").string()?.unwrap()[0])?.to_string();
    let sv_type = match &sv_type[..] {
        "DUP:TANDEM" => "DUP".to_string(),
        "INS" => return Ok(false), // don't care about INS
        _ => sv_type,
    };

    // Extract reference and alternative allele.
    let ref_allele = std::str::from_utf8(src.alleles()[0])?[0..1].to_owned();
    let alt_allele = std::str::from_utf8(src.alleles()[1])?.to_owned();

    // Compute value of INFO/CHR2 and END2.
    let (chr2, end2) = match &sv_type[..] {
        "BND" => parse_bnd_pos(&alt_allele)?,
        "INS" => return Ok(false), // don't care about INS
        _ => (
            src_header.rid2name(src.rid().expect("No REF?"))?,
            src.info(b"END")
                .integer()?
                .expect("Could not read INFO/END")[0],
        ),
    };

    // Set REF and ALT alleles.
    if sv_type == "BND" {
        dst.set_alleles(&[&ref_allele.as_bytes(), &alt_allele.as_bytes()])?;
    } else {
        dst.set_alleles(&[&ref_allele.as_bytes(), format!("<{}>", &sv_type).as_bytes()])?;
    }

    // INFO/END2
    dst.push_info_integer(b"END2", &[end2])?;
    // INFO/CHR2
    dst.push_info_string(b"CHR2", &[chr2])?;
    // INFO/SVTYPE
    dst.push_info_string(b"SVTYPE", &[sv_type.as_bytes()])?;

    // INFO/STRANDS
    let strands: String = match &sv_type[..] {
        "INV" => match src.info(b"INV3").flag()? {
            true => "++".to_string(),
            false => "--".to_string(),
        },
        "BND" => parse_bnd_strands(&alt_allele),
        "DEL" => "+-".to_string(),
        "DUP" => "-+".to_string(),
        "INS" => "+-".to_string(),
        _ => panic!("Unexpected sv_type: {}", &sv_type),
    };
    dst.push_info_string(b"STRANDS", &[strands.as_bytes()])?;
    // INFO/SVLEN
    let svtype = std::str::from_utf8(src.info(b"SVTYPE").string()?.unwrap()[0])?;
    if svtype == "BND" {
        dst.push_info_integer(b"SVLEN", &[-1])?;
    } else {
        dst.push_info_integer(b"SVLEN", &[end2 - (src.pos() as i32)])?;
    }
    // INFO/ALGORITHMS
    dst.push_info_string(b"ALGORITHMS", &[b"manta"])?;

    // FORMAT/GT
    let sample_count = dst_header.sample_count() as usize;
    let mut gts: Vec<i32> = Vec::new();
    src.format(b"GT").integer()?[..sample_count]
        .iter()
        .for_each(|xs| gts.extend_from_slice(xs));
    dst.push_format_integer(b"GT", &gts)?;

    // FORMAT/manta
    // For each sample, collect wehther there is any called allele.
    let called: Vec<i32> = (0..sample_count)
        .map(|i| {
            let ith_genotype: bcf::record::Genotype = src.genotypes().unwrap().get(i as usize);
            let any_var: bool = ith_genotype.iter().any(|gt_allele| match gt_allele {
                bcf::record::GenotypeAllele::Unphased(i)
                | bcf::record::GenotypeAllele::Phased(i) => (*i != 0),
                _ => false,
            });
            if any_var {
                1
            } else {
                0
            }
        })
        .collect();
    dst.push_format_integer(b"manta", &called)?;

    Ok(true)
}

/// Summarize a single BCF record.
fn summarize_record(
    src: &mut bcf::Record,
    src_header: &bcf::header::HeaderView,
    dst: &mut bcf::Record,
    dst_header: &bcf::header::HeaderView,
    config: &Config,
    tool: &Algorithm,
) -> Result<bool, Error> {
    if config.stdvcf_apply_filters
        && src
            .filters()
            .map(|id| String::from_utf8(src.header().id_to_name(id)).unwrap())
            .any(|f| f != "PASS")
    {
        return Ok(false);
    }

    dst.set_rid(Some(
        dst_header.name2rid(src_header.rid2name(src.rid().unwrap())?)?,
    ));
    dst.set_pos(src.pos());

    match tool {
        Algorithm::CNMOPS => summarize_record_cnmops(src, src_header, dst, dst_header),
        Algorithm::Delly => summarize_record_delly(src, src_header, dst, dst_header),
        Algorithm::Manta => summarize_record_manta(src, src_header, dst, dst_header),
    }
}

/// Main entry point after parsing command line and loading options.
fn perform_extraction(options: &Options, config: &Config) -> Result<(), Error> {
    info!("Starting to extract BCF file...");
    let mut reader = bcf::IndexedReader::from_path(&options.path_input)?;
    let header = build_vcf_header(reader.header())?;
    let guessed = guess_bcf_format(&options.path_output);
    let mut writer = bcf::Writer::from_path(
        &options.path_output,
        &header,
        guessed.uncompressed,
        guessed.format,
    )?;

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

    let regions = if let Some(regions) = &options.regions {
        regions.clone()
    } else {
        collect_contigs(&reader)?
            .iter()
            .map(|name| Interval::new(name.clone(), 0..10_000_000_000))
            .collect()
    };

    let mut counter: usize = 0;
    let mut buffer_read = reader.empty_record();
    for region in &regions {
        let rid = reader.header().name2rid(region.contig().as_bytes())?;
        if reader
            .fetch(rid, region.range().start, region.range().end)
            .is_err()
        {
            warn!("No data for {:?}", &region);
        }

        while reader.read(&mut buffer_read)? {
            let mut buffer_write = writer.empty_record();
            if summarize_record(
                &mut buffer_read,
                reader.header(),
                &mut buffer_write,
                writer.header(),
                &config,
                &options.tool,
            )? {
                writer.write(&buffer_write)?;
            }

            counter += 1;
            if counter % 1_000 == 0 {
                if let Some(spinner) = &spinner {
                    spinner.set_message(&format!(
                        "currently at {}:{}",
                        &std::str::from_utf8(
                            reader.header().rid2name(buffer_read.rid().unwrap())?
                        )?,
                        buffer_read.pos() + 1
                    ));
                }
            }
        }
    }
    if let Some(spinner) = &spinner {
        spinner.finish();
    }

    info!("Done extracting BCF file...");
    Ok(())
}

fn main() -> Result<(), Error> {
    // Setup command line parser and parse options.
    let matches = App::new("maelstrom-stdvcf")
        .version(git_version!())
        .author("Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>")
        .about("Extract and standardize records from tool VCF files")
        .args(&[
            Arg::from_usage("-v... 'Increase verbosity'"),
            Arg::from_usage("-t, --tool=<cnmops|delly|manta> 'Name of SV calling tool'"),
            Arg::from_usage("--overwrite 'Allow overwriting of output file'"),
            Arg::from_usage("-c, --config=[FILE] 'Sets a custom config file'"),
            Arg::from_usage("-r, --regions=[REGIONS] 'comma-separated list of regions'"),
            Arg::from_usage("<input> 'input file to read from'"),
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
        .level(if matches.is_present("v") {
            LevelFilter::Debug
        } else {
            LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();
    info!("Starting maelstrom-stdvcf");
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
    perform_extraction(&options, &config)?;

    info!("All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{Algorithm, Interval};
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempdir::TempDir;

    /// Helper that runs `perform_extraction()` and compares the result.
    fn _perform_extraction_and_test(
        tmp_dir: &TempDir,
        path_input: &str,
        path_expected: &str,
        config_str: &str,
        regions: &Option<Vec<Interval>>,
        algorithm: &Algorithm,
    ) -> Result<(), super::Error> {
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());
        let options = super::Options {
            verbosity: 1, // disable progress bar
            tool: algorithm.clone(),
            regions: regions.clone(),
            path_config: None,
            path_input: String::from(path_input),
            path_output: path_output.clone(),
            overwrite: false,
        };
        let config: super::Config = toml::from_str(config_str).unwrap();

        super::perform_extraction(&options, &config)?;

        assert_eq!(
            fs::read_to_string(path_expected).unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }

    #[test]
    fn test_stdvcf_delly2_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-delly-filter.vcf.gz",
            "./src/tests/data/ex-delly-filter.expected.vcf",
            "stdvcf_apply_filters = true",
            &None,
            &Algorithm::Delly,
        )?;
        Ok(())
    }

    #[test]
    fn test_stdvcf_delly2_nofilter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-delly-nofilter.vcf.gz",
            "./src/tests/data/ex-delly-nofilter.expected.vcf",
            "stdvcf_apply_filters = false",
            &None,
            &Algorithm::Delly,
        )?;
        Ok(())
    }
    #[test]
    fn test_stdvcf_manta_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-manta-filter.vcf.gz",
            "./src/tests/data/ex-manta-filter.expected.vcf",
            "stdvcf_apply_filters = true",
            &None,
            &Algorithm::Manta,
        )?;
        Ok(())
    }

    #[test]
    fn test_stdvcf_manta_nofilter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-manta-nofilter.vcf.gz",
            "./src/tests/data/ex-manta-nofilter.expected.vcf",
            "stdvcf_apply_filters = false",
            &None,
            &Algorithm::Manta,
        )?;
        Ok(())
    }
    #[test]
    fn test_stdvcf_cnmops_filter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-cnmops-filter.vcf.gz",
            "./src/tests/data/ex-cnmops-filter.expected.vcf",
            "stdvcf_apply_filters = true",
            &None,
            &Algorithm::CNMOPS,
        )?;
        Ok(())
    }

    #[test]
    fn test_stdvcf_cnmops_nofilter() -> Result<(), super::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_extraction_and_test(
            &tmp_dir,
            "./src/tests/data/ex-cnmops-nofilter.vcf.gz",
            "./src/tests/data/ex-cnmops-nofilter.expected.vcf",
            "stdvcf_apply_filters = false",
            &None,
            &Algorithm::CNMOPS,
        )?;
        Ok(())
    }
}
