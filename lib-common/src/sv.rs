use std::cmp::Ordering;
use std::collections::HashMap;

use rust_htslib::bcf;

use super::error::Error;

/// Representation of a structural variant record suitable for clustering.
#[derive(Debug, Clone)]
pub struct StandardizedRecord {
    pub chrom: String,
    pub pos: i64,
    pub reference: String,
    pub alt: String,
    pub filters: Vec<String>,
    pub chrom2: String,
    pub end2: i64,
    pub sv_type: String,
    pub strands: String,
    pub sv_len: i64,
    pub algorithms: Vec<String>,
    pub samples: Vec<String>,
    pub gts: Vec<String>,
    pub called_by: Vec<Vec<String>>,
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

impl Default for StandardizedRecord {
    fn default() -> Self {
        Self::new()
    }
}

impl StandardizedRecord {
    pub fn new() -> Self {
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

    pub fn from_bcf_record(record: &mut bcf::Record) -> Result<Self, Error> {
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

        let strands = if record.info(b"STRANDS").string()?.is_some() {
            String::from_utf8(record.info(b"STRANDS").string()?.unwrap()[0].to_vec())?
        } else {
            ".".to_string()
        };

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
            strands,
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

    pub fn update_bcf_record(&self, record: &mut bcf::Record) -> Result<(), Error> {
        let sample_count = record.header().sample_count() as usize;

        record.set_rid(Some(record.header().name2rid(self.chrom.as_bytes())?));
        record.set_pos(self.pos);
        let alleles: Vec<&[u8]> = vec![&self.reference.as_bytes(), &self.alt.as_bytes()];
        record.set_alleles(&alleles)?;
        for filter in &self.filters {
            record.push_filter(record.header().name_to_id(&filter.as_bytes())?);
        }
        record.push_info_integer(b"END2", vec![self.end2 as i32].as_slice())?;
        record.push_info_string(b"CHR2", vec![self.chrom2.as_bytes()].as_slice())?;
        record.push_info_string(b"SVTYPE", vec![self.sv_type.as_bytes()].as_slice())?;
        record.push_info_string(b"STRANDS", vec![self.strands.as_bytes()].as_slice())?;
        record.push_info_integer(b"SVLEN", vec![self.sv_len as i32].as_slice())?;
        let algorithms: Vec<&[u8]> = self.algorithms.iter().map(|a| a.as_bytes()).collect();
        record.push_info_string(b"ALGORITHMS", algorithms.as_slice())?;

        // Prepare buffers for FORMAT fields.
        let mut called_by: HashMap<String, Vec<i32>> = HashMap::new();
        for algorithm in &self.algorithms {
            called_by.insert(algorithm.clone(), vec![0; sample_count as usize]);
        }

        // Write information from "sparse" StandardizedRecord into BCF record.
        let mut genotypes = vec![0; 2 * sample_count as usize];
        for (i, sample) in self.samples.iter().enumerate() {
            let sample_idx = *record.header().sample_to_id(sample.as_bytes())? as usize;

            genotypes[2 * sample_idx] = (self.gts[i].starts_with('1') as i32 + 1) << 1;
            genotypes[2 * sample_idx + 1] = (self.gts[i].ends_with('1') as i32 + 1) << 1;
            for algorithm in &self.called_by[i] {
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

        Ok(())
    }

    pub fn interval(&self) -> std::ops::Range<i64> {
        self.extended_interval(0)
    }

    pub fn extended_interval(&self, delta: i64) -> std::ops::Range<i64> {
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

pub fn overlap(lhs: &std::ops::Range<i64>, rhs: &std::ops::Range<i64>) -> std::ops::Range<i64> {
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rust_htslib::{bcf, bcf::Read};
    use std::fs;
    use tempdir::TempDir;

    #[test]
    fn test_from_bcf_record() -> Result<(), super::Error> {
        let mut reader = bcf::Reader::from_path("./src/tests/data/ex-delly-1.vcf")?;
        let mut bcf_record = reader.empty_record();
        assert!(reader.read(&mut bcf_record)?);

        let record = super::StandardizedRecord::from_bcf_record(&mut bcf_record)?;

        assert_eq!(
            record,
            super::StandardizedRecord {
                chrom: "1".to_owned(),
                pos: 10_002,
                reference: "A".to_owned(),
                alt: "<DEL>".to_owned(),
                filters: Vec::new(),
                chrom2: "1".to_owned(),
                end2: 20_000,
                sv_type: "DEL".to_owned(),
                strands: "-+".to_owned(),
                sv_len: 10_003,
                algorithms: vec!["delly".to_owned()],
                samples: vec!["sample-1".to_owned()],
                gts: vec!["0/1".to_owned()],
                called_by: vec![vec!["delly".to_owned()]],
            }
        );

        Ok(())
    }

    #[test]
    fn test_update_bcf_record() -> Result<(), super::Error> {
        let reader = bcf::Reader::from_path("./src/tests/data/ex-delly-1.vcf")?;

        let header = bcf::Header::from_template(reader.header());
        let tmp_dir = TempDir::new("tests")?;
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());

        let record = super::StandardizedRecord {
            chrom: "1".to_owned(),
            pos: 10_002,
            reference: "A".to_owned(),
            alt: "<DEL>".to_owned(),
            filters: Vec::new(),
            chrom2: "1".to_owned(),
            end2: 20_000,
            sv_type: "DEL".to_owned(),
            strands: "-+".to_owned(),
            sv_len: 10_003,
            algorithms: vec!["delly".to_owned()],
            samples: vec!["sample-1".to_owned()],
            gts: vec!["0/1".to_owned()],
            called_by: vec![vec!["delly".to_owned()]],
        };

        {
            let mut writer = bcf::Writer::from_path(&path_output, &header, true, bcf::Format::VCF)?;
            let mut bcf_record: bcf::Record = writer.empty_record();
            record.update_bcf_record(&mut bcf_record)?;
            writer.write(&bcf_record)?;
        }

        assert_eq!(
            fs::read_to_string("./src/tests/data/ex-delly-1.vcf").unwrap(),
            fs::read_to_string(&path_output).unwrap()
        );

        Ok(())
    }
}
