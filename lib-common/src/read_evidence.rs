use std::io::prelude::*;

use rust_htslib::{tbx, tbx::Read};
use serde::{Deserialize, Serialize};

use super::error;

/// Strand.
#[derive(Debug, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Sides.
#[derive(Debug, Serialize, Deserialize)]
pub enum Sides {
    Left,
    Right,
    Both,
    Neither,
}

/// Read pair/split read annotation from one read alignment.
#[derive(Debug, Serialize, Deserialize)]
pub enum Record {
    /// Paired read based evidence.
    PairedRead {
        /// Name of the read pair.
        read_name: String,
        /// Whether this alignment is first in pair.
        is_first1: bool,
        /// Contig of this alignment.
        contig1: String,
        /// Start position on this contig.
        start1: i64,
        /// End position on this contig.
        end1: i64,
        /// Read orientation in this contig.
        strand1: Strand,
        /// Other alignment's contig.
        contig2: Option<String>,
        /// Start position on other contig.
        start2: Option<i64>,
        /// Strand on other contig.
        strand2: Option<Strand>,
        /// Template size.
        tlen: Option<i64>,
    },
    /// Split read based evidence.
    SplitRead {
        /// Name of the read pair.
        read_name: String,
        /// Whether is first read in this alignment.
        is_first: bool,
        /// The alignment's contig.
        contig: String,
        /// The alignment's start position.
        start: i64,
        /// The alignment's end position.
        end: i64,
        /// The side that the read has been clipped on.
        clipped_sides: Sides,
    },
}

impl Record {
    pub fn interval(&self) -> std::ops::Range<i64> {
        match self {
            Self::PairedRead { start1, end1, .. } => *start1..*end1,
            Self::SplitRead { start, end, .. } => *start..*end,
        }
    }
}

pub struct Writer {
    handle: std::fs::File,
}

impl Writer {
    pub fn from_path(path: &str) -> Result<Self, std::io::Error> {
        let mut handle = std::fs::File::create(path)?;
        handle.write_all(b"#contig\tstart\tend\tsignal\n")?;
        Ok(Self { handle })
    }

    pub fn write(&mut self, e: &Record) -> Result<(), std::io::Error> {
        match e {
            Record::PairedRead {
                contig1,
                start1,
                end1,
                ..
            } => {
                self.handle
                    .write_all(format!("{}\t{}\t{}\t", &contig1, start1, end1).as_bytes())?;
            }
            Record::SplitRead {
                contig, start, end, ..
            } => {
                self.handle
                    .write_all(format!("{}\t{}\t{}\t", &contig, start, end).as_bytes())?;
            }
        }
        self.handle
            .write_all(serde_json::to_string(&e)?.as_bytes())?;
        self.handle.write_all(b"\n")?;

        Ok(())
    }
}

pub struct IndexedReader {
    inner: tbx::Reader,
    buffer: Vec<u8>,
}

impl IndexedReader {
    pub fn from_path(path: &str) -> Result<Self, error::Error> {
        Ok(Self {
            inner: tbx::Reader::from_path(&path)?,
            buffer: Vec::new(),
        })
    }

    pub fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<(), error::Error> {
        Ok(self.inner.fetch(self.inner.tid(&contig)?, start, end)?)
    }

    pub fn read_record(&mut self) -> Result<Option<Record>, error::Error> {
        if self.inner.read(&mut self.buffer)? {
            let arr: Vec<&[u8]> = self.buffer.split(|&c| c == b'\t').collect();
            Ok(serde_json::de::from_slice(arr[3])?)
        } else {
            Ok(None)
        }
    }
}
