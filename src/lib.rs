//! # bcf_reader
//! This is an attempt to create a small, lightweight, pure-Rust library to allow
//! efficient, cross-platform access to genotype data in BCF files.
//!
//! Currently, the `rust_htslib` crate works only on Linux and macOS (not Windows?).
//! The noodles crate is a pure Rust library for many bioinformatic file formats and
//! works across Windows, Linux, and macOS. However, the `noodles` API for reading
//! genotype data from BCF files can be slow due to its memory allocation patterns.
//! Additionally, both crates have a large number of dependencies, as they provide
//! many features and support a wide range of file formats.
//!
//! One way to address the memory allocation and dependency issues is to manually
//! parse BCF records according to its specification
//! (`<https://samtools.github.io/hts-specs/VCFv4.2.pdf>`) and use iterators whenever
//! possible, especially for the per-sample fields, like GT and AD.
//!
//! Note: This crate is in its early stages of development.
//!
//! ## Usage
//! ```
//! use bcf_reader::*;
//! let mut reader = smart_reader("testdata/test2.bcf").unwrap();
//! let header = Header::from_string(&read_header(&mut reader).unwrap()).unwrap();
//! // find key for a field in INFO or FORMAT or FILTER
//! let key = header.get_idx_from_dictionary_str("FORMAT", "GT").unwrap();
//! // access header dictionary
//! let d = &header.dict_strings()[&key];
//! assert_eq!(d["ID"], "GT");
//! assert_eq!(d["Dictionary"], "FORMAT");
//! /// get chromosome name
//! assert_eq!(header.get_chrname(0), "Pf3D7_01_v3");
//! let fmt_ad_key = header
//!     .get_idx_from_dictionary_str("FORMAT", "AD")
//!     .expect("FORMAT/AD not found");
//! let info_af_key = header
//!     .get_idx_from_dictionary_str("INFO", "AF")
//!     .expect("INFO/AF not found");
//!
//! // this can be and should be reused to reduce allocation
//! let mut record = Record::default();
//! while let Ok(_) = record.read(&mut reader) {
//!     let pos = record.pos();
//!
//!     // use byte ranges and shared buffer to get allele string values
//!     let allele_byte_ranges = record.alleles();
//!     let share_buf = record.buf_shared();
//!     let ref_rng = &allele_byte_ranges[0];
//!     let ref_allele_str = std::str::from_utf8(&share_buf[ref_rng.clone()]).unwrap();
//!     let alt1_rng = &allele_byte_ranges[1];
//!     let alt1_allele_str = std::str::from_utf8(&share_buf[alt1_rng.clone()]).unwrap();
//!     // ...
//!
//!     // access FORMAT/GT via iterator
//!     for nv in record.fmt_gt(&header) {
//!         let nv = nv.unwrap();
//!         let (has_no_ploidy, is_missing, is_phased, allele_idx) = nv.gt_val();
//!         // ...
//!     }
//!
//!     // access FORMAT/AD via iterator
//!     for nv in record.fmt_field(fmt_ad_key) {
//!         let nv = nv.unwrap();
//!         match nv.int_val() {
//!             None => {}
//!             Some(ad) => {
//!                 // ...
//!             }
//!         }
//!         // ...
//!     }
//!
//!     // access FILTERS via itertor
//!     record.filters().for_each(|nv| {
//!         let nv = nv.unwrap();
//!         let filter_key = nv.int_val().unwrap() as usize;
//!         let dict_string_map = &header.dict_strings()[&filter_key];
//!         let filter_name = &dict_string_map["ID"];
//!         // ...
//!     });
//!
//!     // access INFO/AF via itertor
//!     record.info_field_numeric(info_af_key).for_each(|nv| {
//!         let nv = nv.unwrap();
//!         let af = nv.float_val().unwrap();
//!         // ...
//!     });
//! }
//! ```
//!
//! More examples to access each field/column are available in docs of [`Record`] and [`Header`].
//!
//! # Reader types
//! - For parallelized decompression reader, see [`BcfReader`].
//! - For parallelized indexed reader, see [ `IndexedBcfReader`].
//! - For the Lower-level reader underlying `BcfReader` and `IndexedBcfReader`,
//!   see [`ParMultiGzipReader`].
//!
//! # `flate2` backends
//!
//! By default, a rust backend is used. Other `flate2` backends `zlib` and
//! `zlib-ng-compat` has been exported as the corresponding features (`zlib` and
//! `zlib-ng-compat`). See <https://docs.rs/flate2/latest/flate2/> for more details.
//!
#![warn(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro)]
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::bufread::DeflateDecoder;
use rayon::prelude::*;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::Read;
use std::ops::Range;
use std::path::Path;
use std::{collections::HashMap, io::Seek};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    HeaderNotParsed,
    IndexBcfReaderMissingGenomeInterval,
    CsBinIndexNotFound,
    FromUtf8Error(std::string::FromUtf8Error),
    Utf8Error(std::str::Utf8Error),
    ParseHeaderError(ParseHeaderError),
    NumericaValueEmptyInt,
    NumericaValueAsF32Error,
    Other(String),
}

#[derive(Debug)]
pub struct ParseGzipHeaderError {
    pub key: &'static str,
    pub expected: usize,
    pub found: usize,
}

#[derive(Debug)]
pub enum ParseHeaderError {
    HeaderCommentCharError,
    MissingDictionaryname,
    FormatError(&'static str),
}

trait AddContext {
    fn add_context(self, context: &'static str) -> Error;
}
impl AddContext for std::io::Error {
    fn add_context(self, context: &'static str) -> Error {
        Error::from(std::io::Error::new(self.kind(), context))
    }
}

// --- Display
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", *self)
    }
}
impl std::fmt::Display for ParseHeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", *self)
    }
}
impl std::fmt::Display for ParseGzipHeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", *self)
    }
}

// --- Error

impl std::error::Error for Error {}

impl std::error::Error for ParseGzipHeaderError {}

// --- From
impl From<ParseGzipHeaderError> for std::io::Error {
    fn from(value: ParseGzipHeaderError) -> Self {
        std::io::Error::new(std::io::ErrorKind::InvalidData, Box::new(value))
    }
}
impl From<std::io::Error> for Error {
    fn from(value: std::io::Error) -> Self {
        Error::Io(value)
    }
}

/// An iterator used to split a `str` by a separator with separators within pairs
/// of quotes ignored.
pub struct QuotedSplitter<'a> {
    data: &'a str,
    in_quotes: bool,
    sep: char,
    quote: char,
}

impl<'a> QuotedSplitter<'a> {
    /// Creates a new `QuotedSplitter` iterator.
    ///
    /// # Arguments
    ///
    /// * `buffer` - The input string to be split.
    /// * `separator` - The separator character used to split the string.
    /// * `quote` - The quote character used to ignore separators within quotes.
    ///
    /// # Examples
    ///
    /// ```
    /// use bcf_reader::QuotedSplitter;
    /// let input_string = "hello,\"world, this is fun\",test";
    /// let result: Vec<_> = QuotedSplitter::new(input_string, ',', '"').collect();
    /// assert_eq!(result, vec!["hello", "\"world, this is fun\"", "test"]);
    /// ```
    pub fn new(buffer: &'a str, separator: char, quote: char) -> Self {
        Self {
            data: buffer,
            in_quotes: false,
            sep: separator,
            quote,
        }
    }
}

impl<'a> Iterator for QuotedSplitter<'a> {
    type Item = &'a str;

    /// Advances the iterator and returns the next split substring.
    ///
    /// # Returns
    ///
    /// * `Some(&str)` - The next split substring.
    /// * `None` - If there are no more substrings to split.
    fn next(&mut self) -> Option<Self::Item> {
        for (idx, ch) in self.data.char_indices() {
            if ch == self.quote {
                self.in_quotes = !self.in_quotes;
            }
            if (!self.in_quotes) && ch == self.sep {
                let (out, remain) = self.data.split_at(idx);
                self.data = remain.strip_prefix(self.sep).unwrap_or(remain);
                return Some(out);
            }
        }
        if !self.data.is_empty() {
            let out = self.data;
            self.data = "";
            Some(out)
        } else {
            None
        }
    }
}

/// Represents a header of a BCF file.
///
/// The `Header` struct contains information about the dictionar of strings and
/// contigs, samples, and the index of the FORMAT field "GT".  It provides
/// methods to parse a header from a string, retrieve the index of a dictionary
/// string, retrieve the chromosome name by index, retrieve the index of the
/// FORMAT field "GT", and access the dictionary strings, contigs, and samples.
///
/// # Examples
///
/// ```
/// use bcf_reader::Header;
///
/// let header_text = concat!(
///     r#"##fileformat=VCFv4.3"#,
///     "\n",
///     r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
///     "\n",
///     r#"##FILTER=<ID=FAILED1,Description="failed due to something">"#,
///     "\n",
///     r#"##contig=<ID=chr1,length=123123123>"#,
///     "\n",
///     r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#,
///     "\n",
///     r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
///     "\n",
///     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2",
///     "\n",
/// );
///
/// let header = Header::from_string(&header_text).unwrap();
///
/// assert_eq!(header.get_idx_from_dictionary_str("INFO", "DP"), Some(2));
/// assert_eq!(header.get_chrname(0), "chr1");
/// assert_eq!(header.get_fmt_gt_id(), Some(3));
/// assert_eq!(header.dict_contigs().len(), 1);
/// assert_eq!(header.dict_strings().len(), 4);
/// assert_eq!(header.get_samples().len(), 2);
/// ```
#[derive(Debug)]
pub struct Header {
    dict_strings: HashMap<usize, HashMap<String, String>>,
    dict_contigs: HashMap<usize, HashMap<String, String>>,
    samples: Vec<String>,
    fmt_gt_idx: Option<usize>,
}
impl Header {
    /// parse header lines to structured data `Header`
    pub fn from_string(text: &str) -> std::result::Result<Self, ParseHeaderError> {
        use ParseHeaderError::*;
        let mut dict_strings = HashMap::<usize, HashMap<String, String>>::new();
        let mut dict_contigs = HashMap::<usize, HashMap<String, String>>::new();
        let mut samples = Vec::<String>::new();

        // implicit FILTER/PASS header line
        let mut m = HashMap::<String, String>::new();
        m.insert("Dictionary".into(), "FILTER".into());
        m.insert("ID".into(), "PASS".into());
        m.insert("Description".into(), r#""All filters passed""#.into());
        dict_strings.insert(0, m);
        //
        let mut dict_str_idx_counter = 1;
        let mut dict_contig_idx_counter = 0;
        for line in QuotedSplitter::new(text.trim_end_matches('\0').trim(), '\n', '"') {
            if line.starts_with("#CHROM") {
                line.split("\t")
                    .skip(9)
                    .for_each(|s| samples.push(s.into()));
                continue;
            }
            if line.trim().is_empty() {
                continue;
            }
            let mut it = QuotedSplitter::new(
                line.strip_prefix("##").ok_or(HeaderCommentCharError)?,
                '=',
                '"',
            );
            let dict_name = it.next().ok_or(MissingDictionaryname)?;
            let valid_dict = matches!(it.next(), Some(x) if x.starts_with("<"));
            if !valid_dict {
                continue;
            }
            let l = line.find('<').ok_or(FormatError("'<' not found"))?;
            let s = line.split_at(l + 1).1;
            let r = s.rfind('>').ok_or(FormatError("> not found"))?;
            let s = s.split_at(r).0;
            let mut m = HashMap::<String, String>::new();
            for kv_str in QuotedSplitter::new(s, ',', '"') {
                let kv_str = kv_str.trim();

                let mut it = QuotedSplitter::new(kv_str, '=', '"');
                let k = it.next().ok_or(FormatError("key not found"))?;
                let v = it
                    .next()
                    .ok_or(FormatError("value not found"))?
                    .trim_end_matches('"')
                    .trim_start_matches('"');
                m.insert(k.into(), v.into());
            }
            match dict_name {
                "contig" => {
                    if m.contains_key("IDX") {
                        let idx: usize = m["IDX"]
                            .parse()
                            .map_err(|_| FormatError("IDX value parsing error"))?;
                        if dict_contig_idx_counter != 0 {
                            // if one dict string has IDX all of them should have IDX in the dictionary
                            return Err(FormatError("not all CONTIG lines have key IDX"));
                        }
                        dict_contigs.insert(idx, m);
                    } else {
                        dict_contigs.insert(dict_contig_idx_counter, m);
                        dict_contig_idx_counter += 1;
                    }
                }
                _ => {
                    if ((dict_name != "FILTER") || (&m["ID"] != "PASS"))
                        && ["INFO", "FILTER", "FORMAT"].iter().any(|x| *x == dict_name)
                    {
                        m.insert("Dictionary".into(), dict_name.into());
                        // dbg!(&m, dict_str_idx_counter);
                        if m.contains_key("IDX") {
                            let idx: usize = m["IDX"]
                                .parse()
                                .map_err(|_| FormatError("IDX value parsing error"))?;
                            if dict_str_idx_counter != 1 {
                                // if one dict string has IDX all of them should have IDX in the dictionary
                                return Err(FormatError(
                                    "not all INFO/FILTER/FORMAT lines have key IDX",
                                ));
                            }
                            dict_strings.insert(idx, m);
                        } else {
                            dict_strings.insert(dict_str_idx_counter, m);
                            dict_str_idx_counter += 1;
                        }
                    }
                }
            };
        }

        // find fmt_key for FORMAT/GT for convenience
        let mut fmt_gt_idx = None;
        for (k, m) in dict_strings.iter() {
            if (&m["Dictionary"] == "FORMAT") && (&m["ID"] == "GT") {
                fmt_gt_idx = Some(*k);
            }
        }

        Ok(Self {
            dict_strings,
            dict_contigs,
            samples,
            fmt_gt_idx,
        })
    }

    /// Find the key (offset in header line) for a given INFO/xx or FILTER/xx or FORMAT/xx field.
    ///
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let key_found = header.get_idx_from_dictionary_str("FORMAT", "GT").unwrap();
    /// assert_eq!(key_found, header.get_fmt_gt_id().unwrap());
    /// ```
    pub fn get_idx_from_dictionary_str(&self, dictionary: &str, field: &str) -> Option<usize> {
        for (k, m) in self.dict_strings.iter() {
            if (m["Dictionary"] == dictionary) && (m["ID"] == field) {
                return Some(*k);
            }
        }
        None
    }

    /// Get chromosome name from the contig index
    pub fn get_chrname(&self, idx: usize) -> &str {
        &self.dict_contigs[&idx]["ID"]
    }

    /// Get key for FORMAT/GT field.
    pub fn get_fmt_gt_id(&self) -> Option<usize> {
        self.fmt_gt_idx
    }

    /// Get hashmap of hashmap of dictionary of contigs
    /// outer key: contig_idx
    /// inner key: is the key of the dictionary of contig, such as 'ID', 'Description'
    pub fn dict_contigs(&self) -> &HashMap<usize, HashMap<String, String>> {
        &self.dict_contigs
    }

    /// Get hashmap of hashmap of dictionary of strings
    /// outer key: item_idx, for FILTER/xx, FORMAT/xx, INFO/xx,
    /// inner key: is the key of the dictionary of string, such as 'ID', 'Description'
    pub fn dict_strings(&self) -> &HashMap<usize, HashMap<String, String>> {
        &self.dict_strings
    }

    /// Get samples names from sample idx
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// // read data generated by bcftools
    /// // bcftools query -l test.bcf | bgzip -c > test_samples.gz
    /// let mut samples_str = String::new();
    /// smart_reader("./testdata/test_samples.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut samples_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let samples_str2 = header.get_samples().join("\n");
    /// // compare bcftools results and bcf-reader results
    /// assert_eq!(samples_str.trim(), samples_str2.trim());
    /// ```
    pub fn get_samples(&self) -> &Vec<String> {
        &self.samples
    }
}

/// map bcf2 type to width in bytes
///
/// `typ`:
/// - 0: MISSING
/// - 1: u8 (1 byte)
/// - 2: u16 (2 bytes)
/// - 3: u32 (3 bytes)
/// - 5: f32 (4 bytes)
/// - 7: c-char (u8, 1 byte)
pub fn bcf2_typ_width(typ: u8) -> usize {
    match typ {
        0x0 => 0,
        0x1 => 1,
        0x2 => 2,
        0x3 => 4,
        0x5 => 4,
        0x7 => 1,
        _ => panic!(),
    }
}

#[derive(Debug, PartialEq)]
/// Represents a numeric value in the context of the bcf-reader.
pub enum NumericValue {
    /// Represents an unsigned 8-bit integer value.
    U8(u8),
    /// Represents an unsigned 16-bit integer value.
    U16(u16),
    /// Represents an unsigned 32-bit integer value.
    U32(u32),
    /// Represents a 32-bit floating-point value. (Note that a u32 is used to
    /// hold the bits for the f32 value)
    F32(u32),
}

impl Default for NumericValue {
    fn default() -> Self {
        NumericValue::U8(0)
    }
}

impl From<u8> for NumericValue {
    fn from(value: u8) -> Self {
        Self::U8(value)
    }
}
impl From<u16> for NumericValue {
    fn from(value: u16) -> Self {
        Self::U16(value)
    }
}
impl From<u32> for NumericValue {
    fn from(value: u32) -> Self {
        Self::U32(value)
    }
}

impl NumericValue {
    fn is_missing(&self) -> bool {
        match *self {
            NumericValue::U8(x) => x == 0x80,
            NumericValue::U16(x) => x == 0x8000,
            NumericValue::U32(x) => x == 0x80000000,
            NumericValue::F32(x) => x == 0x7F800001,
        }
    }

    fn is_end_of_vector(&self) -> bool {
        match *self {
            NumericValue::U8(x) => x == 0x81,
            NumericValue::U16(x) => x == 0x8001,
            NumericValue::U32(x) => x == 0x80000001,
            NumericValue::F32(x) => x == 0x7F800002,
        }
    }

    fn as_f32(&self) -> Result<Self> {
        match *self {
            NumericValue::U32(x) => Ok(NumericValue::F32(x)),
            _ => Err(Error::NumericaValueAsF32Error),
        }
    }

    /// Returns the integer value if the NumericValue is an unsigned integer and not missing.
    ///
    /// # Examples
    ///
    /// ```
    /// use bcf_reader::NumericValue;
    ///
    /// let value = NumericValue::U8(42);
    /// assert_eq!(value.int_val(), Some(42u32));
    ///
    /// let missing_value = NumericValue::U8(0x80u8);
    /// assert_eq!(missing_value.int_val(), None);
    /// ```
    pub fn int_val(&self) -> Option<u32> {
        if self.is_end_of_vector() || self.is_missing() {
            None
        } else {
            match *self {
                Self::U8(x) => Some(x as u32),
                Self::U16(x) => Some(x as u32),
                Self::U32(x) => Some(x),
                _ => None,
            }
        }
    }

    /// Returns the floating-point value if the NumericValue is a 32-bit float and not missing.
    ///
    /// # Examples
    ///
    /// ```
    /// use bcf_reader::NumericValue;
    ///
    /// let value = NumericValue::F32(3.14f32.to_bits());
    /// assert_eq!(value.float_val(), Some(3.14f32));
    ///
    /// let missing_value = NumericValue::F32(0x7F800001);
    /// let missing_value2 = NumericValue::F32(0x7F800001);
    /// assert_eq!(missing_value.float_val(), missing_value2.float_val());
    /// dbg!(&missing_value);
    /// assert_eq!(missing_value.float_val(), None);
    /// ```
    pub fn float_val(&self) -> Option<f32> {
        if self.is_end_of_vector() || self.is_missing() {
            None
        } else {
            match *self {
                Self::F32(x) => Some(f32::from_bits(x)),
                _ => panic!(),
            }
        }
    }

    /// Returns a tuple representing the GT value.
    ///
    /// The tuple contains the following elements:
    /// - `noploidy`: A boolean indicating if the ploidy is missing (for individuals with fewer ploids compared to individuals with the maximum ploidy).
    /// - `dot`: A boolean indicating if the genotype call is a dot.
    /// - `phased`: A boolean indicating if the allele is phased (the first allele of a call is always unphased).
    /// - `allele`: The allele value (index).
    ///
    /// # Examples
    ///
    /// ```
    /// use bcf_reader::NumericValue;
    ///
    /// let value = NumericValue::U8(3);
    /// assert_eq!(value.gt_val(), (false, false, true, 0));
    ///
    /// let value = NumericValue::U8(5);
    /// assert_eq!(value.gt_val(), (false, false, true, 1));
    ///
    /// let missing_value = NumericValue::U8(0);
    /// assert_eq!(missing_value.gt_val(), (false, true, false, u32::MAX));
    /// ```
    pub fn gt_val(&self) -> (bool, bool, bool, u32) {
        let mut noploidy = false;
        let mut dot = false;
        let mut phased = false;
        let mut allele = u32::MAX;

        match self.int_val() {
            None => {
                noploidy = true;
            }
            Some(int_val) => {
                phased = (int_val & 0x1) != 0;

                let int_val = int_val >> 1;
                if int_val == 0 {
                    dot = true;
                } else {
                    allele = int_val - 1;
                }
            }
        };

        (noploidy, dot, phased, allele)
    }
}

/// Read typed descriptor from the reader (of decompressed BCF buffer)
///
/// Return `typ` for type and `n` for count of elements of the type.
pub fn read_typed_descriptor_bytes<R>(reader: &mut R) -> std::io::Result<(u8, usize)>
where
    R: std::io::Read + ReadBytesExt,
{
    let tdb = reader.read_u8()?;
    let typ = tdb & 0xf;
    let mut n = (tdb >> 4) as usize;
    if n == 15 {
        n = read_single_typed_integer(reader)? as usize;
    }
    Ok((typ, n))
}

/// Read a single typed integer from the reader (of decompressed BCF buffer)
pub fn read_single_typed_integer<R>(reader: &mut R) -> std::io::Result<u32>
where
    R: std::io::Read + ReadBytesExt,
{
    let (typ, n) = read_typed_descriptor_bytes(reader)?;
    assert_eq!(n, 1);
    Ok(match typ {
        1 => reader.read_u8()? as u32,
        2 => reader.read_u16::<LittleEndian>()? as u32,
        3 => reader.read_u32::<LittleEndian>()?,
        _ => panic!(),
    })
}

/// Iterator for accessing arrays of numeric values (integers or floats)
/// directly from the buffer bytes without building Vec<_> or Vec<Vec<_>>
/// for each site.
#[derive(Default, Debug)]
pub struct NumericValueIter<'r> {
    reader: std::io::Cursor<&'r [u8]>,
    typ: u8,
    len: usize,
    cur: usize,
}

impl<'r> Iterator for NumericValueIter<'r> {
    type Item = std::io::Result<NumericValue>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.len {
            None
        } else {
            match self.typ {
                0 => None,
                1 => {
                    self.cur += 1;
                    Some(self.reader.read_u8().map(NumericValue::from))
                }
                2 => {
                    self.cur += 1;
                    Some(
                        self.reader
                            .read_u16::<LittleEndian>()
                            .map(NumericValue::from),
                    )
                }
                3 => {
                    self.cur += 1;
                    Some(
                        self.reader
                            .read_u32::<LittleEndian>()
                            .map(NumericValue::from),
                    )
                }
                5 => {
                    self.cur += 1;
                    let val_res = match self
                        .reader
                        .read_u32::<LittleEndian>()
                        .map(NumericValue::from)
                    {
                        Ok(nv) => match nv.as_f32() {
                            Ok(nv) => Ok(nv),
                            Err(_e) => Err(std::io::Error::new(
                                std::io::ErrorKind::Other,
                                Error::NumericaValueAsF32Error,
                            )),
                        },
                        Err(e) => Err(e),
                    };

                    Some(val_res)
                }
                _ => panic!(),
            }
        }
    }
}

/// Generate an iterator of numbers from a continuous bytes buffer
/// - typ: data type byte
/// - n: total number of elements to iterate
/// - buffer: the bytes buffer
pub fn iter_typed_integers(typ: u8, n: usize, buffer: &[u8]) -> NumericValueIter {
    NumericValueIter {
        reader: std::io::Cursor::new(buffer),
        typ,
        len: n,
        cur: 0,
    }
}

/// Read a typed string from the reader to a Rust String
pub fn read_typed_string<R>(reader: &mut R, buffer: &mut Vec<u8>) -> std::io::Result<usize>
where
    R: std::io::Read + ReadBytesExt,
{
    let (typ, n) = read_typed_descriptor_bytes(reader)?;
    assert_eq!(typ, 0x7);
    let s = buffer.len();
    buffer.resize(s + n, b'\0');
    reader.read_exact(&mut buffer.as_mut_slice()[s..s + n])?;
    Ok(n)
}

/// read the header lines to a String
/// use Header::from_string(text) to convert the string into structured data
pub fn read_header<R>(reader: &mut R) -> Result<String>
where
    R: std::io::Read + ReadBytesExt,
{
    // read magic
    let mut magic = [0u8; 3];
    reader.read_exact(&mut magic)?;
    assert_eq!(&magic, b"BCF");

    // read major verion and minor version
    let major = reader.read_u8()?;
    let minor = reader.read_u8()?;
    assert_eq!(major, 2);
    assert_eq!(minor, 2);

    // read text length
    let l_length = reader.read_u32::<LittleEndian>()?;
    let mut text = vec![0u8; l_length as usize];
    reader.read_exact(&mut text)?;

    String::from_utf8(text).map_err(Error::FromUtf8Error)
}

/// Represents a record (a line or a site) in BCF file
#[derive(Default, Debug)]
pub struct Record {
    buf_shared: Vec<u8>,
    buf_indiv: Vec<u8>,
    chrom: i32,
    pos: i32,
    rlen: i32,
    qual: NumericValue,
    n_info: u16,
    n_allele: u16,
    n_sample: u32,
    n_fmt: u8,
    id: Range<usize>,
    alleles: Vec<Range<usize>>,
    /// (typ, n, byte_range)
    filters: (u8, usize, Range<usize>),
    /// (info_key, typ, n, byte_range)
    info: Vec<(usize, u8, usize, Range<usize>)>,
    /// (fmt_key, typ, n, byte_range)
    gt: Vec<(usize, u8, usize, Range<usize>)>,
}
impl Record {
    /// read a record (copy bytes from the reader to the record's interval
    /// buffers), and separate fields
    pub fn read<R>(&mut self, reader: &mut R) -> Result<()>
    where
        R: std::io::Read + ReadBytesExt,
    {
        let l_shared = match reader.read_u32::<LittleEndian>() {
            Ok(x) => x,
            Err(_x) => Err(_x)?,
        };
        let l_indv = reader.read_u32::<LittleEndian>()?;
        // dbg!(l_shared, l_indv);
        self.buf_shared.resize(l_shared as usize, 0u8);
        self.buf_indiv.resize(l_indv as usize, 0u8);
        reader.read_exact(self.buf_shared.as_mut_slice())?;
        reader.read_exact(self.buf_indiv.as_mut_slice())?;
        self.parse_shared()?;
        self.parse_indv()?;
        // dbg!(self.pos);
        Ok(())
    }
    /// parse shared fields
    fn parse_shared(&mut self) -> std::io::Result<()> {
        let mut reader = std::io::Cursor::new(self.buf_shared.as_slice());
        self.chrom = reader.read_i32::<LittleEndian>()?;
        self.pos = reader.read_i32::<LittleEndian>()?;
        self.rlen = reader.read_i32::<LittleEndian>()?;
        let qual_u32 = reader.read_u32::<LittleEndian>()?;
        self.qual = NumericValue::from(qual_u32)
            .as_f32()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        self.n_info = reader.read_u16::<LittleEndian>()?;
        self.n_allele = reader.read_u16::<LittleEndian>()?;
        let combined = reader.read_u32::<LittleEndian>()?;
        self.n_sample = combined & 0xffffff;
        self.n_fmt = (combined >> 24) as u8;
        // id
        let (typ, n) = read_typed_descriptor_bytes(&mut reader)?;
        assert_eq!(typ, 0x7);
        let cur = reader.position() as usize;
        self.id = cur..cur + n;
        reader.seek(std::io::SeekFrom::Current(n as i64))?;
        // alleles
        self.alleles.clear();
        for _ in 0..self.n_allele {
            let (typ, n) = read_typed_descriptor_bytes(&mut reader)?;
            assert_eq!(typ, 0x7);
            let cur = reader.position() as usize;
            self.alleles.push(cur..cur + n);
            reader.seek(std::io::SeekFrom::Current(n as i64))?;
        }
        //filters
        let (typ, n) = read_typed_descriptor_bytes(&mut reader)?;
        let width: usize = bcf2_typ_width(typ);
        let s = reader.position() as usize;
        let e = s + width * n;
        reader.seek(std::io::SeekFrom::Current((e - s) as i64))?;
        self.filters = (typ, n, s..e);
        // infos
        self.info.clear();
        for _idx in 0..(self.n_info as usize) {
            let info_key = read_single_typed_integer(&mut reader)?;
            let (typ, n) = read_typed_descriptor_bytes(&mut reader)?;
            let width = bcf2_typ_width(typ);
            let s = reader.position() as usize;
            let e = s + width * n;
            reader.seek(std::io::SeekFrom::Current((e - s) as i64))?;
            self.info.push((info_key as usize, typ, n, s..e));
        }
        Ok(())
    }
    /// parse indiv fields, complicated field will need further processing
    fn parse_indv(&mut self) -> std::io::Result<()> {
        let mut reader = std::io::Cursor::new(self.buf_indiv.as_slice());
        self.gt.clear();
        for _idx in 0..(self.n_fmt as usize) {
            let fmt_key = read_single_typed_integer(&mut reader)?;
            let (typ, n) = read_typed_descriptor_bytes(&mut reader)?;
            let width = bcf2_typ_width(typ);
            let s = reader.position() as usize;
            let e = s + width * self.n_sample as usize * n;
            reader.seek(std::io::SeekFrom::Current((e - s) as i64))?;
            self.gt.push((fmt_key as usize, typ, n, s..e));
        }
        Ok(())
    }

    /// get chromosome offset
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// use std::io::Write;
    /// // read data generated by bcftools
    /// // bcftools query -f '%CHROM\n' test.bcf | bgzip -c > test_chrom.gz
    /// let mut chrom_str = String::new();
    /// smart_reader("testdata/test_chrom.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut chrom_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let mut record = Record::default();
    /// let mut chrom_str2 = Vec::<u8>::new();
    /// while let Ok(_) = record.read(&mut f) {
    ///     write!(
    ///         chrom_str2,
    ///         "{}\n",
    ///         header.get_chrname(record.chrom() as usize)
    ///     )
    ///     .unwrap();
    /// }
    /// let chrom_str2 = String::from_utf8(chrom_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// assert_eq!(chrom_str, chrom_str2);
    /// ```
    pub fn chrom(&self) -> i32 {
        self.chrom
    }

    /// Returns the reference length of the record.
    pub fn rlen(&self) -> i32 {
        self.rlen
    }

    /// Returns the quality score of the record, if available
    pub fn qual(&self) -> Option<f32> {
        self.qual.float_val()
    }

    pub fn n_allele(&self) -> u16 {
        self.n_allele
    }

    /// Returns an iterator over the genotype values in the record's FORMAT field.
    /// If no FORMAT/GT field available, the returned iterator will have items.
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// use std::io::Write;
    /// // read data generated by bcftools
    /// // bcftools query -f '[\t%GT]\n' test.bcf | bgzip -c > test_gt.gz
    /// let mut gt_str = String::new();
    /// smart_reader("testdata/test_gt.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut gt_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let mut record = Record::default();
    /// let mut gt_str2 = Vec::<u8>::new();
    /// while let Ok(_) = record.read(&mut f) {
    ///     for (i, bn) in record.fmt_gt(&header).enumerate() {
    ///         let (noploidy, dot, phased, allele) = bn.unwrap().gt_val();
    ///         assert_eq!(noploidy, false); // missing ploidy
    ///         let mut sep = '\t';
    ///         if i % 2 == 1 {
    ///             if phased {
    ///                 sep = '|';
    ///             } else {
    ///                 sep = '/';
    ///             }
    ///         }
    ///         if dot {
    ///             write!(gt_str2, "{sep}.").unwrap();
    ///         } else {
    ///             write!(gt_str2, "{sep}{allele}").unwrap();
    ///         }
    ///     }
    ///     write!(gt_str2, "\n").unwrap();
    /// }
    /// let gt_str2 = String::from_utf8(gt_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// for (a, b) in gt_str
    ///     .split(|c| (c == '\n') || (c == '\t'))
    ///     .zip(gt_str2.split(|c| (c == '\n') || (c == '\t')))
    /// {
    ///     assert_eq!(a, b);
    /// }
    /// ```
    pub fn fmt_gt(&self, header: &Header) -> NumericValueIter<'_> {
        let mut it = NumericValueIter::default();
        match header.get_fmt_gt_id() {
            None => it,
            Some(fmt_gt_id) => {
                // find the right field for gt
                self.gt.iter().for_each(|e| {
                    if e.0 == fmt_gt_id {
                        it = iter_typed_integers(
                            e.1,
                            e.2 * self.n_sample as usize,
                            &self.buf_indiv[e.3.start..e.3.end],
                        );
                    }
                });
                it
            }
        }
    }

    /// Returns an iterator over all values for a field in the record's FORMATs (indiv).
    ///
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// use std::io::Write;
    /// // read data generated by bcftools
    /// //  bcftools query -f '[\t%AD]\n' test.bcf | bgzip -c > test_ad.gz
    /// let mut ad_str = String::new();
    /// smart_reader("testdata/test_ad.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut ad_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let mut record = Record::default();
    /// let mut ad_str2 = Vec::<u8>::new();
    /// let ad_filed_key = header.get_idx_from_dictionary_str("FORMAT", "AD").unwrap();
    /// while let Ok(_) = record.read(&mut f) {
    ///     for (i, val) in record.fmt_field(ad_filed_key).enumerate() {
    ///         let val = val.unwrap();
    ///         if i % record.n_allele() as usize == 0 {
    ///             if ad_str2.last().map(|c| *c == b',') == Some(true) {
    ///                 ad_str2.pop(); // trim last allele separator
    ///             }
    ///             ad_str2.push(b'\t'); // sample separator
    ///         }
    ///         match val.int_val() {
    ///             None => {}
    ///             Some(ad) => {
    ///                 write!(ad_str2, "{ad},").unwrap(); // allele separator
    ///             }
    ///         }
    ///     }
    ///     // site separator
    ///     *ad_str2.last_mut().unwrap() = b'\n'; // sample separator
    /// }
    /// let ad_str2 = String::from_utf8(ad_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// for (a, b) in ad_str
    ///     .split(|c| (c == '\n') || (c == '\t'))
    ///     .zip(ad_str2.split(|c| (c == '\n') || (c == '\t')))
    /// {
    ///     assert_eq!(a, b);
    /// }
    /// ```
    pub fn fmt_field(&self, fmt_key: usize) -> NumericValueIter<'_> {
        // default iterator
        let mut it = NumericValueIter::default();

        // find the right field for gt
        self.gt.iter().for_each(|e| {
            if e.0 == fmt_key {
                it = iter_typed_integers(
                    e.1,
                    e.2 * self.n_sample as usize,
                    &self.buf_indiv[e.3.start..e.3.end],
                );
            }
        });
        it
    }

    /// get 0-based position (bp) value
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// // read data generated by bcftools
    /// // bcftools query -f '%POS\n' test.bcf | bgzip -c > test_pos.gz
    /// let mut pos_str = String::new();
    /// smart_reader("testdata/test_pos.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut pos_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let _s = read_header(&mut f).unwrap();
    /// let mut record = Record::default();
    /// let mut pos_str2 = Vec::<u8>::new();
    /// use std::io::Write;
    /// while let Ok(_) = record.read(&mut f) {
    ///     write!(pos_str2, "{}\n", record.pos() + 1).unwrap();
    /// }
    /// let pos_str2 = String::from_utf8(pos_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// assert_eq!(pos_str, pos_str2);
    /// ```
    pub fn pos(&self) -> i32 {
        self.pos
    }

    /// get variant ID as &str
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// // read data generated by bcftools
    /// // bcftools query -f '%ID\n' test.bcf | bgzip -c > test_id.gz
    /// let mut id_str = String::new();
    ///
    /// smart_reader("testdata/test_id.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut id_str)
    ///     .unwrap();
    ///
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let _s = read_header(&mut f).unwrap();
    /// let mut record = Record::default();
    /// let mut id_str2 = Vec::<u8>::new();
    /// use std::io::Write;
    /// while let Ok(_) = record.read(&mut f) {
    ///     let record_id = record.id().unwrap();
    ///     let id = if record_id.is_empty() { "." } else { record_id };
    ///     write!(id_str2, "{}\n", id).unwrap();
    /// }
    /// let id_str2 = String::from_utf8(id_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// assert_eq!(id_str, id_str2);
    /// ```
    pub fn id(&self) -> Result<&str> {
        std::str::from_utf8(&self.buf_shared[self.id.start..self.id.end]).map_err(Error::Utf8Error)
    }

    /// Returns the ranges of bytes in buf_shared for all alleles in the record.
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// // read data generated by bcftools
    /// // bcftools query -f '%REF,%ALT\n' test.bcf | bgzip -c > test_allele.gz
    /// let mut allele_str = String::new();
    /// smart_reader("testdata/test_allele.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut allele_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test.bcf").unwrap();
    /// let _s = read_header(&mut f).unwrap();
    /// let mut record = Record::default();
    /// let mut allele_str2 = Vec::<u8>::new();
    /// while let Ok(_) = record.read(&mut f) {
    ///     for rng in record.alleles().iter() {
    ///         let slice = &record.buf_shared()[rng.clone()];
    ///         allele_str2.extend(slice);
    ///         allele_str2.push(b',');
    ///     }
    ///     *allele_str2.last_mut().unwrap() = b'\n';
    /// }
    /// let allele_str2 = String::from_utf8(allele_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// assert_eq!(allele_str, allele_str2);
    /// ```
    pub fn alleles(&self) -> &[Range<usize>] {
        &self.alleles[..]
    }

    /// Return an iterator of numeric values for an INFO/xxx field.
    /// If the key is not found, the returned iterator will have a zero length.
    ///
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// use std::io::Write;
    /// // read data generated by bcftools
    /// // bcftools query -f '%INFO/AF\n' testdata/test2.bcf | bgzip -c > testdata/test2_info_af.gz
    /// let mut info_af_str = String::new();
    /// smart_reader("testdata/test2_info_af.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut info_af_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test2.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let mut record = Record::default();
    /// let mut info_af_str2 = Vec::<u8>::new();
    /// let info_af_key = header.get_idx_from_dictionary_str("INFO", "AF").unwrap();
    /// while let Ok(_) = record.read(&mut f) {
    ///     record.info_field_numeric(info_af_key).for_each(|nv| {
    ///         let af = nv.unwrap().float_val().unwrap();
    ///         write!(info_af_str2, "{af},").unwrap();
    ///     });
    ///     *info_af_str2.last_mut().unwrap() = b'\n'; // line separators
    /// }
    /// let filter_str2 = String::from_utf8(info_af_str2).unwrap();
    /// assert_eq!(info_af_str, filter_str2);
    /// ```
    pub fn info_field_numeric(&self, info_key: usize) -> NumericValueIter {
        // default
        let mut it = NumericValueIter {
            reader: std::io::Cursor::new(&[0u8; 0]),
            typ: 0,
            len: 0,
            cur: 0,
        };
        for (key, typ, n, rng) in self.info.iter() {
            if *key == info_key {
                it = NumericValueIter {
                    reader: std::io::Cursor::new(&self.buf_shared[rng.start..rng.end]),
                    typ: *typ,
                    len: *n,
                    cur: 0,
                };
                break;
            }
        }
        it
    }

    /// Return str value for an INFO/xxx field.
    /// If the key is not found or data type is not string, then return None.
    pub fn info_field_str(&self, info_key: usize) -> Result<Option<&str>> {
        let mut res = None;
        for (key, typ, _n, rng) in self.info.iter() {
            if *key == info_key {
                if *typ != 0x7 {
                    return Ok(None);
                }
                let s = std::str::from_utf8(&self.buf_shared[rng.start..rng.end])
                    .map_err(Error::Utf8Error)?;
                res = Some(s);
                break;
            }
        }
        Ok(res)
    }

    /// iterate an integer for each filter key.
    /// If the length of the iterator is 0, it means no filter label is set.
    /// Example:
    /// ```
    /// use bcf_reader::*;
    /// use std::io::Write;
    /// // read data generated by bcftools
    /// // bcftools query -f '%FILTER\n' testdata/test2.bcf | bgzip -c > testdata/test2_filters.gz
    /// let mut filter_str = String::new();
    /// smart_reader("testdata/test2_filters.gz")
    ///     .unwrap()
    ///     .read_to_string(&mut filter_str)
    ///     .unwrap();
    /// // read data via bcf-reader
    /// let mut f = smart_reader("testdata/test2.bcf").unwrap();
    /// let s = read_header(&mut f).unwrap();
    /// let header = Header::from_string(&s).unwrap();
    /// let mut record = Record::default();
    /// let mut filter_str2 = Vec::<u8>::new();
    /// let d = header.dict_strings();
    /// while let Ok(_) = record.read(&mut f) {
    ///     record.filters().for_each(|nv| {
    ///         let nv = nv.unwrap();
    ///         let filter_key = nv.int_val().unwrap() as usize;
    ///         let dict_string_map = &d[&filter_key];
    ///         let filter_name = &dict_string_map["ID"];
    ///         write!(filter_str2, "{filter_name};").unwrap();
    ///     });
    ///     *filter_str2.last_mut().unwrap() = b'\n'; // line separators
    /// }
    /// let filter_str2 = String::from_utf8(filter_str2).unwrap();
    /// // compare bcftools results and bcf-reader results
    /// assert_eq!(filter_str, filter_str2);
    /// ```
    pub fn filters(&self) -> NumericValueIter {
        let (typ, n, rng) = &self.filters;
        NumericValueIter {
            reader: std::io::Cursor::new(&self.buf_shared[rng.start..rng.end]),
            typ: *typ,
            len: *n,
            cur: 0,
        }
    }

    /// Returns the buffer containing indv (sample-level) information
    pub fn buf_indiv(&self) -> &[u8] {
        &self.buf_indiv[..]
    }

    /// Returns the buffer containing the shared (site-level) information
    pub fn buf_shared(&self) -> &[u8] {
        &self.buf_shared[..]
    }
}

/// Open a file from a path as a MultiGzDecoder or a BufReader depending on
/// whether the file has the magic number for gzip (0x1f and 0x8b)
pub fn smart_reader(p: impl AsRef<std::path::Path>) -> std::io::Result<Box<dyn std::io::Read>> {
    let mut f = std::fs::File::open(p.as_ref())?;
    if (f.read_u8()? == 0x1fu8) && (f.read_u8()? == 0x8bu8)
    //.expect("can not read second byte")
    {
        // gzip format
        f.rewind()?;
        Ok(Box::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        // not gzip format
        f.rewind()?;
        Ok(Box::new(std::io::BufReader::new(f)))
    }
}

/// This reader facilitates parallel decompression of BCF data compressed in
/// the BGZF format—a specialized version of the multi-member gzip file format.
/// It utilizes internal buffers to sequentially ingest compressed data from
/// various gzip blocks, leveraging the `rayon` crate to achieve concurrent
/// decompression. This design addresses the potential bottleneck in data
/// processing speed that occurs when decompression is not executed in parallel,
/// ensuring more efficient handling of compressed data streams.
/// Example:
/// ```
/// use bcf_reader::*;
/// use std::fs::File;
/// use std::io::BufReader;
/// use std::io::Write;
/// // read data generated by bcftools
/// // bcftools query -f '[\t%GT]\n' test.bcf | bgzip -c > test_gt.gz
/// let mut gt_str = String::new();
/// smart_reader("testdata/test_gt.gz")
///     .unwrap()
///     .read_to_string(&mut gt_str)
///     .unwrap();
/// // read data via bcf-reader
/// let max_gzip_block_in_buffer = 10;
/// let reader = File::open("testdata/test.bcf").map(BufReader::new).unwrap();
/// let mut f =
///     ParMultiGzipReader::from_reader(reader, max_gzip_block_in_buffer, None, None).unwrap();
/// let s = read_header(&mut f).unwrap();
/// let header = Header::from_string(&s).unwrap();
/// let mut record = Record::default();
/// let mut gt_str2 = Vec::<u8>::new();
/// while let Ok(_) = record.read(&mut f) {
///     for (i, bn) in record.fmt_gt(&header).enumerate() {
///         let bn = bn.unwrap();
///         let (noploidy, dot, phased, allele) = bn.gt_val();
///         assert_eq!(noploidy, false); // missing ploidy
///         let mut sep = '\t';
///         if i % 2 == 1 {
///             if phased {
///                 sep = '|';
///             } else {
///                 sep = '/';
///             }
///         }
///         if dot {
///             write!(gt_str2, "{sep}.").unwrap();
///         } else {
///             write!(gt_str2, "{sep}{allele}").unwrap();
///         }
///     }
///     write!(gt_str2, "\n").unwrap();
/// }
/// let gt_str2 = String::from_utf8(gt_str2).unwrap();
/// // compare bcftools results and bcf-reader results
/// for (a, b) in gt_str
///     .split(|c| (c == '\n') || (c == '\t'))
///     .zip(gt_str2.split(|c| (c == '\n') || (c == '\t')))
/// {
///     assert_eq!(a, b);
/// }
/// ```
///
/// See [`ParMultiGzipReader::from_reader`] for an example to jump to a target
/// genome interval.
pub struct ParMultiGzipReader<R>
where
    R: Read,
{
    inner: R,
    // vector of pairs of vector; For each pair one is for compressed data, the other for uncompressed
    buffer: Vec<BgzfBuffer>,
    ngzip: usize, // number gzip blocks read into buffer
    igzip: usize, // the current block to be consumed
    ibyte: usize, // the current byte to be consumed
    coffset: u64,
    inner_eof: bool,
}

#[derive(Default, Clone)]
struct BgzfBuffer {
    compressed: Vec<u8>,
    uncompressed: Vec<u8>,
    coffset: u64, // inclusive
    gzip_size: u16,
    uncompressed_data_size: u32,
}

impl<R> ParMultiGzipReader<R>
where
    R: Read,
{
    /// Constructs a new `ParMultiGzipReader` by specifying the `ngzip_max` parameter,
    /// which defines the maximum number of gzip blocks that the internal buffers can
    /// handle simultaneously. This parameter should ideally be set to the number of
    /// CPU cores available to optimize parallel decompression performance, thereby
    /// leveraging the hardware's concurrency capabilities.
    ///
    /// The `coffset` parameter indicates the offset to the first byte of a
    /// target gzip block of the input reader.  The input reader should point to
    /// the start byte of a gzip block as indicated by `coffset` before passing
    /// the reader to `ParMultiGzipReader::from_reader`; otherwise,
    /// [`Seek::seek`] should be used on the input reader to adjust the position
    /// accordingly. Note that `ParMultiGzipReader` does not call [`Seek::seek`]
    /// on the reader.
    ///
    /// The `uoffset` parameter specifies the number of bytes to skip within the
    /// first decompressed gzip data. Since skipping within uncompressed data
    /// requires decompression, this offset is applied within the
    /// `ParMultiGzipReader::from_reader` method.
    ///
    /// # Examples
    /// ```
    /// use bcf_reader::*;
    /// use std::{
    ///     fs::File,
    ///     io::{BufReader, Seek},
    /// };
    /// // index file
    /// let csi = Csi::from_path("testdata/test3.bcf.csi").unwrap();
    /// // reader
    /// let mut reader = File::open("testdata/test3.bcf")
    ///     .map(BufReader::new)
    ///     .unwrap();
    ///
    /// // calculate first offsets of the target postion
    /// let start = 1495403 - 1;
    /// let end = 1495746 - 1;
    /// let chrom_id = 0;
    /// let bin_id = csi.get_bin_id(start, start + 1 as i64);
    /// let bin_details = csi.get_bin_details(chrom_id, bin_id).unwrap();
    /// let (coffset, uoffset) = bin_details.chunks()[0].chunk_beg.get_coffset_uoffset();
    ///
    /// // seek to the target bgzip block
    /// reader.seek(std::io::SeekFrom::Start(coffset)).unwrap();
    /// // create the parallelizable reader by wraping around the existing reader
    /// // and specifing offsets
    /// let mut reader =
    ///     ParMultiGzipReader::from_reader(reader, 1, Some(coffset), Some(uoffset)).unwrap();
    ///
    /// let mut record = Record::default();
    /// let mut pos_found = vec![];
    /// while let Ok(_) = record.read(&mut reader) {
    ///     let pos = record.pos() as i64;
    ///     // the bin containing the start position of target interval may have records
    ///     // before the start position, so skip them.
    ///     if pos < start {
    ///         continue;
    ///     }
    ///     // read the record until out of the target interval
    ///     else if pos >= end {
    ///         break;
    ///     }
    ///     pos_found.push(pos);
    /// }
    /// assert_eq!(pos_found, vec![start]);
    /// ```
    ///
    /// # Parameters
    /// - `reader`: The input reader from which gzip blocks will be read.
    /// - `ngzip_max`: The maximum number of gzip blocks that can be processed in parallel.
    /// - `coffset`: An optional offset to the start of a gzip block in the input reader.
    /// - `uoffset`: An optional offset within the first block of decompressed data.
    ///
    /// # Returns
    /// Returns a new instance of `ParMultiGzipReader`.
    pub fn from_reader(
        reader: R,
        ngzip_max: usize,
        coffset: Option<u64>,
        uoffset: Option<u64>,
    ) -> Result<Self> {
        let mut this = Self {
            inner: reader,
            buffer: vec![BgzfBuffer::default(); ngzip_max],
            ngzip: 0,
            igzip: 0,
            ibyte: 0,
            coffset: coffset.unwrap_or(0),
            inner_eof: false,
        };
        this.clear_and_fill_buffers()?;
        this.decomp_all()?;
        this.ibyte = uoffset.unwrap_or(0) as usize;
        Ok(this)
    }
    pub fn get_coffset_uoffset(&self) -> (u64, u64) {
        let coffset = self.buffer[self.igzip].coffset;
        let uoffset = self.ibyte as u64;
        (coffset, uoffset)
    }
    fn read_single_gzip(&mut self) -> io::Result<()> {
        let this_buffer_offset = if self.ngzip == 0 {
            self.coffset
        } else {
            let prev_buffer = &self.buffer[self.ngzip - 1];
            prev_buffer.coffset + prev_buffer.gzip_size as u64
        };
        let this_buffer = &mut self.buffer[self.ngzip];

        // dbg!(this_buffer_offset);

        // let coffset_beg = self.inner.stream_position().unwrap();
        // dbg!(coffset_beg);
        let id1 = match self.inner.read_u8() {
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                self.inner_eof = true;
                return Ok(());
            }
            Err(e) => {
                return Err(e);
            }
            Ok(id1) => id1,
        };
        if id1 != 31 {
            Err(ParseGzipHeaderError {
                key: "id1",
                expected: 31,
                found: id1 as usize,
            })?;
        }

        let id2 = self.inner.read_u8()?;
        if id2 != 139 {
            Err(ParseGzipHeaderError {
                key: "id2",
                expected: 139,
                found: id2 as usize,
            })?;
        }
        // assert_eq!(id2, 139);
        let cm = self.inner.read_u8()?;
        if cm != 8 {
            Err(ParseGzipHeaderError {
                key: "cm",
                expected: 8,
                found: cm as usize,
            })?;
        }
        // assert_eq!(cm, 8);
        let flg = self.inner.read_u8()?;
        if flg != 4 {
            Err(ParseGzipHeaderError {
                key: "flg",
                expected: 4,
                found: flg as usize,
            })?;
        }
        // assert_eq!(flg, 4);
        let _mtime = self.inner.read_u32::<LittleEndian>()?;
        let _xfl = self.inner.read_u8()?;
        let _os = self.inner.read_u8()?;
        let xlen = self.inner.read_u16::<LittleEndian>()?;
        let si1 = self.inner.read_u8()?;
        assert_eq!(si1, 66);
        let si2 = self.inner.read_u8()?;
        assert_eq!(si2, 67);
        let slen = self.inner.read_u16::<LittleEndian>()?;
        assert_eq!(slen, 2);
        let bsize = self.inner.read_u16::<LittleEndian>()?;

        let buffer_compressed = &mut this_buffer.compressed;
        let cdata_sz = bsize - xlen - 19;

        buffer_compressed.resize(cdata_sz as usize, 0u8);
        self.inner.read_exact(buffer_compressed.as_mut_slice())?;

        let _crc32 = self.inner.read_u32::<LittleEndian>()?;
        let isize = self.inner.read_u32::<LittleEndian>()?;

        let buffer_uncompressed = &mut this_buffer.uncompressed;
        buffer_uncompressed.resize(isize as usize, 0u8);
        this_buffer.coffset = this_buffer_offset;
        this_buffer.gzip_size = bsize + 1;
        this_buffer.uncompressed_data_size = isize;

        // increment counter
        self.ngzip += 1;

        Ok(())
    }
    // clear buffer and refill (sequential)
    fn clear_and_fill_buffers(&mut self) -> std::io::Result<()> {
        let Self {
            inner: _,
            buffer,
            coffset,
            ngzip,
            igzip,
            ibyte,
            inner_eof: _,
        } = self;

        // update coffset for the buffer vector based on last used buffer
        if *ngzip > 0 {
            let last_buffer = &buffer[*ngzip - 1];
            *coffset = last_buffer.coffset + last_buffer.gzip_size as u64;
        }

        buffer.iter_mut().for_each(|bgzf_buffer| {
            bgzf_buffer.compressed.clear();
            bgzf_buffer.uncompressed.clear();
            bgzf_buffer.coffset = 0;
            bgzf_buffer.gzip_size = 0;
            bgzf_buffer.uncompressed_data_size = 0;
        });
        *ngzip = 0;
        *igzip = 0;
        *ibyte = 0;
        for _i in 0..self.buffer.len() {
            self.read_single_gzip()?;
            if self.inner_eof {
                break;
            }
        }
        Ok(())
    }

    /// decompress all read gzip file in memory (parallel)
    fn decomp_all(&mut self) -> std::io::Result<()> {
        self.buffer
            .par_iter_mut()
            .try_for_each(|buffer| -> std::io::Result<()> {
                let compressed = buffer.compressed.as_slice();
                let uncompressed = &mut buffer.uncompressed.as_mut_slice();
                let mut deflater = DeflateDecoder::new(compressed);
                deflater.read_exact(uncompressed)?;
                //.map_err(|e| e.add_context("deflater.read_exact"))?;
                Ok(())
            })?;
        Ok(())
    }
}

impl<R> Read for ParMultiGzipReader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // no more data in the buffer
        if self.ngzip == self.igzip {
            // no more data to read from file
            if self.inner_eof {
                return Ok(0);
            }
            self.clear_and_fill_buffers()?;
            self.decomp_all()?;
        }
        // read from the buffer
        //  check current
        let uncompressed = self.buffer[self.igzip].uncompressed.as_slice();
        let mut n = uncompressed.len() - self.ibyte;
        if n > buf.len() {
            n = buf.len();
        }
        // copy data
        buf[0..n]
            .iter_mut()
            .zip(uncompressed[self.ibyte..].iter())
            .for_each(|(d, s)| *d = *s);
        self.ibyte += n;
        if self.ibyte == uncompressed.len() {
            self.igzip += 1;
            self.ibyte = 0;
            if self.ngzip == self.igzip {
                // no more data to read from file
                if self.inner_eof {
                    return Ok(0);
                }
                self.clear_and_fill_buffers()?;
                self.decomp_all()?;
            }
        }
        Ok(n)
    }
}

/// Virutal File offset used to jump to specific indexed bin within BCF-format
/// genotype data separated into BGZF blocks
#[derive(Default)]
pub struct VirtualFileOffsets(u64);

impl VirtualFileOffsets {
    /// Get the `coffset` and `uoffset` tuple from the virutalfileoffset
    pub fn get_coffset_uoffset(&self) -> (u64, u64) {
        (self.0 >> 16, self.0 & 0xffff)
    }
}

impl From<u64> for VirtualFileOffsets {
    /// Convert u64 into `VirtualFileOffsets`
    fn from(value: u64) -> Self {
        VirtualFileOffsets(value)
    }
}
impl Debug for VirtualFileOffsets {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (coffset, uoffset) = self.get_coffset_uoffset();
        write!(f, "coffset: {}, uoffset: {}", coffset, uoffset)
    }
}

#[derive(Default, Debug)]
struct CsiIndex {
    n_bin: i32,
    bins: Vec<CsiBin>,
}

/// A chunk within a bin in the  CSI data structure
#[derive(Default, Debug)]
pub struct CsiChunk {
    pub chunk_beg: VirtualFileOffsets,
    pub chunk_end: VirtualFileOffsets,
}

/// A bin in the CSI data structure
#[derive(Default, Debug)]
pub struct CsiBin {
    bin: u32,
    _loffset: VirtualFileOffsets,
    n_chunk: i32,
    chunks: Vec<CsiChunk>,
}

impl CsiBin {
    /// return a slice of chunks with a bin in the CSI data structure
    pub fn chunks(&self) -> &[CsiChunk] {
        &self.chunks[..]
    }
}

/// A struct representing CSI index file content
#[derive(Default, Debug)]
pub struct Csi {
    magic: [u8; 4],
    min_shift: i32,
    depth: i32,
    l_aux: i32,
    aux: Vec<u8>,
    n_ref: i32,
    indices: Vec<CsiIndex>,
    n_no_coor: Option<u64>,
}

impl Csi {
    /// Create Csi from a path to a `*.csi` file
    pub fn from_path(p: impl AsRef<Path>) -> Result<Self> {
        let mut csi = Csi::default();
        let mut file = smart_reader(p.as_ref())?;
        // magic
        file.read_exact(csi.magic.as_mut())
            .map_err(|e| e.add_context("error in reading csi magic bytes"))?;
        if csi.magic != [b'C', b'S', b'I', 1] {
            return Err(Error::Io(std::io::Error::new(
                io::ErrorKind::Other,
                "csi magic is not 'CSI1'".to_owned(),
            )));
        }
        // min_shift
        csi.min_shift = file
            .read_i32::<LittleEndian>()
            .map_err(|e| e.add_context("error in reading csi min_shift field"))?;
        // dbg!(csi.min_shift);
        // depth
        csi.depth = file
            .read_i32::<LittleEndian>()
            .map_err(|e| e.add_context("error in reading csi depth field"))?;
        // dbg!(csi.depth);
        // l_aux
        csi.l_aux = file
            .read_i32::<LittleEndian>()
            .map_err(|e| e.add_context("error in reading csi l_aux field"))?;
        // dbg!(csi.l_aux);
        // aux
        csi.aux.resize(csi.l_aux as usize, 0u8);
        file.read_exact(csi.aux.as_mut())
            .map_err(|e| e.add_context("error in reading csi aux field"))?;
        // n_ref
        csi.n_ref = file
            .read_i32::<LittleEndian>()
            .map_err(|e| e.add_context("error in reading csi n_ref field"))?;

        // iterate over chromosomes
        for _ in 0..csi.n_ref {
            let mut idx = CsiIndex {
                n_bin: {
                    file.read_i32::<LittleEndian>()
                        .map_err(|e| e.add_context("error in reading csi index n_bin field"))?
                },
                ..Default::default()
            };
            for _ in 0..idx.n_bin {
                let mut bin = CsiBin {
                    // bin
                    bin: file
                        .read_u32::<LittleEndian>()
                        .map_err(|e| e.add_context("error in reading csi bin bin field"))?,
                    _loffset: file
                        .read_u64::<LittleEndian>()
                        .map_err(|e| {
                            Error::from(std::io::Error::new(
                                e.kind(),
                                "error in reading csi bin loffset field",
                            ))
                        })?
                        .into(),
                    // n_chunk
                    n_chunk: file.read_i32::<LittleEndian>().map_err(|e| {
                        Error::from(std::io::Error::new(
                            e.kind(),
                            "error in reading csi bin n_chunk field",
                        ))
                    })?,
                    ..Default::default()
                };

                for _ in 0..bin.n_chunk {
                    let chunk = CsiChunk {
                        chunk_beg: file
                            .read_u64::<LittleEndian>()
                            .map_err(|e| {
                                Error::from(std::io::Error::new(
                                    e.kind(),
                                    "error in reading csi chunk chunk_beg",
                                ))
                            })?
                            .into(),
                        chunk_end: file
                            .read_u64::<LittleEndian>()
                            .map_err(|e| {
                                Error::from(std::io::Error::new(
                                    e.kind(),
                                    "error in reading csi chunk chunk_end",
                                ))
                            })?
                            .into(),
                    };
                    bin.chunks.push(chunk);
                }
                idx.bins.push(bin);
            }
            idx.bins.sort_by_key(|x| x.bin);
            csi.indices.push(idx);
        }
        // n_no_coor
        csi.n_no_coor = file.read_u64::<LittleEndian>().ok();

        Ok(csi)
    }

    /// Convert positional coordinate range to a bin number
    ///
    /// `beg`, `end`` coordinates are 0-based. It is exclusive for end.
    /// For some reason, not all length bin are searchable.
    /// It is seems that setting `end` to `beg` + 1 works well.
    // no sure how it works, but it is from CSIv1.pdf c code
    pub fn get_bin_id(&self, beg: i64, end: i64) -> u32 {
        let mut l = self.depth;
        let end = end - 1;
        let mut s = self.min_shift;
        let mut t = ((1 << (self.depth * 3)) - 1) / 7;

        while l > 0 {
            // dbg!(s);
            if (beg >> s) == (end >> s) {
                return (t + (beg >> s)) as u32;
            }
            s += 3;
            t += 1 << (l * 3);
            l += 1;
        }

        0
    }

    /// Get CsiBin based the chromosome id and bin number.
    ///
    /// The return CsiBin can provide details of the included chunks.
    pub fn get_bin_details(&self, seqid: usize, bin_id: u32) -> Result<&CsiBin> {
        // assert!(bin_id <= self.get_bin_limit());
        let bins = &self.indices[seqid].bins;
        // dbg!(self.indices.len());
        // dbg!(bins.len());
        let i = bins
            .as_slice()
            .binary_search_by(|x| x.bin.cmp(&bin_id))
            .map_err(|_| Error::CsBinIndexNotFound)?;
        // .expect("bin index not found\n");
        Ok(&bins[i])
    }

    /// Get the max possible bin number in theory. Note, the maximum bin may not
    /// be present in the Csi index file.
    pub fn get_bin_limit(&self) -> u32 {
        (1 << (((self.depth + 1) * 3) - 1)) / 7
    }
}

/// BcfReader suitable for read through the BCF file.
///
/// This assumes that the source reader is in BCF format and stored as BGZF blocks.
///
/// # Example
/// ```
/// use bcf_reader::*;
/// use std::{
///     fs::File,
///     io::{BufReader, Seek},
/// };
/// // underlying reader can be stdin or a file
/// let reader = std::fs::File::open("testdata/test3.bcf")
///     .map(BufReader::new)
///     .unwrap();
/// // 1. for sequential decompression (works for data from stdin or a file)
/// let reader = flate2::bufread::MultiGzDecoder::new(reader);
/// // 2. for parallel decompression (works for data from stdin or a file)
/// // however, when data is from stdin, jump is not supported for now.
/// // let reader = ParMultiGzipReader::from_reader(reader, 3, None, None);
///
/// // create bcf reader
/// let mut reader = BcfReader::from_reader(reader);
/// // read though header
/// let _header = reader.read_header();
/// // create a reusable record
/// let mut record = Record::default();
///
/// let mut pos_found = vec![];
/// // iterate records from the targeted genome interval
/// while let Ok(_) = reader.read_record(&mut record) {
///     pos_found.push(record.pos() + 1);
/// }
///
/// assert_eq!(pos_found[0], 72);
/// assert_eq!(*pos_found.last().unwrap(), 1498841);
/// ```
pub struct BcfReader<R>
where
    R: Read,
{
    inner: R,
    header_parsed: bool,
}

impl<R> BcfReader<R>
where
    R: Read,
{
    /// `max_gzip`, the number of gzip blocks to read before batch
    /// decompression.  See `ParMultiGzipReader::from_reader` (by default or
    /// None, 3 gzip buffers will be used.)
    pub fn from_reader(reader: R) -> Self {
        Self {
            inner: reader,
            header_parsed: false,
        }
    }

    /// Read the header
    pub fn read_header(&mut self) -> Result<Header> {
        let header =
            Header::from_string(&read_header(&mut self.inner)?).map_err(Error::ParseHeaderError)?;
        self.header_parsed = true;
        Ok(header)
    }

    /// Read one record. This should be called after the header is read and parsed.
    /// Otherwise, it will panic.
    pub fn read_record(&mut self, record: &mut Record) -> Result<()> {
        assert!(
            self.header_parsed,
            "header should be parsed before reading records"
        );
        record.read(&mut self.inner)
    }
}

/// A genome interval defined by chromosome id, start, and end positions
pub struct GenomeInterval {
    pub chrom_id: usize,
    pub start: i64,
    pub end: Option<i64>,
}

/// IndexedBcfReader allows random access to a specific genome interval of the
/// BCF file using a CSI index file. It is an wrapper around
/// [`ParMultiGzipReader<BufReader<File>>`] to allow parallelizable bgzip
/// decompression.
///
/// The source should be BCF file consisting of small BGZF blocks.
///
/// # Example
///
/// ```
/// use bcf_reader::*;
/// // create indexed bcf reader
/// let mut reader =
///     IndexedBcfReader::from_path("testdata/test3.bcf", "testdata/test3.bcf.csi", None).unwrap();
/// // read though header
/// let _header = reader.read_header();
/// // define targeted genome interval
/// let interval = GenomeInterval {
///     chrom_id: 0,
///     start: 1489230 - 1,
///     end: Some(1498509 - 1),
/// };
/// // set interval
/// reader.set_interval(interval);
/// // create a reusable record
/// let mut record = Record::default();
///
/// let mut pos_found = vec![];
/// // iterate records from the targeted genome interval
/// while let Ok(_) = reader.read_record(&mut record) {
///     pos_found.push(record.pos() + 1);
/// }
///
/// assert_eq!(
///     pos_found,
///     vec![
///         1489230, 1489979, 1490069, 1490311, 1492233, 1492337, 1493505, 1494178, 1495207,
///         1495403, 1495746, 1496047, 1497964, 1498188,
///     ]
/// )
/// ```
pub struct IndexedBcfReader {
    inner: ParMultiGzipReader<BufReader<File>>,
    csi: Csi,
    header_parsed: bool,
    genome_interval: Option<GenomeInterval>,
}

impl IndexedBcfReader {
    /// Create an IndexedBcfReader from paths to a bcf file and a corresponding
    /// csi index file.
    ///
    ///  - `max_gzip`, the number of gzip blocks to read before each batch
    ///    parallelized decompression. See [`ParMultiGzipReader::from_reader`]
    ///    (by default (None) use 3); this construct will automaticall read and
    ///    parse the header
    pub fn from_path(
        path_bcf: impl AsRef<Path>,
        path_csi: impl AsRef<Path>,
        max_gzip: Option<usize>,
    ) -> Result<Self> {
        let reader = File::open(path_bcf.as_ref()).map(BufReader::new)?;
        let csi = Csi::from_path(path_csi.as_ref())?;
        let reader = ParMultiGzipReader::from_reader(reader, max_gzip.unwrap_or(3), None, None)?;
        Ok(Self {
            inner: reader,
            csi,
            header_parsed: false,
            genome_interval: None,
        })
    }
    /// Read the header bytes, parse them and return a `Header`
    pub fn read_header(&mut self) -> Result<Header> {
        let header =
            Header::from_string(&read_header(&mut self.inner)?).map_err(Error::ParseHeaderError)?;
        self.header_parsed = true;
        Ok(header)
    }

    /// Jump the file pointer to the begining to the targeted genome interval
    ///
    /// If no site within the genome interval, read_record will return Err(_)
    pub fn set_interval(&mut self, genome_interval: GenomeInterval) -> Result<()> {
        // find the target based on csi
        let bin_id = self
            .csi
            .get_bin_id(genome_interval.start, genome_interval.start + 1);

        let (coffset, uoffset) = self
            .csi
            .get_bin_details(genome_interval.chrom_id, bin_id)?
            .chunks[0]
            .chunk_beg
            .get_coffset_uoffset();

        let par_reader = &mut self.inner;
        par_reader.inner.seek(io::SeekFrom::Start(coffset))?;

        // clear buffer, especially things related to coffset
        par_reader.buffer.iter_mut().for_each(|bgzf_buffer| {
            bgzf_buffer.compressed.clear();
            bgzf_buffer.uncompressed.clear();
            bgzf_buffer.coffset = 0;
            bgzf_buffer.gzip_size = 0;
            bgzf_buffer.uncompressed_data_size = 0;
        });
        par_reader.ngzip = 0;
        par_reader.igzip = 0;
        par_reader.ibyte = 0;

        // fill buffer
        par_reader.clear_and_fill_buffers()?;
        par_reader.decomp_all()?;

        // jump for uoffset
        par_reader.ibyte = uoffset as usize;

        self.genome_interval = Some(genome_interval);
        Ok(())
    }

    /// Read one record. Should be called after header is parsed.
    ///
    /// If `set_interval` has been called, only records within the given interval
    /// will be read.
    pub fn read_record(&mut self, record: &mut Record) -> Result<()> {
        if !self.header_parsed {
            return Err(Error::HeaderNotParsed);
        }

        let interval = self
            .genome_interval
            .as_ref()
            .ok_or(Error::IndexBcfReaderMissingGenomeInterval)?;
        let start = interval.start;
        let end = interval.end;
        loop {
            match record.read(&mut self.inner) {
                Ok(_) => {
                    if let Some(end) = end {
                        if record.pos as i64 >= end {
                            let e =
                                std::io::Error::new(std::io::ErrorKind::NotFound, "out of range");
                            Err(e)?;
                        }
                    }
                    if record.pos as i64 >= start {
                        return Ok(());
                    }
                }
                Err(e) => return Err(e),
            }
        }
    }
}
