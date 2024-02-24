use byteorder::{LittleEndian, ReadBytesExt};
use std::ops::Range;
use std::{collections::HashMap, io::Seek};

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
                self.data = remain.strip_prefix(self.sep).unwrap();
                return Some(out);
            }
        }
        if self.data.len() > 0 {
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
/// let header_text =concat!(
///     r#"##fileformat=VCFv4.3"#, "\n",
///     r#"##FILTER=<ID=PASS,Description="All filters passed">"#,"\n",
///     r#"##FILTER=<ID=FAILED1,Description="failed due to something">"#,"\n",
///     r#"##contig=<ID=chr1,length=123123123>"#,"\n",
///     r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#,"\n",
///     r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,"\n",
///     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2","\n",
/// );
///
/// let header = Header::from_string(&header_text);
///
/// assert_eq!(header.get_idx_from_dictionary_str("INFO", "DP"), Some(2));
/// assert_eq!(header.get_chrname(0), "chr1");
/// assert_eq!(header.get_fmt_gt_id(), Some(3));
/// assert_eq!(header.get_contigs().len(), 1);
/// assert_eq!(header.get_dict_strings().len(), 4);
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
    pub fn from_string(text: &str) -> Self {
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
            if line.trim().len() == 0 {
                continue;
            }
            let mut it = QuotedSplitter::new(line.strip_prefix("##").unwrap(), '=', '"');
            let dict_name = it.next().unwrap();
            let valid_dict = match it.next() {
                Some(x) if x.starts_with("<") => true,
                _ => false,
            };
            if !valid_dict {
                continue;
            }
            let l = line.find('<').unwrap();
            let s = line.split_at(l + 1).1;
            let r = s.rfind('>').unwrap();
            let s = s.split_at(r).0;
            let mut m = HashMap::<String, String>::new();
            for kv_str in QuotedSplitter::new(s, ',', '"') {
                let kv_str = kv_str.trim();

                let mut it = QuotedSplitter::new(kv_str, '=', '"');
                let k = it.next().unwrap();
                let v = it
                    .next()
                    .unwrap()
                    .trim_end_matches('"')
                    .trim_start_matches('"');
                m.insert(k.into(), v.into());
            }
            match dict_name {
                "contig" => {
                    if m.contains_key("IDX") {
                        assert_eq!(dict_contig_idx_counter, 0, "if one dict string has IDX all of them should have IDX in the dictionary");
                        let idx: usize = m["IDX"].parse().unwrap();
                        dict_contigs.insert(idx, m);
                    } else {
                        dict_contigs.insert(dict_contig_idx_counter, m);
                        dict_contig_idx_counter += 1;
                    }
                }
                _ => {
                    if (dict_name == "FILTER") && (&m["ID"] == "PASS") {
                        // skip FILTER/PASS already added
                    } else {
                        if ["INFO", "FILTER", "FORMAT"].iter().any(|x| *x == dict_name) {
                            m.insert("Dictionary".into(), dict_name.into());
                            if m.contains_key("IDX") {
                                assert_eq!(dict_str_idx_counter, 1, "if one dict string has IDX all of them should have IDX in the dictionary");
                                let idx: usize = m["IDX"].parse().unwrap();
                                dict_strings.insert(idx, m);
                            } else {
                                dict_strings.insert(dict_str_idx_counter, m);
                                dict_str_idx_counter += 1;
                            }
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

        Self {
            dict_strings,
            dict_contigs,
            samples,
            fmt_gt_idx,
        }
    }

    /// Find the key (offset in header line) for a given INFO/xx or FILTER/xx or FORMAT/xx field.
    pub fn get_idx_from_dictionary_str(&self, dictionary: &str, field: &str) -> Option<usize> {
        for (k, m) in self.dict_strings.iter() {
            if (&m["Dictionary"] == dictionary) && (&m["ID"] == field) {
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
    pub fn get_contigs(&self) -> &HashMap<usize, HashMap<String, String>> {
        &self.dict_contigs
    }

    /// Get hashmap of hashmap of dictionary of strings
    /// outer key: item_idx, for FILTER/xx, FORMAT/xx, INFO/xx,
    /// inner key: is the key of the dictionary of string, such as 'ID', 'Description'
    pub fn get_dict_strings(&self) -> &HashMap<usize, HashMap<String, String>> {
        &self.dict_strings
    }

    /// Get samples names from sample idx
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
    /// Represents a 32-bit floating-point value. (Note use a u32 to
    /// hold the bits for the f32 value
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

    fn as_f32(&self) -> Self {
        match *self {
            NumericValue::U32(x) => NumericValue::F32(x),
            _ => panic!(),
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
                Self::U32(x) => Some(x as u32),
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
    /// assert_eq!(missing_value.float_val(), missing_value2.float_val()) ;
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

pub fn read_typed_descriptor_bytes<R>(reader: &mut R) -> (u8, usize)
where
    R: std::io::Read + ReadBytesExt,
{
    let tdb = reader.read_u8().unwrap();
    let typ = tdb & 0xf;
    let mut n = (tdb >> 4) as usize;
    if n == 15 {
        n = read_single_typed_integer(reader) as usize;
    }
    (typ, n)
}

pub fn read_single_typed_integer<R>(reader: &mut R) -> u32
where
    R: std::io::Read + ReadBytesExt,
{
    let (typ, n) = read_typed_descriptor_bytes(reader);
    assert_eq!(n, 1);
    match typ {
        1 => reader.read_u8().unwrap() as u32,
        2 => reader.read_u16::<LittleEndian>().unwrap() as u32,
        3 => reader.read_u32::<LittleEndian>().unwrap(),
        _ => panic!(),
    }
}

#[derive(Default, Debug)]
pub struct NumericValueIter<'r> {
    reader: std::io::Cursor<&'r [u8]>,
    typ: u8,
    len: usize,
    cur: usize,
}

impl<'r> Iterator for NumericValueIter<'r> {
    type Item = NumericValue;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.len {
            None
        } else {
            match self.typ {
                0 => None,
                1 => {
                    self.cur += 1;
                    Some(self.reader.read_u8().unwrap().into())
                }
                2 => {
                    self.cur += 1;
                    Some(self.reader.read_u16::<LittleEndian>().unwrap().into())
                }
                3 => {
                    self.cur += 1;
                    Some(self.reader.read_u32::<LittleEndian>().unwrap().into())
                }
                5 => {
                    self.cur += 1;
                    let val = self.reader.read_u32::<LittleEndian>().unwrap();
                    Some(NumericValue::from(val).as_f32())
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
pub fn read_typed_string<R>(reader: &mut R, buffer: &mut Vec<u8>) -> usize
where
    R: std::io::Read + ReadBytesExt,
{
    let (typ, n) = read_typed_descriptor_bytes(reader);
    assert_eq!(typ, 0x7);
    let s = buffer.len();
    buffer.resize(s + n, b'\0');
    reader.read(&mut buffer.as_mut_slice()[s..s + n]).unwrap();
    n
}

/// read the header lines to a String
/// use Header::from_string(text) to convert the string into structured data
pub fn read_header<R>(reader: &mut R) -> String
where
    R: std::io::Read + ReadBytesExt,
{
    // read magic
    let mut magic = [0u8; 3];
    reader.read(&mut magic).unwrap();
    assert_eq!(&magic, b"BCF");

    // read major verion and minor version
    let major = reader.read_u8().unwrap();
    let minor = reader.read_u8().unwrap();
    assert_eq!(major, 2);
    assert_eq!(minor, 2);

    // read text length
    let l_length = reader.read_u32::<LittleEndian>().unwrap();
    let mut text = vec![0u8; l_length as usize];
    reader.read_exact(&mut text).unwrap();

    String::from_utf8(text).unwrap()
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
    pub fn read<R>(&mut self, reader: &mut R) -> Result<(), Box<dyn std::error::Error>>
    where
        R: std::io::Read + ReadBytesExt,
    {
        let l_shared;
        let l_indv;
        l_shared = match reader.read_u32::<LittleEndian>() {
            Ok(x) => x,
            Err(_x) => Err(_x)?,
        };
        l_indv = reader.read_u32::<LittleEndian>()?;
        self.buf_shared.resize(l_shared as usize, 0u8);
        self.buf_indiv.resize(l_indv as usize, 0u8);
        reader.read_exact(self.buf_shared.as_mut_slice()).unwrap();
        reader.read_exact(self.buf_indiv.as_mut_slice()).unwrap();
        self.parse_shared();
        self.parse_indv();
        Ok(())
    }
    /// parse shared fields
    fn parse_shared(&mut self) {
        let mut reader = std::io::Cursor::new(self.buf_shared.as_slice());
        self.chrom = reader.read_i32::<LittleEndian>().unwrap();
        self.pos = reader.read_i32::<LittleEndian>().unwrap();
        self.rlen = reader.read_i32::<LittleEndian>().unwrap();
        let qual_u32 = reader.read_u32::<LittleEndian>().unwrap();
        self.qual = NumericValue::from(qual_u32).as_f32();
        self.n_info = reader.read_u16::<LittleEndian>().unwrap();
        self.n_allele = reader.read_u16::<LittleEndian>().unwrap();
        let combined = reader.read_u32::<LittleEndian>().unwrap();
        self.n_sample = combined & 0xffffff;
        self.n_fmt = (combined >> 24) as u8;
        // id
        let (typ, n) = read_typed_descriptor_bytes(&mut reader);
        assert_eq!(typ, 0x7);
        let cur = reader.position() as usize;
        self.id = cur..cur + n as usize;
        reader.seek(std::io::SeekFrom::Current(n as i64)).unwrap();
        // alleles
        self.alleles.clear();
        for _ in 0..self.n_allele {
            let (typ, n) = read_typed_descriptor_bytes(&mut reader);
            assert_eq!(typ, 0x7);
            let cur = reader.position() as usize;
            self.alleles.push(cur..cur + n as usize);
            reader.seek(std::io::SeekFrom::Current(n as i64)).unwrap();
        }
        //filters
        let (typ, n) = read_typed_descriptor_bytes(&mut reader);
        let width: usize = bcf2_typ_width(typ);
        let s = reader.position() as usize;
        let e = s + width * n as usize;
        reader
            .seek(std::io::SeekFrom::Current((e - s) as i64))
            .unwrap();
        self.filters = (typ, n as usize, s..e);
        // infos
        self.info.clear();
        for _idx in 0..(self.n_info as usize) {
            let info_key = read_single_typed_integer(&mut reader);
            let (typ, n) = read_typed_descriptor_bytes(&mut reader);
            let width = bcf2_typ_width(typ);
            let s = reader.position() as usize;
            let e = s + width * n as usize;
            reader
                .seek(std::io::SeekFrom::Current((e - s) as i64))
                .unwrap();
            self.info.push((info_key as usize, typ, n as usize, s..e));
        }
    }
    /// parse indiv fields, complicated field will need further processing
    fn parse_indv(&mut self) {
        let mut reader = std::io::Cursor::new(self.buf_indiv.as_slice());
        self.gt.clear();
        for _idx in 0..(self.n_fmt as usize) {
            let fmt_key = read_single_typed_integer(&mut reader);
            let (typ, n) = read_typed_descriptor_bytes(&mut reader);
            let width = bcf2_typ_width(typ);
            let s = reader.position() as usize;
            let e = s + width * self.n_sample as usize * n as usize;
            reader
                .seek(std::io::SeekFrom::Current((e - s) as i64))
                .unwrap();
            self.gt.push((fmt_key as usize, typ, n as usize, s..e));
        }
    }

    /// get chromosome offset
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

    /// Returns an iterator over the genotype values in the record's FORMAT field.
    /// See example in `test_read_fmt_gt`
    /// If no FORMAT/GT field available, the returned iterator will have items.
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
                            e.2 as usize * self.n_sample as usize,
                            &self.buf_indiv[e.3.start..e.3.end],
                        );
                    }
                });
                it
            }
        }
    }

    /// Returns an iterator over all values for a field in the record's FORMATs (indiv).
    /// See example in `test_read_fmt_ad`
    pub fn fmt_field(&self, fmt_key: usize) -> NumericValueIter<'_> {
        // default iterator
        let mut it = NumericValueIter::default();

        // find the right field for gt
        self.gt.iter().for_each(|e| {
            if e.0 == fmt_key {
                it = iter_typed_integers(
                    e.1,
                    e.2 as usize * self.n_sample as usize,
                    &self.buf_indiv[e.3.start..e.3.end],
                );
            }
        });
        it
    }

    /// get 0-based position (bp) value
    pub fn pos(&self) -> i32 {
        self.pos
    }

    /// Returns the ranges of bytes in buf_shared for all alleles in the record.
    ///
    /// Examples:
    ///
    /// // let ref_rng = &record.alleles()[0];
    /// // let alt1_rng = &record.alleles()[1];
    /// // let ref_str = std::str::from_utf8(&record.buf_shared()[ref_rng.start..ref_rng.end]).unwrap();
    /// // let alt1_str =
    /// //     std::str::from_utf8(&record.buf_shared()[alt1_rng.start..alt1_rng.end]).unwrap();
    ///
    pub fn alleles(&self) -> &[Range<usize>] {
        &self.alleles[..]
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
/// whether the file has the magic number for gzip (1f 8b)
pub fn smart_reader(p: impl AsRef<std::path::Path>) -> Box<dyn std::io::Read> {
    let mut f = std::fs::File::open(p.as_ref()).expect("can not open file");
    if (f.read_u8().expect("can not read first byte") == 0x1fu8)
        && (f.read_u8().expect("can not read second byte") == 0x8bu8)
    {
        // gzip format
        f.rewind().unwrap();
        Box::new(flate2::read::MultiGzDecoder::new(f))
    } else {
        // not gzip format
        f.rewind().unwrap();
        Box::new(std::io::BufReader::new(f))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_find_fmt_idx() {
        let mut f = smart_reader("testdata/test.bcf");
        let s = read_header(&mut f);
        let header = Header::from_string(&s);
        let key_found = header.get_idx_from_dictionary_str("FORMAT", "GT").unwrap();
        assert_eq!(key_found, header.get_fmt_gt_id().unwrap());
    }

    #[test]
    fn test_read_samples() {
        // read data generated by bcftools
        // bcftools query -l test.bcf | bgzip -c > test_samples.gz
        let mut samples_str = String::new();
        smart_reader("./testdata/test_samples.gz")
            .read_to_string(&mut samples_str)
            .unwrap();

        // read data via bcf-reader
        let mut f = smart_reader("testdata/test.bcf");
        let s = read_header(&mut f);
        let header = Header::from_string(&s);
        let samples_str2 = header.get_samples().join("\n");

        // compare bcftools results and bcf-reader results
        assert_eq!(samples_str.trim(), samples_str2.trim());
    }

    #[test]
    fn test_read_site_chrom() {
        // read data generated by bcftools
        // bcftools query -f '%CHROM\n' test.bcf | bgzip -c > test_chrom.gz
        let mut chrom_str = String::new();
        smart_reader("testdata/test_chrom.gz")
            .read_to_string(&mut chrom_str)
            .unwrap();

        // read data via bcf-reader
        let mut f = smart_reader("testdata/test.bcf");
        let s = read_header(&mut f);
        let header = Header::from_string(&s);
        let mut record = Record::default();
        let mut chrom_str2 = Vec::<u8>::new();
        use std::io::Write;
        while let Ok(_) = record.read(&mut f) {
            write!(
                chrom_str2,
                "{}\n",
                header.get_chrname(record.chrom as usize)
            )
            .unwrap();
        }
        let chrom_str2 = String::from_utf8(chrom_str2).unwrap();

        // compare bcftools results and bcf-reader results
        assert_eq!(chrom_str, chrom_str2);
    }

    #[test]
    fn test_read_fmt_gt() {
        // read data generated by bcftools
        // bcftools query -f '[\t%GT]\n' test.bcf | bgzip -c > test_gt.gz
        let mut gt_str = String::new();
        smart_reader("testdata/test_gt.gz")
            .read_to_string(&mut gt_str)
            .unwrap();

        // read data via bcf-reader
        let mut f = smart_reader("testdata/test.bcf");
        let s = read_header(&mut f);
        let header = Header::from_string(&s);
        let mut record = Record::default();
        let mut gt_str2 = Vec::<u8>::new();

        use std::io::Write;
        while let Ok(_) = record.read(&mut f) {
            for (i, bn) in record.fmt_gt(&header).enumerate() {
                let (noploidy, dot, phased, allele) = bn.gt_val();
                assert_eq!(noploidy, false); // missing ploidy
                let mut sep = '\t';
                if i % 2 == 1 {
                    if phased {
                        sep = '|';
                    } else {
                        sep = '/';
                    }
                }
                if dot {
                    write!(gt_str2, "{sep}.").unwrap();
                } else {
                    write!(gt_str2, "{sep}{allele}").unwrap();
                }
            }
            write!(gt_str2, "\n").unwrap();
        }

        let gt_str2 = String::from_utf8(gt_str2).unwrap();

        // compare bcftools results and bcf-reader results
        // assert_eq!(gt_str, gt_str2);
        for (a, b) in gt_str
            .split(|c| (c == '\n') || (c == '\t'))
            .zip(gt_str2.split(|c| (c == '\n') || (c == '\t')))
        {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_read_fmt_ad() {
        // read data generated by bcftools
        //  bcftools query -f '[\t%AD]\n' test.bcf | bgzip -c > test_ad.gz
        let mut ad_str = String::new();
        smart_reader("testdata/test_ad.gz")
            .read_to_string(&mut ad_str)
            .unwrap();

        // read data via bcf-reader
        let mut f = smart_reader("testdata/test.bcf");
        let s = read_header(&mut f);
        let header = Header::from_string(&s);
        let mut record = Record::default();
        let mut ad_str2 = Vec::<u8>::new();

        use std::io::Write;
        let ad_filed_key = header.get_idx_from_dictionary_str("FORMAT", "AD").unwrap();
        while let Ok(_) = record.read(&mut f) {
            for (i, val) in record.fmt_field(ad_filed_key).enumerate() {
                if i % record.n_allele as usize == 0 {
                    if ad_str2.last().map(|c| *c == b',') == Some(true) {
                        ad_str2.pop(); // trim last allele separator
                    }
                    ad_str2.push(b'\t'); // sample separator
                }
                match val.int_val() {
                    None => {}
                    Some(ad) => {
                        write!(ad_str2, "{ad},").unwrap(); // allele separator
                    }
                }
            }
            // site separator
            *ad_str2.last_mut().unwrap() = b'\n'; // sample separator
        }

        let ad_str2 = String::from_utf8(ad_str2).unwrap();

        // compare bcftools results and bcf-reader results
        for (a, b) in ad_str
            .split(|c| (c == '\n') || (c == '\t'))
            .zip(ad_str2.split(|c| (c == '\n') || (c == '\t')))
        {
            assert_eq!(a, b);
        }
    }
    #[test]
    fn test_read_site_pos() {
        // read data generated by bcftools
        // bcftools query -f '%POS\n' test.bcf | bgzip -c > test_pos.gz
        let mut pos_str = String::new();
        smart_reader("testdata/test_pos.gz")
            .read_to_string(&mut pos_str)
            .unwrap();

        // read data via bcf-reader
        let mut f = smart_reader("testdata/test.bcf");
        let _s = read_header(&mut f);
        let mut record = Record::default();
        let mut pos_str2 = Vec::<u8>::new();

        use std::io::Write;
        while let Ok(_) = record.read(&mut f) {
            write!(pos_str2, "{}\n", record.pos + 1).unwrap();
        }

        let pos_str2 = String::from_utf8(pos_str2).unwrap();
        // compare bcftools results and bcf-reader results
        assert_eq!(pos_str, pos_str2);
    }

    #[test]
    fn test_read_site_alleles() {
        // read data generated by bcftools
        // bcftools query -f '%REF,%ALT\n' test.bcf | bgzip -c > test_allele.gz
        let mut allele_str = String::new();
        smart_reader("testdata/test_allele.gz")
            .read_to_string(&mut allele_str)
            .unwrap();

        // read data via bcf-reader
        let mut f = smart_reader("testdata/test.bcf");
        let _s = read_header(&mut f);
        let mut record = Record::default();
        let mut allele_str2 = Vec::<u8>::new();

        while let Ok(_) = record.read(&mut f) {
            for rng in record.alleles.iter() {
                let slice = &record.buf_shared[rng.start..rng.end];
                allele_str2.extend(slice);
                allele_str2.push(b',');
            }
            *allele_str2.last_mut().unwrap() = b'\n';
        }

        let allele_str2 = String::from_utf8(allele_str2).unwrap();
        // compare bcftools results and bcf-reader results
        assert_eq!(allele_str, allele_str2);
    }
}
