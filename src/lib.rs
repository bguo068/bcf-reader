use byteorder::{LittleEndian, ReadBytesExt};
use std::ops::Range;
use std::{collections::HashMap, io::Seek};

pub struct QuotedSplitter<'a> {
    data: &'a str,
    in_quotes: bool,
    sep: char,
    quote: char,
}

impl<'a> QuotedSplitter<'a> {
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

#[test]
fn test_quoted_splitter() {
    let input_string = "hello,\"world, this is fun\",test";
    let result: Vec<_> = QuotedSplitter::new(input_string, ',', '"').collect();
    assert_eq!(result, vec!["hello", "\"world, this is fun\"", "test"]);
}

#[derive(Debug)]
pub struct Header {
    dict_strings: HashMap<usize, HashMap<String, String>>,
    dict_contigs: HashMap<usize, HashMap<String, String>>,
    samples: Vec<String>,
    fmt_gt_idx: usize,
}
impl Header {
    pub fn from_string(text: &str) -> Self {
        let mut dict_strings = HashMap::<usize, HashMap<String, String>>::new();
        let mut dict_contigs = HashMap::<usize, HashMap<String, String>>::new();
        let mut samples = Vec::<String>::new();

        // implicit FILTER/PASS header lines
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
            } else if line.trim().len() == 0 {
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

        // reorder items if the header line has IDX key
        let mut fmt_gt_idx = 0;
        for (k, m) in dict_strings.iter() {
            if (&m["Dictionary"] == "FORMAT") && (&m["ID"] == "GT") {
                fmt_gt_idx = *k;
            }
        }

        Self {
            dict_strings,
            dict_contigs,
            samples,
            fmt_gt_idx,
        }
    }

    pub fn get_idx_from_dictionary_str(&self, dictionary: &str, field: &str) -> Option<usize> {
        for (k, m) in self.dict_strings.iter() {
            if (&m["Dictionary"] == dictionary) && (&m["ID"] == field) {
                return Some(*k);
            }
        }
        None
    }

    pub fn get_chrname(&self, idx: usize) -> &str {
        &self.dict_contigs[&idx]["ID"]
    }
    pub fn get_fmt_gt_id(&self) -> usize {
        self.fmt_gt_idx
    }
    pub fn get_contigs(&self) -> &HashMap<usize, HashMap<String, String>> {
        &self.dict_contigs
    }
    pub fn get_dict_strings(&self) -> &HashMap<usize, HashMap<String, String>> {
        &self.dict_strings
    }
    pub fn get_samples(&self) -> &Vec<String> {
        &self.samples
    }
}

pub trait Bcf2Number {
    fn is_missing(&self) -> bool;
    fn is_end_of_vector(&self) -> bool;
    fn is_reserved_value(&self) -> bool;
}

impl Bcf2Number for u8 {
    fn is_missing(&self) -> bool {
        *self == 0x80
    }
    fn is_end_of_vector(&self) -> bool {
        *self == 0x81
    }
    fn is_reserved_value(&self) -> bool {
        (*self >= 0x80) && (*self <= 0x87)
    }
}

impl Bcf2Number for u16 {
    fn is_missing(&self) -> bool {
        *self == 0x8000
    }
    fn is_end_of_vector(&self) -> bool {
        *self == 0x8001
    }
    fn is_reserved_value(&self) -> bool {
        (*self >= 0x8000) && (*self <= 0x8007)
    }
}

impl Bcf2Number for u32 {
    fn is_missing(&self) -> bool {
        *self == 0x80000000
    }
    fn is_end_of_vector(&self) -> bool {
        *self == 0x80000001
    }
    fn is_reserved_value(&self) -> bool {
        (*self >= 0x80000000) && (*self <= 0x80000007)
    }
}
impl Bcf2Number for f32 {
    fn is_missing(&self) -> bool {
        (*self) as u32 == 0x7FC00000
    }
    fn is_end_of_vector(&self) -> bool {
        (*self) as u32 == 0x7FC00001
    }
    fn is_reserved_value(&self) -> bool {
        ((*self) as u32 >= 0x7FC00001) && ((*self) as u32 <= 0x7FC00007)
    }
}

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

#[derive(Debug)]
pub enum NumbericValue {
    U8(u8),
    U16(u16),
    U32(u32),
    F32(f32),
}

impl From<u8> for NumbericValue {
    fn from(value: u8) -> Self {
        Self::U8(value)
    }
}
impl From<u16> for NumbericValue {
    fn from(value: u16) -> Self {
        Self::U16(value)
    }
}
impl From<u32> for NumbericValue {
    fn from(value: u32) -> Self {
        Self::U32(value)
    }
}
impl From<f32> for NumbericValue {
    fn from(value: f32) -> Self {
        Self::F32(value)
    }
}

impl NumbericValue {
    pub fn int_val(&self) -> Option<u32> {
        match *self {
            Self::U8(x) if !x.is_missing() => Some(x as u32),
            Self::U16(x) if !x.is_missing() => Some(x as u32),
            Self::U32(x) if !x.is_missing() => Some(x as u32),
            _ => None,
        }
    }
    pub fn float_val(&self) -> Option<f32> {
        match *self {
            Self::F32(x) if !x.is_missing() => Some(x),
            _ => None,
        }
    }

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
pub struct NumberIter<'r> {
    reader: std::io::Cursor<&'r [u8]>,
    typ: u8,
    len: usize,
    cur: usize,
}

impl<'r> Iterator for NumberIter<'r> {
    type Item = NumbericValue;
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
                    Some(self.reader.read_f32::<LittleEndian>().unwrap().into())
                }
                _ => panic!(),
            }
        }
    }
}

pub fn iter_typed_integers(typ: u8, n: usize, buffer: &[u8]) -> NumberIter {
    NumberIter {
        reader: std::io::Cursor::new(buffer),
        typ,
        len: n,
        cur: 0,
    }
}

/// if 0 is return, it means the string is missing
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

#[derive(Default, Debug)]
pub struct Record {
    buf_shared: Vec<u8>,
    buf_indiv: Vec<u8>,
    chrom: i32,
    pos: i32,
    rlen: i32,
    qual: f32,
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
    /// read a record, copy bytes and separate fields
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
        self.qual = reader.read_f32::<LittleEndian>().unwrap();
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
    pub fn rlen(&self) -> i32 {
        self.rlen
    }
    pub fn qual(&self) -> Option<f32> {
        match self.qual.is_missing() {
            true => None,
            false => Some(self.qual),
        }
    }
    pub fn fmt_gt(&self, header: &Header) -> NumberIter<'_> {
        let fmt_gt_id = header.get_fmt_gt_id();
        // default iterator
        let mut it = NumberIter::default();

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
    pub fn fmt_field(&self, fmt_key: usize) -> NumberIter<'_> {
        // default iterator
        let mut it = NumberIter::default();

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
    pub fn pos(&self) -> i32 {
        self.pos
    }
    pub fn alleles(&self) -> &[Range<usize>] {
        &self.alleles[..]
    }
    pub fn buf_gt(&self) -> &[u8] {
        &self.buf_indiv[..]
    }
    pub fn buf_site(&self) -> &[u8] {
        &self.buf_shared[..]
    }
}

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
        assert_eq!(key_found, header.fmt_gt_idx);
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
