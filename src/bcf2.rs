use byteorder::{LittleEndian, ReadBytesExt};

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

pub enum IntValue {
    U8(u8),
    U16(u16),
    U32(u32),
}

impl From<u8> for IntValue {
    fn from(value: u8) -> Self {
        Self::U8(value)
    }
}
impl From<u16> for IntValue {
    fn from(value: u16) -> Self {
        Self::U16(value)
    }
}
impl From<u32> for IntValue {
    fn from(value: u32) -> Self {
        Self::U32(value)
    }
}

impl IntValue {
    pub fn val(&self) -> Option<u32> {
        match *self {
            Self::U8(x) if !x.is_missing() => Some(x as u32),
            Self::U16(x) if !x.is_missing() => Some(x as u32),
            Self::U32(x) if !x.is_missing() => Some(x as u32),
            _ => None,
        }
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

pub struct IntValueIter<'r, R>
where
    R: std::io::Read + ReadBytesExt,
{
    reader: &'r mut R,
    typ: u8,
    len: usize,
    cur: usize,
}

impl<'r, R> IntValueIter<'r, R>
where
    R: std::io::Read + ReadBytesExt,
{
    /// skip over the END_OF_VECTOR value
    pub fn finish(&mut self) {
        let mut buffer = [0u8; 1024];
        let nbytes = match self.typ {
            1 => 1,
            2 => 2,
            3 => 4,
            _ => panic!(),
        };
        let a = self.cur * nbytes;
        let b = (self.len + 1) * nbytes;
        let mut s = a;
        let mut e;
        while s < b {
            e = s + 1024;
            if e > b {
                e = b;
            }
            self.reader.read(&mut buffer[s..e]).unwrap();
            s = e;
        }
    }
}

impl<'r, R> Iterator for IntValueIter<'r, R>
where
    R: std::io::Read + ReadBytesExt,
{
    type Item = IntValue;
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
                _ => panic!(),
            }
        }
    }
}

pub fn iter_typed_integers<R>(reader: &mut R) -> IntValueIter<R>
where
    R: std::io::Read + ReadBytesExt,
{
    let (typ, n) = read_typed_descriptor_bytes(reader);
    IntValueIter {
        reader,
        typ,
        len: n,
        cur: 0,
    }
}

pub struct FloatValueIter<'r, R>
where
    R: std::io::Read + ReadBytesExt,
{
    reader: &'r mut R,
    len: usize,
    cur: usize,
}

impl<'r, R> FloatValueIter<'r, R>
where
    R: std::io::Read + ReadBytesExt,
{
    /// skip over the END_OF_VECTOR value
    pub fn finish(&mut self) {
        let mut buffer = [0u8; 1024];
        let nbytes = 4;
        let a = self.cur * nbytes;
        let b = (self.len + 1) * nbytes;
        let mut s = a;
        let mut e;
        while s < b {
            e = s + 1024;
            if e > b {
                e = b;
            }
            self.reader.read(&mut buffer[s..e]).unwrap();
            s = e;
        }
    }
}

impl<'r, R> Iterator for FloatValueIter<'r, R>
where
    R: std::io::Read + ReadBytesExt,
{
    type Item = f32;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.len {
            None
        } else {
            self.cur += 1;
            Some(self.reader.read_f32::<LittleEndian>().unwrap())
        }
    }
}

pub fn iter_typed_floats<R>(reader: &mut R) -> FloatValueIter<R>
where
    R: std::io::Read + ReadBytesExt,
{
    let (typ, n) = read_typed_descriptor_bytes(reader);
    assert_eq!(typ, 0x7);
    FloatValueIter {
        reader,
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
    reader.read(&mut text).unwrap();

    String::from_utf8(text).unwrap()
}

#[test]
fn test_read_header() {
    let mut f = std::fs::File::open("test_flat.bcf").unwrap();
    let s = read_header(&mut f);
    println!("{s}");
}
