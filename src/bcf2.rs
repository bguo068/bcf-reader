use byteorder::{LittleEndian, ReadBytesExt};

pub trait TypeDescriptorByte {
    fn typ(&self) -> Bcf2Type;
    fn n(&self) -> usize;
}

pub trait Bcf2Number {
    fn is_missing(&self) -> bool;
    fn is_end_of_vector(&self) -> bool;
    fn is_reserved_value(&self) -> bool;
}

#[derive(Debug, PartialEq, Eq)]
pub enum Bcf2Type {
    MISSING,
    I8,
    I16,
    I32,
    F32,
    CHR,
}

impl TypeDescriptorByte for u8 {
    fn typ(&self) -> Bcf2Type {
        match *self & 0xff {
            0u8 => Bcf2Type::MISSING,
            1u8 => Bcf2Type::I8,
            2u8 => Bcf2Type::I16,
            3u8 => Bcf2Type::I32,
            5u8 => Bcf2Type::F32,
            7u8 => Bcf2Type::CHR,
            _ => panic!(),
        }
    }
    fn n(&self) -> usize {
        (*self >> 4) as usize
    }
}

#[test]
fn test_bcft2type() {
    assert_eq!(1u8.typ(), Bcf2Type::I8);
}
