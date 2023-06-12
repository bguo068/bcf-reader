use noodles::bcf::{self, lazy::Record};
use noodles::vcf::record::genotypes::keys::{key::Standard, Key};

fn main() {
    let filename = "test.bcf";
    let mut reader = std::fs::File::open(filename).map(bcf::Reader::new).unwrap();
    let header = reader.read_header().unwrap();
    let gt_key = Key::Standard(Standard::Genotype);
    let gt_fmt_key = header.formats()[&gt_key].idx().unwrap();
    let samples = header.sample_names();
    let nsam = samples.len();

    let mut gt = vec![0u8; nsam * 2];
    let mut record = Record::default();

    use std::io::{stdout, Write};
    let mut lock = stdout().lock();
    while reader.read_lazy_record(&mut record).unwrap() > 0 {
        let mut r = record.genotypes().as_ref().iter();

        // loop to find the gt field and parse gt (skipping non-gt field)
        loop {
            let (field_idx, num_element, element_width, bcf2_type) =
                read_bcf_gtbuf_field_info(&mut r);
            let total_value_byte = num_element * element_width * (nsam as usize);

            // if GT:
            if field_idx == gt_fmt_key {
                assert_eq!(element_width, 1);
                assert_eq!(num_element, 2);
                assert_eq!(bcf2_type, Bcf2Type::Int8);
                for i in 0..total_value_byte {
                    // +48 to covert int value to ascii char value
                    // this is only used for printing purpose
                    gt[i] = ((*r.next().unwrap()) >> 1) - 1 + 48;
                }
                break;
            } else {
                r.nth(total_value_byte);
            }
        }
        lock.write_all(&gt[..]).unwrap();
        writeln!(lock).unwrap();
    }
}

pub fn read_bcf_gtbuf_integer(it: &mut std::slice::Iter<u8>) -> usize {
    let mut fmt_type = *it.next().unwrap();
    // integer type and width
    let int_len = fmt_type >> 4;
    assert_eq!(int_len, 1);
    // integer type
    fmt_type &= 0xf;
    let int_width = match fmt_type {
        1 => 1usize,
        2 => 2,
        3 => 4,
        _ => panic!("not valid integer"),
    };
    // read bytes to integer
    let mut value = 0usize;
    for ibyte in 0..int_width {
        value |= (*it.next().unwrap() as usize) << (8 * ibyte);
    }
    value
}

#[derive(PartialEq, Eq, Debug)]
pub enum Bcf2Type {
    MISSING,
    Int8,
    Int16,
    Int32,
    Float32,
    Char8,
}

pub fn read_bcf_gtbuf_field_info(it: &mut std::slice::Iter<u8>) -> (usize, usize, usize, Bcf2Type) {
    // field index: used for identify which field it is
    let field_idx = read_bcf_gtbuf_integer(it);

    let fmt_type = *it.next().unwrap();

    // len is how man elements are there for each sample
    let num_element = match fmt_type >> 4 {
        15 => read_bcf_gtbuf_integer(it),
        x => x as usize,
    };

    // how many bytes is used by each element
    use Bcf2Type::*;
    let (element_width, bcf2_type) = match fmt_type & 0xf {
        0 => (0, MISSING), // MISSING
        1 => (1, Int8),    // Int8
        2 => (2, Int16),   // Int16
        3 => (4, Int32),   // Int32
        5 => (5, Float32), // Float32
        7 => (1, Char8),   // Int8 (char ascii)
        _ => panic!("invalid BCF2 type"),
    };

    (field_idx, num_element, element_width, bcf2_type)
}
