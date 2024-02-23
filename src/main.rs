use byteorder::{LittleEndian, ReadBytesExt};
use noodles::bcf;
use noodles::vcf::record::genotypes::keys::{key::Standard, Key};
use std::io::{stdout, Read, Write};
pub mod bcf2;
use bcf2::*;


fn main () {
    // let mut f = std::fs::File::open("test_flat.bcf").map(std::io::BufReader::new).unwrap();
    let mut f = std::fs::File::open("test_flat.bcf").map(std::io::BufReader::new).unwrap();
    let s = read_header(&mut f);
    let header = Header::from_string(&s);
    let mut record = Record::default();

    // let mut buf = vec![0u8; 0];

    let mut cnt0 = 0;
    let mut cnt1 = 0;
    while let Ok(_) = record.read(&mut f) {
        // use std::io::Write;
        
        for bn in record.gt(&header) {
            // write!(buf, "{}",bn.gt_val().3 ).unwrap();
            let allele = bn.gt_val().3;
            // let allele = 0;
            if allele == 0 {
                cnt0 += 1;
            } else {
                cnt1 += 1;
            }
        }
        // write!(buf, "\n").unwrap();
    }
    // let buf = String::from_utf8(buf).unwrap();
    // let buf2 = std::fs::read_to_string("test_gt.txt").unwrap();
    // assert_eq!(buf, buf2);
    eprintln!("cnt0= {cnt0}, cnt1={cnt1}");
}
fn main2() {
    let filename = "test.bcf";

    // read bcf header
    let mut reader = std::fs::File::open(filename).map(bcf::Reader::new).unwrap();
    let header = reader.read_header().unwrap();

    // get fmt/gt field idx
    let gt_key = Key::Standard(Standard::Genotype);
    let gt_fmt_key = header.formats()[&gt_key].idx().unwrap();

    // get samples
    let samples = header.sample_names();
    let nsam = samples.len();

    // genotype buffer
    let mut gt = vec![0u8; nsam * 2];

    // record reading buffer
    let mut var_bytes = Vec::<u8>::new();
    let mut fmt_bytes = Vec::<u8>::new();
    let mut s_bytes = Vec::<u8>::new();
    let mut r = reader.into_inner();

    let mut lock = stdout().lock();

    let mut irecord = 0usize;

    // loop over bcf records
    loop {
        // shared block size
        let l_shared = match r.read_u32::<LittleEndian>() {
            Ok(x) => x as usize,
            Err(_) => {
                break; // break when reaching error due to EOF
            }
        };

        // indv block size
        let l_indv = r.read_u32::<LittleEndian>().unwrap() as usize;

        // set buffer sizes
        if l_shared > var_bytes.capacity() {
            let additional = l_shared - var_bytes.capacity();
            var_bytes.reserve(additional);
        }
        if l_indv > fmt_bytes.capacity() {
            let additional = l_indv - fmt_bytes.capacity();
            fmt_bytes.reserve(additional);
        }
        unsafe {
            var_bytes.set_len(l_shared);
            fmt_bytes.set_len(l_indv)
        };

        // read bytes to shared/indv blocks
        r.read_exact(&mut var_bytes[..]).unwrap();
        r.read_exact(&mut fmt_bytes[..]).unwrap();

        irecord += 1;

        let mut r = var_bytes.as_slice();

        // read shared block info
        let chrom = r.read_i32::<LittleEndian>().unwrap();
        let pos = r.read_i32::<LittleEndian>().unwrap();
        let rlen = r.read_i32::<LittleEndian>().unwrap();
        let qual = r.read_f32::<LittleEndian>().unwrap();
        print!("irecord={irecord}, chrom={chrom}, pos={pos}, rlen={rlen}, qual={qual}, ");
        let _n_info = r.read_u16::<LittleEndian>().unwrap();
        let n_allele = r.read_u16::<LittleEndian>().unwrap();
        let _n_sample = r.read_u24::<LittleEndian>().unwrap();
        let _n_fmt = r.read_u8().unwrap();
        // print!("n_info={n_info}, n_allele={n_allele}, n_sample={n_sample}, n_fmt={n_fmt}");

        // make iterator of shared block bytes
        let mut r = r.iter();

        // read typed string -- ids
        read_bcf2_typed_string(&mut r, &mut s_bytes);

        // read a list of typed string -- list of alleles
        for _ in 0..n_allele {
            read_bcf2_typed_string(&mut r, &mut s_bytes);
            print!("{}", std::str::from_utf8(&s_bytes).unwrap());
        }
        print!(",\n");

        // skip Filter/INFO parsing and go to indv block

        // make iterator of indv block bytes
        let mut r = fmt_bytes.as_slice().iter();

        // loop to find the gt field and parse gt (skipping non-gt field)
        loop {
            // values each fmt field are array of arrays
            // (nsam, num_element)
            let (field_idx, num_element, element_width, bcf2_type) = read_bcf_gt_field_info(&mut r);
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

pub fn read_bcf2_typed_int(it: &mut std::slice::Iter<u8>) -> usize {
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

pub fn read_bcf2_typed_string(it: &mut std::slice::Iter<u8>, s: &mut Vec<u8>) {
    let mut fmt_type = *it.next().unwrap();
    // integer type and width
    let mut int_len = (fmt_type >> 4).into();
    // integer type
    fmt_type &= 0xf;
    assert_eq!(fmt_type, 7);
    if int_len == 15 {
        // overflow
        int_len = read_bcf2_typed_int(it);
    }

    s.clear();
    s.extend(it.as_slice()[..int_len].iter());
    if int_len > 0 {
        it.nth(int_len - 1);
    }
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

pub fn read_bcf_gt_field_info(it: &mut std::slice::Iter<u8>) -> (usize, usize, usize, Bcf2Type) {
    // field index: used for identify which field it is
    let field_idx = read_bcf2_typed_int(it);

    let fmt_type = *it.next().unwrap();

    // len is how man elements are there for each sample
    let num_element = match fmt_type >> 4 {
        15 => read_bcf2_typed_int(it),
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
