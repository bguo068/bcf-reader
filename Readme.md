[![Crates.io](https://img.shields.io/crates/d/bcf_reader.svg)](https://crates.io/crates/bcf_reader)
[![Crates.io](https://img.shields.io/crates/v/bcf_reader.svg)](https://crates.io/crates/bcf_reader)
[![Crates.io](https://img.shields.io/crates/l/bcf_reader.svg)](https://crates.io/crates/bcf_reader)
[![docs.rs](https://docs.rs/bcf_reader/badge.svg)](https://docs.rs/bcf_reader)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/bguo068/bcf-reader/rust.yml?branch=main&label=tests)


# bcf_reader
This is an attempt to create a small, lightweight, pure-Rust library to allow
efficient, cross-platform access to genotype data in BCF files.

Currently, the `rust_htslib` crate works only on Linux and macOS (not Windows?).
The noodles crate is a pure Rust library for many bioinformatic file formats and
works across Windows, Linux, and macOS. However, the `noodles` API for reading
genotype data from BCF files can be slow due to its memory allocation patterns.
Additionally, both crates have a large number of dependencies, as they provide
many features and support a wide range of file formats.

One way to address the memory allocation and dependency issues is to manually
parse BCF records according to its specification
(https://samtools.github.io/hts-specs/VCFv4.2.pdf) and use iterators whenever
possible, especially for the per-sample fields, like GT and AD.

Note: This crate is in its early stages of development.

## Usage

```rust
use bcf_reader::*;
let mut reader = smart_reader("testdata/test2.bcf");
let header = Header::from_string(&read_header(&mut reader));
// find key for a field in INFO or FORMAT or FILTER
let key = header.get_idx_from_dictionary_str("FORMAT", "GT").unwrap();
// access header dictionary
let d = &header.dict_strings()[&key];
assert_eq!(d["ID"], "GT");
assert_eq!(d["Dictionary"], "FORMAT");
/// get chromosome name
assert_eq!(header.get_chrname(0), "Pf3D7_01_v3");
let fmt_ad_key = header
    .get_idx_from_dictionary_str("FORMAT", "AD")
    .expect("FORMAT/AD not found");
let info_af_key = header
    .get_idx_from_dictionary_str("INFO", "AF")
    .expect("INFO/AF not found");

// this can be and should be reused to reduce allocation
let mut record = Record::default();
while let Ok(_) = record.read(&mut reader) {
    let pos = record.pos();

    // use byte ranges and shared buffer to get allele string values
    let allele_byte_ranges = record.alleles();
    let share_buf = record.buf_shared();
    let ref_rng = &allele_byte_ranges[0];
    let ref_allele_str =
        std::str::from_utf8(&share_buf[ref_rng.start..ref_rng.end]).unwrap();
    let alt1_rng = &allele_byte_ranges[1];
    let alt1_allele_str =
        std::str::from_utf8(&share_buf[alt1_rng.start..alt1_rng.end]).unwrap();
    // ...

    // access FORMAT/GT via iterator
    for nv in record.fmt_gt(&header) {
        let (has_ploidy, is_missing, is_phased, allele_idx) = nv.gt_val();
        // ...
    }

    // access FORMAT/AD via iterator
    for nv in record.fmt_field(fmt_ad_key) {
        match nv.int_val() {
            None => {}
            Some(ad) => {
                // ...
            }
        }
        // ...
    }

    // access FILTERS via itertor
    record.filters().for_each(|nv| {
        let filter_key = nv.int_val().unwrap() as usize;
        let dict_string_map = &header.dict_strings()[&filter_key];
        let filter_name = &dict_string_map["ID"];
        // ...
    });

    // access INFO/AF via itertor
    record.info_field_numeric(info_af_key).for_each(|nv| {
        let af = nv.float_val().unwrap();
        // ...
    });
}
```
