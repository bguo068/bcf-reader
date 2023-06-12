# noodles_bcf_gt_reading_example

The `rust_htslib` crate only works on linux and macos, not windows (?).
The `noodles` crate is a pure rust library for many bioinformatic file format
and works for windows, linux and macos.
However, currently, the `noodles::bcf` api of reading genotype data from bcf is
quite slow which involes a lot memory allocation.

One way to get around memory allocation is get the raw bytes of genotype data
per record and parse the FORMAT/GT field according to BCF specification
https://samtools.github.io/hts-specs/VCFv4.2.pdf

Steps:
1. get the raw bytes can be done by `noodles::bcf::Reader::read_lazy_record` method
and `bcf::lazy::Record::genotypes().as_ref()`
2. get genotype field information and skip fields until the GT field is found
3. parse GT field: `allele = raw >> 1 - 1`