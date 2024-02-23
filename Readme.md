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