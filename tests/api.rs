use std::path::Path;

use sketchlib::api::{self, DistOptions, SketchOptions};
use sketchlib::cli::{DEFAULT_MINCOUNT, DEFAULT_MINQUAL};
use sketchlib::hashing::{HashType, DEFAULT_LEVEL};

pub mod common;
use crate::common::*;

#[test]
fn api_sketch_dist_and_repeated_thread_setup() {
    let sandbox = TestSetup::setup();
    let rfile = sandbox.create_named_rfile(
        "api_file_list.txt",
        &["R6.fa.gz", "TIGR4.fa.gz", "14412_3#82.contigs_velvet.fa.gz"],
    );
    let rfile = sandbox.file_string(&rfile, TestDir::Output);

    for output_name in ["api_db1", "api_db2"] {
        let output = sandbox.file_string(output_name, TestDir::Output);
        let paths = api::sketch(SketchOptions {
            seq_files: None,
            file_list: Some(rfile.clone()),
            output,
            kmers: vec![21],
            sketch_size: 1000,
            seq_type: HashType::DNA,
            level: DEFAULT_LEVEL,
            concat_fasta: false,
            single_strand: false,
            min_count: DEFAULT_MINCOUNT,
            min_qual: DEFAULT_MINQUAL,
            threads: 2,
            quiet: true,
        })
        .expect("API sketch should succeed");
        assert!(Path::new(&paths.skm).exists());
        assert!(Path::new(&paths.skd).exists());
    }

    let db_prefix = sandbox.file_string("api_db1", TestDir::Output);
    let info = api::db_info(&db_prefix, false).expect("API db_info should succeed");
    assert!(!info.inverted);
    assert_eq!(info.kmers, vec![21]);
    assert_eq!(info.n_samples, 3);

    let dist_path = sandbox.file_string("api.dists", TestDir::Output);
    api::dist(DistOptions {
        ref_db: db_prefix,
        query_db: None,
        output: Some(dist_path.clone()),
        knn: Some(1),
        subset: None,
        kmer: Some(21),
        ani: true,
        threads: 2,
        ref_completeness_file: None,
        query_completeness_file: None,
        completeness_cutoff: 0.64,
        quiet: true,
    })
    .expect("API dist should succeed");
    assert!(Path::new(&dist_path).exists());
}
