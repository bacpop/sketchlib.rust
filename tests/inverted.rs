use sketchlib::hashing::HashType;
// use snapbox::cmd::{cargo_bin, Command};
pub mod common;
use crate::common::*;
use sketchlib::multisketch::MultiSketch;
// use sketchlib::sketch::Sketch;
use sketchlib::sketch::sketch_files;

#[test]
fn test_identical_sequences() {
    let sandbox = TestSetup::setup();

    let input_files = vec![(
        "seq1".to_string(),
        sandbox.file_string("inverted_test_genome1.fa", TestDir::Input),
        None,
    )];

    let sketch_size = 3;
    let seq_type = HashType::DNA;
    let kmers: &[usize] = &[2];

    let mut sketches = sketch_files(
        format!("{}/test_identical", sandbox.get_wd()).as_str(),
        &input_files,
        false,       // concat_fasta
        &kmers,      // k-mer sizes
        sketch_size, // sketch_size
        &seq_type,
        false, // rc
        1,     // min_count
        0,     // min_qual
        true,  // inverted
        false, // quiet
    );

    let mut multisketches = MultiSketch::new(&mut sketches, sketch_size, &kmers, seq_type, true);

    // multisketches.update_sketches();

    let inverted_index1 = MultiSketch::invert_index(&multisketches);
    let inverted_index2 = MultiSketch::invert_index(&multisketches);

    // test: same size
    assert_eq!(
        inverted_index1.len(),
        inverted_index2.len(),
        "Both inverted indices should have the same number of entries"
    );

    // test: key and value for both
    for (sig, sketch_indices1) in &inverted_index1 {
        assert!(
            inverted_index2.contains_key(sig),
            "Signature {} should exist in both indices",
            sig
        );

        let sketch_indices2 = inverted_index2.get(sig).unwrap();
        assert_eq!(
            sketch_indices1, sketch_indices2,
            "Sets of sketch indices for signature {} should be identical",
            sig
        );
    }
}

// #[test]
// fn test_different_sequences() {
//     let sandbox = TestSetup::setup();

//     let input_files = vec![
//         (
//             "seq1".to_string(),
//             sandbox.file_string(
//                 "inverted_test_genome1.fa",
//                 TestDir::Input,
//             ),
//             None,
//         ),
//         (
//             "seq2".to_string(),
//             sandbox.file_string("inverted_test_genome2.fa", TestDir::Input),
//             None,
//         ),
//     ];

//     assert_eq!(sketches.len(), 2, "Should have two sketches");
//     assert_ne!(
//         sketches[0].get_usigs(),
//         sketches[1].get_usigs(),
//         "Different sequences should produce different sketches"
//     );
// }
