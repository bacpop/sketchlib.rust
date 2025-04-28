use predicates::prelude::*;
use snapbox::cmd::{cargo_bin, Command};
use std::path::Path;

pub mod common;
use crate::common::*;

use sketchlib::sketch::multisketch::MultiSketch;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_delete_sketches() {
        let sandbox = TestSetup::setup();

        let rfile1 = sandbox.create_rfile(&vec![
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ]);
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile1])
            .arg("-v")
            .args(&["-o", "full_db"])
            .assert()
            .success();

        let rfile2 = sandbox.create_rfile(&vec![
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
        ]);
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile2])
            .arg("-v")
            .args(&["-o", "result_db"])
            .assert()
            .success();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("delete")
            .arg("full_db")
            // .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            // .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("delete_test.txt", TestDir::Input))
            .arg("deleted_db")
            .assert()
            .success();

        let predicate_file = predicate::path::eq_file(Path::new(
            &sandbox.file_string("deleted_db.skd", TestDir::Output),
        ));
        assert_eq!(
            true,
            predicate_file.eval(Path::new(
                &sandbox.file_string("result_db.skd", TestDir::Output)
            )),
            "Merged sketch data does not match"
        );

        // Check .skm the same
        let merged_sketch: MultiSketch =
            MultiSketch::load_metadata(&sandbox.file_string("deleted_db", TestDir::Output))
                .expect("Failed to load output merged sketch");
        let expected_sketch =
            MultiSketch::load_metadata(&sandbox.file_string("result_db", TestDir::Output))
                .expect("Failed to load expected merged sketch");
        assert_eq!(
            merged_sketch, expected_sketch,
            "Merged sketch metadata does not match"
        );
    }
}
