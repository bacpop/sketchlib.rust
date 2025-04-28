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
    fn test_is_compatible_with() {
        let sandbox = TestSetup::setup();
        let ref_db1 = sandbox.file_string("sketches1", TestDir::Input);
        let ref_db2 = sandbox.file_string("sketches2", TestDir::Input);
        let ref_db3 = sandbox.file_string("sketches3", TestDir::Input);

        let sketches1: MultiSketch = MultiSketch::load_metadata(&ref_db1).expect(&format!(
            "Could not read sketch metadata from {}.skm",
            ref_db1
        ));

        let sketches2: MultiSketch = MultiSketch::load_metadata(&ref_db2).expect(&format!(
            "Could not read sketch metadata from {}.skm",
            ref_db2
        ));

        let sketches3: MultiSketch = MultiSketch::load_metadata(&ref_db3).expect(&format!(
            "Could not read sketch metadata from {}.skm",
            ref_db3
        ));

        // Test case 1: Sketches are compatible
        assert!(
            sketches1.is_compatible_with(&sketches1),
            "Sketches should be compatible"
        );

        // Test case 2: Sketches are not compatible
        // Because of differeing sketch size
        assert!(
            !sketches1.is_compatible_with(&sketches2),
            "Sketches should not be compatible"
        );
        // Because of difference in number of k-mers
        assert!(
            !sketches1.is_compatible_with(&sketches3),
            "Sketches should not be compatible"
        );
    }

    #[test]
    fn test_merge_sketches() {
        let sandbox = TestSetup::setup();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "part1"])
            .assert()
            .success();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "part2"])
            .assert()
            .success();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "merged_ref"])
            .assert()
            .success();

        // Overlapping labels fails
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("merge")
            .arg("part1")
            .arg("part1")
            .arg("-v")
            .args(&["-o", "merged_test"])
            .assert()
            .failure();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("merge")
            .arg("part1")
            .arg("part2")
            .arg("-v")
            .args(&["-o", "merged_test"])
            .assert()
            .success();

        // Check .skd the same
        let predicate_file = predicate::path::eq_file(Path::new(
            &sandbox.file_string("merged_test.skd", TestDir::Output),
        ));
        assert_eq!(
            true,
            predicate_file.eval(Path::new(
                &sandbox.file_string("merged_ref.skd", TestDir::Output)
            )),
            "Merged sketch data does not match"
        );

        // Check .skm the same
        let merged_sketch: MultiSketch =
            MultiSketch::load_metadata(&sandbox.file_string("merged_test", TestDir::Output))
                .expect("Failed to load output merged sketch");
        let expected_sketch =
            MultiSketch::load_metadata(&sandbox.file_string("merged_ref", TestDir::Output))
                .expect("Failed to load expected merged sketch");
        assert_eq!(
            merged_sketch, expected_sketch,
            "Merged sketch metadata does not match"
        );
    }
}
