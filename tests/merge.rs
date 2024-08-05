use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

use sketchlib::multisketch::MultiSketch;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_is_compatible_with() {
        let sandbox = TestSetup::setup();
        let ref_db1 = sandbox.file_string("sketches1", TestDir::Input);
        let ref_db2 = sandbox.file_string("sketches2", TestDir::Input);
        let ref_db3 = sandbox.file_string("sketches3", TestDir::Input);

        let sketches1: MultiSketch = MultiSketch::load(&ref_db1).expect(&format!(
            "Could not read sketch metadata from {}.skm",
            ref_db1
        ));

        let sketches2: MultiSketch = MultiSketch::load(&ref_db2).expect(&format!(
            "Could not read sketch metadata from {}.skm",
            ref_db2
        ));

        let sketches3: MultiSketch = MultiSketch::load(&ref_db3).expect(&format!(
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
            .arg("merge")
            .arg(sandbox.file_string("test_merge_sketches1", TestDir::Input))
            .arg(sandbox.file_string("test_merge_sketches2", TestDir::Input))
            .arg("-v")
            .args(&["-o", "merge_test"])
            .assert()
            .success();

        let merged_sketch: MultiSketch =
            MultiSketch::load(&sandbox.file_string("merge_test", TestDir::Output))
                .expect("Failed to load output merged sketch");
        let expected_sketch =
            MultiSketch::load(&sandbox.file_string("sketches_all", TestDir::Correct))
                .expect("Failed to load expected merged sketch");

        // assess if sketches are the same
        assert_eq!(
            merged_sketch, expected_sketch,
            "Merged sketch does not match expected sketch"
        );
    }
}
