use predicates::prelude::*;
use snapbox::cmd::{cargo_bin, Command};
use std::path::Path;

pub mod common;
use crate::common::*;

use sketchlib::multisketch::MultiSketch;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_delete_sketches() {
        let sandbox = TestSetup::setup();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg("--seq-files")
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "full_db"])
            .assert()
            .success();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg("--seq-files")
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
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
            MultiSketch::load(&sandbox.file_string("deleted_db", TestDir::Output))
                .expect("Failed to load output merged sketch");
        let expected_sketch = MultiSketch::load(&sandbox.file_string("result_db", TestDir::Output))
            .expect("Failed to load expected merged sketch");
        assert_eq!(
            merged_sketch, expected_sketch,
            "Merged sketch metadata does not match"
        );
    }
}

// use predicates::prelude::*;
// use snapbox::cmd::{cargo_bin, Command};
// use std::path::Path;

// pub mod common;
// use crate::common::*;
// use sketchlib::io::*;

// use sketchlib::multisketch::MultiSketch;

// #[cfg(test)]

// mod tests {
//     use super::*;
//     use sketchlib::io;

//     #[test]
//     fn test_delete_sketches() {
//         let sandbox = TestSetup::setup();

//         Command::new(cargo_bin("sketchlib"))
//             .current_dir(sandbox.get_wd())
//             .arg("sketch")
//             .args(&["--k-vals", "17"])
//             .arg("--seq-files")
//             .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
//             .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
//             .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
//             .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
//             .arg("-v")
//             .args(&["-o", "full_db"])
//             .assert()
//             .success();

//         Command::new(cargo_bin("sketchlib"))
//             .current_dir(sandbox.get_wd())
//             .arg("sketch")
//             .args(&["--k-vals", "17"])
//             .arg("--seq-files")
//             .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
//             .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
//             .arg("-v")
//             .args(&["-o", "result_db"])
//             .assert()
//             .success();

//         // samples in delete_test.txt are TIGR4, and R6
//         Command::new(cargo_bin("sketchlib"))
//             .current_dir(sandbox.get_wd())
//             .arg("delete")
//             .arg("full_db")
//             .arg(sandbox.file_string("delete_test.txt", TestDir::Input))
//             .arg("deleted_db")
//             .assert()
//             .success();

//         let predicate_file = predicate::path::eq_file(Path::new(
//             &sandbox.file_string("deleted_db.skd", TestDir::Output),
//         ));
//         assert_eq!(
//             true,
//             predicate_file.eval(Path::new(
//                 &sandbox.file_string("result_db.skd", TestDir::Output)
//             )),
//             "Sketch data after deletion incorrect"
//         );

//         // Check .skm the same
//         let merged_sketch: MultiSketch =
//             MultiSketch::load(&sandbox.file_string("deleted_db", TestDir::Output))
//                 .expect("Failed to load output merged sketch");
//         let expected_sketch = MultiSketch::load(&sandbox.file_string("result_db", TestDir::Output))
//             .expect("Failed to load expected merged sketch");
//         assert_eq!(
//             merged_sketch, expected_sketch,
//             "Metadata data after deletion incorrect"
//         );
//     }
// }
