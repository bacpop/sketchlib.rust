use predicates::prelude::*;
use snapbox::cmd::{cargo_bin, Command};
use std::path::Path;

pub mod common;
use crate::common::*;
use sketchlib::io::*;

use sketchlib::multisketch::MultiSketch;

#[cfg(test)]

mod tests {
    use sketchlib::io;
    use super::*;
   

    #[test]
    fn test_concat_sketches() {
        let sandbox = TestSetup::setup();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg("--seq-files")
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
            .arg("--seq-files")
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
            .arg("--seq-files")
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "concat_ref"])
            .assert()
            .success();

        // Overlapping labels fails
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("concat")
            .arg("part1")
            .arg("--seq-files")
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "concat_test"])
            .assert()
            .failure();

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("append")
            .arg("part1")
            .arg("--seq-files")
            // .arg(sandbox.file_string("fasta_part2.txt", TestDir::Input))
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "concat_test"])
            .assert()
            .success();

        // Check .skm the same
        let concat_sketch: MultiSketch =
            MultiSketch::load(&sandbox.file_string("concat_test", TestDir::Output))
                .expect("Failed to load output merged sketch");
        let expected_sketch =
            MultiSketch::load(&sandbox.file_string("concat_ref", TestDir::Output))
                .expect("Failed to load expected merged sketch");
        println!("{}", concat_sketch);
        println!("{}", expected_sketch);
        assert_eq!(
            concat_sketch, expected_sketch,
            "Concat sketch metadata does not match"
        );

        // Check .skd the same
        let predicate_file = predicate::path::eq_file(Path::new(
            &sandbox.file_string("concat_test.skd", TestDir::Output),
        ));
        assert_eq!(
            true,
            predicate_file.eval(Path::new(
                &sandbox.file_string("concat_ref.skd", TestDir::Output)
            )),
            "Concat sketch data does not match"
        );
    }
}
