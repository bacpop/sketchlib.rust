use predicates::prelude::*;
use snapbox::cmd::{cargo_bin, Command};
use std::path::Path;

pub mod common;
use crate::common::*;

#[path = "../src/io.rs"]
pub mod io;

use sketchlib::multisketch::MultiSketch;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn concat_competibility_test() {
        let sandbox = TestSetup::setup();
        let ref_db1 = sandbox.file_string("sketches1", TestDir::Input);
        let ref_db2 = sandbox.file_string("sketches2", TestDir::Input);
        let file_list_name = sandbox.file_string("fasta.txt", TestDir::Input);
        let non_seq_files: Option<Vec<String>> = None;


        let file_list: Option<String> = Some(file_list_name);
        //get input files with file_list
        log::info!("Getting input files for test");
        let input_file_list = io::get_input_list(&file_list, &non_seq_files);
        log::info!("Parsed {} samples in input list", input_file_list.len());

        //check if any of the new files are already existant in the db
        let db1_metadata: MultiSketch = MultiSketch::load(&ref_db1)
            .unwrap_or_else(|_| panic!("Could not read sketch metadata from {}.skm", ref_db1));
        println!("{:?}", db1_metadata);

        //check if any of the new files are already existant in the db
        let db2_metadata: MultiSketch = MultiSketch::load(&ref_db2)
            .unwrap_or_else(|_| panic!("Could not read sketch metadata from {}.skm", ref_db2));
        println!("{:?}", db2_metadata);

        // Test case 1:
        assert!(
            !db1_metadata.concat_competibility(&input_file_list),
            "Sketches should not be compatible"
        );
    }

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
        log::info!("Part 1 Sketched");

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
        log::info!("Part 2 Sketched");

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
        log::info!("concat_ref Sketched");

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
            .arg("concat")
            .arg("part1")
            .arg("--seq-files")
            // .arg(sandbox.file_string("fasta_part2.txt", TestDir::Input))
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "concat_test"])
            .assert()
            .success();
        log::info!("concat_test Sketched");

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
