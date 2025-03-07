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
    fn sketch_fasta() {
        let sandbox = TestSetup::setup();

        // Create a fasta rfile in the tmp dir
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("assembly")
            .args(["-v", "-k", "31"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("assembly.skm"));
        assert_eq!(true, sandbox.file_exists("assembly.skd"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("assembly")
            .assert()
            .stdout_matches_path(sandbox.file_string("assembly_sketch_info.stdout", TestDir::Correct));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("assembly")
            .arg("-v")
            .assert()
            .stdout_matches_path(sandbox.file_string("assembly_sketch_full_info.stdout", TestDir::Correct));
    }

    #[test]
    fn sketch_fastq() {
        let sandbox = TestSetup::setup();

        // Create a fastq rfile in the tmp dir
        let rfile_name = sandbox.create_fastq_rfile("test");
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-f")
            .arg(rfile_name)
            .arg("-o")
            .arg("reads")
            .args(["--min-count", "2", "-v", "-k", "9", "--min-qual", "2"])
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("reads.skm"));
        assert_eq!(true, sandbox.file_exists("reads.skd"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("reads")
            .assert()
            .stdout_matches_path(sandbox.file_string("read_sketch_info.stdout", TestDir::Correct));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("reads")
            .arg("-v")
            .assert()
            .stdout_matches_path(sandbox.file_string("read_sketch_full_info.stdout", TestDir::Correct));
    }
}
