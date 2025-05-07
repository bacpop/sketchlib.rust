use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn sketch_fasta() {
        let sandbox = TestSetup::setup();

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
            .stdout_eq(sandbox.snapbox_file("assembly_sketch_info.stdout", TestDir::Correct));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("assembly")
            .arg("-v")
            .assert()
            .stdout_eq(sandbox.snapbox_file("assembly_sketch_full_info.stdout", TestDir::Correct));
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
            .stdout_eq(sandbox.snapbox_file("read_sketch_info.stdout", TestDir::Correct));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("reads")
            .arg("-v")
            .assert()
            .stdout_eq(sandbox.snapbox_file("read_sketch_full_info.stdout", TestDir::Correct));
    }

    #[test]
    /// Check that older databases can still be read correctly
    /// (some fields added on .skm in v0.2.0)
    fn legacy_databases() {
        use sketchlib::sketch::multisketch::MultiSketch;
        let sandbox = TestSetup::setup();
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");

        // Old command:
        // sketchlib sketch -v -o legacy_db --k-vals 17,21,25 -s 100 R6.fa.gz TIGR4.fa.gz

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("new_db")
            .args(["-v", "--k-vals", "17,21,25", "-s", "100"])
            .arg("R6.fa.gz")
            .arg("TIGR4.fa.gz")
            .assert()
            .success();

        // Check .skm the same
        let new_sketch: MultiSketch =
            MultiSketch::load_metadata(&sandbox.file_string("new_db", TestDir::Output))
                .expect("Failed to load new sketch");
        let expected_sketch =
            MultiSketch::load_metadata(&sandbox.file_string("legacy_db", TestDir::Input))
                .expect("Failed to load legacy sketch");

        assert_eq!(
            new_sketch, expected_sketch,
            "Sketch metadata does not match legacy and new version"
        );
    }
}
