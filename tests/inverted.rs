use snapbox::cmd::{self, Command};

pub mod common;
use crate::common::*;

#[cfg(test)]

mod tests {
    use snapbox::assert_data_eq;

    use super::*;

    #[test]
    fn inverted_sketch() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("inverted.ski"));

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("inverted.ski")
            .arg("-v")
            .assert()
            .stdout_eq(sandbox.snapbox_file("inverted_sketch_info.stdout", TestDir::Correct));

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("inverted.ski")
            .arg("-v")
            .assert()
            .stdout_eq(sandbox.snapbox_file("inverted_sketch_full_info.stdout", TestDir::Correct));
    }

    #[test]
    fn inverted_sketch_with_metadata() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");
        sandbox.copy_input_file_to_wd("metadata.txt");
        sandbox.copy_input_file_to_wd("species_names.txt");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31"])
            .args(["-f", "rfile.txt"])
            .args(["--species-names", "species_names.txt"])
            .args(["--metadata", "metadata.txt"])
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("inverted.ski"));
    }

    #[test]
    fn inverted_sketch_errors() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");
        sandbox.copy_input_file_to_wd("species_names_empty.txt");

        // This shouldn't give an error, but it is an unfortunate situation with respect to the species-names file the user provided
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31"])
            .args(["-f", "rfile.txt"])
            .args(["--species-names", "species_names_empty.txt"])
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("inverted.ski"));
    }

    #[test]
    fn inverted_sketch_multifile() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile_multi.txt");

        // This shouldn't give an error, but it is an unfortunate situation with respect to the species-names file the user provided
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31", "-s", "1000"])
            .args(["-f", "rfile_multi.txt"])
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("inverted.ski"));
    }

    #[test]
    fn inverted_reorder() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "61", "-s", "63"])
            .arg("--species-names")
            .arg(sandbox.file_string("species_names.txt", TestDir::Input))
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .arg("R6.fa.gz")
            .arg("TIGR4.fa.gz")
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("inverted.ski"));

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("inverted.ski")
            .arg("--sample-info")
            .arg("-v")
            .assert()
            .stdout_eq(
                sandbox.snapbox_file("inverted_sketch_info_reorder.stdout", TestDir::Correct),
            );
    }

    #[test]
    fn inverted_query() {
        let sandbox = TestSetup::setup();

        // Move files to test dir
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Build an inverted index
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        // Query with the same samples, default is 'match-count'
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("query")
            .arg("-v")
            .arg("-f")
            .arg("rfile.txt")
            .arg("inverted.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_query_count.stdout", TestDir::Correct)
                    .unordered(),
            );

        // Any bin matching gives two pairs
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("query")
            .arg("-v")
            .arg("-f")
            .arg("rfile.txt")
            .arg("inverted.ski")
            .args(&["--query-type", "any-bins"])
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_query_any.stdout", TestDir::Correct)
                    .unordered(),
            );

        // All matching gives self only
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("query")
            .arg("-v")
            .arg("-f")
            .arg("rfile.txt")
            .arg("inverted.ski")
            .args(&["--query-type", "all-bins"])
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_query_all.stdout", TestDir::Correct)
                    .unordered(),
            );
    }

    #[test]
    fn inverted_precluster() {
        let sandbox = TestSetup::setup();

        // Move files to test dir
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Build an inverted index .ski and .skq
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg("rfile.txt")
            .arg("--write-skq")
            .assert()
            .success();

        // Test that the written .skq is correct (order should be right)
        assert_data_eq!(
            sandbox.snapbox_file("inverted.skq", TestDir::Output),
            sandbox.snapbox_file("inverted.skq", TestDir::Correct)
        );

        // See the match-count results in inverted_query_count.stdout
        // 14412_3#82.contigs_velvet.fa.gz and 14412_3#84.contigs_velvet.fa.gz match
        // R6.fa.gz and TIGR4.fa.gz match
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .arg("-v")
            .arg("--count")
            .arg("inverted.ski")
            .assert()
            .stdout_eq("Identified 2 prefilter pairs from a max of 6\n");

        // Build a standard .skd
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("standard")
            .args(["-v", "--k-vals", "21", "-s", "1000"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        // Run preclustering mode
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("inverted.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster.stdout", TestDir::Correct)
                    .unordered(),
            );

        // Same with ANI
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--ani")
            .arg("--skd")
            .arg("standard")
            .arg("inverted.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster_ani.stdout", TestDir::Correct)
                    .unordered(),
            );

        // Default knn, check no junk distances printed
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "50"])
            .arg("--skd")
            .arg("standard")
            .arg("inverted.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster.stdout", TestDir::Correct)
                    .unordered(),
            );
    }
}
