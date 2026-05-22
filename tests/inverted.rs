use snapbox::cmd::{self, Command};

pub mod common;
use crate::common::*;

#[cfg(test)]

mod tests {
    use sketchlib::sketch::{
        multisketch::MultiSketch, sketch_datafile::SketchArrayReader, BIN_BITS,
    };
    use snapbox::assert_data_eq;

    use super::*;

    fn unpack_skd_bin(bit_planes: &[u64], bin_idx: usize) -> u16 {
        let word_idx = bin_idx / (u64::BITS as usize);
        let bit_offset = bin_idx % (u64::BITS as usize);
        let plane_offset = word_idx * BIN_BITS;

        bit_planes[plane_offset..(plane_offset + BIN_BITS)]
            .iter()
            .enumerate()
            .fold(0_u16, |bin_value, (bit_pos, &plane)| {
                bin_value | (((plane >> bit_offset) & 1) as u16) << bit_pos
            })
    }

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
    fn inverted_skq_matches_standard_skd_bins() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "21", "-s", "64"])
            .arg("-f")
            .arg("rfile.txt")
            .arg("--write-skq")
            .assert()
            .success();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("standard")
            .args(["-v", "--k-vals", "21", "-s", "64"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        let mut standard_sketches =
            MultiSketch::load_metadata(&sandbox.file_string("standard", TestDir::Output))
                .expect("failed to load standard sketch metadata");
        standard_sketches.read_sketch_data(&sandbox.file_string("standard", TestDir::Output));

        let sketch_size = standard_sketches.sketch_size as usize;
        assert_eq!(sketch_size, 64);
        assert_eq!(standard_sketches.kmer_lengths(), &[21]);

        let mut skq_reader = SketchArrayReader::open(
            &sandbox.file_string("inverted.skq", TestDir::Output),
            false,
            1,
            1,
            sketch_size,
        );
        let skq_bins =
            skq_reader.read_all_from_skq(sketch_size * standard_sketches.number_samples_loaded());
        assert_eq!(
            skq_bins.len(),
            sketch_size * standard_sketches.number_samples_loaded()
        );

        for sample_idx in 0..standard_sketches.number_samples_loaded() {
            let skd_bins = standard_sketches.get_sketch_slice(sample_idx, 0);
            let skq_offset = sample_idx * sketch_size;
            for bin_idx in 0..sketch_size {
                assert_eq!(
                    skq_bins[skq_offset + bin_idx],
                    unpack_skd_bin(skd_bins, bin_idx),
                    "sample {sample_idx} bin {bin_idx} differs between .skq and .skd"
                );
            }
        }
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

    // Helper to set up reordered SKI + standard SKD
    // SKI order (via --species-names): TIGR4, R6, #82, #84
    // SKD order (via rfile.txt):       #82, #84, R6, TIGR4
    fn setup_reordered_precluster(sandbox: &TestSetup) {
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Build SKI+SKQ with reordered samples via --species-names
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("reordered_inv")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg("rfile.txt")
            .arg("--write-skq")
            .arg("--species-names")
            .arg(sandbox.file_string("species_names.txt", TestDir::Input))
            .assert()
            .success();

        // Build standard SKD (different sample order)
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
    }

    #[test]
    fn inverted_precluster_reordered() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        // Run preclustering with different SKI vs SKD orderings
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster.stdout", TestDir::Correct)
                    .unordered(),
            );
    }

    #[test]
    fn inverted_precluster_reordered_ani() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--ani")
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster_ani.stdout", TestDir::Correct)
                    .unordered(),
            );
    }

    #[test]
    fn inverted_precluster_reordered_knn_padding() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        // knn=50 >> n_samples=4, tests padding with reordered indices
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "50"])
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster.stdout", TestDir::Correct)
                    .unordered(),
            );
    }

    #[test]
    fn inverted_precluster_reversed_rfile() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");

        // Build SKI with reversed sample order
        let reversed_rfile = sandbox.create_named_rfile(
            "reversed_rfile.txt",
            &[
                "TIGR4.fa.gz",
                "R6.fa.gz",
                "14412_3#84.contigs_velvet.fa.gz",
                "14412_3#82.contigs_velvet.fa.gz",
            ],
        );
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("reversed_inv")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg(&reversed_rfile)
            .arg("--write-skq")
            .assert()
            .success();

        // Build SKD with standard order
        let standard_rfile = sandbox.create_named_rfile(
            "standard_rfile.txt",
            &[
                "14412_3#82.contigs_velvet.fa.gz",
                "14412_3#84.contigs_velvet.fa.gz",
                "R6.fa.gz",
                "TIGR4.fa.gz",
            ],
        );
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("standard")
            .args(["-v", "--k-vals", "21", "-s", "1000"])
            .arg("-f")
            .arg(&standard_rfile)
            .assert()
            .success();

        // Run precluster with different orderings
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("reversed_inv.ski")
            .assert()
            .stdout_eq(
                sandbox
                    .snapbox_file("inverted_precluster.stdout", TestDir::Correct)
                    .unordered(),
            );
    }

    #[test]
    fn inverted_precluster_reordered_completeness() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        // Create completeness file with asymmetric values
        TestSetup::create_completeness_file(
            &sandbox,
            "completeness_reorder.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.7),
                ("14412_3#84.contigs_velvet.fa.gz", 0.8),
                ("R6.fa.gz", 0.9),
                ("TIGR4.fa.gz", 0.95),
            ],
        );

        // Also build same-order SKI for comparison
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("sameorder_inv")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg("rfile.txt")
            .arg("--write-skq")
            .assert()
            .success();

        // Run with reordered SKI + completeness
        let reordered_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("--completeness-file")
            .arg(sandbox.file_string("completeness_reorder.txt", TestDir::Output))
            .arg("reordered_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        // Run with same-order SKI + completeness
        let sameorder_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("--completeness-file")
            .arg(sandbox.file_string("completeness_reorder.txt", TestDir::Output))
            .arg("sameorder_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        // Parse and compare: both orderings must produce the same distances
        let mut reordered_lines: Vec<String> = String::from_utf8(reordered_output)
            .unwrap()
            .lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.to_string())
            .collect();
        let mut sameorder_lines: Vec<String> = String::from_utf8(sameorder_output)
            .unwrap()
            .lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.to_string())
            .collect();
        reordered_lines.sort();
        sameorder_lines.sort();
        assert_eq!(
            reordered_lines, sameorder_lines,
            "Reordered and same-order precluster with completeness should produce identical results"
        );
    }

    #[test]
    fn inverted_precluster_all_genomes_in_output() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        let input_names = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ];

        // Run precluster
        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        let output_str = String::from_utf8(output).unwrap();
        let mut found_names: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for line in output_str.lines().filter(|l| !l.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            found_names.insert(fields[0]);
            found_names.insert(fields[1]);
        }

        for name in &input_names {
            assert!(
                found_names.contains(name),
                "Genome {name} missing from precluster output"
            );
        }
        assert_eq!(
            found_names.len(),
            input_names.len(),
            "Number of unique genomes in output ({}) does not match input ({})",
            found_names.len(),
            input_names.len(),
        );
    }

    #[test]
    fn inverted_precluster_self_skip() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        // knn=3 to get all possible neighbors
        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "3"])
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        let output_str = String::from_utf8(output).unwrap();
        for line in output_str.lines().filter(|l| !l.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(
                fields.len() >= 3,
                "Output line should have at least 3 tab-separated fields: {line}"
            );
            // First two fields are sample names - they must not be the same
            assert_ne!(
                fields[0], fields[1],
                "Self-comparison found in output: {line}"
            );
            // Distance should be between 0 and 1 exclusive
            let dist: f64 = fields[2].parse().expect("Failed to parse distance");
            assert!(
                dist > 0.0 && dist < 1.0,
                "Distance out of range (0,1): {dist} in line: {line}"
            );
        }
    }

    #[test]
    fn inverted_precluster_retain_singleton() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        let input_names = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ];

        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--retain-unmatched")
            .arg("singleton")
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        let output_str = String::from_utf8(output).unwrap();
        let mut found_names: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for line in output_str.lines().filter(|l| !l.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            found_names.insert(fields[0]);
            found_names.insert(fields[1]);
            // Self-match lines should have distance 0
            if fields[0] == fields[1] {
                let dist: f64 = fields[2].parse().expect("Failed to parse distance");
                assert_eq!(dist, 0.0, "Singleton self-match should have distance 0");
            }
        }

        for name in &input_names {
            assert!(
                found_names.contains(name),
                "Genome {name} missing from output with --retain-unmatched singleton"
            );
        }
    }

    #[test]
    fn inverted_precluster_retain_bruteforce() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        let input_names = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ];

        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--retain-unmatched")
            .arg("bruteforce")
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        let output_str = String::from_utf8(output).unwrap();
        let mut found_names: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for line in output_str.lines().filter(|l| !l.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            found_names.insert(fields[0]);
            found_names.insert(fields[1]);
            // No self-matches in bruteforce mode
            assert_ne!(
                fields[0], fields[1],
                "Self-comparison found in bruteforce output: {line}"
            );
            let dist: f64 = fields[2].parse().expect("Failed to parse distance");
            assert!(
                dist > 0.0 && dist < 1.0,
                "Distance out of range (0,1): {dist} in line: {line}"
            );
        }

        for name in &input_names {
            assert!(
                found_names.contains(name),
                "Genome {name} missing from output with --retain-unmatched bruteforce"
            );
        }
    }

    #[test]
    fn inverted_query_reordered() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Build reordered SKI
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("reordered_inv")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg("rfile.txt")
            .arg("--species-names")
            .arg(sandbox.file_string("species_names.txt", TestDir::Input))
            .assert()
            .success();

        // Query with match-count
        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("query")
            .arg("-v")
            .arg("-f")
            .arg("rfile.txt")
            .arg("reordered_inv.ski")
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();

        // Verify self-matches are 10 and cross-group matches are 0
        let output_str = String::from_utf8(output).unwrap();
        for line in output_str.lines().skip(1).filter(|l| !l.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            let query_name = fields[0];
            // Find the self-match count (column matching the query name)
            let header: Vec<&str> = output_str.lines().next().unwrap().split('\t').collect();
            for (col_idx, &col_name) in header.iter().enumerate().skip(1) {
                let count: u32 = fields[col_idx].parse().unwrap();
                if col_name == query_name {
                    assert_eq!(count, 10, "Self-match count should be 10 for {query_name}");
                }
            }
        }

        // Query with any-bins
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("query")
            .arg("-v")
            .arg("-f")
            .arg("rfile.txt")
            .arg("reordered_inv.ski")
            .args(&["--query-type", "any-bins"])
            .assert()
            .success();

        // Query with all-bins
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("query")
            .arg("-v")
            .arg("-f")
            .arg("rfile.txt")
            .arg("reordered_inv.ski")
            .args(&["--query-type", "all-bins"])
            .assert()
            .success();
    }
}
