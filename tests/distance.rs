use snapbox::cmd::{cargo_bin, Command};
use std::collections::{HashMap, HashSet};
use std::path::Path;

pub mod common;
use crate::common::*;

use std::fs::File;
use std::io::{BufRead, BufReader};

#[cfg(test)]

mod tests {
    use super::*;
    // assert function with tolerance for distances
    fn assert_with_tolerance(actual: f64, expected: f64) {
        let rounded_actual = (actual * 1000.0).round() / 1000.0;
        let rounded_expected = (expected * 1000.0).round() / 1000.0;

        let abs_tolerance = 0.05;
        let diff = (rounded_actual - rounded_expected).abs();

        assert!(
            diff <= abs_tolerance,
            "Absolute difference exceeds tolerance of {}. Actual: {}, Expected: {}",
            abs_tolerance,
            rounded_actual,
            rounded_expected
        );
    }

    fn read_expected_distances(
        true_output: &str,
        sketchlib_true_dict: &mut HashMap<String, Vec<f64>>,
    ) {
        let sandbox = TestSetup::setup();

        let file = File::open(sandbox.file_string(true_output, TestDir::Correct))
            .expect("Failed to open file");
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            let parts: Vec<&str> = line.splitn(2, ": ").collect();
            if parts.len() == 2 {
                let key = parts[0].to_string();
                let value_str = parts[1].trim();

                if value_str.starts_with('[') && value_str.ends_with(']') {
                    // Handle list case
                    let values: Vec<f64> = value_str[1..value_str.len() - 1]
                        .split(',')
                        .map(|s| s.trim().parse().expect("Failed to parse float in list"))
                        .collect();
                    sketchlib_true_dict.insert(key, values);
                } else {
                    // Handle single value case
                    let value = value_str.parse().expect("Failed to parse float");
                    sketchlib_true_dict.insert(key, vec![value]);
                }
            }
        }
    }

    fn read_calculated_distances(full_file_name: &str) -> f64 {
        //Read in the rust sketchlib results from file
        BufReader::new(File::open(full_file_name).expect("Failed to open file"))
            .lines()
            .next()
            .expect("File is empty")
            .expect("Failed to read line")
            .split_whitespace()
            .last()
            .expect("Line has incorrect format")
            .parse()
            .expect("Failed to parse number")
    }

    #[test]
    fn dense_distances() {
        let sandbox = TestSetup::setup();

        //Test 1: One short sequence vs one short seqeunce (1 SNP apart)

        //Test 1 begin -------------
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "3"])
            .arg(sandbox.file_string("short_sequence.fa", TestDir::Input))
            .arg("-v")
            .args(&["-o", "test1_part1"])
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test1_part1.skd"));
        assert_eq!(true, sandbox.file_exists("test1_part1.skm"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "3"])
            .arg(sandbox.file_string("short_sequence_SNP.fa", TestDir::Input))
            .arg("-v")
            .args(&["-o", "test1_part2"])
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test1_part2.skd"));
        assert_eq!(true, sandbox.file_exists("test1_part2.skm"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test1_part1")
            .arg("test1_part2")
            .args(&["-k", "3"])
            .args(&["-o", "short_sequence_dist_3"])
            .arg("-v")
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("short_sequence_dist_3"));

        //Test 2 begin -------------
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "test2_part1"])
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test2_part1.skd"));
        assert_eq!(true, sandbox.file_exists("test2_part1.skm"));

        // removed second to last contigs which is 3610 bps
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .arg(sandbox.file_string(
                "14412_3#82.contigs_velvet_removed_block.fa.gz",
                TestDir::Input,
            ))
            .arg("-v")
            .args(&["-o", "test2_part2"])
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test2_part2.skd"));
        assert_eq!(true, sandbox.file_exists("test2_part2.skm"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test2_part1")
            .arg("test2_part2")
            .args(&["-k", "17"])
            .args(&["-o", "test2_rust_results"])
            .arg("-v")
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test2_rust_results"));

        //Test 3 begin -------------
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "31"])
            .args(&["-s", "10000"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg("-v")
            .args(&["-o", "test3_part1"])
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test3_part1.skd"));
        assert_eq!(true, sandbox.file_exists("test3_part1.skm"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test3_part1")
            .args(&["-k", "31"])
            .args(&["-o", "test3_rust_results"])
            .arg("-v")
            .assert()
            .success();
        assert_eq!(true, sandbox.file_exists("test3_rust_results"));

        // -------------------------------------------------------------------------------------------------------------------------------------------
        // Read in sketchlib C++ true distance results
        let file_path = sandbox.file_string("sketchlib_output_true.txt", TestDir::Correct);
        let _sketchlib_file = match File::open(Path::new(&file_path)) {
            Ok(file) => file,
            Err(error) => {
                eprintln!("Failed to open file: {}", file_path);
                eprintln!("Error: {}", error);
                panic!("Test failed due to file open error");
            }
        };

        let mut sketchlib_true_dict: HashMap<String, Vec<f64>> = HashMap::new();
        read_expected_distances("sketchlib_output_true.txt", &mut sketchlib_true_dict);

        let sketchlib_dist_3 = sketchlib_true_dict
            .get("short_sequence_jaccard_dists_3")
            .expect("Key not found");

        let rust_dist_3 = read_calculated_distances(
            &sandbox.file_string("short_sequence_dist_3", TestDir::Output),
        );

        // Assert or use the value as needed
        let sketchlib_dist_3 = sketchlib_dist_3[0];

        assert_with_tolerance(sketchlib_dist_3, rust_dist_3);

        // TEST 2:
        let rust_whole_genome: f64 = BufReader::new(
            File::open(sandbox.file_string("test2_rust_results", TestDir::Output))
                .expect("Failed to open file"),
        )
        .lines()
        .next()
        .expect("File is empty")
        .expect("Failed to read line")
        .split_whitespace()
        .last()
        .expect("Line has incorrect format")
        .parse()
        .expect("Failed to parse number");

        let whole_genome_block_removed = sketchlib_true_dict
            .get("whole_genome_block_removed")
            .expect("Key not found");

        assert_with_tolerance(rust_whole_genome, whole_genome_block_removed[0]);

        // TEST 3:
        let multiple_genome_rust: Vec<f64> = BufReader::new(
            File::open(sandbox.file_string("test3_rust_results", TestDir::Output))
                .expect("Failed to open file"),
        )
        .lines()
        .map(|line| {
            line.expect("Failed to read line")
                .split_whitespace()
                .last()
                .expect("Line has incorrect format")
                .parse()
                .expect("Failed to parse number")
        })
        .collect();

        let multiple_genome_c = sketchlib_true_dict
            .get("multiple_genomes")
            .expect("Key not found");

        for (_i, (v1, v2)) in multiple_genome_c
            .iter()
            .zip(multiple_genome_rust.iter())
            .enumerate()
        {
            assert_with_tolerance(*v1, *v2);
        }
    }

    #[test]
    fn knn_dists() {
        let sandbox = TestSetup::setup();

        // Move files to test dir
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Sketch the files
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("sketch_db")
            .args(["-v", "--k-seq", "17,31,4", "-s", "10000"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        // C-a dists at knn=1
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--knn")
            .arg("1")
            .assert()
            .stdout_eq(sandbox.snapbox_file("dists_knn_ca.stdout", TestDir::Correct));

        // Jaccard dists at knn=1
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--knn")
            .arg("1")
            .arg("-k")
            .arg("21")
            .assert()
            .stdout_eq(sandbox.snapbox_file("dists_knn_jaccard.stdout", TestDir::Correct));

        // ANI dists at knn=1
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--knn")
            .arg("1")
            .arg("-k")
            .arg("21")
            .arg("--ani")
            .assert()
            .stdout_eq(sandbox.snapbox_file("dists_knn_ani.stdout", TestDir::Correct));
    }

    /// Helper: parse ANI distance output lines into (query, reference, ani) triples.
    fn parse_dist_output(stdout: &str) -> Vec<(String, String, f64)> {
        stdout
            .lines()
            .filter(|l| !l.is_empty())
            .map(|line| {
                let parts: Vec<&str> = line.split_whitespace().collect();
                assert!(parts.len() >= 3, "Unexpected dist output line: {line}");
                (
                    parts[0].to_string(),
                    parts[1].to_string(),
                    parts[2].parse::<f64>().expect("Could not parse ANI"),
                )
            })
            .collect()
    }

    /// Sketch databases into the sandbox:
    /// - `bact_db`: 14412_3#82 + 14412_3#84 (disjoint from query genomes)
    /// - `query_db`: R6 + TIGR4
    /// - `ref_db`: all 4 genomes (used for self-kNN consistency test)
    fn sketch_ref_and_query(sandbox: &TestSetup) {
        for f in &[
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
            "rfile.txt",
            "rfile_ref.txt",
            "qfile.txt",
        ] {
            sandbox.copy_input_file_to_wd(f);
        }

        // Disjoint reference: 14412_3#82 + 14412_3#84 only
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["sketch", "-f", "rfile_ref.txt", "--k-seq", "17,31,4", "-s", "1000", "-o", "bact_db"])
            .assert()
            .success();

        // Query: R6 + TIGR4
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["sketch", "-f", "qfile.txt", "--k-seq", "17,31,4", "-s", "1000", "-o", "query_db"])
            .assert()
            .success();

        // All 4 genomes — used for self-kNN consistency test only
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["sketch", "-f", "rfile.txt", "--k-seq", "17,31,4", "-s", "1000", "-o", "ref_db"])
            .assert()
            .success();
    }

    /// Test 1: output has exactly nq × knn rows.
    ///
    /// 2 query genomes × knn=2 → 4 rows.
    #[test]
    fn knn_cross_query_row_count() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        let output = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run sketchlib");

        let stdout = String::from_utf8(output.stdout).unwrap();
        let n_lines = stdout.lines().filter(|l| !l.is_empty()).count();
        assert_eq!(n_lines, 2, "Expected 2 queries × 1 neighbour = 2 rows, got {n_lines}");
    }

    /// Test 2: kNN output contains the same top-k neighbours as the dense output sorted by ANI.
    ///
    /// Runs both dense and kNN cross-query, then verifies that for each query genome
    /// the kNN output matches the top-2 hits from the dense output ranked by ANI.
    #[test]
    fn knn_cross_query_matches_dense_top_k() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        // Dense: all bact_ref × query pairs (format: ref \t query \t ani)
        let dense_out = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run dense dist");
        let dense_stdout = String::from_utf8(dense_out.stdout).unwrap();

        // Dense output columns: ref(0) \t query(1) \t ani(2) — group by query genome
        let mut dense_by_query: HashMap<String, Vec<(String, f64)>> = HashMap::new();
        for line in dense_stdout.lines().filter(|l| !l.is_empty()) {
            let parts: Vec<&str> = line.split_whitespace().collect();
            assert!(parts.len() >= 3, "Unexpected dense output line: {line}");
            let reference = parts[0].to_string();
            let query = parts[1].to_string();
            let ani: f64 = parts[2].parse().expect("Could not parse ANI");
            dense_by_query.entry(query).or_default().push((reference, ani));
        }

        // For each query, sort by ANI descending and keep top-1
        let mut dense_top1: HashMap<String, (String, f64)> = HashMap::new();
        for (query, mut hits) in dense_by_query {
            hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            dense_top1.insert(query, hits.into_iter().next().unwrap());
        }

        // kNN output columns: query(0) \t ref(1) \t ani(2)
        let knn_out = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run kNN dist");
        let knn_triples = parse_dist_output(&String::from_utf8(knn_out.stdout).unwrap());

        // Assert kNN top-1 matches dense top-1 for every query genome
        for (query, reference, knn_ani) in &knn_triples {
            let (dense_ref, dense_ani) = dense_top1.get(query)
                .unwrap_or_else(|| panic!("Query {query} not found in dense output"));
            assert_eq!(reference, dense_ref,
                "Top neighbour mismatch for {query}: knn={reference}, dense={dense_ref}");
            assert!((knn_ani - dense_ani).abs() < 1e-5,
                "ANI mismatch for {query}/{reference}: knn={knn_ani}, dense={dense_ani}");
        }
    }

    /// Test 3: every row name is a query genome and every column name is a reference genome.
    ///
    /// Verifies that the Display impl uses row_names (query) for rows and ref_names
    /// (reference) for the neighbour column — not ref_names for both.
    #[test]
    fn knn_cross_query_correct_name_columns() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        let output = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run sketchlib");

        let stdout = String::from_utf8(output.stdout).unwrap();
        let triples = parse_dist_output(&stdout);

        let query_names: HashSet<&str> = ["R6.fa.gz", "TIGR4.fa.gz"].iter().cloned().collect();
        let ref_names: HashSet<&str> = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
        ].iter().cloned().collect();

        for (query, reference, _) in &triples {
            assert!(query_names.contains(query.as_str()),
                "Row '{query}' is not a query genome");
            assert!(ref_names.contains(reference.as_str()),
                "Column '{reference}' is not a reference genome");
        }
    }

    /// Test 4: cross-query kNN distances are consistent with self-query kNN.
    ///
    /// Runs self-kNN on all 4 genomes. For R6 and TIGR4, extracts their nearest
    /// neighbours from the self-kNN output restricted to the 2 bacterial reference
    /// genomes (14412_3#82, 14412_3#84). Then runs cross-query kNN with those 2
    /// genomes as reference and R6/TIGR4 as query. Both should give the same distances.
    #[test]
    fn knn_cross_query_consistent_with_self_knn() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        // Self-kNN on all 4 genomes with knn=3
        let self_out = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "ref_db", "--knn", "3", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run self kNN dist");
        let self_triples = parse_dist_output(&String::from_utf8(self_out.stdout).unwrap());

        // Cross-query kNN=1: query=R6+TIGR4 against ref=14412_3#82+14412_3#84
        // kNN output columns: query(0) \t ref(1) \t ani(2)
        let cross_out = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run cross-query kNN dist");
        let cross_triples = parse_dist_output(&String::from_utf8(cross_out.stdout).unwrap());

        // From self-kNN output, extract the best bacterial hit for R6 and TIGR4.
        // Self-kNN output columns: query(0) \t neighbour(1) \t ani(2)
        let bact_names: HashSet<&str> = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
        ].iter().cloned().collect();

        let mut self_best_bact: HashMap<String, (String, f64)> = HashMap::new();
        for query in &["R6.fa.gz", "TIGR4.fa.gz"] {
            let best = self_triples.iter()
                .filter(|(q, r, _)| q == query && bact_names.contains(r.as_str()))
                .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or_else(|| panic!("No bacterial hits for {query} in self-kNN"));
            self_best_bact.insert(query.to_string(), (best.1.clone(), best.2));
        }

        // Compare: cross-query kNN top-1 should match self-kNN top-1 from bacterial genomes
        for (query, cross_ref, cross_ani) in &cross_triples {
            let (self_ref, self_ani) = self_best_bact.get(query)
                .unwrap_or_else(|| panic!("Query {query} not found in self-kNN"));
            assert_eq!(cross_ref, self_ref,
                "Nearest bacterial genome mismatch for {query}: cross={cross_ref}, self={self_ref}");
            assert!((cross_ani - self_ani).abs() < 1e-5,
                "ANI mismatch for {query}/{cross_ref}: cross={cross_ani}, self={self_ani}");
        }
    }

    /// Test 5: cross-query kNN with ref and query completeness files.
    ///
    /// Runs cross-query kNN twice — without and with completeness correction.
    /// Verifies the command succeeds with both completeness flags, output has
    /// the correct row count, ANI values are in [0.0, 1.0], and that out-of-range
    /// completeness values (percentages instead of fractions) are rejected.
    #[test]
    fn knn_cross_query_completeness() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        TestSetup::create_completeness_file(
            &sandbox,
            "ref_completeness.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.8),
                ("14412_3#84.contigs_velvet.fa.gz", 0.85),
            ],
        );
        TestSetup::create_completeness_file(
            &sandbox,
            "query_completeness.txt",
            &[("R6.fa.gz", 0.9), ("TIGR4.fa.gz", 0.75)],
        );

        // Both completeness flags accepted; output has correct row count and valid ANI range
        let comp_out = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db",
                "--knn", "1", "-k", "21", "--ani",
                "--ref-completeness-file", "ref_completeness.txt",
                "--query-completeness-file", "query_completeness.txt",
            ])
            .output()
            .expect("Failed to run with completeness");
        assert!(
            comp_out.status.success(),
            "Command failed: {}",
            String::from_utf8_lossy(&comp_out.stderr)
        );
        let comp_triples = parse_dist_output(&String::from_utf8(comp_out.stdout).unwrap());
        assert_eq!(comp_triples.len(), 2, "Expected 2 rows (2 queries × knn=1)");
        for (query, ref_name, ani) in &comp_triples {
            assert!(
                (0.0..=1.0).contains(ani),
                "ANI out of range for query={query} ref={ref_name}: {ani}"
            );
        }

        // Out-of-range completeness values (percentages) must be rejected
        TestSetup::create_completeness_file(
            &sandbox,
            "bad_completeness.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 80.0),
                ("14412_3#84.contigs_velvet.fa.gz", 85.0),
            ],
        );
        let bad_out = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db",
                "--knn", "1", "-k", "21", "--ani",
                "--ref-completeness-file", "bad_completeness.txt",
            ])
            .output()
            .expect("Failed to run with bad completeness");
        assert!(
            !bad_out.status.success(),
            "Expected failure for out-of-range completeness values"
        );
        let stderr = String::from_utf8_lossy(&bad_out.stderr);
        assert!(
            stderr.contains("[0.0, 1.0]"),
            "Error message should mention [0.0, 1.0] range, got: {stderr}"
        );
    }

    /// Test 6: cross-query kNN in CoreAcc mode (no -k flag).
    ///
    /// Verifies that cross-query kNN works without a k-mer length, producing
    /// 4-column output (query, ref, core, acc).
    #[test]
    fn knn_cross_query_core_acc() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        let output = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "--knn", "1"])
            .output()
            .expect("Failed to run CoreAcc cross-query kNN");

        assert!(
            output.status.success(),
            "CoreAcc cross-query failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8(output.stdout).unwrap();
        let lines: Vec<&str> = stdout.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(lines.len(), 2, "Expected 2 rows (2 queries × knn=1), got {}", lines.len());

        for line in &lines {
            let parts: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(parts.len(), 4, "Expected 4 columns (query, ref, core, acc): {line}");
            parts[2].parse::<f64>().expect("Core distance not a float");
            parts[3].parse::<f64>().expect("Acc distance not a float");
        }
    }

    /// Test 7: cross-query kNN with knn equal to the number of reference genomes.
    ///
    /// bact_db has n=2 reference genomes. With knn=2 every query genome should
    /// get all 2 reference genomes as neighbours (2 queries × 2 = 4 rows).
    /// Previously Bug 2 silently clamped knn=n to knn=n-1, giving only 2 rows.
    #[test]
    fn knn_cross_query_knn_equals_n_ref() {
        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        let output = std::process::Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "bact_db", "query_db", "--knn", "2", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run knn=n_ref cross-query");

        assert!(
            output.status.success(),
            "Command failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8(output.stdout).unwrap();
        let n_lines = stdout.lines().filter(|l| !l.is_empty()).count();
        assert_eq!(
            n_lines, 4,
            "Expected 2 queries × 2 neighbours = 4 rows, got {n_lines}"
        );
    }

    #[test]
    fn subset_dists() {
        let sandbox = TestSetup::setup();

        // Move files to test dir
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Sketch the files
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("sketch_db")
            .args(["-v", "--k-seq", "17,31,4", "-s", "10000"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        // Subset three samples and calc dists
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--subset")
            .arg(sandbox.file_string("subset.txt", TestDir::Input))
            .assert()
            .stdout_eq(sandbox.snapbox_file("dists_subset.stdout", TestDir::Correct));
    }
}
