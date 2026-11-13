use snapbox::cmd::{self, Command};
use std::collections::{HashMap, HashSet};
use std::path::Path;

pub mod common;
use crate::common::*;

use std::fs::File;
use std::io::{BufRead, BufReader};

#[cfg(test)]

mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[cfg(not(feature = "gpu"))]
    #[test]
    fn gpu_request_requires_gpu_build_before_opening_files() {
        let sandbox = TestSetup::setup();
        for query in [None, Some("missing_query")] {
            let mut command = Command::new(cmd::cargo_bin!("sketchlib"))
                .current_dir(sandbox.get_wd())
                .args([
                    "dist",
                    "--gpu",
                    "0,1",
                    "-o",
                    "should-not-exist.tsv",
                    "missing_reference",
                ]);
            if let Some(query) = query {
                command = command.arg(query);
            }
            command.assert().failure().stderr_eq(
                "Error: GPU support is not enabled; rebuild with `cargo oxide build --features gpu`.\n",
            );
            assert!(!sandbox.file_exists("should-not-exist.tsv"));
        }
    }

    // Assert function with tolerance for distances computed against the external
    // C++ sketchlib reference values. Some of these (e.g. very short sequences)
    // are known to diverge more than floating-point noise between implementations,
    // so this stays a fairly loose absolute tolerance rather than a tight epsilon.
    fn assert_with_tolerance(actual: f64, expected: f64) {
        assert_abs_diff_eq!(actual, expected, epsilon = 0.05);
    }

    // Compares tab-separated distance output against a golden file, allowing
    // numeric fields to differ by a small epsilon (SIMD reordering can shift
    // the last few bits of a float) while requiring an exact match on names.
    //
    // Rows are sorted by their (name1, name2) key columns before comparing, since
    // streamed dense output order is not guaranteed to match the golden file's
    // row-major order. Sorting by the name columns (not the whole line) avoids the
    // sort itself being perturbed by the numeric noise this function already
    // tolerates.
    fn assert_dist_stdout_with_tolerance(actual: &str, expected: &snapbox::Data) {
        let expected_str = expected
            .render()
            .expect("Failed to render expected snapshot data");
        let name_key = |line: &&str| {
            let mut fields = line.splitn(3, '\t');
            (
                fields.next().unwrap_or("").to_string(),
                fields.next().unwrap_or("").to_string(),
            )
        };
        let mut actual_lines: Vec<&str> = actual.lines().filter(|l| !l.is_empty()).collect();
        let mut expected_lines: Vec<&str> =
            expected_str.lines().filter(|l| !l.is_empty()).collect();
        actual_lines.sort_by_key(name_key);
        expected_lines.sort_by_key(name_key);
        assert_eq!(
            actual_lines.len(),
            expected_lines.len(),
            "Line count mismatch.\nActual:\n{actual}\nExpected:\n{expected_str}"
        );
        for (actual_line, expected_line) in actual_lines.iter().zip(expected_lines.iter()) {
            let actual_fields: Vec<&str> = actual_line.split('\t').collect();
            let expected_fields: Vec<&str> = expected_line.split('\t').collect();
            assert_eq!(
                actual_fields.len(),
                expected_fields.len(),
                "Field count mismatch. Actual: {actual_line}, Expected: {expected_line}"
            );
            for (actual_field, expected_field) in actual_fields.iter().zip(expected_fields.iter()) {
                match (actual_field.parse::<f64>(), expected_field.parse::<f64>()) {
                    (Ok(actual_val), Ok(expected_val)) => {
                        assert_abs_diff_eq!(actual_val, expected_val, epsilon = 1e-4);
                    }
                    _ => assert_eq!(
                        actual_field, expected_field,
                        "Field mismatch. Actual line: {actual_line}, Expected line: {expected_line}"
                    ),
                }
            }
        }
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

    #[test]
    fn dense_distances() {
        let sandbox = TestSetup::setup();

        //Test 2 begin -------------
        Command::new(cmd::cargo_bin!("sketchlib"))
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
        Command::new(cmd::cargo_bin!("sketchlib"))
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

        Command::new(cmd::cargo_bin!("sketchlib"))
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
        Command::new(cmd::cargo_bin!("sketchlib"))
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

        Command::new(cmd::cargo_bin!("sketchlib"))
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
        Command::new(cmd::cargo_bin!("sketchlib"))
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
        let ca_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--knn")
            .arg("1")
            .output()
            .expect("Failed to run C-a dist");
        assert!(ca_output.status.success());
        assert_dist_stdout_with_tolerance(
            &String::from_utf8(ca_output.stdout).unwrap(),
            &sandbox.snapbox_file("dists_knn_ca.stdout", TestDir::Correct),
        );

        // Jaccard dists at knn=1
        let jaccard_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--knn")
            .arg("1")
            .arg("-k")
            .arg("21")
            .output()
            .expect("Failed to run Jaccard dist");
        assert!(jaccard_output.status.success());
        assert_dist_stdout_with_tolerance(
            &String::from_utf8(jaccard_output.stdout).unwrap(),
            &sandbox.snapbox_file("dists_knn_jaccard.stdout", TestDir::Correct),
        );

        // ANI dists at knn=1
        let ani_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--knn")
            .arg("1")
            .arg("-k")
            .arg("21")
            .arg("--ani")
            .output()
            .expect("Failed to run ANI dist");
        assert!(ani_output.status.success());
        assert_dist_stdout_with_tolerance(
            &String::from_utf8(ani_output.stdout).unwrap(),
            &sandbox.snapbox_file("dists_knn_ani.stdout", TestDir::Correct),
        );
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
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "sketch",
                "-f",
                "rfile_ref.txt",
                "--k-seq",
                "17,31,4",
                "-s",
                "1000",
                "-o",
                "bact_db",
            ])
            .assert()
            .success();

        // Query: R6 + TIGR4
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "sketch",
                "-f",
                "qfile.txt",
                "--k-seq",
                "17,31,4",
                "-s",
                "1000",
                "-o",
                "query_db",
            ])
            .assert()
            .success();

        // All 4 genomes — used for self-kNN consistency test only
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "sketch",
                "-f",
                "rfile.txt",
                "--k-seq",
                "17,31,4",
                "-s",
                "1000",
                "-o",
                "ref_db",
            ])
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

        let output = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani",
            ])
            .output()
            .expect("Failed to run sketchlib");

        let stdout = String::from_utf8(output.stdout).unwrap();
        let n_lines = stdout.lines().filter(|l| !l.is_empty()).count();
        assert_eq!(
            n_lines, 2,
            "Expected 2 queries × 1 neighbour = 2 rows, got {n_lines}"
        );
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
        let dense_out = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
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
            dense_by_query
                .entry(query)
                .or_default()
                .push((reference, ani));
        }

        // For each query, sort by ANI descending and keep top-1
        let mut dense_top1: HashMap<String, (String, f64)> = HashMap::new();
        for (query, mut hits) in dense_by_query {
            hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            dense_top1.insert(query, hits.into_iter().next().unwrap());
        }

        // kNN output columns: query(0) \t ref(1) \t ani(2)
        let knn_out = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani",
            ])
            .output()
            .expect("Failed to run kNN dist");
        let knn_triples = parse_dist_output(&String::from_utf8(knn_out.stdout).unwrap());

        // Assert kNN top-1 matches dense top-1 for every query genome
        for (query, reference, knn_ani) in &knn_triples {
            let (dense_ref, dense_ani) = dense_top1
                .get(query)
                .unwrap_or_else(|| panic!("Query {query} not found in dense output"));
            assert_eq!(
                reference, dense_ref,
                "Top neighbour mismatch for {query}: knn={reference}, dense={dense_ref}"
            );
            assert!(
                (knn_ani - dense_ani).abs() < 1e-4,
                "ANI mismatch for {query}/{reference}: knn={knn_ani}, dense={dense_ani}"
            );
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

        let output = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani",
            ])
            .output()
            .expect("Failed to run sketchlib");

        let stdout = String::from_utf8(output.stdout).unwrap();
        let triples = parse_dist_output(&stdout);

        let query_names: HashSet<&str> = ["R6.fa.gz", "TIGR4.fa.gz"].iter().cloned().collect();
        let ref_names: HashSet<&str> = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
        ]
        .iter()
        .cloned()
        .collect();

        for (query, reference, _) in &triples {
            assert!(
                query_names.contains(query.as_str()),
                "Row '{query}' is not a query genome"
            );
            assert!(
                ref_names.contains(reference.as_str()),
                "Column '{reference}' is not a reference genome"
            );
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
        let self_out = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args(["dist", "ref_db", "--knn", "3", "-k", "21", "--ani"])
            .output()
            .expect("Failed to run self kNN dist");
        let self_triples = parse_dist_output(&String::from_utf8(self_out.stdout).unwrap());

        // Cross-query kNN=1: query=R6+TIGR4 against ref=14412_3#82+14412_3#84
        // kNN output columns: query(0) \t ref(1) \t ani(2)
        let cross_out = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db", "--knn", "1", "-k", "21", "--ani",
            ])
            .output()
            .expect("Failed to run cross-query kNN dist");
        let cross_triples = parse_dist_output(&String::from_utf8(cross_out.stdout).unwrap());

        // From self-kNN output, extract the best bacterial hit for R6 and TIGR4.
        // Self-kNN output columns: query(0) \t neighbour(1) \t ani(2)
        let bact_names: HashSet<&str> = [
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
        ]
        .iter()
        .cloned()
        .collect();

        let mut self_best_bact: HashMap<String, (String, f64)> = HashMap::new();
        for query in &["R6.fa.gz", "TIGR4.fa.gz"] {
            let best = self_triples
                .iter()
                .filter(|(q, r, _)| q == query && bact_names.contains(r.as_str()))
                .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or_else(|| panic!("No bacterial hits for {query} in self-kNN"));
            self_best_bact.insert(query.to_string(), (best.1.clone(), best.2));
        }

        // Compare: cross-query kNN top-1 should match self-kNN top-1 from bacterial genomes
        for (query, cross_ref, cross_ani) in &cross_triples {
            let (self_ref, self_ani) = self_best_bact
                .get(query)
                .unwrap_or_else(|| panic!("Query {query} not found in self-kNN"));
            assert_eq!(
                cross_ref, self_ref,
                "Nearest bacterial genome mismatch for {query}: cross={cross_ref}, self={self_ref}"
            );
            assert!(
                (cross_ani - self_ani).abs() < 1e-4,
                "ANI mismatch for {query}/{cross_ref}: cross={cross_ani}, self={self_ani}"
            );
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
        let comp_out = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist",
                "bact_db",
                "query_db",
                "--knn",
                "1",
                "-k",
                "21",
                "--ani",
                "--ref-completeness-file",
                "ref_completeness.txt",
                "--query-completeness-file",
                "query_completeness.txt",
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
        let bad_out = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist",
                "bact_db",
                "query_db",
                "--knn",
                "1",
                "-k",
                "21",
                "--ani",
                "--ref-completeness-file",
                "bad_completeness.txt",
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

        let output = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
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
        assert_eq!(
            lines.len(),
            2,
            "Expected 2 rows (2 queries × knn=1), got {}",
            lines.len()
        );

        for line in &lines {
            let parts: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(
                parts.len(),
                4,
                "Expected 4 columns (query, ref, core, acc): {line}"
            );
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

        let output = std::process::Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .args([
                "dist", "bact_db", "query_db", "--knn", "2", "-k", "21", "--ani",
            ])
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
        Command::new(cmd::cargo_bin!("sketchlib"))
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
        let subset_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("sketch_db")
            .arg("-v")
            .arg("--subset")
            .arg(sandbox.file_string("subset.txt", TestDir::Input))
            .output()
            .expect("Failed to run subset dist");
        assert!(subset_output.status.success());
        assert_dist_stdout_with_tolerance(
            &String::from_utf8(subset_output.stdout).unwrap(),
            &sandbox.snapbox_file("dists_subset.stdout", TestDir::Correct),
        );
    }

    /// Helper: parse `DistanceMatrix`'s `Display`/TSV output into (primary, accessory) pairs.
    fn parse_display_pairs(display: &str, core_acc: bool) -> Vec<(f64, Option<f64>)> {
        display
            .lines()
            .filter(|l| !l.is_empty())
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                let primary = fields[2].parse().expect("Could not parse distance");
                if core_acc {
                    let accessory = fields[3].parse().expect("Could not parse accessory dist");
                    (primary, Some(accessory))
                } else {
                    (primary, None)
                }
            })
            .collect()
    }

    #[test]
    fn dists_iter_matches_display() {
        use sketchlib::distances::{self_dists_all, set_k};
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("iter_db")
            .args(["-v", "--k-seq", "17,31,4", "-s", "1000"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        let sketches = MultiSketch::load(&sandbox.file_string("iter_db", TestDir::Output))
            .expect("failed to load sketches");
        let n = sketches.number_samples_loaded();

        // CoreAcc matrix (multi-k, no -k given): dists_iter should yield Some(accessory)
        let core_acc_type = set_k(&sketches, None, false).expect("set_k failed");
        let core_acc_matrix = self_dists_all(&sketches, n, core_acc_type, true, None, 0.0);
        let display_pairs = parse_display_pairs(&core_acc_matrix.to_string(), true);
        let iter_pairs: Vec<(f32, Option<f32>)> = core_acc_matrix.dists_iter().collect();
        assert_eq!(display_pairs.len(), iter_pairs.len());
        for ((core, accessory), (iter_core, iter_accessory)) in
            display_pairs.iter().zip(iter_pairs.iter())
        {
            assert_abs_diff_eq!(*core, *iter_core as f64, epsilon = 1e-4);
            let iter_accessory = iter_accessory.expect("dists_iter should yield Some for CoreAcc");
            assert_abs_diff_eq!(
                accessory.expect("Some in CoreAcc"),
                iter_accessory as f64,
                epsilon = 1e-4
            );
        }

        // Jaccard matrix (-k given): dists_iter should yield None
        let jaccard_type = set_k(&sketches, Some(17), false).expect("set_k failed");
        let jaccard_matrix = self_dists_all(&sketches, n, jaccard_type, true, None, 0.0);
        let display_pairs = parse_display_pairs(&jaccard_matrix.to_string(), false);
        let iter_pairs: Vec<(f32, Option<f32>)> = jaccard_matrix.dists_iter().collect();
        assert_eq!(display_pairs.len(), iter_pairs.len());
        for ((dist, _), (iter_dist, iter_accessory)) in display_pairs.iter().zip(iter_pairs.iter())
        {
            assert_abs_diff_eq!(*dist, *iter_dist as f64, epsilon = 1e-4);
            assert!(iter_accessory.is_none());
        }
    }

    #[test]
    fn distance_matrix_api_metadata_and_storage() {
        use sketchlib::distances::distance_matrix::{DistVec, Distances};
        use sketchlib::distances::{
            cross_dists_all, cross_dists_knn, self_dists_all, self_dists_knn, set_k,
        };
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        let references = MultiSketch::load(&sandbox.file_string("bact_db", TestDir::Output))
            .expect("failed to load reference sketches");
        let queries = MultiSketch::load(&sandbox.file_string("query_db", TestDir::Output))
            .expect("failed to load query sketches");
        let all_references = MultiSketch::load(&sandbox.file_string("ref_db", TestDir::Output))
            .expect("failed to load all reference sketches");

        let n_ref = references.number_samples_loaded();
        let n_query = queries.number_samples_loaded();
        let n_all = all_references.number_samples_loaded();

        // Dense self-query Core/Accessory matrices are condensed to n * (n - 1) / 2 rows
        // with one core and one accessory value per row.
        let dense_self = self_dists_all(
            &all_references,
            n_all,
            set_k(&all_references, None, false).expect("set_k failed"),
            true,
            None,
            0.0,
        );
        assert_eq!(dense_self.n_samples(), (n_all, None));
        assert_eq!(dense_self.shape(), (n_all * (n_all - 1) / 2, 2));
        assert_eq!(dense_self.dists_as_ref().len(), n_all * (n_all - 1));

        // Dense cross-query Jaccard matrices are rectangular and have one value per pair.
        let dense_cross = cross_dists_all(
            &references,
            &queries,
            n_ref,
            n_query,
            set_k(&references, Some(17), false).expect("set_k failed"),
            true,
            None,
            None,
            0.0,
        );
        assert_eq!(dense_cross.n_samples(), (n_ref, Some(n_query)));
        assert_eq!(dense_cross.shape(), (n_ref * n_query, 1));
        assert_eq!(dense_cross.dists_as_ref().len(), n_ref * n_query);

        // Sparse self-query Jaccard matrices retain knn values per reference row.
        let self_knn = 1;
        let sparse_self = self_dists_knn(
            &all_references,
            n_all,
            self_knn,
            set_k(&all_references, Some(17), false).expect("set_k failed"),
            true,
            None,
            0.0,
        );
        assert_eq!(sparse_self.n_samples(), (n_all, None));
        assert_eq!(sparse_self.shape(), (n_all * self_knn, self_knn));
        match sparse_self.dists_as_ref() {
            DistVec::Jaccard(distances) => assert_eq!(distances.len(), n_all * self_knn),
            DistVec::CoreAcc(_) => panic!("expected sparse Jaccard storage"),
        }

        // Sparse cross-query Core/Accessory matrices use query rows and reference columns.
        let cross_knn = 1;
        let sparse_cross = cross_dists_knn(
            &references,
            &queries,
            n_ref,
            n_query,
            cross_knn,
            set_k(&references, None, false).expect("set_k failed"),
            true,
            None,
            None,
            0.0,
        );
        assert_eq!(sparse_cross.n_samples(), (n_ref, Some(n_query)));
        assert_eq!(sparse_cross.shape(), (n_query * cross_knn, cross_knn));
        match sparse_cross.dists_as_ref() {
            DistVec::CoreAcc(distances) => assert_eq!(distances.len(), n_query * cross_knn),
            DistVec::Jaccard(_) => panic!("expected sparse CoreAcc storage"),
        }
    }

    /// Streaming self-mode dense output (`self_dists_all_stream`) must match the
    /// non-streaming `self_dists_all`/`Display` output exactly, as a sorted line
    /// multiset (chunk order across threads is not guaranteed, but both paths call
    /// the same per-pair math, so results should be bit-identical once sorted).
    /// Covers CoreAcc, Jaccard, and ANI.
    #[test]
    fn self_dists_stream_matches_non_streaming() {
        use sketchlib::distances::{self_dists_all, self_dists_all_stream, set_k};
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("stream_db")
            .args(["-v", "--k-seq", "17,31,4", "-s", "1000"])
            .arg("-f")
            .arg("rfile.txt")
            .assert()
            .success();

        let sketches = MultiSketch::load(&sandbox.file_string("stream_db", TestDir::Output))
            .expect("failed to load sketches");
        let n = sketches.number_samples_loaded();

        for (dist_type_a, dist_type_b) in [
            (
                set_k(&sketches, None, false).expect("set_k failed"),
                set_k(&sketches, None, false).expect("set_k failed"),
            ),
            (
                set_k(&sketches, Some(17), false).expect("set_k failed"),
                set_k(&sketches, Some(17), false).expect("set_k failed"),
            ),
            (
                set_k(&sketches, Some(17), true).expect("set_k failed"),
                set_k(&sketches, Some(17), true).expect("set_k failed"),
            ),
        ] {
            let matrix = self_dists_all(&sketches, n, dist_type_a, true, None, 0.0);
            let non_streamed = matrix.to_string();

            let mut streamed_bytes: Vec<u8> = Vec::new();
            self_dists_all_stream(
                &mut streamed_bytes,
                &sketches,
                n,
                dist_type_b,
                true,
                None,
                0.0,
                2,
            )
            .expect("self_dists_all_stream failed");
            let streamed = String::from_utf8(streamed_bytes).expect("non-utf8 stream output");

            let mut non_streamed_lines: Vec<&str> =
                non_streamed.lines().filter(|l| !l.is_empty()).collect();
            let mut streamed_lines: Vec<&str> =
                streamed.lines().filter(|l| !l.is_empty()).collect();
            non_streamed_lines.sort_unstable();
            streamed_lines.sort_unstable();
            assert_eq!(
                non_streamed_lines, streamed_lines,
                "Streamed and non-streamed self-mode output differ"
            );
        }
    }

    /// Cross-query analogue of `self_dists_stream_matches_non_streaming`.
    #[test]
    fn cross_dists_stream_matches_non_streaming() {
        use sketchlib::distances::{cross_dists_all, cross_dists_all_stream, set_k};
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        sketch_ref_and_query(&sandbox);

        let references = MultiSketch::load(&sandbox.file_string("bact_db", TestDir::Output))
            .expect("failed to load reference sketches");
        let queries = MultiSketch::load(&sandbox.file_string("query_db", TestDir::Output))
            .expect("failed to load query sketches");
        let n = references.number_samples_loaded();
        let n_query = queries.number_samples_loaded();

        for (dist_type_a, dist_type_b) in [
            (
                set_k(&references, None, false).expect("set_k failed"),
                set_k(&references, None, false).expect("set_k failed"),
            ),
            (
                set_k(&references, Some(17), false).expect("set_k failed"),
                set_k(&references, Some(17), false).expect("set_k failed"),
            ),
            (
                set_k(&references, Some(17), true).expect("set_k failed"),
                set_k(&references, Some(17), true).expect("set_k failed"),
            ),
        ] {
            let matrix = cross_dists_all(
                &references,
                &queries,
                n,
                n_query,
                dist_type_a,
                true,
                None,
                None,
                0.0,
            );
            let non_streamed = matrix.to_string();

            let mut streamed_bytes: Vec<u8> = Vec::new();
            cross_dists_all_stream(
                &mut streamed_bytes,
                &references,
                &queries,
                n,
                n_query,
                dist_type_b,
                true,
                None,
                None,
                0.0,
                2,
            )
            .expect("cross_dists_all_stream failed");
            let streamed = String::from_utf8(streamed_bytes).expect("non-utf8 stream output");

            let mut non_streamed_lines: Vec<&str> =
                non_streamed.lines().filter(|l| !l.is_empty()).collect();
            let mut streamed_lines: Vec<&str> =
                streamed.lines().filter(|l| !l.is_empty()).collect();
            non_streamed_lines.sort_unstable();
            streamed_lines.sort_unstable();
            assert_eq!(
                non_streamed_lines, streamed_lines,
                "Streamed and non-streamed cross-mode output differ"
            );
        }
    }

    /// Loads `legacy_db` (v0.1.3, legacy 14-bit format, R6 vs TIGR4, k=[17,21,25],
    /// sketch_size=128) as the legacy side of a mismatched-generation pair, and
    /// freshly sketches R6.fa.gz as the new-format side. Used only by the
    /// rejection tests below — the two sides deliberately don't need matching
    /// k-mer lengths/sketch size, since the generation guard fires first.
    fn legacy_and_fresh_dbs(
        sandbox: &TestSetup,
    ) -> (
        sketchlib::sketch::multisketch::MultiSketch,
        sketchlib::sketch::multisketch::MultiSketch,
    ) {
        use sketchlib::sketch::multisketch::MultiSketch;

        let legacy = MultiSketch::load(&sandbox.file_string("legacy_db", TestDir::Input))
            .expect("failed to load legacy_db");

        sandbox.copy_input_file_to_wd("R6.fa.gz");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "fresh_for_mismatch"])
            .args(["-v", "-k", "21"])
            .arg("R6.fa.gz")
            .assert()
            .success();
        let fresh = MultiSketch::load(&sandbox.file_string("fresh_for_mismatch", TestDir::Output))
            .expect("failed to load fresh_for_mismatch");

        (legacy, fresh)
    }

    /// Jaccard distance still computes correctly in legacy mode: `legacy_db`'s
    /// R6-vs-TIGR4 pair at k=17, via the public `self_dists_all` dispatch (which
    /// resolves to the legacy code path via `MultiSketch::is_legacy_format`).
    /// The expected value is a deterministic, previously-measured baseline for
    /// this fixture — this is a regression check, not a tolerance-based one.
    #[test]
    fn legacy_db_jaccard_distance_via_self_dists_all() {
        use sketchlib::distances::{self_dists_all, set_k};
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        let sketches = MultiSketch::load(&sandbox.file_string("legacy_db", TestDir::Input))
            .expect("failed to load legacy_db");
        assert!(sketches.is_legacy_format());
        let n = sketches.number_samples_loaded();

        let dist_type = set_k(&sketches, Some(17), false).expect("set_k failed");
        let matrix = self_dists_all(&sketches, n, dist_type, true, None, 0.0);
        let (dist, accessory) = matrix.dists_iter().next().expect("expected one pair");
        assert!(accessory.is_none());
        assert_abs_diff_eq!(dist as f64, 0.2343893, epsilon = 1e-4);
    }

    /// Core/accessory distance still computes correctly in legacy mode, using
    /// all three of `legacy_db`'s k-mer lengths (17, 21, 25).
    #[test]
    fn legacy_db_core_accessory_distance_via_self_dists_all() {
        use sketchlib::distances::{self_dists_all, set_k};
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        let sketches = MultiSketch::load(&sandbox.file_string("legacy_db", TestDir::Input))
            .expect("failed to load legacy_db");
        let n = sketches.number_samples_loaded();

        let dist_type = set_k(&sketches, None, false).expect("set_k failed");
        let matrix = self_dists_all(&sketches, n, dist_type, true, None, 0.0);
        let (core, accessory) = matrix.dists_iter().next().expect("expected one pair");
        let accessory = accessory.expect("CoreAcc pair should have an accessory value");
        assert_abs_diff_eq!(core as f64, 0.022036541, epsilon = 1e-4);
        assert_abs_diff_eq!(accessory as f64, 0.0, epsilon = 1e-4);
    }

    /// Cross-checks legacy-mode distance calculation against an independently
    /// re-sketched (new-format) version of the same two genomes. `legacy_db`
    /// (R6 vs TIGR4, v0.1.3, k=[17,21,25], sketch_size=128) is compared in
    /// legacy mode; R6.fa.gz/TIGR4.fa.gz are freshly sketched with the current
    /// version at matching k/sketch_size and compared in new mode. The two
    /// computed distances estimate the same true genomic distance via
    /// different (incompatible) binning schemes, so they should agree within a
    /// generous tolerance but are not expected to be bit-identical.
    #[test]
    fn legacy_and_new_format_distances_agree_for_same_genomes() {
        use sketchlib::distances::{self_dists_all, set_k};
        use sketchlib::sketch::multisketch::MultiSketch;

        let sandbox = TestSetup::setup();
        let legacy = MultiSketch::load(&sandbox.file_string("legacy_db", TestDir::Input))
            .expect("failed to load legacy_db");
        assert!(legacy.is_legacy_format());

        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "fresh_db"])
            .args(["-v", "--k-vals", "17,21,25", "-s", "128"])
            .arg("R6.fa.gz")
            .arg("TIGR4.fa.gz")
            .assert()
            .success();
        let fresh = MultiSketch::load(&sandbox.file_string("fresh_db", TestDir::Output))
            .expect("failed to load fresh_db");
        assert!(!fresh.is_legacy_format());

        let legacy_n = legacy.number_samples_loaded();
        let fresh_n = fresh.number_samples_loaded();

        // Jaccard at k=17
        let legacy_jaccard = self_dists_all(
            &legacy,
            legacy_n,
            set_k(&legacy, Some(17), false).expect("set_k failed"),
            true,
            None,
            0.0,
        )
        .dists_iter()
        .next()
        .expect("expected one pair")
        .0;
        let fresh_jaccard = self_dists_all(
            &fresh,
            fresh_n,
            set_k(&fresh, Some(17), false).expect("set_k failed"),
            true,
            None,
            0.0,
        )
        .dists_iter()
        .next()
        .expect("expected one pair")
        .0;
        assert_abs_diff_eq!(legacy_jaccard as f64, fresh_jaccard as f64, epsilon = 0.05);

        // Core/accessory
        let (legacy_core, legacy_acc) = self_dists_all(
            &legacy,
            legacy_n,
            set_k(&legacy, None, false).expect("set_k failed"),
            true,
            None,
            0.0,
        )
        .dists_iter()
        .next()
        .map(|(c, a)| (c, a.expect("Some in CoreAcc")))
        .expect("expected one pair");
        let (fresh_core, fresh_acc) = self_dists_all(
            &fresh,
            fresh_n,
            set_k(&fresh, None, false).expect("set_k failed"),
            true,
            None,
            0.0,
        )
        .dists_iter()
        .next()
        .map(|(c, a)| (c, a.expect("Some in CoreAcc")))
        .expect("expected one pair");
        assert_abs_diff_eq!(legacy_core as f64, fresh_core as f64, epsilon = 0.1);
        assert_abs_diff_eq!(legacy_acc as f64, fresh_acc as f64, epsilon = 0.1);
    }

    #[test]
    #[should_panic(
        expected = "Cannot compare reference and query databases with different sketch generations"
    )]
    fn cross_dists_all_rejects_mismatched_generations() {
        use sketchlib::distances::{cross_dists_all, set_k};

        let sandbox = TestSetup::setup();
        let (legacy, fresh) = legacy_and_fresh_dbs(&sandbox);
        let n = legacy.number_samples_loaded();
        let n_query = fresh.number_samples_loaded();
        let dist_type = set_k(&legacy, Some(17), false).expect("set_k failed");
        cross_dists_all(
            &legacy, &fresh, n, n_query, dist_type, true, None, None, 0.0,
        );
    }

    #[test]
    #[should_panic(
        expected = "Cannot compare reference and query databases with different sketch generations"
    )]
    fn cross_dists_knn_rejects_mismatched_generations() {
        use sketchlib::distances::{cross_dists_knn, set_k};

        let sandbox = TestSetup::setup();
        let (legacy, fresh) = legacy_and_fresh_dbs(&sandbox);
        let n = legacy.number_samples_loaded();
        let n_query = fresh.number_samples_loaded();
        let dist_type = set_k(&legacy, Some(17), false).expect("set_k failed");
        cross_dists_knn(
            &legacy, &fresh, n, n_query, 1, dist_type, true, None, None, 0.0,
        );
    }

    #[test]
    fn cross_dists_all_stream_rejects_mismatched_generations() {
        use sketchlib::distances::{cross_dists_all_stream, set_k};

        let sandbox = TestSetup::setup();
        let (legacy, fresh) = legacy_and_fresh_dbs(&sandbox);
        let n = legacy.number_samples_loaded();
        let n_query = fresh.number_samples_loaded();
        let dist_type = set_k(&legacy, Some(17), false).expect("set_k failed");
        let mut buf: Vec<u8> = Vec::new();
        let result = cross_dists_all_stream(
            &mut buf, &legacy, &fresh, n, n_query, dist_type, true, None, None, 0.0, 1,
        );
        let err = result.expect_err("expected mismatched-generation rejection");
        assert!(
            err.to_string().contains(
                "Cannot compare reference and query databases with different sketch generations"
            ),
            "Unexpected error: {err}"
        );
    }

    #[test]
    fn cli_dist_rejects_mismatched_generations() {
        let sandbox = TestSetup::setup();
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "cli_fresh"])
            .args(["-v", "-k", "21"])
            .arg("R6.fa.gz")
            .assert()
            .success();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg(sandbox.file_string("legacy_db", TestDir::Input))
            .arg("cli_fresh")
            .args(["-k", "21"])
            .assert()
            .failure();
    }

    /// Once the hard version error becomes a warning, `Merge` gains a real risk
    /// it didn't need to guard against before: silently byte-concatenating a
    /// legacy (14-bit-packed) and new-format (16-bit-packed) `.skd` into a
    /// corrupted file. `MultiSketch::is_compatible_with` now also requires
    /// matching `is_legacy_format()`, so this must fail loudly instead. The
    /// fresh database is sketched with the same k-mer lengths/sketch_size as
    /// `legacy_db` so the *only* mismatch is generation, isolating this check
    /// from the pre-existing kmer_lengths/sketch_size/hash_type checks.
    #[test]
    fn merge_rejects_mismatched_generations() {
        let sandbox = TestSetup::setup();
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "fresh_matching_shape"])
            .args(["-v", "--k-vals", "17,21,25", "-s", "128"])
            .arg("R6.fa.gz")
            .assert()
            .success();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("merge")
            .arg(sandbox.file_string("legacy_db", TestDir::Input))
            .arg("fresh_matching_shape")
            .args(["-o", "bad_merge"])
            .assert()
            .failure();
    }
}
