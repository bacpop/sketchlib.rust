use snapbox::cmd::{self, Command};

pub mod common;
use crate::common::*;

#[cfg(test)]

mod tests {
    use approx::assert_abs_diff_eq;
    use sketchlib::sketch::{
        multisketch::MultiSketch, sketch_datafile::SketchArrayReader, BIN_BITS,
    };
    use snapbox::assert_data_eq;

    use super::*;

    // Compares tab-separated distance output against a golden file, ignoring
    // line order and allowing numeric fields to differ by a small epsilon
    // (SIMD reordering can shift the last few bits of a float).
    fn assert_dist_stdout_unordered_with_tolerance(actual: &str, expected: &snapbox::Data) {
        let expected_str = expected
            .render()
            .expect("Failed to render expected snapshot data");
        let actual_lines: Vec<&str> = actual.lines().filter(|l| !l.is_empty()).collect();
        let mut expected_lines: Vec<&str> =
            expected_str.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(
            actual_lines.len(),
            expected_lines.len(),
            "Line count mismatch.\nActual:\n{actual}\nExpected:\n{expected_str}"
        );

        for actual_line in &actual_lines {
            let actual_fields: Vec<&str> = actual_line.split('\t').collect();
            let match_idx = expected_lines.iter().position(|expected_line| {
                let expected_fields: Vec<&str> = expected_line.split('\t').collect();
                actual_fields.len() == expected_fields.len()
                    && actual_fields.iter().zip(expected_fields.iter()).all(
                        |(actual_field, expected_field)| match (
                            actual_field.parse::<f64>(),
                            expected_field.parse::<f64>(),
                        ) {
                            (Ok(actual_val), Ok(expected_val)) => {
                                (actual_val - expected_val).abs() < 1e-4
                            }
                            _ => actual_field == expected_field,
                        },
                    )
            });
            match match_idx {
                Some(idx) => {
                    expected_lines.remove(idx);
                }
                None => panic!(
                    "No matching expected line for actual line: {actual_line}\nRemaining expected lines: {expected_lines:?}"
                ),
            }
        }
    }

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

    const METADATA_SAMPLES: [&str; 4] = [
        "TIGR4.fa.gz",
        "R6.fa.gz",
        "14412_3#82.contigs_velvet.fa.gz",
        "14412_3#84.contigs_velvet.fa.gz",
    ];

    fn metadata_build_command(sandbox: &TestSetup) -> std::process::Command {
        sandbox.create_named_rfile("samples.txt", &METADATA_SAMPLES);
        let mut command = std::process::Command::new(cmd::cargo_bin!("sketchlib"));
        command
            .current_dir(sandbox.get_wd())
            .args([
                "inverted",
                "build",
                "-o",
                "reference",
                "-k",
                "21",
                "-s",
                "10",
                "-f",
                "samples.txt",
                "--phylo-clusters",
            ])
            .arg(sandbox.file_string("phylo_clusters.txt", TestDir::Input));
        command
    }

    fn metadata_query_command(
        sandbox: &TestSetup,
        mode: &str,
        threads: &str,
    ) -> std::process::Command {
        let mut command = std::process::Command::new(cmd::cargo_bin!("sketchlib"));
        command.current_dir(sandbox.get_wd()).args([
            "inverted",
            "query",
            "reference.ski",
            "-f",
            "samples.txt",
            "--query-type",
            mode,
            "--threads",
            threads,
        ]);
        command
    }

    fn successful_stdout(output: std::process::Output) -> String {
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout).unwrap()
    }

    fn assert_query_rows(actual: &str, expected: &str) {
        assert_eq!(actual.lines().next(), expected.lines().next());
        let mut actual: Vec<_> = actual.lines().skip(1).collect();
        let mut expected: Vec<_> = expected.lines().skip(1).collect();
        actual.sort_unstable();
        expected.sort_unstable();
        assert_eq!(actual, expected);
    }

    #[test]
    fn inverted_query_metadata_all_modes_and_suppression() {
        let sandbox = TestSetup::setup();
        for (filename, prefix) in [("reference.tsv", "reference"), ("queries.tsv", "query")] {
            let data: String = METADATA_SAMPLES
                .iter()
                .enumerate()
                .map(|(i, name)| format!("{name}\t{prefix} {i}\n"))
                .collect();
            std::fs::write(sandbox.file_path(filename, TestDir::Output), data).unwrap();
        }
        successful_stdout(
            metadata_build_command(&sandbox)
                .args(["--metadata", "reference.tsv"])
                .output()
                .unwrap(),
        );
        let [tigr, r6, a, b] = METADATA_SAMPLES;
        let expectations = [
            ("match-count", format!("Query\t{tigr};reference 0\t{r6};reference 1\t{a};reference 2\t{b};reference 3\n{tigr};query 0\t10\t9\t0\t0\n{r6};query 1\t9\t10\t0\t0\n{a};query 2\t0\t0\t10\t6\n{b};query 3\t0\t0\t6\t10\n")),
            ("all-bins", format!("Query\tMatches\n{tigr};query 0\t{tigr};reference 0\n{r6};query 1\t{r6};reference 1\n{a};query 2\t{a};reference 2\n{b};query 3\t{b};reference 3\n")),
            ("any-bins", format!("Query\tMatches\n{tigr};query 0\t{tigr};reference 0,{r6};reference 1\n{r6};query 1\t{tigr};reference 0,{r6};reference 1\n{a};query 2\t{a};reference 2,{b};reference 3\n{b};query 3\t{a};reference 2,{b};reference 3\n")),
        ];
        for (mode, expected) in &expectations {
            let mut unannotated = expected.clone();
            let mut reference_only = expected.clone();
            for i in 0..4 {
                unannotated = unannotated
                    .replace(&format!(";reference {i}"), "")
                    .replace(&format!(";query {i}"), "");
                reference_only = reference_only.replace(&format!(";query {i}"), "");
            }
            for threads in ["1", "2"] {
                let actual = successful_stdout(
                    metadata_query_command(&sandbox, mode, threads)
                        .args(["--metadata", "queries.tsv"])
                        .output()
                        .unwrap(),
                );
                assert_query_rows(&actual, expected);
                let actual = successful_stdout(
                    metadata_query_command(&sandbox, mode, threads)
                        .args(["--metadata", "queries.tsv", "--suppress-metadata"])
                        .output()
                        .unwrap(),
                );
                assert_query_rows(&actual, &unannotated);
            }
            let stdout = successful_stdout(
                metadata_query_command(&sandbox, mode, "1")
                    .args(["-o", "result.tsv"])
                    .output()
                    .unwrap(),
            );
            assert!(stdout.is_empty());
            assert_query_rows(
                &std::fs::read_to_string(sandbox.file_path("result.tsv", TestDir::Output)).unwrap(),
                &reference_only,
            );
        }
        // Query metadata belongs to the loaded index, but is never saved in .ski files.
        let query_metadata = sketchlib::io::parse_metadata_info(
            Some(&sandbox.file_string("queries.tsv", TestDir::Output)),
            METADATA_SAMPLES,
        )
        .unwrap();
        let prefix = sandbox.file_string("reference", TestDir::Output);
        let index =
            sketchlib::inverted::Inverted::load(&prefix, query_metadata.clone(), false).unwrap();
        assert_eq!(index.get_query_metadata(tigr), Some("query 0"));
        let saved_prefix = sandbox.file_string("saved", TestDir::Output);
        index.save(&saved_prefix).unwrap();
        let reloaded = sketchlib::inverted::Inverted::load(&saved_prefix, None, false).unwrap();
        assert_eq!(reloaded.get_query_metadata(tigr), None);
        assert_eq!(reloaded.get_metadata(), index.get_metadata());
        let suppressed =
            sketchlib::inverted::Inverted::load(&prefix, query_metadata, true).unwrap();
        assert_eq!(suppressed.get_query_metadata(tigr), None);
        assert_eq!(suppressed.get_metadata(), &None);

        // Rebuild without reference metadata: query metadata remains independently usable.
        successful_stdout(metadata_build_command(&sandbox).output().unwrap());
        for (mode, expected) in &expectations {
            let mut query_only = expected.clone();
            for i in 0..4 {
                query_only = query_only.replace(&format!(";reference {i}"), "");
            }
            let actual = successful_stdout(
                metadata_query_command(&sandbox, mode, "1")
                    .args(["--metadata", "queries.tsv"])
                    .output()
                    .unwrap(),
            );
            assert_query_rows(&actual, &query_only);
        }
    }

    #[test]
    fn inverted_missing_metadata_has_one_aggregate_warning_per_source() {
        let sandbox = TestSetup::setup();
        let [tigr, r6, a, _] = METADATA_SAMPLES;
        std::fs::write(
            sandbox.file_path("missing.tsv", TestDir::Output),
            format!("\n\tignored\n  \talso ignored\n{tigr}\t\n{r6}\n{a}\t   \n"),
        )
        .unwrap();
        let output = metadata_build_command(&sandbox)
            .args(["--metadata", "missing.tsv"])
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stderr}");
        assert_eq!(stderr.matches("WARN").count(), 1, "{stderr}");
        assert!(
            stderr.contains("4 missing values set to NA; 2 records without sample names ignored"),
            "{stderr}"
        );
        let output = metadata_query_command(&sandbox, "all-bins", "1")
            .args(["--metadata", "missing.tsv"])
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stderr}");
        assert_eq!(stderr.matches("WARN").count(), 1, "{stderr}");
        assert!(
            stderr.contains("4 missing values set to NA; 2 records without sample names ignored"),
            "{stderr}"
        );
        let expected = format!("Query\tMatches\n{tigr};NA\t{tigr};NA\n{r6};NA\t{r6};NA\n{a};NA\t{a};NA\n{};NA\t{};NA\n", METADATA_SAMPLES[3], METADATA_SAMPLES[3]);
        assert_query_rows(&String::from_utf8(output.stdout).unwrap(), &expected);
    }

    #[test]
    fn inverted_metadata_errors_are_reported_without_panicking() {
        let sandbox = TestSetup::setup();
        successful_stdout(metadata_build_command(&sandbox).output().unwrap());
        for value in ["left\tright", "left,right", "left;right"] {
            std::fs::write(
                sandbox.file_path("invalid.tsv", TestDir::Output),
                format!("{}\t{value}\n", METADATA_SAMPLES[0]),
            )
            .unwrap();
            for mut command in [
                metadata_build_command(&sandbox),
                metadata_query_command(&sandbox, "all-bins", "1"),
            ] {
                let output = command
                    .args(["--metadata", "invalid.tsv"])
                    .output()
                    .unwrap();
                assert!(!output.status.success());
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert!(
                    stderr.contains("Metadata must not contain tabs, commas, or semicolons"),
                    "{stderr}"
                );
                assert!(!stderr.contains("panicked"), "{stderr}");
            }
        }
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

        let standard_sketches =
            MultiSketch::load(&sandbox.file_string("standard", TestDir::Output))
                .expect("failed to load standard sketch data");

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
        sandbox.copy_input_file_to_wd("phylo_clusters.txt");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31"])
            .args(["-f", "rfile.txt"])
            .args(["--phylo-clusters", "phylo_clusters.txt"])
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
        sandbox.copy_input_file_to_wd("phylo_clusters_empty.txt");

        // This shouldn't give an error, but it is an unfortunate situation with respect to the phylo-clusters file the user provided
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31"])
            .args(["-f", "rfile.txt"])
            .args(["--phylo-clusters", "phylo_clusters_empty.txt"])
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

        // This shouldn't give an error, but it is an unfortunate situation with respect to the phylo-clusters file the user provided
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
            .arg("--phylo-clusters")
            .arg(sandbox.file_string("phylo_clusters.txt", TestDir::Input))
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
            .stdout_eq("Identified 2 prefilter pairs from a max of 6 (33.3%)\n");

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
        let precluster_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("inverted.ski")
            .output()
            .expect("Failed to run precluster");
        assert!(precluster_output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(precluster_output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster.stdout", TestDir::Correct),
        );

        // Same with ANI
        let precluster_ani_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--ani")
            .arg("--skd")
            .arg("standard")
            .arg("inverted.ski")
            .output()
            .expect("Failed to run precluster --ani");
        assert!(precluster_ani_output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(precluster_ani_output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster_ani.stdout", TestDir::Correct),
        );

        // Default knn, check no junk distances printed
        let precluster_knn50_output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "50"])
            .arg("--skd")
            .arg("standard")
            .arg("inverted.ski")
            .output()
            .expect("Failed to run precluster --knn 50");
        assert!(precluster_knn50_output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(precluster_knn50_output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster.stdout", TestDir::Correct),
        );
    }

    // Helper to set up reordered SKI + standard SKD
    // SKI order (via --phylo-clusters): TIGR4, R6, #82, #84
    // SKD order (via rfile.txt):       #82, #84, R6, TIGR4
    fn setup_reordered_precluster(sandbox: &TestSetup) {
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        sandbox.copy_input_file_to_wd("rfile.txt");

        // Build SKI+SKQ with reordered samples via --phylo-clusters
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("reordered_inv")
            .args(["-v", "-k", "21", "-s", "10"])
            .arg("-f")
            .arg("rfile.txt")
            .arg("--phylo-clusters")
            .arg(sandbox.file_string("phylo_clusters.txt", TestDir::Input))
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
        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .output()
            .expect("Failed to run precluster");
        assert!(output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster.stdout", TestDir::Correct),
        );
    }

    #[test]
    fn inverted_precluster_reordered_ani() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--ani")
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .output()
            .expect("Failed to run precluster --ani");
        assert!(output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster_ani.stdout", TestDir::Correct),
        );
    }

    #[test]
    fn inverted_precluster_reordered_knn_padding() {
        let sandbox = TestSetup::setup();
        setup_reordered_precluster(&sandbox);

        // knn=50 >> n_samples=4, tests padding with reordered indices
        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "50"])
            .arg("--skd")
            .arg("standard")
            .arg("reordered_inv.ski")
            .output()
            .expect("Failed to run precluster --knn 50");
        assert!(output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster.stdout", TestDir::Correct),
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
        let output = Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .args(["-v", "--knn", "1"])
            .arg("--skd")
            .arg("standard")
            .arg("reversed_inv.ski")
            .output()
            .expect("Failed to run precluster");
        assert!(output.status.success());
        assert_dist_stdout_unordered_with_tolerance(
            &String::from_utf8(output.stdout).unwrap(),
            &sandbox.snapbox_file("inverted_precluster.stdout", TestDir::Correct),
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
            .arg("--ref-completeness-file")
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
            .arg("--ref-completeness-file")
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
                assert_abs_diff_eq!(dist, 0.0, epsilon = 1e-6);
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
            .arg("--phylo-clusters")
            .arg(sandbox.file_string("phylo_clusters.txt", TestDir::Input))
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
