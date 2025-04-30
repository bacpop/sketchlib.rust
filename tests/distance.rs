use snapbox::cmd::{cargo_bin, Command};
use std::path::Path;

pub mod common;
use crate::common::*;

use std::collections::HashMap;
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
            .stdout_eq(
                sandbox
                    .snapbox_file("dists_knn_ca.stdout", TestDir::Correct)
            );

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
            .stdout_eq(
                sandbox
                    .snapbox_file("dists_knn_jaccard.stdout", TestDir::Correct)
            );

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
            .stdout_eq(
                sandbox
                    .snapbox_file("dists_knn_ani.stdout", TestDir::Correct)
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
            .stdout_eq(
                sandbox
                    .snapbox_file("dists_subset.stdout", TestDir::Correct)
            );

    }

}
