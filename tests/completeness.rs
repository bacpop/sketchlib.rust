use snapbox::cmd::{cargo_bin, Command};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub mod common;
use crate::common::*;

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper function to create completeness files with specified entries
    fn create_completeness_file(sandbox: &TestSetup, filename: &str, entries: &[(&str, f64)]) {
        let file_path = sandbox.file_string(filename, TestDir::Output);
        let mut file = File::create(&file_path)
            .unwrap_or_else(|_| panic!("Failed to create completeness file: {}", file_path));
        
        for (genome, completeness) in entries {
            writeln!(file, "{}\t{}", genome, completeness)
                .expect("Failed to write to completeness file");
        }
    }

    /// Helper function to parse distance output and extract distance values
    fn parse_distances(output_file: &str) -> Vec<f64> {
        let file = File::open(output_file).expect("Failed to open distance output file");
        let reader = BufReader::new(file);
        
        reader
            .lines()
            .map(|line| {
                line.expect("Failed to read line")
                    .split_whitespace()
                    .last()
                    .expect("Line has incorrect format")
                    .parse::<f64>()
                    .expect("Failed to parse distance")
            })
            .collect()
    }

    /// Assert function with tolerance for distance comparisons
    fn assert_distance_with_tolerance(actual: f64, expected: f64, tolerance: f64) {
        let diff = (actual - expected).abs();
        assert!(
            diff <= tolerance,
            "Distance difference exceeds tolerance of {}. Actual: {}, Expected: {}, Diff: {}",
            tolerance,
            actual,
            expected,
            diff
        );
    }

    #[test]
    fn test_completeness_ordering() {
        let sandbox = TestSetup::setup();

        // Copy test files to working directory
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");

        // Create sketch database with genomes in alphabetical order
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("test_genomes")
            .args(["-v", "-k", "31", "-s", "1000"])
            .arg("14412_3#82.contigs_velvet.fa.gz")  // First genome
            .arg("14412_3#84.contigs_velvet.fa.gz")  // Second genome  
            .arg("R6.fa.gz")                         // Third genome
            .assert()
            .success();

        // Create completeness file with REVERSE order and very different completeness values
        create_completeness_file(
            &sandbox,
            "completeness_reverse.txt",
            &[
                ("R6.fa.gz", 0.5),                         // Third genome, very low completeness
                ("14412_3#84.contigs_velvet.fa.gz", 0.9),  // Second genome, high completeness
                ("14412_3#82.contigs_velvet.fa.gz", 0.6),  // First genome, low completeness
            ],
        );

        // Calculate distances WITHOUT completeness correction
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test_genomes")
            .args(["-o", "distances_no_correction"])
            .args(["-k", "31"])
            .arg("-v")
            .assert()
            .success();

        // Calculate distances WITH completeness correction
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test_genomes")
            .args(["-o", "distances_with_correction"])
            .arg("--completeness-file")
            .arg("completeness_reverse.txt")
            .args(["-k", "31"])
            .arg("-v")
            .assert()
            .success();

        // Parse both distance outputs
        let distances_no_correction = parse_distances(&sandbox.file_string("distances_no_correction", TestDir::Output));
        let distances_with_correction = parse_distances(&sandbox.file_string("distances_with_correction", TestDir::Output));

        // Verify that distances are different (correction was applied)
        assert_eq!(distances_no_correction.len(), distances_with_correction.len());
        
        // At least one distance should be different due to completeness correction
        let mut differences_found = false;
        for (uncorrected, corrected) in distances_no_correction.iter().zip(distances_with_correction.iter()) {
            if (uncorrected - corrected).abs() > 0.001 {
                differences_found = true;
                break;
            }
        }
        assert!(differences_found, "Completeness correction should change at least one distance");

        // Verify files exist
        assert_eq!(true, sandbox.file_exists("distances_no_correction"));
        assert_eq!(true, sandbox.file_exists("distances_with_correction"));
    }

    #[test]
    fn test_completeness_missing_genomes() {
        let sandbox = TestSetup::setup();

        // Copy test files
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");

        // Create sketch database
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("test_missing")
            .args(["-v", "-k", "21", "-s", "1000"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .arg("R6.fa.gz")
            .assert()
            .success();

        // Create completeness file with MISSING genome (R6 is missing)
        create_completeness_file(
            &sandbox,
            "completeness_missing.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.8),
                ("14412_3#84.contigs_velvet.fa.gz", 0.9),
                // R6 is intentionally missing - should default to 1.0
            ],
        );

        // Calculate distances with missing genome in completeness file
        // This should succeed and use default completeness of 1.0 for R6
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test_missing")
            .args(["-k", "21", "-o", "distances_missing"])
            .arg("--completeness-file")
            .arg("completeness_missing.txt")
            .arg("-v")
            .assert()
            .success();

        // Verify output file exists
        assert_eq!(true, sandbox.file_exists("distances_missing"));

        // Parse distances to ensure they are valid numbers
        let distances = parse_distances(&sandbox.file_string("distances_missing", TestDir::Output));
        assert!(distances.len() > 0, "Should have calculated some distances");
        
        // All distances should be valid (not NaN or infinite)
        for distance in distances {
            assert!(distance.is_finite(), "Distance should be finite: {}", distance);
            assert!(distance >= 0.0 && distance <= 1.0, "Distance should be between 0 and 1: {}", distance);
        }
    }

    #[test]
    fn test_completeness_extra_genomes() {
        let sandbox = TestSetup::setup();

        // Copy test files
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");

        // Create sketch database with only 2 genomes
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("test_extra")
            .args(["-v", "-k", "21", "-s", "1000"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .assert()
            .success();

        // Create completeness file with EXTRA genomes not in sketch database
        create_completeness_file(
            &sandbox,
            "completeness_extra.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.8),
                ("14412_3#84.contigs_velvet.fa.gz", 0.9),
                ("NonExistentGenome1", 0.5),  // Extra genome - should be ignored
                ("NonExistentGenome2", 0.6),  // Extra genome - should be ignored
                ("AnotherFakeGenome", 0.7),   // Extra genome - should be ignored
            ],
        );

        // Calculate distances with extra genomes in completeness file
        // This should succeed and ignore the extra genomes
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test_extra")
            .args(["-k", "21", "-o", "distances_extra"])
            .arg("--completeness-file")
            .arg("completeness_extra.txt")
            .arg("-v")
            .assert()
            .success();

        // Verify output file exists
        assert_eq!(true, sandbox.file_exists("distances_extra"));

        // Parse distances to ensure they are valid
        let distances = parse_distances(&sandbox.file_string("distances_extra", TestDir::Output));
        assert!(distances.len() > 0, "Should have calculated some distances");
        
        // All distances should be valid
        for distance in distances {
            assert!(distance.is_finite(), "Distance should be finite: {}", distance);
            assert!(distance >= 0.0 && distance <= 1.0, "Distance should be between 0 and 1: {}", distance);
        }
    }


    #[test]
    fn test_completeness_precluster() {
        let sandbox = TestSetup::setup();

        // Copy test files
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");

        // Create inverted index for preclustering
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("build")
            .arg("-o")
            .arg("precluster_index")
            .args(["-v", "-k", "21", "-s", "10", "--write-skq"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .arg("R6.fa.gz")
            .assert()
            .success();

        // Create standard sketch database
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("precluster_sketches")
            .args(["-v", "-k", "21", "-s", "1000"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .arg("R6.fa.gz")
            .assert()
            .success();

        // Create completeness file
        create_completeness_file(
            &sandbox,
            "completeness_precluster.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.8),
                ("14412_3#84.contigs_velvet.fa.gz", 0.9),
                ("R6.fa.gz", 0.7),
            ],
        );

        // Run precluster with completeness correction
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("precluster")
            .arg("precluster_index.ski")
            .arg("--skd")
            .arg("precluster_sketches")
            .args(["-v", "--knn", "2"])
            .arg("--completeness-file")
            .arg("completeness_precluster.txt")
            .args(["-o", "precluster_distances"])
            .assert()
            .success();

        // Verify output file exists
        assert_eq!(true, sandbox.file_exists("precluster_distances"));

        // Parse distances to ensure they are valid
        let distances = parse_distances(&sandbox.file_string("precluster_distances", TestDir::Output));
        assert!(distances.len() > 0, "Should have calculated some distances");
        
        // All distances should be valid
        for distance in distances {
            assert!(distance.is_finite(), "Distance should be finite: {}", distance);
            assert!(distance >= 0.0 && distance <= 1.0, "Distance should be between 0 and 1: {}", distance);
        }
    }
}