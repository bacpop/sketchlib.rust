use approx::assert_abs_diff_eq;
use snapbox::cmd::{cargo_bin, Command};
// use std::fs::File;
// use std::io::{BufRead, BufReader, Write};

pub mod common;
use crate::common::TestSetup;
use crate::common::*;

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert function with tolerance for distance comparisons using approx crate
    fn assert_distance_with_tolerance(actual: f32, expected: f32, tolerance: f32) {
        assert_abs_diff_eq!(actual, expected, epsilon = tolerance);
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
            .arg("14412_3#82.contigs_velvet.fa.gz") // First genome
            .arg("14412_3#84.contigs_velvet.fa.gz") // Second genome
            .arg("R6.fa.gz") // Third genome
            .assert()
            .success();

        // Create completeness file with values that test cutoff behavior
        // Use values where all c1 * c2 > 0.64 to ensure corrections are applied
        TestSetup::create_completeness_file(
            &sandbox,
            "completeness_cutoff_test.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.8), // Good completeness
                ("14412_3#84.contigs_velvet.fa.gz", 0.85), // Good completeness
                ("R6.fa.gz", 0.9),                        // High completeness
                                                          // With default cutoff 0.64:
                                                          // - 0.8*0.85=0.68 >= 0.64 (correction applied)
                                                          // - 0.8*0.9=0.72 >= 0.64 (correction applied)
                                                          // - 0.85*0.9=0.765 >= 0.64 (correction applied)
                                                          // With high cutoff 0.8:
                                                          // - 0.8*0.85=0.68 < 0.8 (no correction)
                                                          // - 0.8*0.9=0.72 < 0.8 (no correction)
                                                          // - 0.85*0.9=0.765 < 0.8 (no correction)
            ],
        );

        // First, let's inspect the sketch database to understand genome ordering
        let sketch_info_output = Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("info")
            .arg("test_genomes")
            .output()
            .expect("Failed to get sketch info");

        println!("Sketch database info:");
        println!("{}", String::from_utf8_lossy(&sketch_info_output.stdout));

        // Also print the completeness file contents
        let completeness_content = std::fs::read_to_string(
            sandbox.file_string("completeness_cutoff_test.txt", TestDir::Output),
        )
        .expect("Failed to read completeness file");
        println!("Completeness file contents:");
        println!("{}", completeness_content);

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

        // Calculate distances WITH completeness correction using default cutoff (0.64)
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test_genomes")
            .args(["-o", "distances_with_correction_default_cutoff"])
            .arg("--completeness-file")
            .arg("completeness_cutoff_test.txt")
            .args(["-k", "31"])
            .arg("-v")
            .assert()
            .success();

        // Calculate distances WITH completeness correction using high cutoff (0.8)
        // This should correct no pairs since all c1*c2 products < 0.8
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("test_genomes")
            .args(["-o", "distances_with_correction_high_cutoff"])
            .arg("--completeness-file")
            .arg("completeness_cutoff_test.txt")
            .arg("--completeness-cutoff")
            .arg("0.8")
            .args(["-k", "31"])
            .arg("-v")
            .assert()
            .success();

        // Parse all distance outputs
        let distances_no_correction = TestSetup::parse_distances(
            &sandbox.file_string("distances_no_correction", TestDir::Output),
        );
        let distances_default_cutoff = TestSetup::parse_distances(
            &sandbox.file_string("distances_with_correction_default_cutoff", TestDir::Output),
        );
        let distances_high_cutoff = TestSetup::parse_distances(
            &sandbox.file_string("distances_with_correction_high_cutoff", TestDir::Output),
        );

        // Debug: Print all distance values
        println!(
            "Distances without correction: {:?}",
            distances_no_correction
        );
        println!(
            "Distances with default cutoff (0.64): {:?}",
            distances_default_cutoff
        );
        println!(
            "Distances with high cutoff (0.8): {:?}",
            distances_high_cutoff
        );

        // Verify that all outputs have the same number of distances
        assert_eq!(
            distances_no_correction.len(),
            distances_default_cutoff.len()
        );
        assert_eq!(distances_no_correction.len(), distances_high_cutoff.len());

        // With default cutoff (0.64), all meaningful pairs should be corrected since all c1*c2 >= 0.64
        let mut default_cutoff_differences = 0;
        for (uncorrected, corrected) in distances_no_correction
            .iter()
            .zip(distances_default_cutoff.iter())
        {
            if (uncorrected - corrected).abs() > 0.001 {
                default_cutoff_differences += 1;
                println!(
                    "Default cutoff difference found: {} vs {}",
                    uncorrected, corrected
                );
            }
        }

        // With high cutoff (0.8), no pairs should be corrected since all c1*c2 < 0.8
        let mut high_cutoff_differences = 0;
        for (uncorrected, corrected) in distances_no_correction
            .iter()
            .zip(distances_high_cutoff.iter())
        {
            if (uncorrected - corrected).abs() > 0.001 {
                high_cutoff_differences += 1;
                println!(
                    "High cutoff difference found: {} vs {}",
                    uncorrected, corrected
                );
            }
        }

        // Debug: Print the actual differences found
        println!(
            "Default cutoff (0.64) differences: {}",
            default_cutoff_differences
        );
        println!("High cutoff (0.8) differences: {}", high_cutoff_differences);

        // Count how many pairs have meaningful distances (not 1.0)
        let meaningful_pairs = distances_no_correction
            .iter()
            .filter(|&&d| d < 0.99)
            .count();
        println!("Number of meaningful distance pairs: {}", meaningful_pairs);

        // With default cutoff (0.64): All meaningful pairs should be corrected
        // - 0.8*0.85=0.68 >= 0.64 (correction applied)
        // - 0.8*0.9=0.72 >= 0.64 (correction applied)
        // - 0.85*0.9=0.765 >= 0.64 (correction applied)
        assert_eq!(
            default_cutoff_differences, meaningful_pairs,
            "Default cutoff (0.64) should result in corrections for all meaningful pairs (all c1*c2 >= 0.64)"
        );

        // With high cutoff (0.8): No pairs should be corrected
        // - 0.8*0.85=0.68 < 0.8 (no correction)
        // - 0.8*0.9=0.72 < 0.8 (no correction)
        // - 0.85*0.9=0.765 < 0.8 (no correction)
        assert_eq!(
            high_cutoff_differences, 0,
            "High cutoff (0.8) should result in 0 corrections (all c1*c2 < 0.8)"
        );

        // At least some correction should happen with default cutoff
        assert!(
            default_cutoff_differences > 0,
            "Default cutoff should result in at least some completeness corrections"
        );

        // Verify files exist
        assert_eq!(true, sandbox.file_exists("distances_no_correction"));
        assert_eq!(
            true,
            sandbox.file_exists("distances_with_correction_default_cutoff")
        );
        assert_eq!(
            true,
            sandbox.file_exists("distances_with_correction_high_cutoff")
        );

        println!("✅ Completeness cutoff test passed!");
        println!(
            "   - Default cutoff (0.64): {} corrections",
            default_cutoff_differences
        );
        println!(
            "   - High cutoff (0.8): {} corrections",
            high_cutoff_differences
        );
        println!("   - Meaningful distance pairs: {}", meaningful_pairs);
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
        TestSetup::create_completeness_file(
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
        let distances =
            TestSetup::parse_distances(&sandbox.file_string("distances_missing", TestDir::Output));
        assert!(distances.len() > 0, "Should have calculated some distances");

        // All distances should be valid (not NaN or infinite)
        for distance in distances {
            assert!(
                distance.is_finite(),
                "Distance should be finite: {}",
                distance
            );
            assert!(
                distance >= 0.0 && distance <= 1.0,
                "Distance should be between 0 and 1: {}",
                distance
            );
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
        TestSetup::create_completeness_file(
            &sandbox,
            "completeness_extra.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", 0.8),
                ("14412_3#84.contigs_velvet.fa.gz", 0.9),
                ("NonExistentGenome1", 0.5), // Extra genome - should be ignored
                ("NonExistentGenome2", 0.6), // Extra genome - should be ignored
                ("AnotherFakeGenome", 0.7),  // Extra genome - should be ignored
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
        let distances =
            TestSetup::parse_distances(&sandbox.file_string("distances_extra", TestDir::Output));
        assert!(distances.len() > 0, "Should have calculated some distances");

        // All distances should be valid
        for distance in distances {
            assert!(
                distance.is_finite(),
                "Distance should be finite: {}",
                distance
            );
            assert!(
                distance >= 0.0 && distance <= 1.0,
                "Distance should be between 0 and 1: {}",
                distance
            );
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
        TestSetup::create_completeness_file(
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
        let distances = TestSetup::parse_distances(
            &sandbox.file_string("precluster_distances", TestDir::Output),
        );
        assert!(distances.len() > 0, "Should have calculated some distances");

        // All distances should be valid
        for distance in distances {
            assert!(
                distance.is_finite(),
                "Distance should be finite: {}",
                distance
            );
            assert!(
                distance >= 0.0 && distance <= 1.0,
                "Distance should be between 0 and 1: {}",
                distance
            );
        }
    }

    #[test]
    fn test_completeness_correction_formula_exact() {
        let sandbox = TestSetup::setup();

        // Copy test files
        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");

        // Create sketch database
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("formula_test")
            .args(["-v", "-k", "31", "-s", "1000"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .assert()
            .success();

        // Create completeness file with known values for exact calculation
        let c1 = 0.8;
        let c2 = 0.9;
        TestSetup::create_completeness_file(
            &sandbox,
            "completeness_formula.txt",
            &[
                ("14412_3#82.contigs_velvet.fa.gz", c1),
                ("14412_3#84.contigs_velvet.fa.gz", c2),
            ],
        );

        // Calculate distances WITHOUT completeness correction
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("formula_test")
            .args(["-o", "distances_uncorrected"])
            .args(["-k", "31"])
            .arg("-v")
            .assert()
            .success();

        // Calculate distances WITH completeness correction
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("dist")
            .arg("formula_test")
            .args(["-o", "distances_corrected"])
            .arg("--completeness-file")
            .arg("completeness_formula.txt")
            .args(["-k", "31"])
            .arg("-v")
            .assert()
            .success();

        // Parse both distance outputs
        let uncorrected_distances = TestSetup::parse_distances(
            &sandbox.file_string("distances_uncorrected", TestDir::Output),
        );
        let corrected_distances = TestSetup::parse_distances(
            &sandbox.file_string("distances_corrected", TestDir::Output),
        );

        // Verify we have exactly one distance pair
        assert_eq!(
            uncorrected_distances.len(),
            1,
            "Should have exactly one distance pair"
        );
        assert_eq!(
            corrected_distances.len(),
            1,
            "Should have exactly one distance pair"
        );

        let uncorrected_jaccard = 1.0 - uncorrected_distances[0] as f32;
        let corrected_jaccard = 1.0 - corrected_distances[0] as f32;

        // Calculate expected corrected Jaccard using the formula:
        // corrected_jaccard = jaccard / (c1 * c2 / (c1 + c2 - c1 * c2))
        let correction_factor = (c1 * c2 / (c1 + c2 - c1 * c2)) as f32;
        let expected_corrected_jaccard = uncorrected_jaccard / correction_factor;
        let expected_corrected_distance = 1.0 - expected_corrected_jaccard;

        println!("Uncorrected Jaccard: {}", uncorrected_jaccard);
        println!("Corrected Jaccard: {}", corrected_jaccard);
        println!("Expected corrected Jaccard: {}", expected_corrected_jaccard);
        println!(
            "Correction factor (c1*c2/(c1+c2-c1*c2)): {}",
            correction_factor
        );
        println!("Uncorrected distance: {}", uncorrected_distances[0]);
        println!("Corrected distance: {}", corrected_distances[0]);
        println!(
            "Expected corrected distance: {}",
            expected_corrected_distance
        );

        // Verify the correction formula produces the expected result within floating-point tolerance
        // Use f32 precision throughout to match the actual computation precision
        assert_distance_with_tolerance(
            corrected_distances[0] as f32,
            expected_corrected_distance,
            1e-6,
        );

        // Also verify that correction was actually applied (distances should be different)
        assert!(
            (uncorrected_distances[0] - corrected_distances[0]).abs() > 1e-6,
            "Correction should have been applied - distances should be different"
        );

        // Verify correction makes distance smaller (Jaccard larger) as expected
        assert!(
            corrected_distances[0] < uncorrected_distances[0],
            "Completeness correction should reduce distance (increase Jaccard similarity)"
        );

        // Verify the corrected Jaccard doesn't exceed 1.0 (which would give negative distance)
        assert!(
            corrected_jaccard <= 1.0,
            "Corrected Jaccard should not exceed 1.0: {}",
            corrected_jaccard
        );
    }
}
