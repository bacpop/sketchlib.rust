use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn inverted_sketch() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
            .arg("-o")
            .arg("inverted")
            .args(["-v", "-k", "31"])
            .arg("14412_3#82.contigs_velvet.fa.gz")
            .arg("14412_3#84.contigs_velvet.fa.gz")
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("inverted.ski"));

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("inverted.ski")
            .arg("-v")
            .assert()
            .stdout_matches_path(
                sandbox.file_string("inverted_sketch_info.stdout", TestDir::Correct),
            );

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("inverted.ski")
            .arg("-v")
            .assert()
            .stdout_matches_path(
                sandbox.file_string("inverted_sketch_full_info.stdout", TestDir::Correct),
            );
    }

    #[test]
    fn inverted_reorder() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("14412_3#82.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("14412_3#84.contigs_velvet.fa.gz");
        sandbox.copy_input_file_to_wd("R6.fa.gz");
        sandbox.copy_input_file_to_wd("TIGR4.fa.gz");
        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("inverted")
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

        Command::new(cargo_bin("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("inverted.ski")
            .arg("--sample-info")
            .arg("-v")
            .assert()
            .stdout_matches_path(
                sandbox.file_string("inverted_sketch_info_reorder.stdout", TestDir::Correct),
            );
    }
}
