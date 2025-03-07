use std::{
    fs::{read_dir, File},
    io::{LineWriter, Write},
    path::{Path, PathBuf},
};

use assert_fs::{prelude::*, TempDir};
use predicates::prelude::*;

#[cfg(test)]

// Creates correct path for input/output files
static FILE_IN: &'static str = "tests/test_files_in";
static FILE_TEST: &'static str = "tests/test_results_correct";
static SYM_IN: &'static str = "input";
static SYM_TEST: &'static str = "correct";
static RFILE_NAME: &'static str = "file_list.txt";

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum TestDir {
    Input,
    Output,
    Correct,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum FxType {
    Fasta,
    Fastq,
}

pub struct TestSetup {
    wd: TempDir,
}

impl TestSetup {
    pub fn setup() -> Self {
        let wd = assert_fs::TempDir::new().unwrap();

        let file_in_path = Path::new(FILE_IN);
        let file_test_path = Path::new(FILE_TEST);

        // Attempt to create symlink for SYM_IN
        match wd
            .child(SYM_IN)
            .symlink_to_dir(file_in_path.canonicalize().unwrap_or_else(|e| {
                panic!(
                    "Failed to canonicalize FILE_IN path: {:?}. Error: {:?}",
                    file_in_path, e
                );
            })) {
            Ok(_) => println!("Successfully created symlink for {}", SYM_IN),
            Err(e) => println!("Failed to create symlink for {}: {:?}", SYM_IN, e),
        }

        // Attempt to create symlink for SYM_TEST
        match wd
            .child(SYM_TEST)
            .symlink_to_dir(file_test_path.canonicalize().unwrap_or_else(|e| {
                panic!(
                    "Failed to canonicalize FILE_TEST path: {:?}. Error: {:?}",
                    file_test_path, e
                );
            })) {
            Ok(_) => println!("Successfully created symlink for {}", SYM_TEST),
            Err(e) => println!("Failed to create symlink for {}: {:?}", SYM_TEST, e),
        }

        Self { wd }
    }

    pub fn get_wd(&self) -> String {
        self.wd.path().display().to_string()
    }

    pub fn file_path(&self, name: &str, file_type: TestDir) -> PathBuf {
        match file_type {
            TestDir::Input => {
                PathBuf::from(&format!("{}/{}/{}", self.wd.path().display(), SYM_IN, name))
            }
            TestDir::Output => PathBuf::from(&format!("{}/{}", self.wd.path().display(), name)),
            TestDir::Correct => PathBuf::from(&format!(
                "{}/{}/{}",
                self.wd.path().display(),
                SYM_TEST,
                name
            )),
        }
    }

    pub fn file_string(&self, name: &str, file_type: TestDir) -> String {
        self.file_path(name, file_type)
            .to_str()
            .expect("Could not unpack file path")
            .to_owned()
    }

    pub fn file_check(&self, name_out: &str, name_correct: &str) -> bool {
        let predicate_file = predicate::path::eq_file(self.wd.child(name_out).path());
        predicate_file.eval(self.file_path(name_correct, TestDir::Correct).as_path())
    }

    pub fn file_exists(&self, name_out: &str) -> bool {
        let predicate_fn = predicate::path::is_file();
        predicate_fn.eval(self.wd.child(name_out).path())
    }

    pub fn create_rfile(&self, filenames: &[&str]) -> &str {
        // Create an rfile in the tmp dir
        let mut rfile = LineWriter::new(
            File::create(format!("{}/{}", self.get_wd(), RFILE_NAME))
                .expect("Could not write rfile"),
        );
        for file in filenames {
            writeln!(
                rfile,
                "{}",
                &format!(
                    "{}\t{}",
                    file,
                    self.file_string(file, TestDir::Input),
                )
            )
            .unwrap();
        }
        RFILE_NAME
    }

    /// Create an rfile wtih two fastqs in the tmp dir
    pub fn create_fastq_rfile(&self, prefix: &str) -> &str {
        let mut rfile = LineWriter::new(
            File::create(format!("{}/{}", self.get_wd(), RFILE_NAME))
                .expect("Could not write rfile"),
        );
        for file_idx in 1..=2 {
            let sample_name = format!("{prefix}_{file_idx}");
            writeln!(
                rfile,
                "{}",
                &format!(
                    "{sample_name}\t{}\t{}",
                    self.file_string(&format!("{sample_name}_fwd.fastq.gz"), TestDir::Input),
                    self.file_string(&format!("{sample_name}_rev.fastq.gz"), TestDir::Input),
                )
            )
            .unwrap();
        }
        RFILE_NAME
    }
}
