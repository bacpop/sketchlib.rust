//! Helper functions for file manipulation

use std::fs::File;
use std::io::copy;

use indicatif::{ProgressBar, ProgressStyle};

/// Removes .skm or .skd, if they exist, at the end of a filename
pub fn strip_sketch_extension(file_name: &str) -> &str {
    if file_name.ends_with(".skm") || file_name.ends_with(".skd") {
        &file_name[..file_name.len() - 4]
    } else {
        file_name
    }
}

/// Concatenates two .skd files using [`std::io::copy`]
pub fn save_sketch_data(
    db1_prefix: &str,
    db2_prefix: &str,
    str_output: &str,
) -> Result<(), anyhow::Error> {
    let mut output_file = File::create(format!("{}.skd", str_output))?;
    // Open and copy the contents of the first input file
    let mut db_sketch1 = File::open(format!("{}.skd", db1_prefix))?;
    copy(&mut db_sketch1, &mut output_file)?;

    // Open and copy the contents of the second input file
    let mut db_sketch2 = File::open(format!("{}.skd", db2_prefix))?;
    copy(&mut db_sketch2, &mut output_file)?;

    Ok(())
}

/// Create a progress bar for use in iterators
pub fn get_progress_bar(length: usize, percent: bool, quiet: bool) -> ProgressBar {
    if quiet {
        ProgressBar::hidden()
    } else {
        let style = if percent {
            ProgressStyle::with_template("{percent}% {bar:80.cyan/blue} eta:{eta}").unwrap()
        } else {
            ProgressStyle::with_template("{human_pos}/{human_len} {bar:80.cyan/blue} eta:{eta}")
                .unwrap()
        };
        ProgressBar::new(length as u64).with_style(style)
    }
}
