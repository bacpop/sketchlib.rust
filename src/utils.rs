use std::error::Error;
use std::fs::File;
use std::io::copy;

/// Some helper functions
pub fn strip_sketch_extension(file_name: &str) -> &str {
    if file_name.ends_with(".skm") || file_name.ends_with(".skd") {
        &file_name[..file_name.len() - 4]
    } else {
        file_name
    }
}

pub fn save_sketch_data(db1: &str, db2: &str, str_output: &str) -> Result<(), anyhow::Error> {
    let mut output_file = File::create(str_output)?;
    // Open and copy the contents of the first input file
    let mut db_sketch1 = File::open(db1)?;
    copy(&mut db_sketch1, &mut output_file)?;

    // Open and copy the contents of the second input file
    let mut db_sketch2 = File::open(db2)?;
    copy(&mut db_sketch2, &mut output_file)?;

    log::info!("Databases merged successfully to {}", str_output);
    Ok(())
}
