//! Support for .pdb files and the 3di alphabet
#[cfg(feature = "3di")]
use crate::io::InputFastx;
#[cfg(feature = "3di")]
use anyhow::Error;

#[cfg(feature = "3di")]
use indicatif::{ProgressIterator, ProgressStyle};
#[cfg(feature = "3di")]
use pyo3::prelude::*;

#[cfg(feature = "3di")]
/// Uses python library to convert pdb files into 1D 3di representation
pub fn pdb_to_3di(input_files: &[InputFastx]) -> Result<Vec<String>, Error> {
    pyo3::prepare_freethreaded_python();
    let py_file = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/python_mini3di/3di_convert.py"
    ));

    let mut struct_strings = Vec::with_capacity(input_files.len());
    let bar_style =
        ProgressStyle::with_template("{human_pos}/{human_len} {bar:80.cyan/blue} eta:{eta}")
            .unwrap();
    Python::with_gil(|py| {
        let converter_fun: Py<PyAny> = PyModule::from_code_bound(py, py_file, "", "")
            .unwrap()
            .getattr("pdb_to_3di")
            .unwrap()
            .into();
        //         println!("n f1 f2: {:?} {:?} {:?}", input_files[0], input_files[1], input_files[2]);

        for ftup in input_files.iter().progress_with_style(bar_style) {
            //             println!("{:?}", ftup);
            if ftup.1.len() == 1 {
                let struct_string: Py<PyAny> = converter_fun
                    .call1(py, (ftup.0.clone(), ftup.1[0].clone()))
                    .unwrap();
                struct_strings.push(struct_string.extract::<String>(py).unwrap());
            } else {
                let mut struct_string = "".to_owned();
                for file in ftup.1.iter() {
                    let struct_string_tmp: Py<PyAny> =
                        converter_fun.call1(py, (ftup.0.clone(), file)).unwrap();
                    struct_string.push_str(",");
                    struct_string
                        .push_str(struct_string_tmp.extract::<String>(py).unwrap().as_str());
                }
                struct_strings.push(struct_string);
            }
        }
    });

    //     println!("Finished!");

    Ok(struct_strings)
}

// This shouldn't be needed, but I am putting it here just in case
/*
#[cfg(not(feature = "3di"))]
pub fn pdb_to_3di(_input_files: &[InputFastx]) -> Result<Vec<String>, Error> {
    unimplemented!("This executable was not compiled with the 3di feature")
}
*/
