use anyhow::Error;
use crate::io::InputFastx;

#[cfg(feature = "3di")]
use pyo3::prelude::*;

#[cfg(feature = "3di")]
pub fn pdb_to_3di(input_files: &[InputFastx]) -> Result<Vec<String>, Error> {
    pyo3::prepare_freethreaded_python();
    let py_file = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/python_mini3di/3di_convert.py"));

    let mut struct_strings = Vec::with_capacity(input_files.len());
    let bar_style =
    ProgressStyle::with_template("{human_pos}/{human_len} {bar:80.cyan/blue} eta:{eta}")
        .unwrap();
    Python::with_gil(|py| {
        let converter_fun: Py<PyAny> = PyModule::from_code_bound(py, py_file, "", "")?
            .getattr("pdb_to_3di")?
            .into();
        input_files
                .iter()
                .progress_with_style(bar_style)
                .map(|(name, fastx1, fastx2)| {
                    let struct_string: Py<PyAny> = converter_fun.call1(py, (name, fastx1));
                    struct_strings.push(struct_string.extract()?)
        });
    });

    Ok((struct_strings))
}

// This shouldn't be needed, but I am putting it here just in case
/*
#[cfg(not(feature = "3di"))]
pub fn pdb_to_3di(_input_files: &[InputFastx]) -> Result<Vec<String>, Error> {
    unimplemented!("This executable was not compiled with the 3di feature")
}
*/