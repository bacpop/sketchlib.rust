use anyhow::Error;

#[cfg(feature = "3di")]
use pyo3::prelude::*;

#[cfg(feature = "3di")]
fn main() -> Result<(), Error> {
    pyo3::prepare_freethreaded_python();

    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/python/3di_convert.py"));
    let from_python = Python::with_gil(|py| -> PyResult<Py<PyAny>> {
        PyModule::import_bound(py, "mini3di");
        PyModule::import_bound(py, "biopython");
        let app: Py<PyAny> = PyModule::from_code_bound(py, py_app, "", "")?
            .getattr("pdb_to_3di")?
            .into();
        app.call0(py)
    });

    println!("py: {}", from_python?);
    Ok(())
}

#[cfg(not(feature = "3di"))]
fn main() -> Result<(), Error> {
    sketchlib::main();
    Ok(())
}
