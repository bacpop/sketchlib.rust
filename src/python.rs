//! Python bindings for sketchlib.

use std::collections::HashMap;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::api::{
    self, DistOptions, InvertedBuildOptions, InvertedPreclusterOptions, SketchOptions,
};
use crate::cli::{RetainUnmatched, DEFAULT_MINCOUNT, DEFAULT_MINQUAL};
use crate::hashing::{HashType, DEFAULT_LEVEL};

fn py_runtime_error(err: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

fn parse_retain_unmatched(value: Option<&str>) -> PyResult<Option<RetainUnmatched>> {
    match value {
        None => Ok(None),
        Some("singleton") => Ok(Some(RetainUnmatched::Singleton)),
        Some("bruteforce") => Ok(Some(RetainUnmatched::Bruteforce)),
        Some(other) => Err(PyValueError::new_err(format!(
            "retain_unmatched must be 'singleton', 'bruteforce', or None, got {other:?}"
        ))),
    }
}

/// Create a sketch database and return `(skm_path, skd_path)`.
#[pyfunction]
#[pyo3(signature = (
    output_prefix,
    kmer_length,
    file_list=None,
    seq_files=None,
    sketch_size=1000,
    threads=1,
    min_count=DEFAULT_MINCOUNT,
    min_qual=DEFAULT_MINQUAL,
    single_strand=false,
    quiet=true
))]
fn sketch(
    output_prefix: String,
    kmer_length: usize,
    file_list: Option<String>,
    seq_files: Option<Vec<String>>,
    sketch_size: u64,
    threads: usize,
    min_count: u16,
    min_qual: u8,
    single_strand: bool,
    quiet: bool,
) -> PyResult<(String, String)> {
    let paths = api::sketch(SketchOptions {
        seq_files,
        file_list,
        output: output_prefix,
        kmers: vec![kmer_length],
        sketch_size,
        seq_type: HashType::DNA,
        level: DEFAULT_LEVEL,
        concat_fasta: false,
        single_strand,
        min_count,
        min_qual,
        threads,
        quiet,
    })
    .map_err(py_runtime_error)?;
    Ok((paths.skm, paths.skd))
}

/// Calculate distances and return the output path, or `None` when writing to stdout.
#[pyfunction]
#[pyo3(signature = (
    ref_db,
    output=None,
    query_db=None,
    knn=None,
    subset=None,
    kmer=None,
    ani=false,
    threads=1,
    ref_completeness_file=None,
    query_completeness_file=None,
    completeness_cutoff=0.64,
    quiet=true
))]
fn dist(
    ref_db: String,
    output: Option<String>,
    query_db: Option<String>,
    knn: Option<usize>,
    subset: Option<String>,
    kmer: Option<usize>,
    ani: bool,
    threads: usize,
    ref_completeness_file: Option<String>,
    query_completeness_file: Option<String>,
    completeness_cutoff: f64,
    quiet: bool,
) -> PyResult<Option<String>> {
    api::dist(DistOptions {
        ref_db,
        query_db,
        output,
        knn,
        subset,
        kmer,
        ani,
        threads,
        ref_completeness_file,
        query_completeness_file,
        completeness_cutoff,
        quiet,
    })
    .map_err(py_runtime_error)
}

/// Build an inverted index and return `(ski_path, skq_path_or_none)`.
#[pyfunction]
#[pyo3(signature = (
    output_prefix,
    kmer_length,
    file_list=None,
    seq_files=None,
    sketch_size=1000,
    threads=1,
    write_skq=true,
    species_names=None,
    metadata=None,
    min_count=DEFAULT_MINCOUNT,
    min_qual=DEFAULT_MINQUAL,
    single_strand=false,
    quiet=true
))]
fn inverted_build(
    output_prefix: String,
    kmer_length: usize,
    file_list: Option<String>,
    seq_files: Option<Vec<String>>,
    sketch_size: u64,
    threads: usize,
    write_skq: bool,
    species_names: Option<String>,
    metadata: Option<String>,
    min_count: u16,
    min_qual: u8,
    single_strand: bool,
    quiet: bool,
) -> PyResult<(String, Option<String>)> {
    let paths = api::inverted_build(InvertedBuildOptions {
        seq_files,
        file_list,
        output: output_prefix,
        write_skq,
        species_names,
        metadata,
        sketch_size,
        kmer_length,
        single_strand,
        min_count,
        min_qual,
        threads,
        quiet,
    })
    .map_err(py_runtime_error)?;
    Ok((paths.ski, paths.skq))
}

/// Run inverted-index preclustering and return the output path, or `None` for stdout.
#[pyfunction]
#[pyo3(signature = (
    ski,
    skd,
    output=None,
    knn=50,
    ani=false,
    threads=1,
    ref_completeness_file=None,
    completeness_cutoff=0.64,
    retain_unmatched=None,
    quiet=true
))]
fn inverted_precluster(
    ski: String,
    skd: String,
    output: Option<String>,
    knn: usize,
    ani: bool,
    threads: usize,
    ref_completeness_file: Option<String>,
    completeness_cutoff: f64,
    retain_unmatched: Option<&str>,
    quiet: bool,
) -> PyResult<Option<String>> {
    api::inverted_precluster(InvertedPreclusterOptions {
        ski,
        skd,
        output,
        knn,
        ani,
        threads,
        ref_completeness_file,
        completeness_cutoff,
        retain_unmatched: parse_retain_unmatched(retain_unmatched)?,
        quiet,
    })
    .map_err(py_runtime_error)
}

/// Compute query-vs-reference distances and return the output path.
#[pyfunction]
#[pyo3(signature = (
    reference_skm,
    query_skm,
    output,
    kmer_length,
    knn,
    threads=1,
    ref_completeness_file=None,
    query_completeness_file=None,
    completeness_cutoff=0.64,
    quiet=true
))]
fn query_dist(
    reference_skm: String,
    query_skm: String,
    output: String,
    kmer_length: usize,
    knn: usize,
    threads: usize,
    ref_completeness_file: Option<String>,
    query_completeness_file: Option<String>,
    completeness_cutoff: f64,
    quiet: bool,
) -> PyResult<String> {
    api::query_dist(
        reference_skm,
        query_skm,
        output,
        kmer_length,
        knn,
        threads,
        ref_completeness_file,
        query_completeness_file,
        completeness_cutoff,
        quiet,
    )
    .map_err(py_runtime_error)
}

/// Return database metadata as a dict.
#[pyfunction]
#[pyo3(signature = (db_file, sample_info=false))]
fn db_info(py: Python<'_>, db_file: String, sample_info: bool) -> PyResult<PyObject> {
    let info = api::db_info(&db_file, sample_info).map_err(py_runtime_error)?;
    let dict = PyDict::new_bound(py);
    dict.set_item("text", info.text)?;
    dict.set_item("inverted", info.inverted)?;
    dict.set_item("sketch_size", info.sketch_size)?;
    dict.set_item("n_samples", info.n_samples)?;
    dict.set_item("kmers", info.kmers)?;
    Ok(dict.into())
}

/// Return available binding functions.
#[pyfunction]
fn functions() -> HashMap<&'static str, &'static str> {
    HashMap::from([
        ("sketch", "create .skm/.skd sketch files"),
        ("dist", "calculate distances from sketch files"),
        ("inverted_build", "create .ski/.skq inverted index files"),
        ("inverted_precluster", "calculate preclustered distances"),
        ("query_dist", "calculate query-vs-reference distances"),
        ("db_info", "read sketch database metadata"),
    ])
}

/// Python module entrypoint.
#[pymodule]
fn sketchlib(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sketch, m)?)?;
    m.add_function(wrap_pyfunction!(dist, m)?)?;
    m.add_function(wrap_pyfunction!(inverted_build, m)?)?;
    m.add_function(wrap_pyfunction!(inverted_precluster, m)?)?;
    m.add_function(wrap_pyfunction!(query_dist, m)?)?;
    m.add_function(wrap_pyfunction!(db_info, m)?)?;
    m.add_function(wrap_pyfunction!(functions, m)?)?;
    Ok(())
}
