#[cfg(feature = "3di")]
fn main() {
    // NB: I ran 'export PYO3_PYTHON=/usr/bin/python3'

    pyo3_build_config::add_extension_module_link_args();
    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,$CONDA_PREFIX/lib"
    );
    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,/Library/Developer/CommandLineTools/Library/Frameworks"
    );
}

#[cfg(not(feature = "3di"))]
fn main() {
}