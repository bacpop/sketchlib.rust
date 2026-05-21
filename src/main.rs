#[cfg(not(target_family = "wasm"))]
use anyhow::Error;

#[cfg(not(target_family = "wasm"))]
fn main() -> Result<(), Error> {
    sketchlib::main()
}

#[cfg(target_family = "wasm")]
fn main() {
    sketchlib::main()
}
