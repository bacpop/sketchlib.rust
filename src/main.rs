#[cfg(not(target_arch = "wasm32"))]
use anyhow::Error;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Error> {
    sketchlib::main()
}

#[cfg(target_arch = "wasm32")]
fn main() {
    sketchlib::main()
}
