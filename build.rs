#[cfg(feature = "gpu")]
extern crate cc;

#[cfg(feature = "gpu")]
fn main() {
    cc::Build::new()
        .cuda(true)
        .cudart("static")
        .std("c++17")
        .flag("-gencode").flag("arch=compute_86,code=sm_86")
        .flag("-ccbin").flag("/home/linuxbrew/.linuxbrew/bin/g++-11")
        .file("cuda/dist.cu")
        .compile("cuda_dist.a");

    /* Link CUDA Runtime (libcudart.so) */

    // Add link directory
    // - This path depends on where you install CUDA (i.e. depends on your Linux distribution)
    // - This should be set by `$LIBRARY_PATH`
    println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    println!("cargo:rustc-link-lib=cudart");

    /* Optional: Link CUDA Driver API (libcuda.so) */

    // println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64/stub");
    // println!("cargo:rustc-link-lib=cuda");
}

#[cfg(not(feature = "gpu"))]
fn main() {
}