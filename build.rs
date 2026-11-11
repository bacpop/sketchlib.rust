#[cfg(feature = "gpu")]
fn main() {
    use std::{env, path::PathBuf};
    if env::var("CARGO_CFG_TARGET_ARCH").as_deref() == Ok("wasm32") {
        return;
    }
    for name in [
        "CUDA_TOOLKIT_PATH",
        "CUDA_HOME",
        "NVCC_CCBIN",
        "CUDAHOSTCXX",
    ] {
        println!("cargo:rerun-if-env-changed={name}");
    }
    println!("cargo:rerun-if-changed=src/gpu/cub_topk.cu");
    let root = PathBuf::from(
        env::var_os("CUDA_TOOLKIT_PATH")
            .or_else(|| env::var_os("CUDA_HOME"))
            .expect(
                "GPU builds require CUDA_TOOLKIT_PATH or CUDA_HOME pointing to CUDA 13.5 or newer",
            ),
    );
    assert!(
        root.join("bin/nvcc").is_file(),
        "CUDA toolkit has no bin/nvcc: {}",
        root.display()
    );
    let include = root.join("include");
    let cccl = if include.join("cccl/cub/version.cuh").is_file() {
        include.join("cccl")
    } else {
        include.clone()
    };
    let version = |path: PathBuf, name: &str| -> u32 {
        println!("cargo:rerun-if-changed={}", path.display());
        std::fs::read_to_string(&path)
            .ok()
            .and_then(|text| {
                text.lines().find_map(|line| {
                    let fields: Vec<_> = line.split_whitespace().collect();
                    (fields.len() >= 3 && fields[0] == "#define" && fields[1] == name)
                        .then(|| fields[2].parse().ok())
                        .flatten()
                })
            })
            .unwrap_or_else(|| {
                panic!(
                    "Cannot read {name} in {}; select CUDA 13.5+ with CCCL 3.5+",
                    path.display()
                )
            })
    };
    assert!(
        version(include.join("cuda.h"), "CUDA_VERSION") >= 13050,
        "GPU kNN requires CUDA Toolkit 13.5 or newer"
    );
    assert!(
        version(cccl.join("cub/version.cuh"), "CUB_VERSION") >= 300500,
        "GPU kNN requires CCCL/CUB 3.5 or newer"
    );
    // Build scripts run before this crate starts; select the same toolkit as
    // cuda-bindings instead of letting cc find a different nvcc on PATH.
    env::set_var("NVCC", root.join("bin/nvcc"));
    let mut build = cc::Build::new();
    build
        .cuda(true)
        .cudart("static")
        .pic(true)
        .ccbin(false)
        .include(cccl)
        .flag("-std=c++17")
        .flag("-gencode=arch=compute_75,code=[sm_75,compute_75]")
        .flag("-gencode=arch=compute_80,code=[sm_80,compute_80]")
        // Native Blackwell code also supports CUDA minor-version compatibility
        // when the installed driver cannot JIT this toolkit's newer PTX.
        .flag("-gencode=arch=compute_120,code=sm_120")
        .file("src/gpu/cub_topk.cu");
    if let Some(host) = env::var_os("NVCC_CCBIN").or_else(|| env::var_os("CUDAHOSTCXX")) {
        build.flag("-ccbin").flag(
            host.to_str()
                .expect("CUDA host compiler path must be UTF-8"),
        );
    }
    build.compile("sketchlib_cub_topk");
}

#[cfg(not(feature = "gpu"))]
fn main() {}
