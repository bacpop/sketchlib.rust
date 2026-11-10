# sketchlib.rust <img src='sketchlib.rust_logo.png' align="right" height="139" />

<!-- badges: start -->
[![Cargo Build & Test](https://github.com/bacpop/sketchlib.rust/actions/workflows/ci.yml/badge.svg)](https://github.com/bacpop/sketchlib.rust/actions/workflows/ci.yml)
[![Clippy check](https://github.com/bacpop/sketchlib.rust/actions/workflows/clippy.yml/badge.svg)](https://github.com/bacpop/sketchlib.rust/actions/workflows/clippy.yml)
[![docs.rs](https://img.shields.io/docsrs/sketchlib)](https://docs.rs/sketchlib)
[![codecov](https://codecov.io/gh/bacpop/sketchlib.rust/graph/badge.svg?token=IBYPTT4J3F)](https://codecov.io/gh/bacpop/sketchlib.rust)
[![Crates.io](https://img.shields.io/crates/v/sketchlib)](https://crates.io/crates/sketchlib)
[![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/bacpop/sketchlib.rust)](https://github.com/bacpop/sketchlib.rust/releases)
<!-- badges: end -->

## Description

This is a reimplementation and extension of [pp-sketchlib](https://github.com/bacpop/pp-sketchlib)
in the rust language. This version is optimised for larger sample numbers, particularly
allowing subsets of samples to be compared.

### News

- **Data format changed in v0.4.0**. If you have an older database (see `sketchlib info <db_prefix> | grep "sketch_version"`)
you should resketch samples with >v0.4 if possible. Otherwise you will need to use an older release available on crates.io.
- v0.2.0 was the first stable release.

## Documentation

See https://docs.rs/sketchlib

## Installation

Choose from:

1. Download a binary from the releases.
2. Use `cargo install sketchlib` or `cargo add sketchlib`.
3. Use conda install -c bioconda sketchlib.
4. Build from source

For 2) or 4) you must have the rust toolchain installed.

### OS X users

If you have an M1-4 (arm64) Mac, we aren't currently automatically building binaries, so would recommend either option 2) or 3) for best performance.

If you get a message saying the binary isn't signed by Apple and can't be run, use the following command to bypass this:

```
xattr -d "com.apple.quarantine" ./sketchlib
```

### Build from source

1. Clone the repository with git clone.
2. Run `cargo install --path .` or `RUSTFLAGS="-C target-cpu=native" cargo install --path .` to optimise for your machine.

GPU builds use `cargo-oxide` and require CUDA Toolkit 13.5 or newer, including
CCCL/CUB 3.5+ and `nvcc`. Select the toolkit explicitly, for example:

```bash
CUDA_HOME=/usr/local/cuda-13.5 cargo oxide build
```

`CUDA_TOOLKIT_PATH` takes precedence over `CUDA_HOME` when both are set. Cargo
compiles the CUB kNN wrapper with `nvcc` and links the CUDA runtime statically;
no separate CUB installation or manual linking step is needed. The `gpu` feature
is currently enabled by default in this branch.

## Citation

Please cite:

von Wachsmann J, Lorenz LJ, Russell MJ, Gurbich TA, Rodríguez-Bouza V, Horsfield ST, Lees JA, Finn RD (2026).\
Rapid and consistent clustering of millions of genomes highlights the diversity of prokaryotic life.\
*bioRxiv*. 

https://doi.org/10.64898/2025.12.30.695181

Lees JA, Tonkin-Hill G, Yang Z, Corander J.\
Mandrake: visualizing microbial population structure by embedding millions of genomes into a low-dimensional representation.\
*Philosophical Transactions of The Royal Society B*. 2022;377: 20210237.

https://doi.org/10.1098/rstb.2021.0237

We rely on algorithms from:

*bindash* (written by XiaoFei Zhao):\
Zhao, X. BinDash, software for fast genome distance estimation on a typical personal laptop.\
*Bioinformatics* **35**:671–673 (2019).\
doi:[10.1093/bioinformatics/bty651](https://dx.doi.org/10.1093/bioinformatics/bty651)

*ntHash* (written by Hamid Mohamadi):\
Mohamadi, H., Chu, J., Vandervalk, B. P. & Birol, I. ntHash: recursive nucleotide hashing.\
*Bioinformatics* **32**:3492–3494 (2016).\
doi:[10.1093/bioinformatics/btw397](https://dx.doi.org/10.1093/bioinformatics/btw397)
