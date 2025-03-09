# sketchlib.rust

<!-- badges: start -->
[![Cargo Build & Test](https://github.com/bacpop/sketchlib.rust/actions/workflows/ci.yml/badge.svg)](https://github.com/bacpop/sketchlib.rust/actions/workflows/ci.yml)
[![Clippy check](https://github.com/bacpop/sketchlib.rust/actions/workflows/clippy.yml/badge.svg)](https://github.com/bacpop/sketchlib.rust/actions/workflows/clippy.yml)
[![codecov](https://codecov.io/gh/bacpop/sketchlib.rust/graph/badge.svg?token=IBYPTT4J3F)](https://codecov.io/gh/bacpop/sketchlib.rust)
<!-- badges: end -->

## Description

This is a reimplementation of [pp-sketchlib](https://github.com/bacpop/pp-sketchlib)
in the rust language. This version is optimised for larger sample numbers, particularly
allowing subsets of samples to be compared.

Sketch databases have two files: `.skm` which is the metadata (samples names, base counts etc)
and `.skd` which is the actual sketch data.

## Usage
With all options we typically recommend using `-v` to see all progress during the run.

### Sketching

Using input fasta/fastq files, create a sketch database. Run `sketchlib sketch -h` to see the help.

- List .fasta files on the command line, or use `-f` to provide a file(s). Inputs can be gzipped or not, this is automatically detected.
From file, these are one line per sample listing:
    - One column (fasta input): file name, which is also used as the sample name
    - Two columns (fasta input): sample name and file name
    - Three columns (fastq input): sample name and two read files
- To set the k-mer size in the sketch database you can either give a list of sizes with `--k-vals`
or a sequence `--k-seq` with start,stop,step. e.g. `--k-seq 17,29,4` would sketch at k=17, 21, 25 and 29.
- Set the sketch size with `-s`. Typically 1000 is enough for species level resolution, 10000 for within-species/strain
resolution and 100000-1000000 for SNP level resolution.
- To sketch amino acid sequences use `--seq-type aa --concat-fasta` if you have the typical case
of each fasta file being a multifasta with many aa sequences. Each one will then be its own sample.
- You can also sketch structures with .pdb input, see 'Enabling PDB->3Di' below. This is experimental.

### Distances

To compute internal all-vs-all core and accessory distances use:
```
sketchlib dist db_name
```
Note the database names can be the prefix, or the full path to the .skm file. The output
is in pairwise 'long' format, which lists the upper triangle of the distance matrix row-by-row.

To calculate distances between two different sample sets, each in their own sketch database, use:
```
sketchlib dist db1 db2
```
For example, if you want to query distances of a new sample against an existing database,
first sketch the new sample with e.g. `sketchlib sketch -o db2 new_sample.fasta`, then
run the above command.

Modifiers:
- Use `-k` to calculate Jaccard distance at the given k. Otherwise the default is to
calculate across multiple k and output core and accessory distances.
- Use `--ani` with `-k` to transform the Jaccard distance into average nucleotide identity.
- Use `--subset` to provide a list of sample names to include in the distance calculations,
only these sample will be loaded from the `.skd` file.
- Use `-o` to write the distances to a file. The default it to write to stdout, so you can also
use `>` to redirect to a file (progress messages are written to stderr).
- Use `--knn` to only keep this many nearest neighbour distances. For very large databases
it may be useful to keep only ~50 distances. This makes the memory use manageable. This sparse output
can be used with e.g. [mandrake](https://github.com/bacpop/mandrake).

### Other operations

- `merge` joins two existing sketch databases.
- `append` sketches new input samples, and adds them to an existing database.
- `delete` removes samples from a sketch database.

## Enabling PDB->3Di
conda doesn't work, so make sure it is deactivated
```
export PYO3_PYTHON=python3
python3 -m venv 3di_venv
source 3di_venv/bin/activate
python3 -m pip install numpy biopython mini3di
cargo run -F 3di
export PYTHONPATH=${PYTHONPATH}:$(realpath ./)/3di_venv/lib/python3.12/site-packages
```