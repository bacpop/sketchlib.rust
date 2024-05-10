# sketchlib.rust

## Description

This is a reimplementation of [pp-sketchlib](https://github.com/bacpop/pp-sketchlib)
in the rust language.

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