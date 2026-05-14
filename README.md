# GENESIS-CG-tool

This project helps you prepare topology and coordinate files for coarse-grained simulations in the MD package [GENESIS](https://github.com/genesis-release-r-ccs/genesis).

## Usage

Please use the following command to show the basic usage of this tool:
```bash
./src/aa_2_cg.jl -h
```

A full documentation can be found from the [wiki-page of this project](https://github.com/genesis-release-r-ccs/genesis_cg_tool/wiki).

## Dependencies

Please run the following command to install dependencies:
```bash
./opt/dependency.jl
```

## Tutorials

A bunch of tutorials on performing CG MD simulations with GENESIS and this project can be found from the [website of GENESIS](https://www.r-ccs.riken.jp/labs/cbrt/tutorials2022/).

## Branches

- The `main` branch is included as a part of official [GENESIS](https://github.com/genesis-release-r-ccs/genesis) release and the code only get updated when there is update in GENESIS.
- The `dev` branch is based on the `main` branch but contains newer features that will be officially released later.

## Citation

If this project helps you in your research, please consider citing this publication:

- Tan C, Jung J, Kobayashi C, Torre DUL, Takada S, Sugita Y (2022) Implementation of residue-level coarse-grained models in GENESIS for large-scale molecular dynamics simulations. PLoS Comput Biol 18(4): e1009578. https://doi.org/10.1371/journal.pcbi.1009578

