# genesis_cg_input
This project helps you prepare topology and coordinate file for CG models in the molecular dynamics package Genesis.

## Usage

```bash
$ ./aa_pdb_2_cg_itp.jl -h
usage: aa_pdb_2_cg_itp.jl [-c RESPAC] [--aicg-scale AICG-SCALE]
                        [--3spn-param] [--pwmcos]
                        [--pwmcos-scale PWMCOS-SCALE]
                        [--pwmcos-shift PWMCOS-SHIFT] [--psf]
                        [--cgpdb] [-p PFM] [--patch PATCH]
                        [--show-sequence] [--debug] [-h] pdb

positional arguments:
  pdb                   PDB file name.

optional arguments:
  -c, --respac RESPAC   RESPAC protein charge distribution data.
                        (default: "")
  --aicg-scale AICG-SCALE
                        Scale AICG2+ local interactions: 0) average;
                        1) general (default). (type: Int64, default:
                        1)
  --3spn-param          Generate 3SPN.2C parameters from x3DNA
                        generated PDB structure.
  --pwmcos              Generate parameters for protein-DNA
                        sequence-specific interactions.
  --pwmcos-scale PWMCOS-SCALE
                        Energy scaling factor for PWMcos. (type:
                        Float64, default: 1.0)
  --pwmcos-shift PWMCOS-SHIFT
                        Energy shifting factor for PWMcos. (type:
                        Float64, default: 0.0)
  --psf                 Prepare PSF file.
  --cgpdb               Prepare CG PDB file.
  -p, --pfm PFM         Position frequency matrix file for protein-DNA
                        sequence-specific interactions. (default: "")
  --patch PATCH         Append (apply patch) to .itp file. (default:
                        "")
  --show-sequence       Show sequence of molecules in PDB.
  --debug               DEBUG.
  -h, --help            show this help message and exit


```

## Dependence

```bash
$ julia
julia> import Pkg
julia> Pkg.add("ProgressMeter")
julia> Pkg.add("ArgParse")
julia> Pkg.add("Formatting")
julia> Pkg.add("TextWrap")
julia> Pkg.add("Compat")
```
