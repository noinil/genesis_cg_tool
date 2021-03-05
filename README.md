# genesis_cg_julia
This project helps you prepare topology and coordinate files for coarse-grained simulations in the MD package Genesis.

## Usage

```bash
./src/aa_2_cg.jl -h
usage: aa_2_cg.jl [--force-field-protein FORCE-FIELD-PROTEIN]
                  [--force-field-DNA FORCE-FIELD-DNA]
                  [--force-field-RNA FORCE-FIELD-RNA] [-f CONFIG]
                  [--CCGO-contact-scale CCGO-CONTACT-SCALE]
                  [-c RESPAC] [--aicg-scale AICG-SCALE]
                  [--use-safe-dihedral USE-SAFE-DIHEDRAL]
                  [--3spn-param _3SPN-PARAM] [--3spn-use-5-phos]
                  [--protein-DNA-Go] [--pwmcos]
                  [--pwmcos-scale PWMCOS-SCALE]
                  [--pwmcos-shift PWMCOS-SHIFT] [--pwmcos-ns]
                  [--pwmcos-ns-ene PWMCOS-NS-ENE] [--psf] [--cgpdb]
                  [--cgconnect] [--cgRNA-phosphate-Go] [-p PFM]
                  [--test-local-only] [--patch PATCH]
                  [--show-sequence] [--output-name OUTPUT-NAME] [-v]
                  [--log] [--debug] [-h] pdb

positional arguments:
  pdb                   PDB file name.

optional arguments:
  --force-field-protein FORCE-FIELD-PROTEIN
                        Force field for protein. (default: "AICG2+")
  --force-field-DNA FORCE-FIELD-DNA
                        Force field for DNA. (default: "3SPN.2C")
  --force-field-RNA FORCE-FIELD-RNA
                        Force field for RNA. (default: "HT")
  -f, --config CONFIG   Force field configuration details. (default:
                        "")
  --CCGO-contact-scale CCGO-CONTACT-SCALE
                        Scaling native contact interaction
                        coefficient. (type: Float64, default: 1.0)
  -c, --respac RESPAC   RESPAC protein charge distribution data.
                        (default: "")
  --aicg-scale AICG-SCALE
                        Scale AICG2+ local interactions: 0) average;
                        1) general (default). (type: Int64, default:
                        1)
  --use-safe-dihedral USE-SAFE-DIHEDRAL
                        Safe dih potential: 0) do nothing (default);
                        1) cos^2(kθ) type; 2) remove dih w/ large
                        angles; 3) sin^3(kθ) type. (type: Int64,
                        default: 0)
  --3spn-param _3SPN-PARAM
                        Generate 3SPN.2C parameters from x3DNA
                        generated PDB structure. (type: Int64,
                        default: 0)
  --3spn-use-5-phos     Generate 3SPN.2C parameters with 5-phosphate.
  --protein-DNA-Go      Generate parameters for protein-DNA Go-like
                        contact interactions.
  --pwmcos              Generate parameters for protein-DNA
                        sequence-specific interactions.
  --pwmcos-scale PWMCOS-SCALE
                        Energy scaling factor for PWMcos. (type:
                        Float64, default: 1.0)
  --pwmcos-shift PWMCOS-SHIFT
                        Energy shifting factor for PWMcos. (type:
                        Float64, default: 0.0)
  --pwmcos-ns           Generate parameters for protein-DNA
                        sequence-NON-specific interactions.
  --pwmcos-ns-ene PWMCOS-NS-ENE
                        Interaction strength for PWMcos-ns
                        (hydrogen-bond). (type: Float64, default:
                        -1.0)
  --psf                 Prepare PSF file.
  --cgpdb               Prepare CG PDB file.
  --cgconnect           Prepare CG PDB file with CONECTed bonds.
  --cgRNA-phosphate-Go  Include phosphate in Go-type contact
                        interactions.
  -p, --pfm PFM         Position frequency matrix file for protein-DNA
                        sequence-specific interactions. (default: "")
  --test-local-only     TEST: only generate local interaction
                        parameters.
  --patch PATCH         Append (apply patch) to .itp file. (default:
                        "")
  --show-sequence       Show sequence of molecules in PDB.
  --output-name OUTPUT-NAME
                        Specify the system name for output. (default:
                        "")
  -v, --verbose         Output more information.
  --log                 Output information to log file.
  --debug               DEBUG.
  -h, --help            show this help message and exit

```

## Dependence

```bash
$ julia
julia> import Pkg
julia> Pkg.add("ArgParse")
```
