#!/usr/bin/env julia

using Printf
using ArgParse
using Formatting

base_pair_parms = Dict(
    #         a-b,  shear, stretch, stagger, buckle, propeller, opening,  
    'A' => ("A-T",   0.07,   -0.19,    0.07,    1.8,     -15.0,     1.5),
    'T' => ("T-A",  -0.07,   -0.19,    0.07,   -1.8,     -15.0,     1.5),
    'C' => ("C-G",   0.16,   -0.17,    0.15,   -4.9,      -8.7,    -0.6),
    'G' => ("G-C",  -0.16,   -0.17,    0.15,    4.9,      -8.7,    -0.6)
)

base_step_parms = Dict(
    #         shift, slide,  rise,  tilt,  roll,  twist
    "AA" => ( -0.05, -0.21,  3.27, -1.84,  0.76,  35.31),
    "AT" => (  0.00, -0.56,  3.39,  0.00, -1.39,  31.21),
    "AC" => (  0.21, -0.54,  3.39, -0.64,  0.91,  31.52),
    "AG" => (  0.12, -0.27,  3.38, -1.48,  3.15,  33.05),
    "TA" => (  0.00,  0.03,  3.34,  0.00,  5.25,  36.20),
    "TT" => (  0.05, -0.21,  3.27,  1.84,  0.91,  35.31),
    "TC" => (  0.27, -0.03,  3.35,  1.52,  3.87,  34.80),
    "TG" => (  0.16,  0.18,  3.38,  0.05,  5.95,  35.02),
    "CA" => ( -0.16,  0.18,  3.38, -0.05,  5.95,  35.02),
    "CT" => ( -0.12, -0.27,  3.38,  1.48,  3.15,  33.05),
    "CC" => (  0.02, -0.47,  3.28,  0.40,  3.86,  33.17),
    "CG" => (  0.00,  0.57,  3.49,  0.00,  4.29,  34.38),
    "GA" => ( -0.27, -0.03,  3.35, -1.52,  3.87,  34.80),
    "GT" => ( -0.21, -0.54,  3.39,  0.64,  0.91,  31.52),
    "GC" => (  0.00, -0.07,  3.38,  0.00,  0.67,  34.38),
    "GG" => ( -0.02, -0.47,  3.28, -0.40,  3.86,  33.17),
    "00" => (  0.00,  0.00,  0.00,  0.00,  0.00,   0.00),
)

base_type = ['A', 'T', 'C', 'G']

# ====================
# Read in DNA sequence
# ====================
function read_DNA_sequence(file_name)
    for line in eachline(file_name)
        if line[1] == '>'
            continue
        end
        seq = strip(line)
        for b in seq
            if ! in(b, "ACGT")
                error("Wrong DNA sequence!")
            end
        end
        return seq
    end
end

# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "fasta"
        help     = "DNA sequence file name."
        required = true
        arg_type = String

    end

    return parse_args(s)
end

# ====
# Main
# ====
function main()

    args = parse_commandline()

    # read in DNA sequence:
    seq_DNA = read_DNA_sequence(args["fasta"])
    len_DNA = length(seq_DNA)
    println(" DNA lenth: ", len_DNA, " bp.")

    # prepare the output file:
    out_file = open("dna2c.curv", "w")
    printfmt(out_file, "{:>4d} # bps \n", len_DNA)
    printfmt(out_file, "   0 # local base-pairing and base-step parameters \n")
    printfmt(out_file, "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist\n")
    parm_line = "{1}  {2:>9.3f} {3:>9.3f} {4:>9.3f} {5:>9.3f} {6:>9.3f} {7:>9.3f} {8:>9.3f} {9:>9.3f} {10:>9.3f} {11:>9.3f} {12:>9.3f} {13:>9.3f}\n"

    # parameter output loop
    for ( i, b ) in enumerate(seq_DNA)
        bp_parm = base_pair_parms[b]
        base_step = i == 1 ? "00" : seq_DNA[i-1:i]
        bs_parm = base_step_parms[base_step]
        printfmt(out_file, parm_line,
                 bp_parm[1],
                 bp_parm[2],
                 bp_parm[3],
                 bp_parm[4],
                 bp_parm[5],
                 bp_parm[6],
                 bp_parm[7],
                 bs_parm[1],
                 bs_parm[2],
                 bs_parm[3],
                 bs_parm[4],
                 bs_parm[5],
                 bs_parm[6]
                 )
    end

    close(out_file)

end

main()
