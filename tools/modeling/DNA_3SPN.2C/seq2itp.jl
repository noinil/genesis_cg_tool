#!/usr/bin/env julia

using Printf
using ArgParse

include("../../../src/lib/constants.jl")
include("../../../src/lib/molecule.jl")
include("../../../src/lib/biomath.jl")
include("../../../src/lib/conformation.jl")
include("../../../src/lib/topology.jl")
include("../../../src/lib/parser_top.jl")
include("../../../src/lib/coarse_graining_subroutines.jl")
include("../../../src/lib/coarse_graining.jl")

base_type = ['A', 'T', 'C', 'G']

# ====================
# Read in DNA sequence
# ====================
function read_DNA_sequence(file_name)
    dna_seq = ""
    for line in eachline(file_name)
        if line[1] == '>'
            if length(dna_seq) > 0
                break
            end
            continue
        end
        seq = strip(line)
        for b in seq
            if ! in(b, "ACGT")
                error("Wrong DNA sequence!")
            end
        end
        dna_seq *= seq
    end
    return dna_seq
end

function get_complementary_seq(seq)
    BP_DICT = Dict(
        'A' => "T", 'T' => "A",
        'C' => "G", 'G' => "C"
    )

    comp_seq = ""
    for b ∈ seq
        comp_seq = BP_DICT[b] * comp_seq
    end

    return comp_seq
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

        "--complementary", "-c"
        help = "Generate complementary strand."
        action = :store_true

        "--5P"
        help = "Starting from Phosphate at the 5-end."
        action = :store_true

        "--circular"
        help = "Build a circular or infinite DNA."
        action = :store_true

        "--output-name", "-o"
        help = "Specify the system name for output."
        arg_type = String
        default = ""

        "--verbose", "-v"
        help = "Output more information."
        action = :store_true

    end

    return parse_args(s)
end

# ====
# Main
# ====
function gen_3spn_itp_from_DNA_seq()

    # -----------------------
    # Parse command line args
    # -----------------------
    args        = parse_commandline()
    verbose     = get(args, "verbose", false)
    comp_strand = get(args, "complementary", false)
    phos_5p     = get(args, "5P", false)
    circ_dna    = get(args, "circular", false)
    mol_name    = get(args, "output-name", "bdna")

    # ===============================
    # Read and make list of sequences
    # ===============================
    sequences = []
    dna_seq_name = args["fasta"]
    seq_DNA_a = read_DNA_sequence(dna_seq_name)
    push!(sequences, seq_DNA_a)

    if comp_strand
        seq_DNA_b = get_complementary_seq(seq_DNA_a)
        push!(sequences, seq_DNA_b)
    end

    # =====================
    # Generate DNA topology
    # =====================

    aa_num_atom  = sum([phos_5p ? length(c) * 4 : length(c) * 4 - 1 for c in sequences])

    aa_atom_name = fill("    ",       aa_num_atom)
    aa_coor      = zeros(Float64, (3, aa_num_atom))
    aa_residues  = []
    aa_chains    = []

    # ---------------------
    # Create fake molecules
    # ---------------------
    i_atom = 0
    i_resi = 0
    for i_chain in 1:length(sequences)
        dnaseq = sequences[i_chain]
        tmp_resi_indices = []
        for (ib, b) ∈ enumerate(dnaseq)
            i_resi += 1
            tmp_atom_indices = []
            # phosphate
            if ib > 1 || phos_5p
                i_atom += 1
                aa_coor[:, i_atom] = [1.0, 0.0, 0.0] # fake coordinates...
                aa_atom_name[i_atom] = "P"
                push!(tmp_atom_indices, i_atom)
            end
            # sugar
            i_atom += 1
            aa_coor[:, i_atom] = [0.0, 1.0, 0.0] # fake coordinates...
            aa_atom_name[i_atom] = "C3'"
            push!(tmp_atom_indices, i_atom)
            # base
            i_atom += 1
            aa_coor[:, i_atom] = [0.0, 0.0, 1.0] # fake coordinates...
            aa_atom_name[i_atom] = "C5"
            push!(tmp_atom_indices, i_atom)
            # O3'
            i_atom += 1
            aa_coor[:, i_atom] = [1.0, 1.0, 1.0] # fake coordinates...
            aa_atom_name[i_atom] = "O3'"
            push!(tmp_atom_indices, i_atom)
            # add new residue
            push!(aa_residues, AAResidue("D"*b, tmp_atom_indices))
            # add residue index to chain
            push!(tmp_resi_indices, i_resi)
        end
        chain_id = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i_chain]
        segment_id = dna_seq_name[1:3] * "_" * chain_id
        push!(aa_chains, AAChain(chain_id, segment_id, MOL_DNA, tmp_resi_indices))
    end

    dna_molecule = AAMolecule(aa_atom_name, aa_coor, aa_residues, aa_chains)

    args["3spn-param"] = 2
    args["3spn-use-5-phos"] = phos_5p
    args["3spn-circular"] = circ_dna
    force_field = ForceFieldCG(1, 1, 1, 0, 0, 0)
    cg_top, cg_conf = coarse_graining(dna_molecule, force_field, args)

    if length(mol_name) == 0
        mol_name = dna_seq_name[1:3] * "_bdna"
    end
    write_grotop(cg_top, mol_name, args)

end

if abspath(PROGRAM_FILE) == @__FILE__
    gen_3spn_itp_from_DNA_seq()
end
