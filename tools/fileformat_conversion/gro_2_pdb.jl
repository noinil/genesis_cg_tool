#!/usr/bin/env julia

using Printf
using ArgParse

include("../../src/lib/constants.jl")
include("../../src/lib/topology.jl")
include("../../src/lib/conformation.jl")
include("../../src/lib/parser_top.jl")
include("../../src/lib/parser_crd.jl")
include("../../src/lib/parser_pdb.jl")
include("../../src/lib/parser_cif.jl")

function main(args)

    verbose = get(args, "verbose", false)

    top_filename = get(args, "top", "")
    crd_filename = get(args, "crd", "")
    pdb_filename = get(args, "output", "")
    out_mmcif    = get(args, "mmCIF", false)

    # ================================
    # Read in topology and coordinates
    # ================================

    mytop = read_grotop(top_filename)
    mycrd = read_grocrd(crd_filename)

    # =================
    # Write to PDB file
    # =================

    if length(pdb_filename) == 0
        system_name = crd_filename[1:end-4] * "_new"
    else
        system_name = pdb_filename[1:end-4]
    end

    if out_mmcif
        write_mmCIF(mytop, mycrd, system_name, args)
    else
        write_pdb(mytop, mycrd, system_name, args)
    end

    if verbose
        println("> converting from *.gro to *.pdb : DONE!")
    end
end

# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "--top", "-t"
        help     = "Topology file name (gromacs style)."
        required = true
        arg_type = String

        "--crd", "-c"
        help     = "Coordinate file name (gromacs style)."
        required = true
        arg_type = String

        "--output", "-o"
        help     = "Output PDB file name."
        arg_type = String
        default  = ""

        "--cgconnect"
        help     = "Output CONECTs in CG PDB."
        action   = :store_true

        "--mmCIF", "-M"
        help     = "Output as mmCIF/PDBx."
        action   = :store_true

        "--verbose"
        help     = "Output more information."
        action   = :store_true

        "--debug"
        help     = "DEBUG."
        action   = :store_true
    end

    return parse_args(s)
end



if abspath(PROGRAM_FILE) == @__FILE__

    args = parse_commandline()

    main(args)

end
