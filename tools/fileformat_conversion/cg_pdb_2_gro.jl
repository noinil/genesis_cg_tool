#!/usr/bin/env julia

using Printf
using ArgParse

include("../../src/lib/constants.jl")
include("../../src/lib/molecule.jl")
include("../../src/lib/topology.jl")
include("../../src/lib/conformation.jl")
include("../../src/lib/file_io.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "pdb"
            help     = "CG PDB file name."
            required = true
            arg_type = String
    end

    return parse_args(s)
end


function pdb_2_gro(args)

    pdb_name = args["pdb"]
    mol_name = pdb_name[1:end-4]

    # ================
    # Read in PDB info
    # ================
    cgmol = read_cgPDB(pdb_name)

    # ==================
    # Output to gro file
    # ==================
    write_cg_grocrd(cgmol, mol_name, args)
    return 1
end


function main()
    args = parse_commandline()
    pdb_2_gro(args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
