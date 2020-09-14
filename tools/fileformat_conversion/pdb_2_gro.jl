#!/usr/bin/env julia

using Printf
using ArgParse

include("../../src/lib/constants.jl")
include("../../src/lib/topology.jl")
include("../../src/lib/molecule.jl")
include("../../src/lib/conformation.jl")
include("../../src/lib/parser_top.jl")
include("../../src/lib/parser_crd.jl")
include("../../src/lib/parser_pdb.jl")

function main(args)

    verbose = get(args, "verbose", false)

    top_filename = get(args, "top", "")
    pdb_filename = get(args, "pdb", "")
    crd_filename = get(args, "output", "")
    pdb_no_ter   = get(args, "pdb-noTER", false)

    # =============
    # Read topology
    # =============
    mytop = read_grotop(top_filename)

    # =========================
    # Read coordinates from PDB
    # =========================
    if pdb_no_ter
        aa_coor = zeros(Float64, (3, mytop.num_atom))
        i_atom  = 0
        for line in eachline(pdb_filename)
            if startswith(line, "ATOM")
                i_atom += 1
                new_atom_data      = parse_PDB_line(rpad(line, 80))
                aa_coor[1, i_atom] = new_atom_data.coor_x
                aa_coor[2, i_atom] = new_atom_data.coor_y
                aa_coor[3, i_atom] = new_atom_data.coor_z
            end
        end
        new_conf = Conformation(mytop.num_atom, aa_coor)
    else
        new_molecule = read_PDB(pdb_filename)
        num_atoms = length(new_molecule.atom_coors[1, :])
        new_conf = Conformation(num_atoms, new_molecule.atom_coors)
    end

    if length(crd_filename) == 0
        system_name = top_filename[1:end-4]
    else
        system_name = crd_filename[1:end-4]
    end

    write_grocrd(mytop, new_conf, system_name, args)

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

        "--pdb", "-p"
        help     = "PDB file name."
        required = true
        arg_type = String

        "--pdb-noTER"
        help     = "PDB file does not have TER lines."
        action   = :store_true

        "--output", "-o"
        help     = "Output file name."
        arg_type = String
        default  = ""

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
