#!/usr/bin/env julia

using ArgParse
using Formatting

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "pdb"
            help     = "PDB file name."
            required = true
            arg_type = String
    end

    return parse_args(s)
end


function pdb_2_gro(pdb_name)
    # HEAD: time in the unit of ps
    GRO_HEAD_STR  = "Gro from CG PDB, t= 0.000 \n"
    # ATOM NUM: free format int
    GRO_ATOM_NUM  = "{1:>12d}\n"
    # XYZ: in the unit of nm!!!
    # "{res_num:>5d}{res_name:>5}{atm_name:>5}{atm_num:>5d}{x:>8.4f}{y:>8.4f}{z:>8.4f}{vx:>8.4f}{vy:>8.4f}{vz:>8.4f} \n"
    GRO_ATOM_LINE = "{1:>5d}{2:>5}{3:>5}{4:>5d}{5:>8.4f}{6:>8.4f}{7:>8.4f}{8:>8.4f}{9:>8.4f}{10:>8.4f} \n"
    GRO_BOX_LINE  = "{1:>15.4f}{2:>15.4f}{3:>15.4f} \n\n"


    # ================
    # Read in PDB info
    # ================
    pdb_lines = []
    for line in eachline(pdb_name)
        if startswith(line, "ATOM")
            push!(pdb_lines, rpad(line, 80))
        end
    end


    # ==================
    # Output to gro file
    # ==================
    gro_name = pdb_name[1 : end - 4] * ".gro"

    open(gro_name, "w") do gro_file
        write(gro_file, GRO_HEAD_STR)
        printfmt(gro_file, GRO_ATOM_NUM, length(pdb_lines))
        for line in pdb_lines
            atom_serial    = parse(Int, line[7:11])
            atom_name      = strip(line[13:16])
            residue_name   = strip(line[18:21])
            chain_id       = line[22]
            residue_serial = parse(Int, line[23:26])
            coor_x         = parse(Float64, line[31:38])
            coor_y         = parse(Float64, line[39:46])
            coor_z         = parse(Float64, line[47:54])
            printfmt(gro_file, GRO_ATOM_LINE,
                     residue_serial,
                     residue_name,
                     atom_name,
                     atom_serial,
                     coor_x / 10.0,
                     coor_y / 10.0,
                     coor_z / 10.0,
                     0.0,
                     0.0,
                     0.0)
        end
        printfmt(gro_file, GRO_BOX_LINE, 0.0, 0.0, 0.0)
    end

    return 1
end


function main()
    args = parse_commandline()
    pdb_2_gro(args["pdb"])
end

main()
