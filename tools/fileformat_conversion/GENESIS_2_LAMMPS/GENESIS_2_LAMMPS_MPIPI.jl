#!/usr/bin/env julia

include("../../../src/lib/gcj.jl")
using ArgParse

function g2l(args)

    # --------------
    # load arguments
    # --------------
    g_top_fname = get(args, "top", "")
    g_crd_fname = get(args, "crd", "")
    l_out_fname = get(args, "outstructure", "")
    κ           = get(args, "kappa", 0.1299)
    lj_cut      = get(args, "LJcutoff", 25.0)
    ele_cut     = get(args, "ELEcutoff", 35.0)
    PBCxlo      = get(args, "PBCxlo", -200.0)
    PBCylo      = get(args, "PBCylo", -200.0)
    PBCzlo      = get(args, "PBCzlo", -200.0)
    PBCxhi      = get(args, "PBCxhi",  200.0)
    PBCyhi      = get(args, "PBCyhi",  200.0)
    PBCzhi      = get(args, "PBCzhi",  200.0)
    verbose     = get(args, "verbose", "")

    if l_out_fname == ""
        L_struct_fname = "OUT_" * g_top_fname[1:end-4] * "_struct.dat"
    else
        L_struct_fname = "OUT_" * l_out_fname * "_struct.dat"
    end
    L_Struct_file = open(L_struct_fname, "w")

    # read in topology
    g_top = read_grotop(g_top_fname)
    g_crd = read_grocrd(g_crd_fname)


    # -----------------
    # make a dictionary
    # -----------------
    atype_2_atypeidx = Dict()
    for (j, atype) in enumerate(g_top.top_default_CGIDR_MPIPI_atomtype)
        atype_2_atypeidx[atype.name] = j
    end

    # =====================
    # generate LAMMPS param
    # =====================
    L_param_fname = "OUT_param_MPIPI.dat"
    L_Param_file = open(L_param_fname, "w")
    # step 1: set type charge
    # println(g_top.top_default_CGIDR_MPIPI_atomtype)
    for (j, a) in enumerate(g_top.top_default_CGIDR_MPIPI_atomtype)
        c = g_top.top_default_CGIDR_MPIPI_atomtype[j].charge
        if abs(c) < 1e-6
            continue
        end
        @printf(L_Param_file, "set type %-2d charge %8.3f\n", j, c)
    end
    @printf(L_Param_file, "\n")
    # step 2: interaction types
    @printf(L_Param_file, "pair_style hybrid/overlay wf/cut %6.1f coul/debye %6.4f %6.1f\n", lj_cut, κ, 0)
    # step 3: bonded coeff
    @printf(L_Param_file, "bond_coeff 1        %6.4f %6.1f \n\n", 40.1664 / 4.184, 3.80)
    # step 4: nonbonded coeff
    for (j, pj) in enumerate(g_top.top_default_CGIDR_MPIPI_pairparams)
        name_a  = pj.name1
        name_b  = pj.name2
        eps = pj.epsilon
        sig = pj.sigma
        mu  = pj.mu
        @printf(L_Param_file,
                "pair_coeff %2d %2d wf/cut %9.6f %18.15f 1 %2d %18.15f\n",
                atype_2_atypeidx[name_a], atype_2_atypeidx[name_b],
                eps, sig * 10, mu, sig * 30)
    end
    for (j, aj) in enumerate(g_top.top_default_CGIDR_MPIPI_atomtype)
        charge_j = g_top.top_default_CGIDR_MPIPI_atomtype[j].charge
        if abs(charge_j) < 1e-6
            continue
        end
        for (k, ak) in enumerate(g_top.top_default_CGIDR_MPIPI_atomtype)
            charge_k = g_top.top_default_CGIDR_MPIPI_atomtype[k].charge
            if abs(charge_k) < 1e-6
                continue
            end
            if k < j
                continue
            end
            @printf(L_Param_file,
                    "pair_coeff %2d %2d coul/debye %8.3f \n",
                    j, k, ele_cut)
        end
    end
    close(L_Param_file)

    # =========================
    # generate LAMMPS structure
    # =========================
    # step 1: output preparation
    @printf(L_Struct_file, "LAMMPS data file via GENESIS-cg-tool \n\n")
    # step 2: atom and bond types
    @printf(L_Struct_file, "%-8d atoms\n", g_top.num_atom)
    @printf(L_Struct_file, "%-8d atom types\n", 40)
    @printf(L_Struct_file, "%-8d bonds\n", length(g_top.top_bonds))
    @printf(L_Struct_file, "%-8d bond types\n", 1)
    @printf(L_Struct_file, "\n")
    # step 3: PBC
    @printf(L_Struct_file, "%-8.1f %-8.1f  xlo  xhi\n", PBCxlo, PBCxhi)
    @printf(L_Struct_file, "%-8.1f %-8.1f  ylo  yhi\n", PBCylo, PBCyhi)
    @printf(L_Struct_file, "%-8.1f %-8.1f  zlo  zhi\n", PBCzlo, PBCzhi)
    @printf(L_Struct_file, "\n")
    # step 4: Masses
    @printf(L_Struct_file, "Masses\n\n")
    for (j, atype) in enumerate(g_top.top_default_CGIDR_MPIPI_atomtype)
        @printf(L_Struct_file, "%-2d  %-8.3f\n", j, atype.mass)
    end
    @printf(L_Struct_file, "\n")
    # step 5: Atoms
    # atom_id; chain_id; atom_type; atom_charge; x; y; z; PBCi; PBCj; PBCk
    @printf(L_Struct_file, "Atoms\n\n")
    for (j, a) in enumerate(g_top.top_atoms)
        x, y, z = g_crd.coors[:, j]
        @printf(L_Struct_file,
                "%-8d %-8d %-3d %8.3f %8.3f %8.3f %8.3f  0 0 0\n",
                a.atom_index, a.chain_id, atype_2_atypeidx[a.atom_type], a.charge,
                x, y, z)
    end
    @printf(L_Struct_file, "\n")
    # step 6: Velocities
    @printf(L_Struct_file, "Velocities\n\n")
    for (j, a) in enumerate(g_top.top_atoms)
        vx, vy, vz = 0.0, 0.0, 0.0
        @printf(L_Struct_file,
                "%-8d  %8.3f %8.3f %8.3f\n",
                a.atom_index,
                vx, vy, vz)
    end
    @printf(L_Struct_file, "\n")
    # step 7: Bonds
    @printf(L_Struct_file, "Bonds\n\n")
    for (j, b) in enumerate(g_top.top_bonds)
        @printf(L_Struct_file,
                "%-8d  %2d %8d %8d\n",
                j, 1, b.i, b.j)
    end
    @printf(L_Struct_file, "\n")
    # closing
    close(L_Struct_file)

end

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "--top", "-t"
        help     = "GENESIS Topology (input)"
        arg_type = String
        required = true

        "--crd", "-c"
        help     = "GENESIS Coordinate (input)"
        arg_type = String
        required = true

        "--outstructure", "-o"
        help     = "LAMMPS Structure (output)"
        arg_type = String
        default  = ""

        "--verbose", "-V"
        help   = "Verbose"
        action = :store_true

        "--kappa", "-k"
        help     = "Kappa (1 / Debye_length)"
        arg_type = Float64
        default  = 0.1299

        "--LJcutoff"
        help     = "Cutoff for LJ"
        arg_type = Float64
        default  = 25.0

        "--ELEcutoff"
        help     = "Cutoff for Debye-Huckel ELE"
        arg_type = Float64
        default  = 35.0

        "--PBCxlo"
        help     = "PBC x low"
        arg_type = Float64
        default  = -200.0

        "--PBCylo"
        help     = "PBC y low"
        arg_type = Float64
        default  = -200.0

        "--PBCzlo"
        help     = "PBC z low"
        arg_type = Float64
        default  = -200.0

        "--PBCxhi"
        help     = "PBC x high"
        arg_type = Float64
        default  = 200.0

        "--PBCyhi"
        help     = "PBC y high"
        arg_type = Float64
        default  = 200.0

        "--PBCzhi"
        help     = "PBC z high"
        arg_type = Float64
        default  = 200.0
    end

    return parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_commandline()
    g2l(args)
end
