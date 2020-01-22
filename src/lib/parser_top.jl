###############################################################################
#                                       _                                     #
#                     __ _  _ __  ___  | |_  ___   _ __                       #
#                    / _` || '__|/ _ \ | __|/ _ \ | '_ \                      #
#                   | (_| || |  | (_) || |_| (_) || |_) |                     #
#                    \__, ||_|   \___/  \__|\___/ | .__/                      #
#                    |___/                        |_|                         #
#                                                                             #
###############################################################################


###############################################################################
#                                function lists                               
# write_grotop(top::GenTopology, system_name::AbstractString, args::Dict{String, Any})
# write_grotop_pwmcos(top::GenTopology, system_name::AbstractString, args::Dict{String, Any})
# read_groitp(itp_filename::AbstractString)
# read_grotop(top_filename::AbstractString)
# write_psf(top::GenTopology, sys_name::AbstractString, args::Dict{String, Any})
###############################################################################

using Printf

function write_grotop(top::GenTopology, system_name::AbstractString, args::Dict{String, Any}=Dict{String, Any}())

    verbose  = get(args, "verbose", false)

    top_name = system_name * ".top"
    itp_name = system_name * ".itp"

    # ========
    # top file
    # ========
    top_file = open(top_name, "w")

    print(top_file, "; common interaction parameters for CG models\n")
    print(top_file, "#include \"./param/atom_types.itp\" \n")
    print(top_file, "; AICG2+ flexible local angle parameters \n")
    print(top_file, "#include \"./param/flexible_local_angle.itp\" \n")
    print(top_file, "; AICG2+ flexible local dihedral parameters \n")
    print(top_file, "#include \"./param/flexible_local_dihedral.itp\" \n")
    print(top_file, "\n")

    print(top_file, "; Molecule topology \n")
    print(top_file, "#include \"./top/", itp_name, "\" \n\n")

    print(top_file, "[ system ] \n")
    print(top_file, system_name, " \n\n")

    print(top_file, "[ molecules ] \n")
    print(top_file, system_name, "  1 \n\n")

    print(top_file, "; [ cg_ele_mol_pairs ] \n")
    print(top_file, "; ON 1 - 2 : 3 - 4 \n")
    print(top_file, "; OFF 1 - 1 : 3 - 3 \n")
    print(top_file, "; OFF 1 - 1 \n\n")

    print(top_file, "; [ pwmcos_mol_pairs ] \n")
    print(top_file, "; ON 1 - 2 : 3 - 4 \n")
    print(top_file, "; OFF 1 - 1 : 3 - 3 \n")
    print(top_file, "; OFF 2 - 3 \n\n")

    close(top_file)
    
    # ========
    # itp file
    # ========
    itp_file = open(itp_name, "w")



    ###########################################################################
    #                         define output functions                         #
    ###########################################################################

    # ----------------
    # [ moleculetype ]
    # ----------------
    wr_itp_mol_head(io::IO) = print(io, "[ moleculetype ]\n")
    wr_itp_mol_comm(io::IO) = @printf(io, ";%15s %6s\n", rpad("name", 15), "nrexcl")
    wr_itp_mol_line(io::IO, name::AbstractString, n::Int) = @printf(io, "%16s %6d\n", rpad(name, 16), n)

    # ---------
    # [ atoms ]
    # ---------
    wr_itp_atm_head(io::IO) = print(io, "[ atoms ]\n")
    wr_itp_atm_comm(io::IO) = @printf(io, ";%9s%5s%10s%5s%5s%5s %8s %8s\n", "nr", "type", "resnr", "res", "atom", "cg", "charge", "mass")
    function wr_itp_atm_chain_info(io::IO, a::GenTopAtom)
        @printf(io, "; +INFO+ CHAIN: %6d     SEGNAME: %5s\n", a.chain_id, a.seg_name)
    end
    function wr_itp_atm_line(io::IO, a::GenTopAtom)
        @printf(io, "%10d%5s%10d%5s%5s%5d %8.3f %8.3f\n",
                a.atom_index, a.atom_type, a.residue_index, a.residue_name, a.atom_name,
                a.function_type, a.charge, a.mass)
    end

    # ---------
    # [ bonds ]
    # ---------
    wr_itp_bnd_head(io::IO) = print(io, "[ bonds ]\n")
    wr_itp_bnd_comm(io::IO) = @printf(io, ";%9s%10s%5s%18s%18s\n", "i", "j", "f", "eq", "coef")
    function wr_itp_bnd_line(io::IO, b::GenTopBond)
        @printf(io, "%10d%10d%5d%18.4E%18.4E\n", b.i, b.j, b.function_type, b.r0, b.coef)
    end

    # ----------
    # [ angles ]
    # ----------
    wr_itp_ang_head(io::IO) = print(io, "[ angles ]\n")
    wr_itp_ang_comm(io::IO) = @printf(io, ";%9s%10s%10s%5s%15s%15s%15s\n", "i", "j", "k", "f", "eq", "coef", "w")
    function wr_itp_ang_line(io::IO, a::GenTopAngle)
        f = a.function_type
        if f == 1
            @printf(io, "%10d%10d%10d%5d%15.4E%15.4E\n", a.i, a.j, a.k, f, a.a0, a.coef)
        elseif f == 21
            @printf(io, "%10d%10d%10d%5d%15.4E%15.4E%15.4E\n", a.i, a.j, a.k, f, a.a0, a.coef, a.w)
        elseif f == 22
            @printf(io, "%10d%10d%10d%5d\n", a.i, a.j, a.k, f)
        end
    end

    # -------------
    # [ dihedrals ]
    # -------------
    wr_itp_dih_head(io::IO) = print(io, "[ dihedrals ]\n")
    wr_itp_dih_comm(io::IO) = @printf(io, ";%9s%10s%10s%10s%5s%15s%15s%15s\n", "i", "j", "k", "l", "f", "eq", "coef", "w/n")
    function wr_itp_dih_line(io::IO, d::GenTopDihedral)
        f = d.function_type
        if f == 1
            @printf(io, "%10d%10d%10d%10d%5d%15.4E%15.4E%15d\n", d.i, d.j, d.k, d.l, f, d.d0, d.coef, d.n)
        elseif f == 21
            @printf(io, "%10d%10d%10d%10d%5d%15.4E%15.4E%15.4E\n", d.i, d.j, d.k, d.l, f, d.d0, d.coef, d.w)
        elseif f == 22
            @printf(io, "%10d%10d%10d%10d%5d\n", d.i, d.j, d.k, d.l, f)
        end
    end

    # ---------
    # [ pairs ]
    # ---------
    wr_itp_pair_head(io::IO) = @printf(io, "[ pairs ]\n")
    wr_itp_pair_comm(io::IO) = @printf(io, ";%9s%10s%10s%15s%15s\n", "i", "j", "f", "eq", "coef")
    function wr_itp_pair_line(io::IO, c::GenTopPair)
        @printf(io, "%10d%10d%10d%15.4E%15.4E\n", c.i, c.j, c.function_type, c.r0, c.coef)
    end

    # -------------
    # [ exclusions]
    # -------------
    wr_itp_exc_head(io::IO) = print(io, "[ exclusions ] ; Genesis exclusion list\n")
    wr_itp_exc_comm(io::IO) = @printf(io, ";%9s%10s\n", "i", "j")
    wr_itp_exc_line(io::IO, e::GenTopExclusion) = @printf(io, "%10d%10d\n", e.i, e.j)

    # ---------------------
    # [ cg_IDR_HPS_region ]
    # ---------------------
    wr_itp_idr_hps_head(io::IO) = print(io, "[ cg_IDR_HPS_region ] ; IDR HPS model \n")
    wr_itp_idr_hps_comm(io::IO) = @printf(io, ";%9s to %10s\n", "i", "j")
    wr_itp_idr_hps_line(io::IO, e::GenTopIDRHPS) = @printf(io, "%10d to %10d\n", e.istart, e.iend)



    ###########################################################################
    #                        Beging writing to file...                        #
    ###########################################################################

    # ----------------
    # [ moleculetype ]
    # ----------------
    wr_itp_mol_head(itp_file)
    wr_itp_mol_comm(itp_file)
    wr_itp_mol_line(itp_file, system_name, MOL_NR_EXCL)
    print(itp_file,"\n")

    # ---------
    # [ atoms ]
    # ---------
    cg_num_particles = top.num_atom
    wr_itp_atm_head(itp_file)
    wr_itp_atm_comm(itp_file)
    tmp_chain_id       = 0
    for atom in top.top_atoms
        i_chain = atom.chain_id
        if i_chain > tmp_chain_id
            wr_itp_atm_chain_info(itp_file, atom)
            tmp_chain_id = i_chain
        end
        wr_itp_atm_line(itp_file, atom)
    end
    print(itp_file,"\n")

    # ---------
    # [ bonds ]
    # ---------
    if length(top.top_bonds) > 0
        wr_itp_bnd_head(itp_file)
        wr_itp_bnd_comm(itp_file)
        for bond in top.top_bonds
            wr_itp_bnd_line(itp_file, bond)
        end
        print(itp_file, "\n")
    end

    # ----------
    # [ angles ]
    # ----------
    if length(top.top_angles) > 0
        wr_itp_ang_head(itp_file)
        wr_itp_ang_comm(itp_file)
        for angle in top.top_angles
            wr_itp_ang_line(itp_file, angle)
        end
        print(itp_file, "\n")
    end

    # -------------
    # [ dihedrals ]
    # -------------
    if length(top.top_dihedrals) > 0
        wr_itp_dih_head(itp_file)
        wr_itp_dih_comm(itp_file)
        for dih in top.top_dihedrals
            wr_itp_dih_line(itp_file, dih)
        end
        print(itp_file, "\n")
    end

    # ---------
    # [ pairs ]
    # ---------
    if length(top.top_pairs) > 0
        wr_itp_pair_head(itp_file)
        wr_itp_pair_comm(itp_file)
        for pair in top.top_pairs
            wr_itp_pair_line(itp_file, pair)
        end
        print(itp_file, "\n")
    end

    # ---------------------
    #        [ exclusions ]
    # ---------------------
    if length(top.top_exclusions) > 0
        wr_itp_exc_head(itp_file)
        wr_itp_exc_comm(itp_file)
        for exclusion in top.top_exclusions
            wr_itp_exc_line(itp_file, exclusion)
        end
        print(itp_file, "\n")
    end

    # ---------------------
    # [ cg_IDR_HPS_region ]
    # ---------------------
    if length(top.top_idr_hps) > 0
        wr_itp_idr_hps_head(itp_file)
        wr_itp_idr_hps_comm(itp_file)
        for idr in top.top_idr_hps
            wr_itp_idr_hps_line(itp_file, idr)
        end
        print(itp_file, "\n")
    end

    close(itp_file)

    if verbose
        println(">           ... .top: DONE!")
    end
end


function write_grotop_pwmcos(top::GenTopology, system_name::AbstractString, args::Dict{String, Any}=Dict{String, Any}())

    verbose = get(args, "verbose", false)
    appendto_filename = get(args, "patch", "")

    if length( appendto_filename ) == 0
        itp_pwmcos_name = system_name * "_pwmcos.itp_patch"
        itp_pwmcos_file = open(itp_pwmcos_name, "w")
    else
        itp_pwmcos_name = appendto_filename
        itp_pwmcos_file = open(itp_pwmcos_name, "a")
    end

    itp_pwmcos_head = "[ pwmcos ]\n"
    itp_pwmcos_comm = @sprintf(";%5s%4s%9s%9s%9s%9s%12s%12s%12s%12s%8s%8s\n",
                               "i", "f", "r0", "theta1", "theta2", "theta3",
                               "ene_A", "ene_C", "ene_G", "ene_T",
                               "gamma", "eps'")
    itp_pwmcos_line = "%6d %3d %8.5f %8.3f %8.3f %8.3f%12.6f%12.6f%12.6f%12.6f%8.3f%8.3f \n"

    print(itp_pwmcos_file, itp_pwmcos_head)
    print(itp_pwmcos_file, itp_pwmcos_comm)
    for p in top.top_pwmcos
        @printf(itp_pwmcos_file,
                "%6d %3d %8.5f %8.3f %8.3f %8.3f%12.6f%12.6f%12.6f%12.6f%8.3f%8.3f \n",
                p.i, p.function_type, p.r0, p.theta1, p.theta2, p.theta3,
                p.ene_A, p.ene_C, p.ene_G, p.ene_T, p.gamma, p.eps)
    end
    print(itp_pwmcos_file, "\n")

    close(itp_pwmcos_file)
    if verbose
        println(">           ... ", itp_pwmcos_name, " pwmcos.itp: DONE!")
    end
end

# ==================================
# General Topology in Gromacs format
# ==================================

function read_groitp(itp_filename::AbstractString)

    mol_name          = ""
    nonlocal_interval = 0
    num_atom          = 0

    top_atoms         = Vector{GenTopAtom}(undef, 0)
    top_bonds         = Vector{GenTopBond}(undef, 0)
    top_angles        = Vector{GenTopAngle}(undef, 0)
    top_dihedrals     = Vector{GenTopDihedral}(undef, 0)
    top_pairs         = Vector{GenTopPair}(undef, 0)
    top_exclusions    = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos        = Vector{GenTopPWMcos}(undef, 0)
    top_idr_hps       = Vector{GenTopIDRHPS}(undef, 0)

    function read_top_atoms(line::AbstractString, c_id::Int, s_name::AbstractString)
        words = split(line)
        a_indx = parse(Int, words[1])
        a_type = words[2]
        r_indx = parse(Int, words[3])
        r_name = words[4]
        a_name = words[5]
        f_type = parse(Int, words[6])
        charge = parse(Float64, words[7])
        mass   = parse(Float64, words[8])
        new_atom = GenTopAtom(a_indx, a_type, r_indx, r_name,
                              a_name, f_type, charge, mass, c_id, s_name)
        push!(top_atoms, new_atom)
    end

    function read_top_bonds(line::AbstractString)
        words  = split(line)
        i      = parse(Int, words[1])
        j      = parse(Int, words[2])
        f_type = parse(Int, words[3])
        r0     = parse(Float64, words[4])
        coef   = parse(Float64, words[5])
        new_bond = GenTopBond(i, j, f_type, r0, coef)
        push!(top_bonds, new_bond)
    end

    function read_top_angles(line::AbstractString)
        words  = split(line)
        i      = parse(Int, words[1])
        j      = parse(Int, words[2])
        k      = parse(Int, words[3])
        f_type = parse(Int, words[4])
        if f_type == 1
            eq   = parse(Float64, words[5])
            coef = parse(Float64, words[6])
            w    = 0.0
        elseif f_type == 21
            eq   = parse(Float64, words[5])
            coef = parse(Float64, words[6])
            w    = parse(Float64, words[7])
        elseif f_type == 22
            eq   = 0.0
            coef = 0.0
            w    = 0.0
        end
        new_angle = GenTopAngle(i, j, k, f_type, eq, coef, w)
        push!(top_angles, new_angle)
    end

    function read_top_dihedrals(line::AbstractString)
        words  = split(line)
        i      = parse(Int, words[1])
        j      = parse(Int, words[2])
        k      = parse(Int, words[3])
        l      = parse(Int, words[4])
        f_type = parse(Int, words[5])
        if f_type == 1
            eq   = parse(Float64, words[6])
            coef = parse(Float64, words[7])
            w    = 0.0
            n    = parse(Int, words[8])
        elseif f_type == 21
            eq   = parse(Float64, words[6])
            coef = parse(Float64, words[7])
            w    = parse(Float64, words[8])
            n    = 0
        elseif f_type == 22
            eq   = 0.0
            coef = 0.0
            w    = 0.0
            n    = 0
        end
        new_dihedral = GenTopDihedral(i, j, k, l, f_type, eq, coef, w, n)
        push!(top_dihedrals, new_dihedral)
    end

    function read_top_pairs(line::AbstractString)
        words  = split(line)
        i      = parse(Int, words[1])
        j      = parse(Int, words[2])
        f_type = parse(Int, words[3])
        r0     = parse(Float64, words[4])
        coef   = parse(Float64, words[5])
        new_pair = GenTopPair(i, j, f_type, r0, coef)
        push!(top_pairs, new_pair)
    end

    function read_top_exclusions(line::AbstractString)
        words  = split(line)
        i      = parse(Int, words[1])
        j      = parse(Int, words[2])
        new_ex = GenTopExclusion(i, j)
        push!(top_exclusions, new_ex)
    end

    function read_top_pwmcos(line::AbstractString)
        words  = split(line)
        i      = parse(Int, words[1])
        f_type = parse(Int, words[2])
        r0 = parse(Float64, words[3])
        t1 = parse(Float64, words[4])
        t2 = parse(Float64, words[5])
        t3 = parse(Float64, words[6])
        eA = parse(Float64, words[7])
        eC = parse(Float64, words[8])
        eG = parse(Float64, words[9])
        eT = parse(Float64, words[10])
        gm = parse(Float64, words[11])
        ep = parse(Float64, words[12])
        new_pwmcos = GenTopPWMcos(i, f_type, r0, t1, t2, t3,
                                  eA, eC, eG, eT, gm, ep)
        push!(top_pwmcos, new_pwmcos)
    end

    function read_top_idr_hps(line::AbstractString)
        words   = split(line)
        words   = split(line, r"\s*to\s*", keepempty =false)
        i       = parse(Int, words[1])
        j       = parse(Int, words[2])
        new_idr = GenTopIDRHPS(i, j)
        push!(top_idr_hps, new_idr)
    end


    # ---------
    # main part
    # ---------
    section_name = ""
    c_id_tmp = 0
    s_name_tmp = ""
    for line in eachline(itp_filename)
        if section_name == "atoms" && startswith(line, "; +INFO+")
            words = split( line[10:end] )
            if words[1] == "CHAIN:"
                c_id_tmp = parse(Int, words[2])
            end
            if words[3] == "SEGNAME:" && length(words) > 3
                s_name_tmp = words[4]
            end
        end

        sep  = findfirst(";", line)
        if sep != nothing
            line = strip(line[1 : sep[1] - 1])
        else
            line = strip(line)
        end
        if length(line) == 0
            continue
        end

        if line[1] == '['
            sep = findfirst("]", line)
            section_name = strip(line[2 : sep[1] - 1])
            continue
        end

        if section_name == "moleculetype"
            words = split(line)
            mol_name = words[1]
            nonlocal_interval = parse(Int, words[2])
        else
            # read_function_name = "read_top_" * section_name * "(line)"
            # read_expression = Meta.parse(read_function_name)
            # eval(read_expression)
            if section_name == "atoms"
                read_top_atoms(line, c_id_tmp, s_name_tmp)
            elseif section_name == "bonds"
                read_top_bonds(line)
            elseif section_name == "angles"
                read_top_angles(line)
            elseif section_name == "dihedrals"
                read_top_dihedrals(line)
            elseif section_name == "pairs"
                read_top_pairs(line)
            elseif section_name == "exclusions"
                read_top_exclusions(line)
            elseif section_name == "pwmcos"
                read_top_pwmcos(line)
            elseif section_name == "cg_IDR_HPS_region"
                read_top_idr_hps(line)
            end
        end
    end

    num_atom = length(top_atoms)

    new_top_mol = GenTopMolecule(mol_name, nonlocal_interval, num_atom,
                                 top_atoms,
                                 top_bonds,
                                 top_angles,
                                 top_dihedrals,
                                 top_pairs,
                                 top_exclusions,
                                 top_pwmcos,
                                 top_idr_hps)

    return new_top_mol
end

function read_grotop(top_filename::AbstractString)

    sys_name = ""
    num_atom = 0
    mol_id   = 0

    top_default_params             = GenTopDefault(0, 0, false, 0.0, 0.0) 
    top_default_atomtype           = Vector{GenTopAtomType}(undef, 0)
    top_default_CGDNA_bp           = Vector{GenTopCGDNABasepairType}(undef, 0)
    top_default_CGDNA_bs           = Vector{GenTopCGDNABasestackType}(undef, 0)
    top_default_CGDNA_cs           = Vector{GenTopCGDNABasecrossType}(undef, 0)
    top_default_CGDNA_exv          = Vector{GenTopCGDNAExvType}(undef, 0)
    top_default_CGPro_flx_angle    = Vector{GenTopCGProAICGFlexAngleType}(undef, 0)
    top_default_CGPro_flx_dihedral = Vector{GenTopCGProAICGFlexDihedralType}(undef, 0)

    global_index_2_local_index = Vector{Int}(undef, 0)
    global_index_2_local_molid = Vector{Int}(undef, 0)
    top_atoms                  = Vector{GenTopAtom}(undef, 0)
    top_bonds                  = Vector{GenTopBond}(undef, 0)
    top_angles                 = Vector{GenTopAngle}(undef, 0)
    top_dihedrals              = Vector{GenTopDihedral}(undef, 0)
    top_pairs                  = Vector{GenTopPair}(undef, 0)
    top_exclusions             = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos                 = Vector{GenTopPWMcos}(undef, 0)
    top_idr_hps                = Vector{GenTopIDRHPS}(undef, 0)
    top_mol_list               = Vector{GenTopMolList}(undef, 0)

    section_name = ""
    mol_topologies = Dict()
    for line in eachline(top_filename)
        sep  = findfirst(";", line)
        if sep != nothing
            line = strip(line[1 : sep[1] - 1])
        else
            line = strip(line)
        end
        if length(line) == 0
            continue
        end

        top_dirname = dirname(top_filename) * "/"
        if startswith(line, "#include")
            mol_file_name = strip(line[9:end], ['\"', '\'', ' '])
            mol_file_basename = basename(mol_file_name)
            if in(mol_file_basename, ["atom_types.itp", "flexible_local_angle.itp", "flexible_local_dihedral.itp"])
                continue
            end
            if !isabspath(mol_file_name)
                mol_file_name = normpath( joinpath( top_dirname, mol_file_name ) ) 
            end
            new_mol = read_groitp(mol_file_name)
            new_mol_name = new_mol.mol_name
            mol_topologies[new_mol_name] = new_mol
        end

        if line[1] == '['
            sep = findfirst("]", line)
            section_name = strip(line[2 : sep[1] - 1])
            continue
        end

        if section_name == "system"
            words = split(line)
            sys_name = words[1]
        elseif section_name == "molecules"
            words = split(line)
            mol_name = words[1]
            mol_count = parse(Int, words[2])

            tmp_mol_list = GenTopMolList(mol_name, mol_count)
            push!(top_mol_list, tmp_mol_list)

            for i = 1 : mol_count
                tmp_mol = mol_topologies[mol_name]
                mol_id  += 1

                # -----------------------
                # add molecules to system
                # -----------------------
                for t in tmp_mol.top_atoms
                    new_index = t.atom_index + num_atom
                    s = GenTopAtom(new_index,
                                   t.atom_type,
                                   t.residue_index,
                                   t.residue_name,
                                   t.atom_name,
                                   t.function_type,
                                   t.charge,
                                   t.mass,
                                   t.chain_id,
                                   t.seg_name)
                    push!(top_atoms, s)
                    push!(global_index_2_local_index, t.atom_index)
                    push!(global_index_2_local_molid, mol_id)
                end
                for t in tmp_mol.top_bonds
                    s = GenTopBond(t.i + num_atom,
                                   t.j + num_atom,
                                   t.function_type,
                                   t.r0,
                                   t.coef
                                   )
                    push!(top_bonds, s)
                end
                for t in tmp_mol.top_angles
                    s = GenTopAngle(t.i + num_atom,
                                    t.j + num_atom,
                                    t.k + num_atom,
                                    t.function_type,
                                    t.a0, t.coef, t.w)
                    push!(top_angles, s)
                end
                for t in tmp_mol.top_dihedrals
                    s = GenTopDihedral(t.i + num_atom,
                                       t.j + num_atom,
                                       t.k + num_atom,
                                       t.l + num_atom,
                                       t.function_type,
                                       t.d0, t.coef, t.w, t.n)
                    push!(top_dihedrals, s)
                end
                for t in tmp_mol.top_pairs
                    s = GenTopPair(t.i + num_atom,
                                   t.j + num_atom,
                                   t.function_type,
                                   t.r0, t.coef)
                    push!(top_pairs, s)
                end
                for t in tmp_mol.top_exclusions
                    s = GenTopExclusion(t.i + num_atom,
                                        t.j + num_atom)
                    push!(top_exclusions, s)
                end
                for t in tmp_mol.top_pwmcos
                    s = GenTopPWMcos(t.i + num_atom,
                                     t.function_type,
                                     t.r0,
                                     t.theta1,
                                     t.theta2,
                                     t.theta3,
                                     t.ene_A, t.ene_C, t.ene_G, t.ene_T,
                                     t.gamma, t.eps)
                    push!(top_pwmcos, s)
                end
                for t in tmp_mol.top_idr_hps
                    s = GenTopIDRHPS(t.istart + num_atom,
                                     t.iend + num_atom)
                    push!(top_idr_hps, s)
                end

                num_atom += tmp_mol.num_atom
            end
        end
    end

    new_top = GenTopology(sys_name, num_atom,
                          top_default_params,
                          top_default_atomtype,
                          top_default_CGDNA_bp,
                          top_default_CGDNA_bs,
                          top_default_CGDNA_cs,
                          top_default_CGDNA_exv,
                          top_default_CGPro_flx_angle,
                          top_default_CGPro_flx_dihedral,
                          global_index_2_local_index,
                          global_index_2_local_molid,
                          top_atoms,
                          top_bonds,
                          top_angles,
                          top_dihedrals,
                          top_pairs,
                          top_exclusions,
                          top_pwmcos,
                          top_idr_hps,
                          top_mol_list)

    return new_top

end

###############################################################################
#                                             __                              #
#                                _ __   ___  / _|                             #
#                               | '_ \ / __|| |_                              #
#                               | |_) |\__ \|  _|                             #
#                               | .__/ |___/|_|                               #
#                               |_|                                           #
#                                                                             #
###############################################################################

function write_psf(top::GenTopology, sys_name::AbstractString="", args::Dict{String, Any}=Dict{String, Any}())

    verbose = get(args, "verbose", false)

    if length(sys_name) > 0
        system_name = sys_name
    else
        system_name = top.system_name
    end
    psf_name = system_name * ".psf"
    psf_file = open(psf_name, "w")

    cg_num_particles = top.num_atom

    @printf(psf_file, "PSF CMAP \n\n")
    @printf(psf_file, "      3 !NTITLE \n")
    @printf(psf_file, "REMARKS PSF file created with Julia. \n")
    @printf(psf_file, "REMARKS System: %s  \n", system_name)
    @printf(psf_file, "REMARKS ======================================== \n")
    @printf(psf_file, "       \n")

    psf_atom_line = " %6d %3s %5d %3s %3s %5s  %10.6f  %10.6f          0 \n"
    chain_id_set = "_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"

    @printf(psf_file, " %6d !NATOM \n", cg_num_particles)
    for i_bead in 1 : cg_num_particles
        @printf(psf_file, " %6d %3s %5d %3s %3s %5s  %10.6f  %10.6f          0 \n",
                i_bead,
                chain_id_set[mod(top.global_index_2_local_molid[i_bead], 63) + 1],
                top.top_atoms[i_bead].residue_index,
                top.top_atoms[i_bead].residue_name,
                top.top_atoms[i_bead].atom_name,
                top.top_atoms[i_bead].atom_type,
                top.top_atoms[i_bead].charge,
                top.top_atoms[i_bead].mass)
    end
    print(psf_file,"\n")

    close(psf_file)

    if verbose
        println(">           ... .psf: DONE!")
    end
end

