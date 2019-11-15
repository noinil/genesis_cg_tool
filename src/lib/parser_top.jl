###############################################################################
#                                       _                                     #
#                     __ _  _ __  ___  | |_  ___   _ __                       #
#                    / _` || '__|/ _ \ | __|/ _ \ | '_ \                      #
#                   | (_| || |  | (_) || |_| (_) || |_) |                     #
#                    \__, ||_|   \___/  \__|\___/ | .__/                      #
#                    |___/                        |_|                         #
#                                                                             #
###############################################################################

using Printf

function write_cg_grotop(top::CGTopology, force_field::ForceFieldCG, system_name::String, args::Dict{String, Any})

    ff_pro     = force_field.ff_protein
    ff_dna     = force_field.ff_DNA
    ff_rna     = force_field.ff_RNA
    ff_pro_dna = force_field.ff_protein_DNA

    top_name = system_name * "_cg.top"
    itp_name = system_name * "_cg.itp"

    gen_3spn_itp       = args["3spn-param"]
    ccgo_contact_scale = args["CCGO-contact-scale"]

    # =====================
    # Prepare modifications
    # =====================
    toml_filename = args["config"]

    has_toml_mod  = false
    if length(toml_filename) > 0
        has_toml_mod    = true
        new_toml_config = read_TOML(toml_filename)

        AICG2p_flexible_local    = []
        AICG2p_flexible_nonlocal = []
        HPS_IDR_region           = []
        if haskey(new_toml_config, "IDR")
            if haskey(new_toml_config["IDR"], "AICG2p_IDR_local")
                index_string = new_toml_config["IDR"]["AICG2p_IDR_local"]
                AICG2p_flexible_local = parse_selection(index_string)
            end
            if haskey(new_toml_config["IDR"], "AICG2p_IDR_nonlocal")
                index_string = new_toml_config["IDR"]["AICG2p_IDR_nonlocal"]
                AICG2p_flexible_nonlocal = parse_selection(index_string)
            end
            if haskey(new_toml_config["IDR"], "HPS_region")
                index_string = new_toml_config["IDR"]["HPS_region"]
                HPS_IDR_region = parse_selection(index_string)
            end
        end
    end


    # ========
    # top file
    # ========
    top_file = open(top_name, "w")

    print(top_file, "; atom types for coarse-grained models\n")
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
    print(top_file, "; OFF 1 - 1 : 3 - 3 \n\n")
    print(top_file, "; OFF 1 - 1 \n\n")

    print(top_file, "; [ pwmcos_mol_pairs ] \n")
    print(top_file, "; ON 1 - 2 : 3 - 4 \n")
    print(top_file, "; OFF 1 - 1 : 3 - 3 \n\n")
    print(top_file, "; OFF 2 - 3 \n\n")

    if has_toml_mod
        if haskey(new_toml_config["IDR"], "HPS_region")
            index_string = new_toml_config["IDR"]["HPS_region"]

            print(top_file, "[ cg_IDR_HPS_region ] \n")
            hps_words = split(index_string, r"\s*,\s*", keepempty=false)
            for w in hps_words
                if occursin("to", w)
                    println(top_file, replace(w, "to" => " - "))
                else
                    println(top_file, w, " - ", w)
                end
            end
        end
    end

    close(top_file)
    
    # ========
    # itp file
    # ========
    itp_file = open(itp_name, "w")

    wr_itp_mol_head(io::IO) = print(io, "[ moleculetype ]\n")
    wr_itp_mol_comm(io::IO) = @printf(io, ";%15s %6s\n", rpad("name", 15), "nrexcl")
    wr_itp_mol_line(io::IO, name::String, n::Int) = @printf(io, "%16s %6d\n", rpad(name, 16), n)

    wr_itp_atm_head(io::IO) = print(io, "[ atoms ]\n")
    wr_itp_atm_comm(io::IO) = @printf(io, ";%9s%5s%10s%5s%5s%5s %8s %8s\n", "nr", "type", "resnr", "res", "atom", "cg", "charge", "mass")
    function wr_itp_atm_line(io::IO, i, beadtype, resid, resname, beadname, n, charge, mass)
        @printf(io, "%10d%5s%10d%5s%5s%5d %8.3f %8.3f\n", i, beadtype, resid, resname, beadname, n, charge, mass)
    end

    wr_itp_bnd_head(io::IO) = print(io, "[ bonds ]\n")
    wr_itp_bnd_comm(io::IO) = @printf(io, ";%9s%10s%5s%18s%18s\n", "i", "j", "f", "eq", "coef")
    function wr_itp_bnd_line(io::IO, b::CGTopBond, f::Int, e::Float64)
        @printf(io, "%10d%10d%5d%18.4E%18.4E\n", b.i, b.j, f, b.r0 * 0.1, e * 100.0 * 2.0 * CAL2JOU)
    end

    wr_itp_13_head(io::IO) = print(io, "[ angles ] ; AICG2+ 1-3 interaction\n")
    wr_itp_13_comm(io::IO) = @printf(io, ";%9s%10s%10s%5s%15s%15s%15s\n", "i", "j", "k", "f", "eq", "coef", "w")
    function wr_itp_13_line(io::IO, a::CGTopAngle, f::Int, e::Float64, sigma::Float64)
        @printf(io, "%10d%10d%10d%5d%15.4E%15.4E%15.4E\n", a.i, a.j, a.k, f, a.a0 * 0.1, e * CAL2JOU, sigma * 0.1)
    end

    wr_itp_ang_f_head(io::IO) = print(io, "[ angles ] ; AICG2+ flexible local interaction\n")
    wr_itp_ang_f_comm(io::IO) = @printf(io, ";%9s%10s%10s%5s\n", "i", "j", "k", "f")
    function wr_itp_ang_f_line(io::IO, a::CGTopAngle, f::Int)
        @printf(io, "%10d%10d%10d%5d\n", a.i, a.j, a.k, f)
    end

    wr_itp_ang_head(io::IO) = print(io, "[ angles ] ; cannonical angle \n")
    wr_itp_ang_comm(io::IO) = @printf(io, ";%9s%10s%10s%5s%18s%18s \n", "i", "j", "k", "f", "eq", "coef")
    function wr_itp_ang_line(io::IO, a::CGTopAngle, f::Int, e::Float64)
        @printf(io, "%10d%10d%10d%5d%18.4E%18.4E\n", a.i, a.j, a.k, f, a.a0, e * 2.0 * CAL2JOU)
    end

    wr_itp_dih_P_head(io::IO) = print(io, "[ dihedrals ] ; periodic dihedrals\n")
    wr_itp_dih_P_comm(io::IO) = @printf(io, ";%9s%10s%10s%10s%5s%18s%18s%5s\n", "i", "j", "k", "l", "f", "eq", "coef", "n")
    function wr_itp_dih_P_line(io::IO, d::CGTopDihedral, f::Int, e::Float64, n::Int)
        if n == 1
            @printf(io, "%10d%10d%10d%10d%5d%18.4E%18.4E%5d\n", d.i, d.j, d.k, d.l, f, d.t0 - 180., e * CAL2JOU, n)
        elseif n== 3
            @printf(io, "%10d%10d%10d%10d%5d%18.4E%18.4E%5d\n", d.i, d.j, d.k, d.l, f, 3*d.t0-180., e * CAL2JOU, n)
        end
    end


    wr_itp_dih_G_head(io::IO) = print(io, "[ dihedrals ] ; Gaussian dihedrals\n")
    wr_itp_dih_G_comm(io::IO) = @printf(io, ";%9s%10s%10s%10s%5s%15s%15s%15s\n", "i", "j", "k", "l", "f", "eq", "coef", "w")
    function wr_itp_dih_G_line(io::IO, d::CGTopDihedral, f::Int, e::Float64, sigma::Float64)
        @printf(io, "%10d%10d%10d%10d%5d%15.4E%15.4E%15.4E\n", d.i, d.j, d.k, d.l, f, d.t0, e * CAL2JOU, sigma)
    end

    wr_itp_dih_F_head(io::IO) = print(io, "[ dihedrals ] ; AICG2+ flexible local interation\n")
    wr_itp_dih_F_comm(io::IO) = @printf(io, ";%9s%10s%10s%10s%5s\n", "i", "j", "k", "l", "f")
    function wr_itp_dih_F_line(io::IO, d::CGTopDihedral, f::Int)
        @printf(io, "%10d%10d%10d%10d%5d\n", d.i, d.j, d.k, d.l, f)
    end

    wr_itp_contact_head(io::IO, s::String) = @printf(io, "[ pairs ] ; %s - Go-type native contact\n", s)
    wr_itp_contact_comm(io::IO) = @printf(io, ";%9s%10s%10s%15s%15s\n", "i", "j", "f", "eq", "coef")
    function wr_itp_contact_line(io::IO, c::CGTopContact, f::Int, e::Float64)
        @printf(io, "%10d%10d%10d%15.4E%15.4E\n", c.i, c.j, f, c.r0 * 0.1, e * CAL2JOU)
    end

    wr_itp_exc_head(io::IO) = print(io, "[ exclusions ] ; Genesis exclusion list\n")
    wr_itp_exc_comm(io::IO) = @printf(io, ";%9s%10s\n", "i", "j")
    wr_itp_exc_line(io::IO, c::CGTopContact) = @printf(io, "%10d%10d\n", c.i, c.j)


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
    cg_num_particles = length(top.cg_bead_name)
    wr_itp_atm_head(itp_file)
    wr_itp_atm_comm(itp_file)
    for i_bead in 1 : cg_num_particles
        wr_itp_atm_line(itp_file, 
                        i_bead,
                        top.cg_bead_type[i_bead],
                        top.cg_resid_index[i_bead],
                        top.cg_resid_name[i_bead],
                        top.cg_bead_name[i_bead],
                        AICG_ATOM_FUNC_NR,
                        top.cg_bead_charge[i_bead],
                        top.cg_bead_mass[i_bead])
    end
    print(itp_file,"\n")


    # ---------
    # [ bonds ]
    # ---------
    if length(top.top_cg_pro_bonds) + length(top.top_cg_DNA_bonds) + length(top.top_cg_RNA_bonds) > 0
        wr_itp_bnd_head(itp_file)
        wr_itp_bnd_comm(itp_file)
        
        # AICG2+ bonds
        if ff_pro == FF_pro_AICG2p
            for bond in top.top_cg_pro_bonds
                wr_itp_bnd_line(itp_file, bond, AICG_BOND_FUNC_TYPE, AICG_BOND_K)
            end
        # Clementi Go bonds
        elseif ff_pro == FF_pro_Clementi_Go
            for bond in top.top_cg_pro_bonds
                wr_itp_bnd_line(itp_file, bond, CCGO_BOND_FUNC_TYPE, CCGO_BOND_K)
            end
        end

        # 3SPN.2C bonds
        if ff_dna == FF_DNA_3SPN2C && gen_3spn_itp
            for bond in top.top_cg_DNA_bonds
                wr_itp_bnd_line(itp_file, bond, DNA3SPN_BOND_FUNC4_TYPE, DNA3SPN_BOND_K_2 / 100.0)
            end
        end

        # Structure-based RNA bonds
        if ff_rna == FF_RNA_HT
            for ( i_bond, bond ) in enumerate( top.top_cg_RNA_bonds )
                wr_itp_bnd_line(itp_file, bond, RNA_BOND_FUNC_TYPE, top.param_cg_RNA_k_bonds[i_bond])
            end
        end

        print(itp_file, "\n")
    end


    # ----------
    # [ angles ]
    # ----------

    if ff_pro == FF_pro_AICG2p
        # AICG2+ 1-3
        if length(top.top_cg_pro_aicg13) > 0
            wr_itp_13_head(itp_file)
            wr_itp_13_comm(itp_file)
            for ( i_13, a13 ) in enumerate( top.top_cg_pro_aicg13 )
                if has_toml_mod
                    if  in(a13.i, AICG2p_flexible_local) ||
                        in(a13.j, AICG2p_flexible_local) ||
                        in(a13.k, AICG2p_flexible_local)
                        continue
                    elseif in(a13.i, HPS_IDR_region) ||
                        in(a13.j, HPS_IDR_region) ||
                        in(a13.k, HPS_IDR_region)
                        continue
                    end
                end
                wr_itp_13_line(itp_file, a13, AICG_ANG_G_FUNC_TYPE, top.param_cg_pro_e_13[i_13], AICG_13_SIGMA)
            end
            print(itp_file, "\n")
        end
        # AICG2+ flexible
        if length(top.top_cg_pro_angles) > 0
            wr_itp_ang_f_head(itp_file)
            wr_itp_ang_f_comm(itp_file)
            for ang in top.top_cg_pro_angles
                if  in(ang.i, HPS_IDR_region) ||
                    in(ang.j, HPS_IDR_region) ||
                    in(ang.k, HPS_IDR_region)
                    continue
                end
                wr_itp_ang_f_line(itp_file, ang, AICG_ANG_F_FUNC_TYPE)
            end
            print(itp_file, "\n")
        end
    # Clementi Go angle
    elseif ff_pro == FF_pro_Clementi_Go
        if length(top.top_cg_pro_angles) > 0
            wr_itp_ang_head(itp_file)
            wr_itp_ang_comm(itp_file)
            for ang in top.top_cg_pro_angles
                wr_itp_ang_line(itp_file, ang, CCGO_ANG_FUNC_TYPE, CCGO_ANGL_K)
            end
            print(itp_file, "\n")
        end
    end

    # 3SPN.2C angles
    if ff_dna == FF_DNA_3SPN2C
        if length(top.top_cg_DNA_angles) > 0 && gen_3spn_itp
            wr_itp_ang_head(itp_file)
            wr_itp_ang_comm(itp_file)
            for ( i_ang, ang ) in enumerate( top.top_cg_DNA_angles )
                wr_itp_ang_line(itp_file, ang, DNA3SPN_ANG_FUNC_TYPE, top.param_cg_DNA_k_angles[i_ang])
            end
            print(itp_file, "\n")
        end
    end

    # RNA structure-based angles
    if ff_rna == FF_RNA_HT
        if length(top.top_cg_RNA_angles) > 0
            wr_itp_ang_head(itp_file)
            wr_itp_ang_comm(itp_file)
            for ( i_ang, ang ) in enumerate( top.top_cg_RNA_angles )
                wr_itp_ang_line(itp_file, ang, RNA_ANG_FUNC_TYPE, top.param_cg_RNA_k_angles[i_ang])
            end
            print(itp_file, "\n")
        end
    end


    # -------------
    # [ dihedrals ]
    # -------------

    if ff_pro == FF_pro_AICG2p
        # AICG2+ Gaussian dihedrals
        if length(top.top_cg_pro_aicg14) > 0
            wr_itp_dih_G_head(itp_file)
            wr_itp_dih_G_comm(itp_file)
            for ( i_dih, dih ) in enumerate( top.top_cg_pro_aicg14 )
                if has_toml_mod
                    if  in(dih.i, AICG2p_flexible_local) ||
                        in(dih.j, AICG2p_flexible_local) ||
                        in(dih.k, AICG2p_flexible_local) || 
                        in(dih.l, AICG2p_flexible_local)
                        continue
                    elseif  in(dih.i, HPS_IDR_region) ||
                        in(dih.j, HPS_IDR_region) ||
                        in(dih.k, HPS_IDR_region) ||
                        in(dih.l, HPS_IDR_region)
                        continue
                    end
                end
                wr_itp_dih_G_line(itp_file, dih, AICG_DIH_G_FUNC_TYPE, top.param_cg_pro_e_14[i_dih], AICG_14_SIGMA)
            end
            print(itp_file, "\n")
        end
        # AICG2+ flexible dihedrals
        if length(top.top_cg_pro_dihedrals) > 0
            wr_itp_dih_F_head(itp_file)
            wr_itp_dih_F_comm(itp_file)
            for dih in top.top_cg_pro_dihedrals
                if  in(dih.i, HPS_IDR_region) ||
                    in(dih.j, HPS_IDR_region) ||
                    in(dih.k, HPS_IDR_region) ||
                    in(dih.l, HPS_IDR_region)
                    continue
                end
                wr_itp_dih_F_line(itp_file, dih, AICG_DIH_F_FUNC_TYPE)
            end
            print(itp_file, "\n")
        end
    # Clementi Go dihedral
    elseif ff_pro == FF_pro_Clementi_Go
        if length(top.top_cg_pro_dihedrals) > 0
            wr_itp_dih_P_head(itp_file)
            wr_itp_dih_P_comm(itp_file)
            for dih in top.top_cg_pro_dihedrals
                wr_itp_dih_P_line(itp_file, dih, CCGO_DIH_P_FUNC_TYPE, CCGO_DIHE_K_1, 1)
            end
            for dih in top.top_cg_pro_dihedrals
                wr_itp_dih_P_line(itp_file, dih, CCGO_DIH_P_FUNC_TYPE, CCGO_DIHE_K_3, 3)
            end
            print(itp_file, "\n")
        end
    end

    if ff_dna == FF_DNA_3SPN2C
        # 3SPN.2C Gaussian dihedrals
        if length(top.top_cg_DNA_dih_Gaussian) > 0 && gen_3spn_itp
            wr_itp_dih_G_head(itp_file)
            wr_itp_dih_G_comm(itp_file)
            for dih in top.top_cg_DNA_dih_Gaussian
                wr_itp_dih_G_line(itp_file, dih, DNA3SPN_DIH_G_FUNC_TYPE, DNA3SPN_DIH_G_K, DNA3SPN_DIH_G_SIGMA)
            end
            print(itp_file, "\n")
        end

        # 3SPN.2C Periodic dihedrals
        if length(top.top_cg_DNA_dih_periodic) > 0 && gen_3spn_itp
            wr_itp_dih_P_head(itp_file)
            wr_itp_dih_P_comm(itp_file)
            for dih in top.top_cg_DNA_dih_periodic
                wr_itp_dih_P_line(itp_file, dih, DNA3SPN_DIH_P_FUNC_TYPE, DNA3SPN_DIH_P_K, DNA3SPN_DIH_P_FUNC_PERI)
            end
            print(itp_file, "\n")
        end
    end

    # RNA structure-based Periodic dihedrals
    if ff_rna == FF_RNA_HT
        if length(top.top_cg_RNA_dihedrals) > 0
            wr_itp_dih_P_head(itp_file)
            wr_itp_dih_P_comm(itp_file)
            for ( i_dih, dih ) in enumerate( top.top_cg_RNA_dihedrals )
                wr_itp_dih_P_line(itp_file, dih, RNA_DIH_FUNC_TYPE, top.param_cg_RNA_k_dihedrals[i_dih], 1)
            end
            for ( i_dih, dih ) in enumerate( top.top_cg_RNA_dihedrals )
                wr_itp_dih_P_line(itp_file, dih, RNA_DIH_FUNC_TYPE, top.param_cg_RNA_k_dihedrals[i_dih] / 2, 3)
            end
            print(itp_file, "\n")
        end
    end

    # ---------
    # [ pairs ]
    # ---------

    # print protein Go-type native contacts
    if ff_pro == FF_pro_AICG2p
        if length(top.top_cg_pro_go_contact) > 0
            wr_itp_contact_head(itp_file, "AICG2+")
            wr_itp_contact_comm(itp_file)
            for (i_c, c) in enumerate(top.top_cg_pro_go_contact)
                if has_toml_mod
                    if  in(c.i, AICG2p_flexible_nonlocal) ||
                        in(c.j, AICG2p_flexible_nonlocal)
                        continue
                    elseif in(c.i, HPS_IDR_region) ||
                        in(c.j, HPS_IDR_region)
                        continue
                    end
                end
                wr_itp_contact_line(itp_file, c, AICG_CONTACT_FUNC_TYPE, top.param_cg_pro_e_contact[i_c])
            end
            print(itp_file, "\n")
        end
    # Clementi Go native contacts
    elseif ff_pro == FF_pro_Clementi_Go
        if length(top.top_cg_pro_go_contact) > 0
            wr_itp_contact_head(itp_file, "Clementi-Go")
            wr_itp_contact_comm(itp_file)
            for c in top.top_cg_pro_go_contact
                wr_itp_contact_line(itp_file, c, CCGO_CONTACT_FUNC_TYPE, CCGO_NATIVE_EPSILON * ccgo_contact_scale)
            end
            print(itp_file, "\n")
        end
    end

    # print RNA HT-type native contacts
    if ff_rna == FF_RNA_HT
        if length(top.top_cg_RNA_base_stack) + length(top.top_cg_RNA_base_pair) + length(top.top_cg_RNA_other_contact) > 0
            wr_itp_contact_head(itp_file, "RNA-RNA")
            wr_itp_contact_comm(itp_file)
            for (i_c, c) in enumerate(top.top_cg_RNA_base_stack)
                wr_itp_contact_line(itp_file, c, RNA_CONTACT_FUNC_TYPE, top.param_cg_RNA_e_base_stack[i_c])
            end
            for (i_c, c) in enumerate(top.top_cg_RNA_base_pair)
                wr_itp_contact_line(itp_file, c, RNA_CONTACT_FUNC_TYPE, top.param_cg_RNA_e_base_pair[i_c])
            end
            for (i_c, c) in enumerate(top.top_cg_RNA_other_contact)
                wr_itp_contact_line(itp_file, c, RNA_CONTACT_FUNC_TYPE, top.param_cg_RNA_e_other_contact[i_c])
            end
            print(itp_file, "\n")
        end
    end


    # print protein-RNA native contacts
    if ( ff_pro == FF_pro_AICG2p || ff_pro == FF_pro_Clementi_Go ) && ff_rna == FF_RNA_HT
        if length(top.top_cg_pro_RNA_contact) > 0
            wr_itp_contact_head(itp_file, "Protein-RNA")
            wr_itp_contact_comm(itp_file)
            for (i_c, c) in enumerate(top.top_cg_pro_RNA_contact)
                wr_itp_contact_line(itp_file, c, RNP_CONTACT_FUNC_TYPE, top.param_cg_pro_RNA_e_contact[i_c])
            end
            print(itp_file, "\n")
        end
    end


    # ---------------------
    #        [ exclusions ]
    # ---------------------

    # print Protein exclusion list
    if ff_pro == FF_pro_AICG2p || ff_pro == FF_pro_Clementi_Go
        if length(top.top_cg_pro_go_contact) > 0
            wr_itp_exc_head(itp_file)
            wr_itp_exc_comm(itp_file)
            for c in top.top_cg_pro_go_contact
                if has_toml_mod
                    if  in(c.i, AICG2p_flexible_nonlocal) ||
                        in(c.j, AICG2p_flexible_nonlocal)
                        continue
                    elseif in(c.i, HPS_IDR_region) ||
                        in(c.j, HPS_IDR_region)
                        continue
                    end
                end
                wr_itp_exc_line(itp_file, c)
            end
            print(itp_file, "\n")
        end
    end

    # print RNA exclusion list
    if ff_rna == FF_RNA_HT
        if length(top.top_cg_RNA_base_stack) + length(top.top_cg_RNA_base_pair) + length(top.top_cg_RNA_other_contact) > 0
            wr_itp_exc_head(itp_file)
            wr_itp_exc_comm(itp_file)
            for c in top.top_cg_RNA_base_stack
                wr_itp_exc_line(itp_file, c)
            end
            for c in top.top_cg_RNA_base_pair
                wr_itp_exc_line(itp_file, c)
            end
            for c in top.top_cg_RNA_other_contact
                wr_itp_exc_line(itp_file, c)
            end
            print(itp_file, "\n")
        end
    end

    # print protein-RNA exclusion contacts
    if ( ff_pro == FF_pro_AICG2p || ff_pro == FF_pro_Clementi_Go ) && ff_rna == FF_RNA_HT
        if length(top.top_cg_pro_RNA_contact) > 0
            wr_itp_exc_head(itp_file)
            wr_itp_exc_comm(itp_file)
            for c in top.top_cg_pro_RNA_contact
                wr_itp_exc_line(itp_file, c)
            end
            print(itp_file, "\n")
        end
    end

    close(itp_file)

    println(">           ... .top: DONE!")
end


function write_cg_grotop_pwmcos(top::CGTopology, force_field::ForceFieldCG, system_name::String, args::Dict{String, Any})

    appendto_filename = args["patch"]
    pwmcos_gamma      = args["pwmcos-scale"]
    pwmcos_epsil      = args["pwmcos-shift"]

    if length( appendto_filename ) == 0
        itp_pwmcos_name = system_name * "_cg_pwmcos.itp_patch"
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
    for p in top.top_cg_pro_DNA_pwmcos
        @printf(itp_pwmcos_file,
                "%6d %3d %8.5f %8.3f %8.3f %8.3f%12.6f%12.6f%12.6f%12.6f%8.3f%8.3f \n",
                p.i, PWMCOS_FUNC_TYPE, p.r0 * 0.1, p.t1, p.t2, p.t3,
                p.eA, p.eC, p.eG, p.eT, pwmcos_gamma, pwmcos_epsil)
    end
    print(itp_pwmcos_file, "\n")

    close(itp_pwmcos_file)
    println(">           ... ", itp_pwmcos_name, " pwmcos.itp: DONE!")

end

# ==================================
# General Topology in Gromacs format
# ==================================

function read_groitp(itp_filename::AbstractString)

    mol_name = ""
    nonlocal_interval = 0
    num_atom = 0

    top_atoms = Vector{GenTopAtom}(undef, 0)
    top_bonds = Vector{GenTopBond}(undef, 0)
    top_angles = Vector{GenTopAngle}(undef, 0)
    top_dihedrals = Vector{GenTopDihedral}(undef, 0)
    top_pairs = Vector{GenTopPair}(undef, 0)
    top_exclusions = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos = Vector{GenTopPWMcos}(undef, 0)

    function read_top_atoms(line::AbstractString)
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
                              a_name, f_type, charge, mass)
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


    # ---------
    # main part
    # ---------
    section_name = ""
    for line in eachline(itp_filename)
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
                read_top_atoms(line)
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
                                 top_pwmcos)

    return new_top_mol
end

function read_grotop(top_filename::String)

    sys_name = ""
    num_atom = 0
    mol_id   = 0

    top_default_params = GenTopDefault(0, 0, false, 0.0, 0.0) 
    top_default_atomtype = Vector{GenTopAtomType}(undef, 0)
    top_default_CGDNA_bp = Vector{GenTopCGDNABasepairType}(undef, 0)
    top_default_CGDNA_bs = Vector{GenTopCGDNABasestackType}(undef, 0)
    top_default_CGDNA_cs = Vector{GenTopCGDNABasecrossType}(undef, 0)
    top_default_CGDNA_exv = Vector{GenTopCGDNAExvType}(undef, 0)
    top_default_CGPro_flx_angle = Vector{GenTopCGProAICGFlexAngleType}(undef, 0)
    top_default_CGPro_flx_dihedral = Vector{GenTopCGProAICGFlexDihedralType}(undef, 0)

    global_index_2_local_index = Vector{Int}(undef, 0)
    global_index_2_local_molid = Vector{Int}(undef, 0)
    top_atoms = Vector{GenTopAtom}(undef, 0)
    top_bonds = Vector{GenTopBond}(undef, 0)
    top_angles = Vector{GenTopAngle}(undef, 0)
    top_dihedrals = Vector{GenTopDihedral}(undef, 0)
    top_pairs = Vector{GenTopPair}(undef, 0)
    top_exclusions = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos = Vector{GenTopPWMcos}(undef, 0)
    top_mol_list = Vector{GenTopMolList}(undef, 0)

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

        if startswith(line, "#include")
            mol_file_name = strip(line[9:end], ['\"', '\'', ' '])
            mol_file_basename = basename(mol_file_name)
            if in(mol_file_basename, ["atom_types.itp", "flexible_local_angle.itp", "flexible_local_dihedral.itp"])
                continue
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
                                   t.mass)
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

function write_cg_psf(top::CGTopology, system_name::String, args::Dict{String, Any}=Dict{String, Any}())

    psf_name = system_name * "_cg.psf"
    psf_file = open(psf_name, "w")

    cg_num_particles = length(top.cg_resid_name)

    @printf(psf_file, "PSF CMAP \n\n")
    @printf(psf_file, "      3 !NTITLE \n")
    @printf(psf_file, "REMARKS PSF file created with Julia. \n")
    @printf(psf_file, "REMARKS System: %s  \n", system_name)
    @printf(psf_file, "REMARKS ======================================== \n")
    @printf(psf_file, "       \n")

    psf_atom_line = " %6d %3s %5d %3s %3s %5s  %10.6f  %10.6f          0 \n"
    chain_id_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"

    @printf(psf_file, " %6d !NATOM \n", cg_num_particles)
    for i_bead in 1 : cg_num_particles
        @printf(psf_file, " %6d %3s %5d %3s %3s %5s  %10.6f  %10.6f          0 \n",
                i_bead,
                chain_id_set[top.cg_chain_id[i_bead]],
                top.cg_resid_index[i_bead],
                top.cg_resid_name[i_bead],
                top.cg_bead_name[i_bead],
                top.cg_bead_type[i_bead],
                top.cg_bead_charge[i_bead],
                top.cg_bead_mass[i_bead])
    end
    print(psf_file,"\n")

    close(psf_file)
    println(">           ... .psf: DONE!")

end

function write_psf(top::GenTopology, sys_name::String="", args::Dict{String, Any}=Dict{String, Any}())

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
    chain_id_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"

    @printf(psf_file, " %6d !NATOM \n", cg_num_particles)
    for i_bead in 1 : cg_num_particles
        @printf(psf_file, " %6d %3s %5d %3s %3s %5s  %10.6f  %10.6f          0 \n",
                i_bead,
                chain_id_set[top.global_index_2_local_molid[i_bead]],
                top.top_atoms[i_bead].residue_index,
                top.top_atoms[i_bead].residue_name,
                top.top_atoms[i_bead].atom_name,
                top.top_atoms[i_bead].atom_type,
                top.top_atoms[i_bead].charge,
                top.top_atoms[i_bead].mass)
    end
    print(psf_file,"\n")

    close(psf_file)
    println(">           ... .psf: DONE!")

end

