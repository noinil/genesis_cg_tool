###############################################################################
#                         Coarse Graining Biomolecules                        #
###############################################################################

using Printf

function coarse_graining(aa_molecule::AAMolecule, force_field::ForceFieldCG, args)

    # -----------------
    # Parsing arguments
    # -----------------
    pdb_name                = get(args, "pdb", "random.pdb")
    scale_scheme            = get(args, "aicg-scale", 1)
    protein_charge_filename = get(args, "respac", "")
    pfm_filename            = get(args, "pfm", "")
    do_debug                = get(args, "debug", false)
    do_output_log           = get(args, "log", false)

    # ===============
    # Step 0: numbers
    # ===============
    i_step     = 0

    ff_pro     = force_field.ff_protein
    ff_dna     = force_field.ff_DNA
    ff_rna     = force_field.ff_RNA
    ff_pro_dna = force_field.ff_protein_DNA


    aa_atom_name  = aa_molecule.atom_names
    aa_coor       = aa_molecule.atom_coors
    aa_residues   = aa_molecule.residues
    aa_chains     = aa_molecule.chains

    aa_num_atom    = length(aa_atom_name)
    aa_num_residue = length(aa_residues)
    aa_num_chain   = length(aa_chains)

    num_chain_pro  = 0
    num_chain_DNA  = 0
    num_chain_RNA  = 0

    # ===============================
    # Step 1: find out molecule types
    # ===============================
    i_step += 1
    println("============================================================")
    println("> Step $(i_step): estimate CG particle number for every chain.")

    cg_num_particles = 0
    cg_chain_length  = zeros(Int, aa_num_chain)

    for i_chain = 1 : aa_num_chain
        chain = aa_chains[i_chain]
        mol_type = chain.moltype

        # ----------------------
        # determine cg chain len
        # ----------------------
        n_res = length(chain.residues)
        if mol_type == MOL_DNA
            n_particles = 3 * n_res - 1
            num_chain_DNA += 1
        elseif mol_type == MOL_RNA
            n_particles = 3 * n_res - 1
            num_chain_RNA += 1
        elseif mol_type == MOL_PROTEIN
            n_particles = n_res
            num_chain_pro += 1
        else
            n_particles = 0
        end
        cg_chain_length[i_chain] = n_particles
        cg_num_particles += n_particles
        @printf("          > Chain %3d | %7s \n", i_chain, MOL_TYPE_LIST[ mol_type ])
    end

    println("------------------------------------------------------------")
    @printf("          In total: %5d protein chains,\n", num_chain_pro)
    @printf("                    %5d DNA strands,\n", num_chain_DNA)
    @printf("                    %5d RNA strands.\n", num_chain_RNA)

    # ===========================
    # Step 2: Assign CG particles
    # ===========================
    i_step += 1
    println("============================================================")
    println("> Step $(i_step): assign coarse-grained particles.")

    cg_residues = []
    cg_chains   = []

    i_offset_cg_particle = 0
    i_offset_cg_residue  = 0

    for i_chain in 1:aa_num_chain
        chain = aa_chains[i_chain]
        mol_type = chain.moltype
        seg_name = chain.segname

        i_bead = i_offset_cg_particle
        i_resi = i_offset_cg_residue
        if mol_type == MOL_PROTEIN
            for i_res in chain.residues
                cg_idx = []
                aa_res_name = aa_residues[i_res].name
                res_name = RES_NAME_PROTEIN_DICT[aa_res_name]
                for i_atom in aa_residues[i_res].atoms
                    atom_name = aa_atom_name[i_atom]
                    if atom_name[1] == 'H'
                        continue
                    else
                        push!(cg_idx, i_atom)
                    end
                end
                i_bead += 1
                i_resi += 1
                push!(cg_residues, AACGResidue(i_resi, res_name, "CA", cg_idx))
            end
        elseif mol_type == MOL_DNA
            tmp_atom_index_O3p = 0
            for (i_local_index, i_res) in enumerate( chain.residues )
                aa_res_name = aa_residues[i_res].name
                res_name = RES_NAME_DNA_DICT[aa_res_name]
                cg_DP_idx = [tmp_atom_index_O3p]
                cg_DS_idx = []
                cg_DB_idx = []
                for i_atom in aa_residues[i_res].atoms
                    atom_name = aa_atom_name[i_atom]
                    if atom_name[1] == 'H'
                        continue
                    elseif in(atom_name, ATOM_NAME_LIST_DP)
                        push!(cg_DP_idx, i_atom)
                    elseif in(atom_name, ATOM_NAME_LIST_DS)
                        push!(cg_DS_idx, i_atom)
                    elseif atom_name == "O3'"
                        tmp_atom_index_O3p = i_atom
                    else
                        push!(cg_DB_idx, i_atom)
                    end
                end
                i_resi += 1
                if i_local_index > 1
                    i_bead += 1
                    push!(cg_residues, AACGResidue(i_resi, res_name, "DP", cg_DP_idx))
                end
                i_bead += 1
                push!(cg_residues, AACGResidue(i_resi, res_name, "DS", cg_DS_idx))
                i_bead += 1
                push!(cg_residues, AACGResidue(i_resi, res_name, "DB", cg_DB_idx))
            end
        elseif mol_type == MOL_RNA
            for (i_local_index, i_res) in enumerate( chain.residues )
                aa_res_name = aa_residues[i_res].name
                res_name = RES_NAME_RNA_DICT[aa_res_name]
                cg_RP_idx = []
                cg_RS_idx = []
                cg_RB_idx = []
                for i_atom in aa_residues[i_res].atoms
                    atom_name = aa_atom_name[i_atom]
                    if atom_name[1] == 'H'
                        continue
                    elseif in(atom_name, ATOM_NAME_LIST_RP)
                        push!(cg_RP_idx, i_atom)
                    elseif in(atom_name, ATOM_NAME_LIST_RS)
                        push!(cg_RS_idx, i_atom)
                    else
                        push!(cg_RB_idx, i_atom)
                    end
                end
                i_resi += 1
                if i_local_index > 1
                    i_bead += 1
                    push!(cg_residues, AACGResidue(i_resi, res_name, "RP", cg_RP_idx))
                end
                i_bead += 1
                push!(cg_residues, AACGResidue(i_resi, res_name, "RS", cg_RS_idx))
                i_bead += 1
                push!(cg_residues, AACGResidue(i_resi, res_name, "RB", cg_RB_idx))
            end
        end
        push!(cg_chains, AACGChain(i_offset_cg_particle + 1, i_bead, mol_type, seg_name))
        i_offset_cg_particle += cg_chain_length[i_chain]
        i_offset_cg_residue  += length(chain.residues)
    end

    for i_chain in 1:aa_num_chain
        @printf("          > Chain %3d | # particles: %5d | %5d -- %5d \n",
                i_chain, cg_chain_length[i_chain],
                cg_chains[i_chain].first, cg_chains[i_chain].last)
    end

    println("------------------------------------------------------------")
    println("          In total: $(cg_num_particles) CG particles.")




    # =========================================================================
    #        ____ ____   _____ ___  ____   ___  _     ___   ______   __
    #       / ___/ ___| |_   _/ _ \|  _ \ / _ \| |   / _ \ / ___\ \ / /
    #      | |  | |  _    | || | | | |_) | | | | |  | | | | |  _ \ V /
    #      | |__| |_| |   | || |_| |  __/| |_| | |__| |_| | |_| | | |
    #       \____\____|   |_| \___/|_|    \___/|_____\___/ \____| |_|
    #
    # =========================================================================

    cg_resid_name  = fill("    ", cg_num_particles)
    cg_resid_index = zeros(Int, cg_num_particles)
    cg_bead_name   = fill("    ", cg_num_particles)
    cg_bead_type   = fill("    ", cg_num_particles)
    cg_bead_charge = zeros(Float64, cg_num_particles)
    cg_bead_mass   = zeros(Float64, cg_num_particles)
    cg_bead_coor   = zeros(Float64, (3, cg_num_particles))
    cg_chain_id    = zeros(Int, cg_num_particles)
    cg_seg_name    = fill("    ", cg_num_particles)

    # protein
    top_cg_pro_bonds         = Vector{CGTopBond}(undef, 0)
    top_cg_pro_angles        = Vector{CGTopAngle}(undef, 0)
    top_cg_pro_dihedrals     = Vector{CGTopDihedral}(undef, 0)
    top_cg_pro_aicg13        = Vector{CGTopAngle}(undef, 0)
    top_cg_pro_aicg14        = Vector{CGTopDihedral}(undef, 0)
    top_cg_pro_go_contact    = Vector{CGTopContact}(undef, 0)

    param_cg_pro_e_13        = []
    param_cg_pro_e_14        = []
    param_cg_pro_e_contact   = []

    # DNA
    top_cg_DNA_bonds         = Vector{CGTopBond}(undef, 0)
    top_cg_DNA_angles        = Vector{CGTopAngle}(undef, 0)
    top_cg_DNA_dih_Gaussian  = Vector{CGTopDihedral}(undef, 0)
    top_cg_DNA_dih_periodic  = Vector{CGTopDihedral}(undef, 0)
    param_cg_DNA_k_angles    = []

    # RNA
    top_cg_RNA_bonds             = Vector{CGTopBond}(undef, 0)
    top_cg_RNA_angles            = Vector{CGTopAngle}(undef, 0)
    top_cg_RNA_dihedrals         = Vector{CGTopDihedral}(undef, 0)
    top_cg_RNA_base_stack        = Vector{CGTopContact}(undef, 0)
    top_cg_RNA_base_pair         = Vector{CGTopContact}(undef, 0)
    top_cg_RNA_other_contact     = Vector{CGTopContact}(undef, 0)
    param_cg_RNA_k_bonds         = []
    param_cg_RNA_k_angles        = []
    param_cg_RNA_k_dihedrals     = []
    param_cg_RNA_e_base_stack    = []
    param_cg_RNA_e_base_pair     = []
    param_cg_RNA_e_other_contact = []

    # protein-DNA
    top_cg_pro_DNA_pwmcos    = Vector{CGTopPWMcos}(undef, 0)

    # protein-RNA
    top_cg_pro_RNA_contact   = Vector{CGTopContact}(undef, 0)
    param_cg_pro_RNA_e_contact = []

    # --------------------
    # geometric properties
    # --------------------
    # center of geometry
    geo_centroid               = zeros(Float64, (3, aa_num_chain))
    # radius of gyration
    geo_radius_of_gyration     = zeros(Float64, aa_num_chain)
    # radius of circumsphere
    geo_radius_of_circumsphere = zeros(Float64, aa_num_chain)

    # =============================
    # Step 4: CG model for proteins
    # =============================
    #                  _       _
    #  _ __  _ __ ___ | |_ ___(_)_ __
    # | '_ \| '__/ _ \| __/ _ \ | '_ \
    # | |_) | | | (_) | ||  __/ | | | |
    # | .__/|_|  \___/ \__\___|_|_| |_|
    # |_|
    #
    # =================================

    num_cg_pro_contact_all = 0
    num_cg_pro_contact_intra = 0
    num_cg_pro_contact_inter = 0
    if num_chain_pro > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): processing proteins.")

        # --------------------------------
        # Step 4.1: find out C-alpha atoms
        # --------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine CA mass, charge, and coordinates.")

        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain.first : chain.last
                res_name = cg_residues[i_res].res_name
                for i_atom in cg_residues[i_res].atoms
                    if aa_atom_name[i_atom] == "CA"
                        cg_resid_name[i_res]   = res_name
                        cg_resid_index[i_res]  = cg_residues[i_res].res_idx
                        cg_bead_name[i_res]    = "CA"
                        cg_bead_type[i_res]    = res_name
                        cg_bead_charge[i_res]  = RES_CHARGE_DICT[res_name]
                        cg_bead_mass[i_res]    = RES_MASS_DICT[res_name]
                        cg_bead_coor[:, i_res] = aa_coor[:, i_atom]
                        cg_chain_id[i_res]     = i_chain
                        cg_seg_name[i_res]     = chain.segname
                        break
                    end
                end
            end
        end

        if length(protein_charge_filename) > 0
            try
                for line in eachline(protein_charge_filename)
                    charge_data = split(line)
                    if length(charge_data) < 1
                        continue
                    end
                    i = parse(Int, charge_data[1])
                    c = parse(Float64, charge_data[2])
                    cg_bead_charge[i] = c
                end
            catch e
                println(e)
                error("ERROR in user-defined charge distribution.\n")
            end
        end
        println(">           ... DONE!")

        # ----------------------------------------
        # Step 4.1.1: determine geometric features
        # ----------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1.1: determine Centroid, Rg, Rc...")

        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            # centroid
            coor_centroid = zeros(Float64, 3)
            for i_res in chain.first : chain.last
                coor_centroid += cg_bead_coor[:, i_res]
            end
            coor_centroid /= (chain.last - chain.first + 1)
            geo_centroid[:, i_chain] = coor_centroid

            # Rg
            tmp_dist = 0
            tmp_dist_sq_sum = 0
            for i_res in chain.first : chain.last
                vec_from_center   = cg_bead_coor[:, i_res] - coor_centroid
                vec_norm_tmp      = norm(vec_from_center)
                tmp_dist          = vec_norm_tmp > tmp_dist ? vec_norm_tmp : tmp_dist
                tmp_dist_sq_sum  += vec_norm_tmp * vec_norm_tmp
            end
            rg = sqrt(tmp_dist_sq_sum / (chain.last - chain.first + 1))
            rc = tmp_dist
            geo_radius_of_gyration[i_chain]     = rg
            geo_radius_of_circumsphere[i_chain] = rc
        end

        # -----------------------------
        # Step 4.2: CG protein topology
        # -----------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).2: CG protein topology.")
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.1: CG protein local interactions.")
        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain.first : chain.last - 1
                coor1 = cg_bead_coor[:, i_res]
                coor2 = cg_bead_coor[:, i_res + 1]
                dist12 = compute_distance(coor1, coor2)
                tmp_top_bond = CGTopBond(i_res, i_res + 1, dist12)
                push!(top_cg_pro_bonds, tmp_top_bond)
            end
        end
        println(">           ... Bond: DONE!")

        e_ground_local = 0.0
        e_ground_13    = 0.0
        num_angle      = 0
        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain.first : chain.last - 2
                coor1    = cg_bead_coor[:, i_res]
                coor2    = cg_bead_coor[:, i_res + 1]
                coor3    = cg_bead_coor[:, i_res + 2]
                dist13   = compute_distance(coor1, coor3)
                angle123 = compute_angle(coor1, coor2, coor3)
                tmp_top_angle = CGTopAngle(i_res, i_res + 1, i_res + 2, angle123)
                push!(top_cg_pro_angles, tmp_top_angle)
                tmp_top_angle = CGTopAngle(i_res, i_res + 1, i_res + 2, dist13)
                push!(top_cg_pro_aicg13, tmp_top_angle)
                # count AICG2+ 1-3 interaction atomic contact
                contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ].atoms,
                                                           cg_residues[ i_res + 2 ].atoms,
                                                           cg_resid_name[i_res],
                                                           cg_resid_name[i_res + 2],
                                                           aa_atom_name,
                                                           aa_coor)

                # calculate AICG2+ 1-3 interaction pairwise energy
                e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
                if e_local > AICG_ENE_UPPER_LIM
                    e_local = AICG_ENE_UPPER_LIM
                end
                if e_local < AICG_ENE_LOWER_LIM
                    e_local = AICG_ENE_LOWER_LIM
                end
                e_ground_local += e_local
                e_ground_13    += e_local
                num_angle      += 1
                push!(param_cg_pro_e_13, e_local)
            end
        end
        println(">           ... Angle: DONE!")

        e_ground_14 = 0.0
        num_dih = 0
        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain.first : chain.last - 3
                coor1 = cg_bead_coor[:, i_res]
                coor2 = cg_bead_coor[:, i_res + 1]
                coor3 = cg_bead_coor[:, i_res + 2]
                coor4 = cg_bead_coor[:, i_res + 3]
                dihed = compute_dihedral(coor1, coor2, coor3, coor4)
                tmp_top_dihe = CGTopDihedral(i_res, i_res + 1, i_res + 2, i_res + 3, dihed)
                push!(top_cg_pro_dihedrals, tmp_top_dihe)
                push!(top_cg_pro_aicg14, tmp_top_dihe)

                # count AICG2+ dihedral atomic contact
                contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ].atoms,
                                                           cg_residues[ i_res + 3 ].atoms,
                                                           cg_resid_name[i_res],
                                                           cg_resid_name[i_res + 3],
                                                           aa_atom_name,
                                                           aa_coor)

                # calculate AICG2+ dihedral pairwise energy
                e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
                if e_local > AICG_ENE_UPPER_LIM
                    e_local = AICG_ENE_UPPER_LIM
                end
                if e_local < AICG_ENE_LOWER_LIM
                    e_local = AICG_ENE_LOWER_LIM
                end
                e_ground_local += e_local
                e_ground_14    += e_local
                num_dih      += 1
                push!(param_cg_pro_e_14, e_local)
            end
        end
        println(">           ... Dihedral: DONE!")

        # ------------------------
        # Normalize local energies
        # ------------------------
        e_ground_local /= (num_angle + num_dih)
        e_ground_13    /= num_angle
        e_ground_14    /= num_dih

        if scale_scheme == 0
            for i in 1:length(param_cg_pro_e_13)
                param_cg_pro_e_13[i] *= AICG_13_AVE / e_ground_13
            end
            for i in 1:length(param_cg_pro_e_14)
                param_cg_pro_e_14[i] *= AICG_14_AVE / e_ground_14
            end
        elseif scale_scheme == 1
            for i in 1:length(param_cg_pro_e_13)
                param_cg_pro_e_13[i] *= -AICG_13_GEN
            end
            for i in 1:length(param_cg_pro_e_14)
                param_cg_pro_e_14[i] *= -AICG_14_GEN
            end
        end

        # -----------------------
        # Go type native contacts
        # -----------------------
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.2: Looking for native contacts.")
        e_ground_contact = 0.0
        num_contact = 0

        # intra-molecular contacts
        @printf("%11s Calculating intra-molecular contacts... \n", " ")
        @printf("              ... chain   : %32s", " ")
        i_progress_count = 0
        for i_chain in 1:aa_num_chain

            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            # -----------------
            # show progress bar
            # -----------------
            i_progress_count += 1
            print("\b"^32)
            progress_percent = trunc(Int, i_progress_count / num_chain_pro * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %2d / %2d ", progress_bar, i_progress_count, num_chain_pro)
            # ------------------

            for i_res in chain.first : chain.last - 4
                coor_cai = cg_bead_coor[:, i_res]
                for j_res in i_res + 4 : chain.last
                    coor_caj = cg_bead_coor[:, j_res]
                    if is_protein_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                        native_dist = compute_distance(coor_cai, coor_caj)
                        num_cg_pro_contact_all += 1
                        num_cg_pro_contact_intra += 1
                        tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                        push!(top_cg_pro_go_contact, tmp_top_cnt)

                        # count AICG2+ atomic contact
                        contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ].atoms,
                                                                   cg_residues[ j_res ].atoms,
                                                                   cg_resid_name[i_res],
                                                                   cg_resid_name[j_res],
                                                                   aa_atom_name,
                                                                   aa_coor)

                        # calculate AICG2+ pairwise energy
                        e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
                        if e_local > AICG_ENE_UPPER_LIM
                            e_local = AICG_ENE_UPPER_LIM
                        end
                        if e_local < AICG_ENE_LOWER_LIM
                            e_local = AICG_ENE_LOWER_LIM
                        end
                        e_ground_contact += e_local
                        num_contact      += 1
                        push!(param_cg_pro_e_contact, e_local)
                    end
                end
            end
        end
        print("\n              ... intra-molecular contacts: DONE! \n")

        # inter-molecular ( protein-protein ) contacts
        if num_chain_pro > 1
            @printf("%11s Calculating inter-molecular contacts... \n", " ")
            @printf("              ... progress: %32s", " ")
            for i_chain in 1 : aa_num_chain - 1
                chain1 = cg_chains[i_chain]
    
                # -----------------
                # show progress bar
                # -----------------
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / ( aa_num_chain - 1 ) * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / ( aa_num_chain - 1 ) * 100)
                # ------------------
    
                if chain1.moltype != MOL_PROTEIN
                    continue
                end
    
                centroid_chain1 = geo_centroid[:, i_chain]
                radius_of_circ1 = geo_radius_of_circumsphere[i_chain]

                for j_chain in i_chain + 1 : aa_num_chain
                    chain2 = cg_chains[j_chain]
    
                    if chain2.moltype != MOL_PROTEIN
                        continue
                    end

                    centroid_chain2 = geo_centroid[:, j_chain]
                    radius_of_circ2 = geo_radius_of_circumsphere[j_chain]

                    cent_cent_dist = compute_distance(centroid_chain1, centroid_chain2)
                    if cent_cent_dist > radius_of_circ1 + radius_of_circ2 + CG_CONTACT_CUTOFF
                        continue
                    end
    
                    for i_res in chain1.first : chain1.last
                        coor_cai = cg_bead_coor[:, i_res]

                        cai_centj_dist = compute_distance(coor_cai, centroid_chain2)
                        if cai_centj_dist > radius_of_circ2 + CG_CONTACT_CUTOFF
                            continue
                        end

                        for j_res in chain2.first : chain2.last
                            coor_caj = cg_bead_coor[:, j_res]
                            if is_protein_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                                native_dist = compute_distance(coor_cai, coor_caj)
                                num_cg_pro_contact_all += 1
                                num_cg_pro_contact_inter += 1
                                tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                push!(top_cg_pro_go_contact, tmp_top_cnt)
    
                                # count AICG2+ atomic contact
                                contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ].atoms,
                                                                           cg_residues[ j_res ].atoms,
                                                                           cg_resid_name[i_res],
                                                                           cg_resid_name[j_res],
                                                                           aa_atom_name,
                                                                           aa_coor)
    
                                # calculate AICG2+ pairwise energy
                                e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
                                if e_local > AICG_ENE_UPPER_LIM
                                    e_local = AICG_ENE_UPPER_LIM
                                end
                                if e_local < AICG_ENE_LOWER_LIM
                                    e_local = AICG_ENE_LOWER_LIM
                                end
                                e_ground_contact += e_local
                                num_contact      += 1
                                push!(param_cg_pro_e_contact, e_local)
                            end
                        end
                    end
                end
            end
            print("\n              ... inter-molecular contacts: DONE! \n")
        end

        # normalize
        e_ground_contact /= num_contact

        if scale_scheme == 0
            for i in 1:length(param_cg_pro_e_contact)
                param_cg_pro_e_contact[i] *= AICG_CONTACT_AVE / e_ground_contact
            end
        elseif scale_scheme == 1
            for i in 1:length(param_cg_pro_e_contact)
                param_cg_pro_e_contact[i] *= -AICG_CONTACT_GEN
            end
        end

        println("------------------------------------------------------------")
        @printf("          > Total number of protein contacts: %12d  \n",
                length( top_cg_pro_go_contact ))

    end

    # =============================
    # Step 5: 3SPN.2C model for DNA
    # =============================
    #         _
    #      __| |_ __   __ _
    #     / _` | '_ \ / _` |
    #    | (_| | | | | (_| |
    #     \__,_|_| |_|\__,_|
    #
    # =============================

    if num_chain_DNA > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): processing DNA.")

        # ----------------------------------
        #        Step 5.1: determine P, S, B
        # ----------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine P, S, B mass, charge, and coordinates.")

        for i_chain in 1 : aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_DNA
                continue
            end

            for i_res in chain.first : chain.last
                res_name  = cg_residues[i_res].res_name
                bead_name = cg_residues[i_res].atm_name
                bead_type = bead_name == "DP" || bead_name == "DS" ? bead_name : res_name
                bead_coor = compute_center_of_mass(cg_residues[i_res].atoms, aa_atom_name, aa_coor)
                cg_resid_name[i_res]   = res_name
                cg_resid_index[i_res]  = cg_residues[i_res].res_idx
                cg_bead_name[i_res]    = bead_name
                cg_bead_type[i_res]    = bead_type
                cg_bead_charge[i_res]  = RES_CHARGE_DICT[bead_type]
                cg_bead_mass[i_res]    = RES_MASS_DICT[bead_type]
                cg_bead_coor[:, i_res] = bead_coor[:]
                cg_chain_id[i_res]     = i_chain
                cg_seg_name[i_res]     = chain.segname
            end
        end

        println(">           ... DONE!")

        # ----------------------------------------
        # Step 5.1.1: determine geometric features
        # ----------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1.1: determine Centroid, Rg, Rc...")

        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_DNA
                continue
            end

            # centroid
            coor_centroid = zeros(Float64, 3)
            for i_res in chain.first : chain.last
                coor_centroid += cg_bead_coor[:, i_res]
            end
            coor_centroid /= (chain.last - chain.first + 1)
            geo_centroid[:, i_chain] = coor_centroid

            # Rg
            tmp_dist = 0
            tmp_dist_sq_sum = 0
            for i_res in chain.first : chain.last
                vec_from_center   = cg_bead_coor[:, i_res] - coor_centroid
                vec_norm_tmp      = norm(vec_from_center)
                tmp_dist          = vec_norm_tmp > tmp_dist ? vec_norm_tmp : tmp_dist
                tmp_dist_sq_sum  += vec_norm_tmp * vec_norm_tmp
            end
            rg = sqrt(tmp_dist_sq_sum / (chain.last - chain.first + 1))
            rc = tmp_dist
            geo_radius_of_gyration[i_chain]     = rg
            geo_radius_of_circumsphere[i_chain] = rc
        end

        # ---------------------------------
        #        Step 5.2: 3SPN.2C topology
        # ---------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).2: 3SPN.2C topology.")
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.1: 3SPN.2C local interactions.")
        for i_chain in 1:aa_num_chain

            chain = cg_chains[i_chain]

            if chain.moltype != MOL_DNA
                continue
            end

            @printf("%11s Calculating DNA strand %d ... \n", " ", i_chain)
            @printf("              ... progress: %32s", " ")

            for i_res in chain.first : chain.last

                # -----------------
                # show progress bar
                # -----------------
                print("\b"^32)
                progress_percent = trunc(Int, ( i_res - chain.first ) / ( chain.last - chain.first ) * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, progress_percent * 5)
                # ------------------

                if cg_bead_name[i_res] == "DS"
                    # bond S--B
                    coor_s = cg_bead_coor[:, i_res]
                    coor_b = cg_bead_coor[:, i_res + 1]
                    r_sb   = compute_distance(coor_s, coor_b)
                    tmp_top_bond = CGTopBond(i_res, i_res + 1, r_sb)
                    push!(top_cg_DNA_bonds, tmp_top_bond)
                    if i_res + 3 < chain.last
                        # bond S--P+1
                        coor_p3 = cg_bead_coor[:, i_res + 2]
                        r_sp3   = compute_distance(coor_s, coor_p3)
                        tmp_top_bond = CGTopBond(i_res, i_res + 2, r_sp3)
                        push!(top_cg_DNA_bonds, tmp_top_bond)
                        # Angle S--P+1--S+1
                        resname5  = cg_resid_name[i_res][end]
                        resname3  = cg_resid_name[i_res + 3][end]
                        coor_s3   = cg_bead_coor[:, i_res + 3]
                        ang_sp3s3 = compute_angle(coor_s, coor_p3, coor_s3)
                        k         = get_DNA3SPN_angle_param("SPS", resname5 * resname3)
                        tmp_top_angl = CGTopAngle(i_res, i_res + 2, i_res + 3, ang_sp3s3)
                        push!(top_cg_DNA_angles, tmp_top_angl)
                        push!(param_cg_DNA_k_angles, k)
                        # Dihedral S--P+1--S+1--B+1
                        coor_b3     = cg_bead_coor[:, i_res + 4]
                        dih_sp3s3b3 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_b3)
                        tmp_top_dih = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 4, dih_sp3s3b3)
                        push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                        # Dihedral S--P+1--S+1--P+2
                        if i_res + 6 < chain.last
                            coor_p33     = cg_bead_coor[:, i_res + 5]
                            dih_sp3s3p33 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_p33)
                            tmp_top_dih = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                            tmp_top_dih = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33)
                            push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                        end
                    end
                elseif cg_bead_name[i_res] == "DP"
                    # bond P--S
                    coor_p = cg_bead_coor[:, i_res]
                    coor_s = cg_bead_coor[:, i_res + 1]
                    r_ps   = compute_distance(coor_p, coor_s)
                    tmp_top_bond = CGTopBond(i_res, i_res + 1, r_ps)
                    push!(top_cg_DNA_bonds, tmp_top_bond)
                    # angle P--S--B
                    resname5 = cg_resid_name[i_res - 1][end]
                    resname3 = cg_resid_name[i_res + 2][end]
                    coor_b   = cg_bead_coor[:, i_res + 2]
                    ang_psb  = compute_angle(coor_p, coor_s, coor_b)
                    k        = get_DNA3SPN_angle_param("PSB", resname5 * resname3)
                    tmp_top_angl = CGTopAngle(i_res, i_res + 1, i_res + 2, ang_psb)
                    push!(top_cg_DNA_angles, tmp_top_angl)
                    push!(param_cg_DNA_k_angles, k)
                    if i_res + 4 < chain.last
                        # angle P--S--P+1
                        coor_p3  = cg_bead_coor[:, i_res + 3]
                        ang_psp3 = compute_angle(coor_p, coor_s, coor_p3)
                        k        = get_DNA3SPN_angle_param("PSP", "all")
                        tmp_top_angl = CGTopAngle(i_res, i_res + 1, i_res + 3, ang_psp3)
                        push!(top_cg_DNA_angles, tmp_top_angl)
                        push!(param_cg_DNA_k_angles, k)
                        # Dihedral P--S--P+1--S+1
                        coor_s3    = cg_bead_coor[:, i_res + 4]
                        dih_psp3s3 = compute_dihedral(coor_p, coor_s, coor_p3, coor_s3)
                        tmp_top_dih = CGTopDihedral(i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3)
                        push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                        tmp_top_dih = CGTopDihedral(i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3)
                        push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                    end
                elseif cg_bead_name[i_res] == "DB"
                    if i_res + 2 < chain.last
                        # angle B--S--P+1
                        resname5 = cg_resid_name[i_res][end]
                        resname3 = cg_resid_name[i_res + 1][end]
                        coor_b   = cg_bead_coor[:, i_res]
                        coor_s   = cg_bead_coor[:, i_res - 1]
                        coor_p3  = cg_bead_coor[:, i_res + 1]
                        ang_bsp3 = compute_angle(coor_b, coor_s, coor_p3)
                        k        = get_DNA3SPN_angle_param("BSP", resname5 * resname3)
                        tmp_top_angl = CGTopAngle(i_res, i_res - 1, i_res + 1, ang_bsp3)
                        push!(top_cg_DNA_angles, tmp_top_angl)
                        push!(param_cg_DNA_k_angles, k)
                        # Dihedral B--S--P+1--S+1
                        coor_s3    = cg_bead_coor[:, i_res + 2]
                        dih_bsp3s3 = compute_dihedral(coor_b, coor_s, coor_p3, coor_s3)
                        tmp_top_dih = CGTopDihedral(i_res, i_res - 1, i_res + 1, i_res + 2, dih_bsp3s3)
                        push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                    end
                else
                    errmsg = @sprintf("BUG: Wrong DNA particle type in chain %d, residue %d : %s ",
                                      i_chain,
                                      i_res, cg_bead_name[i_res])
                    error(errmsg)
                end
            end
            print(" \n")
        end
        println(">           ... Bond, Angle, Dihedral: DONE!")
    end

    # =========================
    # RNA structure based model
    # =========================
    #     ____  _   _    _    
    #    |  _ \| \ | |  / \   
    #    | |_) |  \| | / _ \  
    #    |  _ <| |\  |/ ___ \ 
    #    |_| \_\_| \_/_/   \_\
    # 
    # =========================

    if num_chain_RNA > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): processing RNA.")

        # ----------------------------------
        #         determine P, S, B
        # ----------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine P, S, B mass, charge, and coordinates.")

        for i_chain in 1 : aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            for i_res in chain.first : chain.last
                res_name  = cg_residues[i_res].res_name
                bead_name = cg_residues[i_res].atm_name
                bead_type = bead_name == "RP" || bead_name == "RS" ? bead_name : res_name
                # bead_coor = compute_center_of_mass(cg_residues[i_res].atoms, aa_atom_name, aa_coor)
                cg_resid_name[i_res]   = res_name
                cg_resid_index[i_res]  = cg_residues[i_res].res_idx
                cg_bead_name[i_res]    = bead_name
                cg_bead_type[i_res]    = bead_type
                cg_bead_charge[i_res]  = RES_CHARGE_DICT[bead_type]
                cg_bead_mass[i_res]    = RES_MASS_DICT[bead_type]
                cg_chain_id[i_res]     = i_chain
                cg_seg_name[i_res]     = chain.segname
                if bead_name == "RP"
                    for i_atom in cg_residues[i_res].atoms
                        if aa_atom_name[i_atom][1] == 'P'
                            bead_coor = aa_coor[:, i_atom]
                        end
                    end
                elseif bead_name == "RS"
                    total_mass      = 0
                    tmp_coor        = zeros(Float64, 3)
                    for i_atom in cg_residues[i_res].atoms
                        a_name = aa_atom_name[i_atom]
                        if in(a_name, ["C1'", "C2'", "C3'", "C4'", "O4'"] )
                            a_mass      = ATOM_MASS_DICT[a_name[1]]
                            a_coor      = aa_coor[:, i_atom]
                            total_mass += a_mass
                            tmp_coor   += a_coor * a_mass
                        end
                    end
                    bead_coor = tmp_coor / total_mass
                elseif bead_name == "RB"
                    if res_name[end] == 'A' || res_name[end] == 'G'
                        for i_atom in cg_residues[i_res].atoms
                            if aa_atom_name[i_atom] == "N1"
                                bead_coor = aa_coor[:, i_atom]
                            end
                        end
                    else
                        for i_atom in cg_residues[i_res].atoms
                            if aa_atom_name[i_atom] == "N3"
                                bead_coor = aa_coor[:, i_atom]
                            end
                        end
                    end
                end
                cg_bead_coor[:, i_res] = bead_coor[:]
            end
        end

        println(">           ... DONE!")

        # ----------------------------------------
        # Step 6.1.1: determine geometric features
        # ----------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1.1: determine Centroid, Rg, Rc...")

        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            # centroid
            coor_centroid = zeros(Float64, 3)
            for i_res in chain.first : chain.last
                coor_centroid += cg_bead_coor[:, i_res]
            end
            coor_centroid /= (chain.last - chain.first + 1)
            geo_centroid[:, i_chain] = coor_centroid

            # Rg
            tmp_dist = 0
            tmp_dist_sq_sum = 0
            for i_res in chain.first : chain.last
                vec_from_center   = cg_bead_coor[:, i_res] - coor_centroid
                vec_norm_tmp      = norm(vec_from_center)
                tmp_dist          = vec_norm_tmp > tmp_dist ? vec_norm_tmp : tmp_dist
                tmp_dist_sq_sum  += vec_norm_tmp * vec_norm_tmp
            end
            rg = sqrt(tmp_dist_sq_sum / (chain.last - chain.first + 1))
            rc = tmp_dist
            geo_radius_of_gyration[i_chain]     = rg
            geo_radius_of_circumsphere[i_chain] = rc
        end

        # -------------------------
        # Step 6.2: RNA topology
        # -------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).2: RNA topology.")
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.1: RNA local interactions.")

        @printf("%11s Calculating bonded terms... \n", " ")
        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            for i_res in chain.first : chain.last
                if cg_bead_name[i_res] == "RS"
                    # bond S--B
                    coor_s    = cg_bead_coor[:, i_res]
                    coor_b    = cg_bead_coor[:, i_res + 1]
                    r_sb      = compute_distance(coor_s, coor_b)
                    base_type = cg_resid_name[i_res] in ["RA", "RG"] ? "R" : "Y"
                    bond_type = "S" * base_type
                    k         = RNA_BOND_K_LIST[bond_type]
                    tmp_top_bond = CGTopBond(i_res, i_res + 1, r_sb)
                    push!(top_cg_RNA_bonds, tmp_top_bond)
                    push!(param_cg_RNA_k_bonds, k)
                    # bond S--P+1
                    if i_res + 2 < chain.last
                        coor_p3 = cg_bead_coor[:, i_res + 2]
                        r_sp3   = compute_distance(coor_s, coor_p3)
                        k       = RNA_BOND_K_LIST["SP"]
                        tmp_top_bond = CGTopBond(i_res, i_res + 2, r_sp3)
                        push!(top_cg_RNA_bonds, tmp_top_bond)
                        push!(param_cg_RNA_k_bonds, k)
                    end
                    if i_res + 4 <= chain.last
                        # Angle S--P+1--S+1
                        coor_s3   = cg_bead_coor[:, i_res + 3]
                        ang_sp3s3 = compute_angle(coor_s, coor_p3, coor_s3)
                        k         = RNA_ANGLE_K_LIST["SPS"]
                        tmp_top_angl = CGTopAngle(i_res, i_res + 2, i_res + 3, ang_sp3s3)
                        push!(top_cg_RNA_angles, tmp_top_angl)
                        push!(param_cg_RNA_k_angles, k)
                        # Dihedral S--P+1--S+1--B+1
                        coor_b3     = cg_bead_coor[:, i_res + 4]
                        dih_sp3s3b3 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_b3)
                        base_type   = cg_resid_name[i_res + 4] in ["RA", "RG"] ? "R" : "Y"
                        dihe_type   = "SPS" * base_type
                        k           = RNA_DIHEDRAL_K_LIST[dihe_type]
                        tmp_top_dih = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 4, dih_sp3s3b3)
                        push!(top_cg_RNA_dihedrals, tmp_top_dih)
                        push!(param_cg_RNA_k_dihedrals, k)
                    end
                    # Dihedral S--P+1--S+1--P+2
                    if i_res + 5 < chain.last
                        coor_p33     = cg_bead_coor[:, i_res + 5]
                        dih_sp3s3p33 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_p33)
                        k            = RNA_DIHEDRAL_K_LIST["SPSP"]
                        tmp_top_dih = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33)
                        push!(top_cg_RNA_dihedrals, tmp_top_dih)
                        push!(param_cg_RNA_k_dihedrals, k)
                    end
                elseif cg_bead_name[i_res] == "RP"
                    # bond P--S
                    coor_p = cg_bead_coor[:, i_res]
                    coor_s = cg_bead_coor[:, i_res + 1]
                    r_ps   = compute_distance(coor_p, coor_s)
                    k      = RNA_BOND_K_LIST["PS"]
                    tmp_top_bond = CGTopBond(i_res, i_res + 1, r_ps)
                    push!(top_cg_RNA_bonds, tmp_top_bond)
                    push!(param_cg_RNA_k_bonds, k)
                    # angle P--S--B
                    coor_b    = cg_bead_coor[:, i_res + 2]
                    ang_psb   = compute_angle(coor_p, coor_s, coor_b)
                    base_type = cg_resid_name[i_res + 2] in ["RA", "RG"] ? "R" : "Y"
                    angl_type = "PS" * base_type
                    k         = RNA_ANGLE_K_LIST[angl_type]
                    tmp_top_angl = CGTopAngle(i_res, i_res + 1, i_res + 2, ang_psb)
                    push!(top_cg_RNA_angles, tmp_top_angl)
                    push!(param_cg_RNA_k_angles, k)
                    if i_res + 4 < chain.last
                        # angle P--S--P+1
                        coor_p3  = cg_bead_coor[:, i_res + 3]
                        ang_psp3 = compute_angle(coor_p, coor_s, coor_p3)
                        k        = RNA_ANGLE_K_LIST["PSP"]
                        tmp_top_angl = CGTopAngle(i_res, i_res + 1, i_res + 3, ang_psp3)
                        push!(top_cg_RNA_angles, tmp_top_angl)
                        push!(param_cg_RNA_k_angles, k)
                        # Dihedral P--S--P+1--S+1
                        coor_s3    = cg_bead_coor[:, i_res + 4]
                        dih_psp3s3 = compute_dihedral(coor_p, coor_s, coor_p3, coor_s3)
                        k          = RNA_DIHEDRAL_K_LIST["PSPS"]
                        tmp_top_dih = CGTopDihedral(i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3)
                        push!(top_cg_RNA_dihedrals, tmp_top_dih)
                        push!(param_cg_RNA_k_dihedrals, k)
                    end
                elseif cg_bead_name[i_res] == "RB"
                    # do nothing...
                end
            end
        end

        # -----------------------
        # HT type native contacts
        # -----------------------
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.2: RNA HT-type native contacts.")
        @printf("%11s Calculating intra-molecular contacts... \n", " ")
        @printf("              ... progress: %32s", " ")
        i_progress_count = 0
        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            # -----------------
            # show progress bar
            # -----------------
            i_progress_count += 1
            print("\b"^32)
            progress_percent = trunc(Int, i_progress_count / ( num_chain_RNA ) * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %2d / %2d ", progress_bar, i_progress_count, num_chain_RNA)
            # ------------------

            for i_res in chain.first : chain.last - 3

                if cg_bead_name[i_res] == "RP"
                    continue
                end

                coor_i = cg_bead_coor[:, i_res]

                for j_res in i_res + 3 : chain.last

                    if cg_bead_name[j_res] == "RP"
                        continue
                    end

                    if cg_bead_name[i_res] == "RS" || cg_bead_name[j_res] == "RS"
                        if j_res < i_res + 6
                            continue
                        end
                    end

                    coor_j = cg_bead_coor[:, j_res]
                    native_dist = compute_distance(coor_i, coor_j)
                    adist, nhb  = compute_RNA_native_contact(cg_residues[i_res].atoms,
                                                             cg_residues[j_res].atoms,
                                                             aa_atom_name,
                                                             aa_coor)

                    if adist > RNA_GO_ATOMIC_CUTOFF
                        continue
                    end
                    
                    if j_res == i_res + 3 && cg_bead_name[i_res] == "RB"
                        coor_i_sug = cg_bead_coor[:, i_res - 1]
                        coor_j_sug = cg_bead_coor[:, j_res - 1]
                        st_dih = compute_dihedral(coor_i, coor_i_sug, coor_j_sug, coor_j)
                        if abs( st_dih ) < RNA_STACK_DIH_CUTOFF && adist < RNA_STACK_DIST_CUTOFF
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_RNA_base_stack, tmp_top_cnt)
                            push!(param_cg_RNA_e_base_stack, RNA_STACK_EPSILON)
                        else
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_RNA_other_contact, tmp_top_cnt)
                            push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER["BB"])
                        end
                    elseif cg_bead_name[i_res] == "RB" && cg_bead_name[j_res] == "RB"
                        if nhb == 2
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_RNA_base_pair, tmp_top_cnt)
                            push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_2HB)
                        elseif nhb >= 3
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_RNA_base_pair, tmp_top_cnt)
                            push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_3HB)
                        else
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_RNA_other_contact, tmp_top_cnt)
                            push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER["BB"])
                        end
                    else
                        contact_type = cg_bead_name[i_res][end] * cg_bead_name[j_res][end]
                        tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                        push!(top_cg_RNA_other_contact, tmp_top_cnt)
                        push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER[contact_type])
                    end
                end
            end
        end
        print("\n              ... intra-molecular contacts: DONE! \n")
 
        if num_chain_RNA > 1
            @printf("%11s Calculating inter-molecular contacts... \n", " ")
            @printf("              ... progress: %32s", " ")
            for i_chain in 1:aa_num_chain - 1
    
                chain_1 = cg_chains[i_chain]
    
                # -----------------
                # show progress bar
                # -----------------
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / ( aa_num_chain - 1 ) * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / ( aa_num_chain - 1 ) * 100)
                # ------------------
    
                if chain_1.moltype != MOL_RNA
                    continue
                end
    
                centroid_chain1 = geo_centroid[:, i_chain]
                radius_of_circ1 = geo_radius_of_circumsphere[i_chain]

                for j_chain in i_chain + 1 : aa_num_chain
                    chain_2 = cg_chains[j_chain]
                    if chain_2.moltype != MOL_RNA
                        continue
                    end

                    centroid_chain2 = geo_centroid[:, j_chain]
                    radius_of_circ2 = geo_radius_of_circumsphere[j_chain]

                    cent_cent_dist = compute_distance(centroid_chain1, centroid_chain2)
                    if cent_cent_dist > radius_of_circ1 + radius_of_circ2 + CG_CONTACT_CUTOFF
                        continue
                    end

                    for i_res in chain_1.first : chain_1.last
                        if cg_bead_name[i_res] == "RP"
                            continue
                        end
                        coor_i = cg_bead_coor[:, i_res]

                        coori_centj_dist = compute_distance(coor_i, centroid_chain2)
                        if coori_centj_dist > radius_of_circ2 + CG_CONTACT_CUTOFF
                            continue
                        end

                        for j_res in chain_2.first : chain_2.last
                            if cg_bead_name[j_res] == "RP"
                                continue
                            end
                            coor_j = cg_bead_coor[:, j_res]
                            native_dist = compute_distance(coor_i, coor_j)
                            adist, nhb  = compute_RNA_native_contact(cg_residues[i_res].atoms,
                                                                     cg_residues[j_res].atoms,
                                                                     aa_atom_name,
                                                                     aa_coor)
                            if adist > RNA_GO_ATOMIC_CUTOFF
                                continue
                            end
                            if cg_bead_name[i_res] == "RB" && cg_bead_name[j_res] == "RB"
                                if nhb == 2
                                    tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                    push!(top_cg_RNA_base_pair, tmp_top_cnt)
                                    push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_2HB)
                                elseif nhb >= 3
                                    tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                    push!(top_cg_RNA_base_pair, tmp_top_cnt)
                                    push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_3HB)
                                else
                                    tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                    push!(top_cg_RNA_other_contact, tmp_top_cnt)
                                    push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER["BB"])
                                end
                            else
                                contact_type = cg_bead_name[i_res][end] * cg_bead_name[j_res][end]
                                tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                push!(top_cg_RNA_other_contact, tmp_top_cnt)
                                push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER[contact_type])
                            end
                        end
                    end
                end
            end
            print("\n              ... inter-molecular contacts: DONE! \n")
        end
 
        println("------------------------------------------------------------")
        @printf("          > Total number of RNA contacts:     %12d  \n",
                length(top_cg_RNA_base_stack) + length(top_cg_RNA_base_pair) + length(top_cg_RNA_other_contact) )

    end

    # ===========================================================
    # Protein-RNA structure-based interactions: Go-like potential
    # ===========================================================
    #                  _       _             ____  _   _    _    
    #  _ __  _ __ ___ | |_ ___(_)_ __       |  _ \| \ | |  / \   
    # | '_ \| '__/ _ \| __/ _ \ | '_ \ _____| |_) |  \| | / _ \  
    # | |_) | | | (_) | ||  __/ | | | |_____|  _ <| |\  |/ ___ \ 
    # | .__/|_|  \___/ \__\___|_|_| |_|     |_| \_\_| \_/_/   \_\
    # |_|                                                        
    # 
    # ============================================================

    if num_chain_RNA > 0 && num_chain_pro > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): Generating protein-RNA native contacts.")

        @printf("%11s Calculating protein-RNA contacts... \n", " ")
        @printf("              ... progress: %32s", " ")
        for i_chain in 1:aa_num_chain
            # -----------------
            # show progress bar
            # -----------------
            print("\b"^32)
            progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            # ------------------
    
            chain_pro = cg_chains[i_chain]

            if chain_pro.moltype != MOL_PROTEIN
                continue
            end

            centroid_chain_pro = geo_centroid[:, i_chain]
            radius_of_circ_pro = geo_radius_of_circumsphere[i_chain]

            for j_chain in 1 : aa_num_chain
                chain_RNA = cg_chains[j_chain]
                if chain_RNA.moltype != MOL_RNA
                    continue
                end

                centroid_chain_rna = geo_centroid[:, j_chain]
                radius_of_circ_rna = geo_radius_of_circumsphere[j_chain]

                cent_cent_dist = compute_distance(centroid_chain_pro, centroid_chain_rna)
                if cent_cent_dist > radius_of_circ_pro + radius_of_circ_rna + CG_CONTACT_CUTOFF
                    continue
                end

                for i_res in chain_pro.first : chain_pro.last
                    coor_i = cg_bead_coor[:, i_res]

                    cai_centj_dist = compute_distance(coor_i, centroid_chain_rna)
                    if cai_centj_dist > radius_of_circ_rna + CG_CONTACT_CUTOFF
                        continue
                    end

                    for j_res in chain_RNA.first : chain_RNA.last
                        if cg_bead_name[j_res] == "RP"
                            continue
                        end
                        if !is_protein_RNA_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                            continue
                        end
                        coor_j = cg_bead_coor[:, j_res]
                        native_dist = compute_distance(coor_i, coor_j)
                        if cg_bead_name[j_res] == "RS"
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_pro_RNA_contact, tmp_top_cnt)
                            push!(param_cg_pro_RNA_e_contact, PRO_RNA_GO_EPSILON_S)
                        elseif cg_bead_name[j_res] == "RB"
                            tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            push!(top_cg_pro_RNA_contact, tmp_top_cnt)
                            push!(param_cg_pro_RNA_e_contact, PRO_RNA_GO_EPSILON_B)
                        end
                    end
                end
            end
        end

        # println("\n>           ... DONE!")
        println("\n------------------------------------------------------------")
        @printf("          > Total number of protein-RNA contacts: %8d  \n",
                length(top_cg_pro_RNA_contact) )
    end


    # ============================================================
    # PWMcos parameters: protein-DNA sequence-specific interaction
    # ============================================================
    #        ______        ____  __               
    #       |  _ \ \      / /  \/  | ___ ___  ___ 
    #       | |_) \ \ /\ / /| |\/| |/ __/ _ \/ __|
    #       |  __/ \ V  V / | |  | | (_| (_) \__ \
    #       |_|     \_/\_/  |_|  |_|\___\___/|___/
    # 
    # ============================================================
    if ff_pro_dna == FF_PWMcos
        pwmcos_native_contacts = []

        if num_chain_pro == 0
            error("Cannot generate PWMcos parameters without protein...")
        end
        if num_chain_DNA != 2
            error("Cannot generate PWMcos parameters from more or less than two DNA chains...")
        end

        i_step += 1
        println("============================================================")
        println("> Step $(i_step): Generating PWMcos parameters.")

        # ------------------------------------------------
        #        Step 7.1: determine protein-DNA contacts
        # ------------------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine contacts between protein and DNA.")

        i_count_DNA = 0
        for i_chain in 1:aa_num_chain
            chain_pro = cg_chains[i_chain]

            if chain_pro.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain_pro.first : chain_pro.last
                i_res_N = i_res == chain_pro.first ? i_res : i_res - 1
                i_res_C = i_res == chain_pro.last  ? i_res : i_res + 1

                coor_pro_i = cg_bead_coor[:, i_res]
                coor_pro_N = cg_bead_coor[:, i_res_N]
                coor_pro_C = cg_bead_coor[:, i_res_C]

                for j_chain in 1:aa_num_chain
                    chain_DNA = cg_chains[j_chain]
                    
                    if chain_DNA.moltype != MOL_DNA
                        continue
                    end

                    for j_res in chain_DNA.first + 3 : chain_DNA.last - 3
                        if cg_bead_name[j_res] != "DB"
                            continue
                        end
                        if !is_PWMcos_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                            continue
                        end

                        j_res_5, j_res_3 = j_res - 3, j_res + 3
                        coor_dna_j       = cg_bead_coor[:, j_res]
                        coor_dna_5       = cg_bead_coor[:, j_res_5]
                        coor_dna_3       = cg_bead_coor[:, j_res_3]
                        coor_dna_S       = cg_bead_coor[:, j_res - 1]

                        vec0   = coor_pro_i - coor_dna_j
                        vec1   = coor_dna_S - coor_dna_j
                        vec2   = coor_dna_3 - coor_dna_5
                        vec3   = coor_pro_N - coor_pro_C
                        r0     = norm(vec0)
                        theta1 = compute_vec_angle(vec0, vec1)
                        theta2 = compute_vec_angle(vec0, vec2)
                        theta3 = compute_vec_angle(vec0, vec3)

                        push!(pwmcos_native_contacts, (i_res - chain_pro.first + 1,
                                                       cg_resid_index[j_res],
                                                       r0,
                                                       theta1,
                                                       theta2,
                                                       theta3))

                        if do_debug
                            println("PWMcos | pro ===> ", i_res - chain_pro.first + 1,
                                    " DNA ===> ", j_res, " : ", cg_resid_index[j_res],
                                    " r0 = ", r0,
                                    " theta1 = ", theta1,
                                    " theta2 = ", theta2,
                                    " theta3 = ", theta3)
                        end
                    end
                end
            end
        end

        # ------------------------------------------------
        #        Step 7.2: Read in PFM and convert to PWM
        # ------------------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).2: read in position frequency matrix (PFM).")

        pwmcos_pwm, pwmcos_chain_a, pwmcos_chain_b = read_modified_pfm(pfm_filename)
        num_pwmcos_terms = length(pwmcos_chain_a)


        # ------------------------------------------------
        #        Step 7.2: Read in PFM and convert to PWM
        # ------------------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).3: decomposing PWM.")

        ip_count = zeros(Float64, num_pwmcos_terms)

        contact_to_pwm = []
        for nat_contact in pwmcos_native_contacts
            i_dna = nat_contact[2] # cg_resid_index[dna]
            cnt_pwm_idx_a = indexin(i_dna, pwmcos_chain_a)[1]
            cnt_pwm_idx_b = indexin(i_dna, pwmcos_chain_b)[1]
            if cnt_pwm_idx_a isa Int
                push!( contact_to_pwm, (cnt_pwm_idx_a, 1) )
                ip_count[cnt_pwm_idx_a] += 1
            elseif cnt_pwm_idx_b isa Int
                push!( contact_to_pwm, (cnt_pwm_idx_b, -1) )
                ip_count[cnt_pwm_idx_b] += 1
            else
                error("Index error in CHAIN_A or CHAIN_B!")
            end
        end
        pwm_decomposed = pwmcos_pwm ./ ip_count

        for (i_cnt, nat_cnt) in enumerate(pwmcos_native_contacts)
            pwm_i, pwm_order = contact_to_pwm[i_cnt][1], contact_to_pwm[i_cnt][2]
            if pwm_order == 1
                eA, eC, eG, eT = pwm_decomposed[pwm_i, 1:4]
            elseif pwm_order == -1
                eA, eC, eG, eT = pwm_decomposed[pwm_i, 4:-1:1]
            end
            tmp_top_pwmcos = CGTopPWMcos(nat_cnt[1],
                                       nat_cnt[3],
                                       nat_cnt[4],
                                       nat_cnt[5],
                                       nat_cnt[6],
                                       eA, eC, eG, eT)
            push!(top_cg_pro_DNA_pwmcos, tmp_top_pwmcos)
        end

        if do_debug
            println(size( contact_to_pwm ))
            println(pwm_decomposed)
        end

        println(">           ... DONE!")
    end





    # =============================
    # Make a new topology structure
    # =============================

    mytop = CGTopology(cg_resid_name,
                       cg_resid_index,
                       cg_bead_name,
                       cg_bead_type,
                       cg_bead_charge,
                       cg_bead_mass,
                       cg_chain_id,
                       cg_seg_name,
                       top_cg_pro_bonds,
                       top_cg_pro_angles,
                       top_cg_pro_dihedrals,
                       top_cg_pro_aicg13,
                       top_cg_pro_aicg14,
                       top_cg_pro_go_contact,
                       param_cg_pro_e_13,
                       param_cg_pro_e_14,
                       param_cg_pro_e_contact,
                       top_cg_DNA_bonds,
                       top_cg_DNA_angles,
                       top_cg_DNA_dih_Gaussian,
                       top_cg_DNA_dih_periodic,
                       param_cg_DNA_k_angles,
                       top_cg_RNA_bonds,
                       top_cg_RNA_angles,
                       top_cg_RNA_dihedrals,
                       top_cg_RNA_base_stack,
                       top_cg_RNA_base_pair,
                       top_cg_RNA_other_contact,
                       param_cg_RNA_k_bonds,
                       param_cg_RNA_k_angles,
                       param_cg_RNA_k_dihedrals,
                       param_cg_RNA_e_base_stack,
                       param_cg_RNA_e_base_pair,
                       param_cg_RNA_e_other_contact,
                       top_cg_pro_DNA_pwmcos,
                       top_cg_pro_RNA_contact,
                       param_cg_pro_RNA_e_contact)

    myconf = Conformation(cg_num_particles, cg_bead_coor)




    # ----------
    # output log
    # ----------
    if do_output_log
        log_name = pdb_name[1:end-4] * "_cg.log"
        log_file = open(log_name, "w")

        println(log_file, "================================================================================")
        println(log_file, " PDB info (atomic):")
        println(log_file, " - Number of atoms    : $(aa_num_atom)")
        println(log_file, " - Number of residues : $(aa_num_residue)")
        println(log_file, " - Number of chains   : $(aa_num_chain)")

        println(log_file, "================================================================================")
        println(log_file, " Chain info (CG):")
        @printf(log_file, " - Number of protein chains: %5d \n", num_chain_pro)
        @printf(log_file, " - Number of DNA strands:    %5d \n", num_chain_DNA)
        @printf(log_file, " - Number of RNA strands:    %5d \n", num_chain_RNA)

        println(log_file, " |--------------------------------------------------------------------|")
        println(log_file, " | Chain | Mol Type | # bead | start --   end |   Rg (Å) | net charge | ")
        println(log_file, " |-------+----------+--------+----------------+----------+------------|")
        for i_chain = 1:aa_num_chain
            chain = cg_chains[i_chain]
            charge = sum( cg_bead_charge[chain.first : chain.last] )
            @printf(log_file, " |   %3d | %8s | %6d | %5d -- %5d | %8.3f | %+10.3f | \n",
                    i_chain, MOL_TYPE_LIST[ chain.moltype ], cg_chain_length[i_chain],
                    cg_chains[i_chain].first, cg_chains[i_chain].last,
                    geo_radius_of_gyration[i_chain],
                    charge)
        end

        println(log_file, " |--------------------------------------------------------------------|")
        println(log_file, " CG mol info:")
        charge = sum( cg_bead_charge )
        rg_all = radius_of_gyration(myconf)
        rc_all = radius_of_circumshpere(myconf)
        @printf(log_file, " - Number of CG particles: %8d \n", cg_num_particles)
        @printf(log_file, " - Radius of gyration:     %8.3f Å \n", rg_all)
        @printf(log_file, " - Radius of circumsphere: %8.3f Å \n", rc_all)
        @printf(log_file, " - Net Charge:             %+8.3f e \n", charge)


        println(log_file, "================================================================================")
        println(log_file, " Interaction info:")
        if num_chain_pro > 0
            @printf(log_file, " - Number of protein contacts:     %12d  \n", length(top_cg_pro_go_contact))
        end
        if num_chain_RNA > 0
            @printf(log_file, " - Number of RNA contacts:         %12d  \n", length(top_cg_RNA_base_stack) + length(top_cg_RNA_base_pair) + length(top_cg_RNA_other_contact) )
        end
        if num_chain_RNA > 0 && num_chain_pro > 0
            @printf(log_file, " - Number of protein-RNA contacts: %12d  \n", length(top_cg_pro_RNA_contact) )
        end
        println(log_file, "================================================================================")

        close(log_file)
    end



    return ( mytop, myconf )

end
