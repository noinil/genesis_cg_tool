###############################################################################
#                         Coarse Graining Biomolecules                        #
###############################################################################

using Printf

struct CGTopBond
    i::Int
    j::Int
    r0::Float64
end

struct CGTopAngle
    i::Int
    j::Int
    k::Int
    a0::Float64
end

struct CGTopDihedral
    i::Int
    j::Int
    k::Int
    l::Int
    t0::Float64
end

struct CGTopContact
    i::Int
    j::Int
    r0::Float64
end

struct CGTopPWMcos
    i::Int
    r0::Float64
    t1::Float64
    t2::Float64
    t3::Float64
    eA::Float64
    eC::Float64
    eG::Float64
    eT::Float64
end

function coarse_graining(aa_molecule::AAMolecule, force_field::ForceFieldCG, args)

    # =================
    # Parsing arguments
    # =================
    pdb_name                = get(args, "pdb", "random.pdb")
    protein_charge_filename = get(args, "respac", "")
    pfm_filename            = get(args, "pfm", "")
    verbose                 = get(args, "verbose", false)
    do_debug                = get(args, "debug", false)
    do_output_log           = get(args, "log", false)
    do_test_local_only      = get(args, "test-local-only", false)

    use_safe_dihedral       = get(args, "use-safe-dihedral", 1)

    gen_3spn_itp            = get(args, "3spn-param", 0)
    DNA_use_5_phos          = get(args, "3spn-use-5-phos", false)
    DNA_circular            = get(args, "3spn-circular", false)

    ccgo_contact_scale      = get(args, "CCGO-contact-scale", 1.0)
    aicg_scale_scheme       = get(args, "aicg-scale", 1)
    cgRNA_use_phosphate_go  = get(args, "cgRNA-phosphate-Go", false)
    pwmcos_gamma            = get(args, "pwmcos-scale", 1.0)
    pwmcos_epsil            = get(args, "pwmcos-shift", 0.0)
    pwmcosns_epsil          = get(args, "pwmcos-ns-ene", -1.0)

    # ----------------------
    # More details from TOML
    # ----------------------
    has_toml_mod   = false
    if haskey(args, "modeling-options")
        has_toml_mod    = true
        ff_detail_config = args["modeling-options"]

        if haskey(ff_detail_config, "3SPN.2C")
            if haskey(ff_detail_config["3SPN.2C"], "USE_5_PHOSPHATE")
                val_string = ff_detail_config["3SPN.2C"]["USE_5_PHOSPHATE"]
                if val_string == "YES"
                    DNA_use_5_phos = true
                end
            end
        end
    end

    # ===============
    # Step 0: numbers
    # ===============
    i_step     = 0

    ff_pro     = force_field.ff_protein
    ff_dna     = force_field.ff_DNA
    ff_rna     = force_field.ff_RNA
    ff_pro_dna = force_field.ff_protein_DNA
    ff_pro_rna = force_field.ff_protein_RNA
    ff_dna_rna = force_field.ff_DNA_RNA

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
    if verbose
        println("============================================================")
        println("> Step $(i_step): estimate CG particle number for every chain.")
    end

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
            if DNA_use_5_phos
                n_particles = 3 * n_res
            else
                n_particles = 3 * n_res - 1
            end
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
        if verbose
            @printf("          > Chain %3d | %7s \n", i_chain, MOL_TYPE_LIST[ mol_type ])
        end
    end

    if verbose
        println("------------------------------------------------------------")
        @printf("          In total: %5d protein chains,\n", num_chain_pro)
        @printf("                    %5d DNA strands,\n", num_chain_DNA)
        @printf("                    %5d RNA strands.\n", num_chain_RNA)
    end

    # ===========================
    # Step 2: Assign CG particles
    # ===========================
    i_step += 1
    if verbose
        println("============================================================")
        println("> Step $(i_step): assign coarse-grained particles.")
    end

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
                if i_local_index > 1
                    cg_DP_idx = [tmp_atom_index_O3p]
                else
                    cg_DP_idx = []
                end
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
                if i_local_index > 1 || DNA_use_5_phos
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

    if verbose
        for i_chain in 1:aa_num_chain
            @printf("          > Chain %3d | # particles: %5d | %5d -- %5d \n",
                    i_chain, cg_chain_length[i_chain],
                    cg_chains[i_chain].first, cg_chains[i_chain].last)
        end

        println("------------------------------------------------------------")
        println("          In total: $(cg_num_particles) CG particles.")
    end




    # =========================================================================
    #        ____ ____   _____ ___  ____   ___  _     ___   ______   __
    #       / ___/ ___| |_   _/ _ \|  _ \ / _ \| |   / _ \ / ___\ \ / /
    #      | |  | |  _    | || | | | |_) | | | | |  | | | | |  _ \ V /
    #      | |__| |_| |   | || |_| |  __/| |_| | |__| |_| | |_| | | |
    #       \____\____|   |_| \___/|_|    \___/|_____\___/ \____| |_|
    #
    # =========================================================================

    # ========================================
    # Coarse Grained Model Topology Structures
    # ========================================

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
    top_cg_pro_go_contact    = [[] for i in 1:cg_num_particles]

    param_cg_pro_e_13        = []
    param_cg_pro_e_14        = []

    AICG2p_flexible_local    = []
    AICG2p_flexible_nonlocal = []
    HPS_IDR_region           = []
    KH_IDR_region            = []

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
    # top_cg_RNA_base_stack        = Vector{CGTopContact}(undef, 0)
    # top_cg_RNA_base_pair         = Vector{CGTopContact}(undef, 0)
    # top_cg_RNA_other_contact     = Vector{CGTopContact}(undef, 0)
    top_cg_RNA_base_stack        = [[] for i in 1:cg_num_particles]
    top_cg_RNA_base_pair         = [[] for i in 1:cg_num_particles]
    top_cg_RNA_other_contact     = [[] for i in 1:cg_num_particles]
    param_cg_RNA_k_bonds         = []
    param_cg_RNA_k_angles        = []
    param_cg_RNA_k_dihedrals     = []
    # param_cg_RNA_e_base_stack    = []
    # param_cg_RNA_e_base_pair     = []
    # param_cg_RNA_e_other_contact = []

    # protein-DNA
    top_cg_pro_DNA_pwmcos         = Vector{CGTopPWMcos}(undef, 0)
    top_cg_pro_DNA_pwmcosns       = Vector{CGTopPWMcos}(undef, 0)
    # top_cg_pro_DNA_contact        = Vector{CGTopContact}(undef, 0)
    top_cg_pro_DNA_contact        = [[] for i in 1:cg_num_particles]

    # protein-RNA
    # top_cg_pro_RNA_contact   = Vector{CGTopContact}(undef, 0)
    # param_cg_pro_RNA_e_contact = []
    top_cg_pro_RNA_contact     = [[] for i in 1:cg_num_particles]

    # --------------------
    # geometric properties
    # --------------------
    # center of geometry
    # geo_centroid               = zeros(Float64, (3, aa_num_chain))
    # radius of gyration
    geo_radius_of_gyration     = zeros(Float64, aa_num_chain)
    # radius of circumsphere
    geo_radius_of_circumsphere = zeros(Float64, aa_num_chain)

    # =========================================================
    # Step 4: Determine CG particles and geometry of each chain
    # =========================================================
    i_step += 1
    if verbose
        println("============================================================")
        println("> Step $(i_step): determine protein/DNA/RNA chais.")
    end

    # -------
    # protein
    # -------

    if num_chain_pro > 0
        if verbose
            println("------------------------------------------------------------")
            println(">      PROTEIN: determine CA mass, charge, and coordinates.")
        end

        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if verbose
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / aa_num_chain  * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            end

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
    end

    # ---
    # DNA
    # ---

    if num_chain_DNA > 0

        if verbose
            println("\n------------------------------------------------------------")
            println(">      DNA: determine P, S, B mass, charge, and coordinates.")
        end

        for i_chain in 1 : aa_num_chain
            chain = cg_chains[i_chain]

            if verbose
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            end

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
    end

    # ---
    # RNA
    # ---

    if num_chain_RNA > 0

        if  verbose
            println("\n------------------------------------------------------------")
            println(">      RNA: determine P, S, B mass, charge, and coordinates.")
        end

        for i_chain in 1 : aa_num_chain
            chain = cg_chains[i_chain]

            if verbose
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            end

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
    end

    # ===================
    # Geometry properties
    # ===================
    if verbose
        println("\n------------------------------------------------------------")
        println(">      determine Centroid, Rg, Rc...")
    end

    for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if verbose
            print("\b"^32)
            progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
        end

        # centroid
        coor_centroid = zeros(Float64, 3)
        for i_res in chain.first : chain.last
            coor_centroid += cg_bead_coor[:, i_res]
        end
        coor_centroid /= (chain.last - chain.first + 1)
        # geo_centroid[:, i_chain] = coor_centroid

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

    if verbose
        println("\n>           ... geometric properties : DONE!")
    end


    # ====================================
    # Cell lists for contact determination
    # ================================================
    #   ____ _____ _     _       _     ___ ____ _____
    #  / ___| ____| |   | |     | |   |_ _/ ___|_   _|
    # | |   |  _| | |   | |     | |    | |\___ \ | |
    # | |___| |___| |___| |___  | |___ | | ___) || |
    #  \____|_____|_____|_____| |_____|___|____/ |_|
    #
    # ================================================
    i_step += 1
    if verbose
        println("============================================================")
        println("> Step $(i_step): Cell list construction.")
    end

    if verbose
        println("------------------------------------------------------------")
        println(">      $(i_step).1: Cell Division.")
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).1.1: Determine system size.")
    end

    sys_size_lower = minimum(cg_bead_coor, dims=2) .- 1.0
    sys_size_upper = maximum(cg_bead_coor, dims=2) .+ 1.0
    sys_size_3d = sys_size_upper - sys_size_lower

    if verbose
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).1.2: Determine cell size.")
    end

    cell_num_3d = max.(Int.(ceil.(sys_size_3d / ( 18.0 + AICG_GO_ATOMIC_CUTOFF ))) .- 1, 1)
    cell_num_all = prod(cell_num_3d)
    cell_size_3d = sys_size_3d ./ cell_num_3d

    if verbose
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).1.3: Prepare cells.")
    end

    cell_neighbors = [[] for i in 1:cell_num_all]
    for i in 1:cell_num_3d[1]
        for j in 1:cell_num_3d[2]
            for k in 1:cell_num_3d[3]
                box_indx = (i - 1) * cell_num_3d[2] * cell_num_3d[3] + (j -1) * cell_num_3d[3] + k
                for n_i in max(1, i - 1):min(cell_num_3d[1], i + 1)
                    for n_j in max(1, j - 1):min(cell_num_3d[2], j + 1)
                        for n_k in max(1, k - 1):min(cell_num_3d[3], k + 1)
                            nb_box_indx = (n_i - 1) * cell_num_3d[2] * cell_num_3d[3] + (n_j -1) * cell_num_3d[3] + n_k
                            push!(cell_neighbors[box_indx], nb_box_indx)
                        end
                    end
                end
            end
        end
    end

    # ------------------------
    # put particles into cells
    # ------------------------
    if verbose
        println("------------------------------------------------------------")
        println(">      $(i_step).2: Put particles into cells.")
    end

    cell_particles = [[] for i in 1:cell_num_all]
    cell_index_cg_bead = [0 for i in 1:cg_num_particles]
    for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if verbose
            print("\b"^32)
            progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
        end

        for i_res in chain.first : chain.last
            x, y, z = cg_bead_coor[:, i_res] .- sys_size_lower
            bi = Int.(ceil.(x / cell_size_3d[1]))
            bj = Int.(ceil.(y / cell_size_3d[2]))
            bk = Int.(ceil.(z / cell_size_3d[3]))
            box_indx = (bi - 1) * cell_num_3d[2] * cell_num_3d[3] + (bj -1) * cell_num_3d[3] + bk

            push!(cell_particles[box_indx], i_res)
            cell_index_cg_bead[i_res] = box_indx
        end
    end
    if verbose
        println("\n>           ... cell lists : DONE!")
    end



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
        if verbose
            println("============================================================")
            println("> Step $(i_step): processing proteins.")
        end

        # -----------------------------
        # Step 4.2: CG protein topology
        # -----------------------------
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).1: CG protein topology.")
            println(" - - - - - - - - - - - - - - - - - - - - - - - -")
            println(">      $(i_step).1.1: CG protein local interactions.")
        end
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
                if verbose && ( dist12 > 4.0 || dist12 < 3.6 )
                    errmsg = @sprintf("WARNING: abnormal bond length in chain %d, residue %d %s and %d %s, r0 = %8.3f",
                                      i_chain,
                                      i_res, cg_bead_name[i_res],
                                      i_res + 1, cg_bead_name[i_res + 1],
                                      dist12)
                    println(errmsg)
                end
            end
        end
        if verbose
            println(">           ... Bond: DONE!")
        end

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

                if ff_pro == FF_pro_AICG2p
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
        end
        if verbose
            println(">           ... Angle: DONE!")
        end

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

                if ff_pro == FF_pro_AICG2p
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
        end
        if verbose
            println(">           ... Dihedral: DONE!")
        end

        # ------------------------
        # Normalize local energies
        # ------------------------
        if ff_pro == FF_pro_AICG2p
            e_ground_local /= (num_angle + num_dih)
            e_ground_13    /= num_angle
            e_ground_14    /= num_dih

            if aicg_scale_scheme == 0
                for i in 1:length(param_cg_pro_e_13)
                    param_cg_pro_e_13[i] *= AICG_13_AVE / e_ground_13
                end
                for i in 1:length(param_cg_pro_e_14)
                    param_cg_pro_e_14[i] *= AICG_14_AVE / e_ground_14
                end
            elseif aicg_scale_scheme == 1
                for i in 1:length(param_cg_pro_e_13)
                    param_cg_pro_e_13[i] *= -AICG_13_GEN
                end
                for i in 1:length(param_cg_pro_e_14)
                    param_cg_pro_e_14[i] *= -AICG_14_GEN
                end
            end
        end

        # -----------------------
        # Go type native contacts
        # -----------------------
        if verbose
            println(" - - - - - - - - - - - - - - - - - - - - - - - -")
            println(">      $(i_step).1.2: Looking for native contacts.")
        end

        # intra-molecular contacts
        if verbose
            @printf("%11s Calculating intra-molecular contacts... \n", " ")
            @printf("              ... chain   : %32s", " ")
        end
        for i_chain in 1:aa_num_chain

            if do_test_local_only
                continue
            end

            chain = cg_chains[i_chain]

            # -----------------
            # show progress bar
            # -----------------
            if verbose
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            end

            if chain.moltype != MOL_PROTEIN
                continue
            end

            if chain.last - chain.first < 100
                Threads.@threads for i_res in chain.first : chain.last - 4
                    cell_i = cell_index_cg_bead[i_res]
                    coor_cai = cg_bead_coor[:, i_res]
                    neighbor_cell_i = cell_neighbors[cell_i]
                    for j_res in i_res + 4 : chain.last
                        cell_j = cell_index_cg_bead[j_res]
                        if ! (cell_j in neighbor_cell_i)
                            continue
                        end
                        coor_caj = cg_bead_coor[:, j_res]
                        if is_protein_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                            native_dist = compute_distance(coor_cai, coor_caj)
                            num_cg_pro_contact_all += 1
                            num_cg_pro_contact_intra += 1

                            e_local = 0.0
                            if ff_pro == FF_pro_AICG2p
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
                            end
                            push!(top_cg_pro_go_contact[i_res], [j_res, native_dist, e_local])
                        end
                    end
                end
            else
                Threads.@threads for i_res in chain.first : chain.last
                    cell_i = cell_index_cg_bead[i_res]
                    coor_cai = cg_bead_coor[:, i_res]
                    neighbor_cell_i = cell_neighbors[cell_i]
                    for j_cell in neighbor_cell_i
                        for j_res in cell_particles[j_cell]
                            if cg_chain_id[j_res] != i_chain
                                continue
                            end
                            if j_res < i_res + 4
                                continue
                            end
                            coor_caj = cg_bead_coor[:, j_res]
                            if is_protein_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                                native_dist = compute_distance(coor_cai, coor_caj)
                                num_cg_pro_contact_all += 1
                                num_cg_pro_contact_inter += 1

                                e_local = 0.0
                                if ff_pro == FF_pro_AICG2p
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
                                end
                                push!(top_cg_pro_go_contact[i_res], [j_res, native_dist, e_local])
                            end
                        end
                    end
                end
            end
        end
        if verbose
            print("\n              ... intra-molecular contacts: DONE! \n")
        end

        # inter-molecular ( protein-protein ) contacts
        if num_chain_pro > 1 && !do_test_local_only
            if verbose
                @printf("%11s Calculating inter-molecular contacts... \n", " ")
                @printf("              ... progress: %32s", " ")
            end
            for i_chain in 1 : aa_num_chain - 1
                chain1 = cg_chains[i_chain]

                # -----------------
                # show progress bar
                # -----------------
                if verbose
                    print("\b"^32)
                    progress_percent = trunc(Int, i_chain / ( aa_num_chain - 1 ) * 20)
                    progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                    @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / ( aa_num_chain - 1 ) * 100)
                end
                # ------------------

                if chain1.moltype != MOL_PROTEIN
                    continue
                end

                Threads.@threads for i_res in chain1.first : chain1.last
                    cell_i = cell_index_cg_bead[i_res]
                    coor_cai = cg_bead_coor[:, i_res]
                    neighbor_cell_i = cell_neighbors[cell_i]
                    for j_cell in neighbor_cell_i
                        for j_res in cell_particles[j_cell]
                            if j_res < i_res
                                continue
                            end
                            if cg_chains[cg_chain_id[j_res]].moltype != MOL_PROTEIN
                                continue
                            end
                            if cg_chain_id[j_res] == i_chain
                                continue
                            end
                            coor_caj = cg_bead_coor[:, j_res]
                            if is_protein_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                                native_dist = compute_distance(coor_cai, coor_caj)
                                num_cg_pro_contact_all += 1
                                num_cg_pro_contact_inter += 1

                                e_local = 0.0
                                if ff_pro == FF_pro_AICG2p
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
                                end
                                push!(top_cg_pro_go_contact[i_res], [j_res, native_dist, e_local])
                            end
                        end
                    end
                end
            end
            if verbose
                print("\n              ... inter-molecular contacts: DONE! \n")
            end
        end

        if ff_pro == FF_pro_AICG2p
            e_ground_contact = 0.0
            num_contact = 0
            # count num of contacts, sum up e_ground_contact
            for i_res in 1:cg_num_particles
                for cntct_tmp in top_cg_pro_go_contact[i_res]
                    num_contact += 1
                    e_ground_contact += cntct_tmp[3]
                end
            end
            # normalize
            if num_contact > 0
                e_ground_contact /= num_contact
            else
                e_ground_contact = 0
            end

            if aicg_scale_scheme == 0
                for i in 1:cg_num_particles
                    for cntct_tmp in top_cg_pro_go_contact[i]
                        cntct_tmp[3] *= AICG_CONTACT_AVE / e_ground_contact
                    end
                end
            elseif aicg_scale_scheme == 1
                for i in 1:cg_num_particles
                    for cntct_tmp in top_cg_pro_go_contact[i]
                        cntct_tmp[3] *= -AICG_CONTACT_GEN
                    end
                end
            end
        end

        if verbose
            println("------------------------------------------------------------")
            @printf("          > Total number of protein contacts: %12d  \n",
                    sum( length.( top_cg_pro_go_contact ) ))
        end

    end

    # ===============================
    # Intrinsically disordered region
    # ===============================
    if haskey(args, "modeling-options")
        has_toml_mod    = true
        ff_detail_config = args["modeling-options"]

        if haskey(ff_detail_config, "IDR")
            if haskey(ff_detail_config["IDR"], "AICG2p_IDR_local")
                index_string = ff_detail_config["IDR"]["AICG2p_IDR_local"]
                AICG2p_flexible_local = parse_selection(index_string)
            end
            if haskey(ff_detail_config["IDR"], "AICG2p_IDR_nonlocal")
                index_string = ff_detail_config["IDR"]["AICG2p_IDR_nonlocal"]
                AICG2p_flexible_nonlocal = parse_selection(index_string)
            end
            if haskey(ff_detail_config["IDR"], "HPS_region")
                index_string = ff_detail_config["IDR"]["HPS_region"]
                HPS_IDR_region = parse_selection(index_string)
            end
            if haskey(ff_detail_config["IDR"], "KH_region")
                index_string = ff_detail_config["IDR"]["KH_region"]
                KH_IDR_region = parse_selection(index_string)
            end
        end
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
        if verbose
            println("============================================================")
            println("> Step $(i_step): processing DNA.")
        end

        # ---------------------------------
        #        Step 5.1: 3SPN.2C topology
        # ---------------------------------
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).1: 3SPN.2C topology.")
            println(" - - - - - - - - - - - - - - - - - - - - - - - -")
            println(">      $(i_step).1.1: 3SPN.2C local interactions.")
        end
        for i_chain in 1:aa_num_chain

            chain = cg_chains[i_chain]

            if chain.moltype != MOL_DNA
                continue
            end

            if verbose
                @printf("%11s Calculating DNA strand %d ... \n", " ", i_chain)
                @printf("              ... progress: %32s", " ")
            end

            if DNA_circular
                DNA_basetype_pre  = cg_resid_name[chain.last][end]
                DNA_basetype_post = cg_resid_name[chain.first][end]
            else
                DNA_basetype_pre = "X"
            end
            if gen_3spn_itp == 1

                for i_res in chain.first : chain.last

                    # -----------------
                    # show progress bar
                    # -----------------
                    if verbose
                        print("\b"^32)
                        progress_percent = trunc(Int, ( i_res - chain.first ) / ( chain.last - chain.first ) * 20)
                        progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                        @printf(" [%20s] %5.1f %% ", progress_bar, progress_percent * 5)
                    end
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
                        resname5  = i_res > chain.first ? cg_resid_name[i_res - 1][end] : DNA_basetype_pre
                        # resname5 = cg_resid_name[i_res - 1][end]
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
                end                 # for i_res

            end                 # if gen_3spn_itp == 1

            if gen_3spn_itp == 2

                for i_res in chain.first : chain.last

                    # -----------------
                    # show progress bar
                    # -----------------
                    if verbose
                        print("\b"^32)
                        progress_percent = trunc(Int, ( i_res - chain.first ) / ( chain.last - chain.first ) * 20)
                        progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                        @printf(" [%20s] %5.1f %% ", progress_bar, progress_percent * 5)
                    end
                    # ------------------

                    if cg_bead_name[i_res] == "DS"
                        # bond S--B
                        resname5     = cg_resid_name[i_res][end]
                        r_sb         = get_DNA3SPN_bond_length("SB", resname5 * " ")
                        tmp_top_bond = CGTopBond(i_res, i_res + 1, r_sb)
                        push!(top_cg_DNA_bonds, tmp_top_bond)
                        if i_res + 3 < chain.last
                            resname3  = cg_resid_name[i_res + 3][end]
                            base_step = resname5 * resname3
                            # bond S--P+1
                            r_sp3        = get_DNA3SPN_bond_length("SP", base_step)
                            tmp_top_bond = CGTopBond(i_res, i_res + 2, r_sp3)
                            push!(top_cg_DNA_bonds, tmp_top_bond)
                            # Angle S--P+1--S+1
                            ang_sp3s3    = get_DNA3SPN_angle_equilibrium("SPS", base_step)
                            k            = get_DNA3SPN_angle_param("SPS", base_step)
                            tmp_top_angl = CGTopAngle(i_res, i_res + 2, i_res + 3, ang_sp3s3)
                            push!(top_cg_DNA_angles, tmp_top_angl)
                            push!(param_cg_DNA_k_angles, k)
                            # Dihedral S--P+1--S+1--B+1
                            dih_sp3s3b3 = get_DNA3SPN_dihedral_equilibrium("SPSB", base_step)
                            tmp_top_dih = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 4, dih_sp3s3b3)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                            # Dihedral S--P+1--S+1--P+2
                            if i_res + 6 < chain.last
                                base_two_steps = base_step * cg_resid_name[i_res + 5][end]
                                dih_sp3s3p33   = get_DNA3SPN_dihedral_equilibrium("SPSP", base_two_steps)
                                tmp_top_dih    = CGTopDihedral(i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33)
                                push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                                push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                            elseif DNA_circular
                                base_two_steps = base_step * DNA_basetype_post
                                dih_sp3s3p33   = get_DNA3SPN_dihedral_equilibrium("SPSP", base_two_steps)
                                tmp_top_dih    = CGTopDihedral(i_res, i_res + 2, i_res + 3, chain.first, dih_sp3s3p33)
                                push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                                push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                            end
                        elseif DNA_circular
                            resname3  = DNA_basetype_post
                            base_step = resname5 * resname3
                            # bond S--P+1
                            r_sp3        = get_DNA3SPN_bond_length("SP", base_step)
                            tmp_top_bond = CGTopBond(i_res, chain.first, r_sp3)
                            push!(top_cg_DNA_bonds, tmp_top_bond)
                            # Angle S--P+1--S+1
                            ang_sp3s3    = get_DNA3SPN_angle_equilibrium("SPS", base_step)
                            k            = get_DNA3SPN_angle_param("SPS", base_step)
                            tmp_top_angl = CGTopAngle(i_res, chain.first, chain.first + 1, ang_sp3s3)
                            push!(top_cg_DNA_angles, tmp_top_angl)
                            push!(param_cg_DNA_k_angles, k)
                            # Dihedral S--P+1--S+1--B+1
                            dih_sp3s3b3 = get_DNA3SPN_dihedral_equilibrium("SPSB", base_step)
                            tmp_top_dih = CGTopDihedral(i_res, chain.first, chain.first + 1, chain.first + 2, dih_sp3s3b3)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                            # Dihedral S--P+1--S+1--P+2
                            base_two_steps = base_step * cg_resid_name[chain.first + 3][end]
                            dih_sp3s3p33   = get_DNA3SPN_dihedral_equilibrium("SPSP", base_two_steps)
                            tmp_top_dih    = CGTopDihedral(i_res, chain.first, chain.first + 1, chain.first + 3, dih_sp3s3p33)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                            push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                        end
                    elseif cg_bead_name[i_res] == "DP"
                        resname5  = i_res > chain.first ? cg_resid_name[i_res - 1][end] : DNA_basetype_pre
                        resname3  = cg_resid_name[i_res + 2][end]
                        base_step = resname5 * resname3
                        # bond P--S
                        r_ps         = get_DNA3SPN_bond_length("PS", base_step)
                        tmp_top_bond = CGTopBond(i_res, i_res + 1, r_ps)
                        push!(top_cg_DNA_bonds, tmp_top_bond)
                        # angle P--S--B
                        ang_psb      = get_DNA3SPN_angle_equilibrium("PSB", base_step)
                        k            = get_DNA3SPN_angle_param("PSB", base_step)
                        tmp_top_angl = CGTopAngle(i_res, i_res + 1, i_res + 2, ang_psb)
                        push!(top_cg_DNA_angles, tmp_top_angl)
                        push!(param_cg_DNA_k_angles, k)
                        if i_res + 4 < chain.last
                            base_two_steps = base_step * cg_resid_name[i_res + 3][end]
                            # angle P--S--P+1
                            ang_psp3 = get_DNA3SPN_angle_equilibrium("PSP", base_two_steps)
                            k        = get_DNA3SPN_angle_param("PSP", "all")
                            tmp_top_angl = CGTopAngle(i_res, i_res + 1, i_res + 3, ang_psp3)
                            push!(top_cg_DNA_angles, tmp_top_angl)
                            push!(param_cg_DNA_k_angles, k)
                            # Dihedral P--S--P+1--S+1
                            dih_psp3s3 = get_DNA3SPN_dihedral_equilibrium("PSPS", base_two_steps)
                            tmp_top_dih = CGTopDihedral(i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                            push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                        elseif DNA_circular
                            base_two_steps = base_step * DNA_basetype_post
                            # angle P--S--P+1
                            ang_psp3 = get_DNA3SPN_angle_equilibrium("PSP", base_two_steps)
                            k        = get_DNA3SPN_angle_param("PSP", "all")
                            tmp_top_angl = CGTopAngle(i_res, i_res + 1, chain.first, ang_psp3)
                            push!(top_cg_DNA_angles, tmp_top_angl)
                            push!(param_cg_DNA_k_angles, k)
                            # Dihedral P--S--P+1--S+1
                            dih_psp3s3 = get_DNA3SPN_dihedral_equilibrium("PSPS", base_two_steps)
                            tmp_top_dih = CGTopDihedral(i_res, i_res + 1, chain.first, chain.first + 1, dih_psp3s3)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                            push!(top_cg_DNA_dih_Gaussian, tmp_top_dih)
                        end
                    elseif cg_bead_name[i_res] == "DB"
                        resname5  = cg_resid_name[i_res][end]
                        if i_res + 2 < chain.last
                            resname3  = cg_resid_name[i_res + 1][end]
                            base_step = resname5 * resname3
                            # angle B--S--P+1
                            ang_bsp3     = get_DNA3SPN_angle_equilibrium("BSP", base_step)
                            k            = get_DNA3SPN_angle_param("BSP", base_step)
                            tmp_top_angl = CGTopAngle(i_res, i_res - 1, i_res + 1, ang_bsp3)
                            push!(top_cg_DNA_angles, tmp_top_angl)
                            push!(param_cg_DNA_k_angles, k)
                            # Dihedral B--S--P+1--S+1
                            dih_bsp3s3 = get_DNA3SPN_dihedral_equilibrium("BSPS", base_step)
                            tmp_top_dih = CGTopDihedral(i_res, i_res - 1, i_res + 1, i_res + 2, dih_bsp3s3)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                        elseif DNA_circular
                            resname3  = DNA_basetype_post
                            base_step = resname5 * resname3
                            # angle B--S--P+1
                            ang_bsp3     = get_DNA3SPN_angle_equilibrium("BSP", base_step)
                            k            = get_DNA3SPN_angle_param("BSP", base_step)
                            tmp_top_angl = CGTopAngle(i_res, i_res - 1, chain.first, ang_bsp3)
                            push!(top_cg_DNA_angles, tmp_top_angl)
                            push!(param_cg_DNA_k_angles, k)
                            # Dihedral B--S--P+1--S+1
                            dih_bsp3s3 = get_DNA3SPN_dihedral_equilibrium("BSPS", base_step)
                            tmp_top_dih = CGTopDihedral(i_res, i_res - 1, chain.first, chain.first + 1, dih_bsp3s3)
                            push!(top_cg_DNA_dih_periodic, tmp_top_dih)
                        end
                    else
                        errmsg = @sprintf("BUG: Wrong DNA particle type in chain %d, residue %d : %s ",
                                          i_chain,
                                          i_res, cg_bead_name[i_res])
                        error(errmsg)
                    end
                end                 # for i_res

            end                 # if gen_3spn_itp == 2

            if verbose
                print(" \n")
            end

        end                     # for i_chain

        if verbose
            println(">           ... Bond, Angle, Dihedral: DONE!")
        end
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
        if verbose
            println("============================================================")
            println("> Step $(i_step): processing RNA.")
        end

        # -------------------------
        # Step 6.1: RNA topology
        # -------------------------
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).1: RNA topology.")
            println(" - - - - - - - - - - - - - - - - - - - - - - - -")
            println(">      $(i_step).1.1: RNA local interactions.")

            @printf("%11s Calculating bonded terms... \n", " ")
        end
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
                    # add fake angles and dihedrals...
                    coor_b = cg_bead_coor[:, i_res]
                    coor_s = cg_bead_coor[:, i_res - 1]
                    if i_res + 1 < chain.last
                        # angle B--S--P+1
                        coor_p3  = cg_bead_coor[:, i_res + 1]
                        ang_bsp3 = compute_angle(coor_b, coor_s, coor_p3)
                        k        = 0.0
                        tmp_top_angl = CGTopAngle(i_res, i_res - 1, i_res + 1, ang_bsp3)
                        push!(top_cg_RNA_angles, tmp_top_angl)
                        push!(param_cg_RNA_k_angles, k)
                    end
                    if i_res + 2 < chain.last
                        # Dihedral B--S--P+1--S+1
                        coor_p3    = cg_bead_coor[:, i_res + 1]
                        coor_s3    = cg_bead_coor[:, i_res + 2]
                        dih_bsp3s3 = compute_dihedral(coor_b, coor_s, coor_p3, coor_s3)
                        k          = 0.0
                        tmp_top_dih = CGTopDihedral(i_res, i_res - 1, i_res + 1, i_res + 2, dih_bsp3s3)
                        push!(top_cg_RNA_dihedrals, tmp_top_dih)
                        push!(param_cg_RNA_k_dihedrals, k)
                    end
                end
            end
        end

        # -----------------------
        # Go type native contacts
        # -----------------------
        if verbose
            println(" - - - - - - - - - - - - - - - - - - - - - - - -")
            println(">      $(i_step).1.2: RNA HT-type native contacts.")
            @printf("%11s Calculating intra-molecular contacts... \n", " ")
            @printf("              ... progress: %32s", " ")
        end
        for i_chain in 1:aa_num_chain

            if do_test_local_only
                continue
            end

            chain = cg_chains[i_chain]

            # -----------------
            # show progress bar
            # -----------------
            if verbose
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            end
            # ------------------

            if chain.moltype != MOL_RNA
                continue
            end

            Threads.@threads for i_res in chain.first : chain.last - 3

                if cg_bead_name[i_res] == "RP" && ! cgRNA_use_phosphate_go
                    continue
                end

                coor_i = cg_bead_coor[:, i_res]
                cell_i = cell_index_cg_bead[i_res]
                neighbor_cell_i = cell_neighbors[cell_i]

                for j_cell in neighbor_cell_i
                    for j_res in cell_particles[j_cell]
                        if cg_chain_id[j_res] != i_chain
                            continue
                        end
                        if j_res < i_res + 3
                            continue
                        end
                        if cgRNA_use_phosphate_go
                            if cg_bead_name[i_res] == "RP" || cg_bead_name[j_res] == "RP"
                                if j_res < i_res + 4
                                    continue
                                end
                            end
                        elseif cg_bead_name[j_res] == "RP"
                            continue
                        end

                        if cg_bead_name[i_res] == "RS" || cg_bead_name[j_res] == "RS"
                            if j_res < i_res + 6
                                continue
                            end
                        end

                        coor_j = cg_bead_coor[:, j_res]

                        adist, nhb  = compute_RNA_native_contact(cg_residues[i_res].atoms,
                                                                 cg_residues[j_res].atoms,
                                                                 aa_atom_name,
                                                                 aa_coor)
                        if adist > RNA_GO_ATOMIC_CUTOFF
                            continue
                        end

                        native_dist = compute_distance(coor_i, coor_j)

                        if j_res == i_res + 3 && cg_bead_name[i_res] == "RB"
                            coor_i_sug = cg_bead_coor[:, i_res - 1]
                            coor_j_sug = cg_bead_coor[:, j_res - 1]
                            st_dih = compute_dihedral(coor_i, coor_i_sug, coor_j_sug, coor_j)
                            if abs( st_dih ) < RNA_STACK_DIH_CUTOFF && adist < RNA_STACK_DIST_CUTOFF
                                # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                # push!(top_cg_RNA_base_stack, tmp_top_cnt)
                                # push!(param_cg_RNA_e_base_stack, RNA_STACK_EPSILON)
                                push!(top_cg_RNA_base_stack[i_res], [j_res, native_dist, RNA_STACK_EPSILON])
                            else
                                # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                # push!(top_cg_RNA_other_contact, tmp_top_cnt)
                                # push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER["BB"])
                                push!(top_cg_RNA_other_contact[i_res], [j_res, native_dist, RNA_PAIR_EPSILON_OTHER["BB"]])
                            end
                        elseif cg_bead_name[i_res] == "RB" && cg_bead_name[j_res] == "RB"
                            if nhb == 2
                                # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                # push!(top_cg_RNA_base_pair, tmp_top_cnt)
                                # push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_2HB)
                                push!(top_cg_RNA_base_pair[i_res], [j_res, native_dist, RNA_BPAIR_EPSILON_2HB])
                            elseif nhb >= 3
                                # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                # push!(top_cg_RNA_base_pair, tmp_top_cnt)
                                # push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_3HB)
                                push!(top_cg_RNA_base_pair[i_res], [j_res, native_dist, RNA_BPAIR_EPSILON_3HB])
                            else
                                # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                # push!(top_cg_RNA_other_contact, tmp_top_cnt)
                                # push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER["BB"])
                                push!(top_cg_RNA_other_contact[i_res], [j_res, native_dist, RNA_PAIR_EPSILON_OTHER["BB"]])
                            end
                        else
                            contact_type = cg_bead_name[i_res][end] * cg_bead_name[j_res][end]
                            # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            # push!(top_cg_RNA_other_contact, tmp_top_cnt)
                            # push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER[contact_type])
                            push!(top_cg_RNA_other_contact[i_res], [j_res, native_dist, RNA_PAIR_EPSILON_OTHER[contact_type]])
                        end
                    end
                end
            end
        end
        if verbose
            print("\n              ... intra-molecular contacts: DONE! \n")
        end

        if num_chain_RNA > 1 && !do_test_local_only
            if verbose
                @printf("%11s Calculating inter-molecular contacts... \n", " ")
                @printf("              ... progress: %32s", " ")
            end
            for i_chain in 1:aa_num_chain - 1

                chain_1 = cg_chains[i_chain]

                # -----------------
                # show progress bar
                # -----------------
                if verbose
                    print("\b"^32)
                    progress_percent = trunc(Int, i_chain / ( aa_num_chain - 1 ) * 20)
                    progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                    @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / ( aa_num_chain - 1 ) * 100)
                end
                # ------------------

                if chain_1.moltype != MOL_RNA
                    continue
                end

                Threads.@threads for i_res in chain_1.first : chain_1.last
                    if cg_bead_name[i_res] == "RP" && ! cgRNA_use_phosphate_go
                        continue
                    end
                    coor_i = cg_bead_coor[:, i_res]
                    cell_i = cell_index_cg_bead[i_res]
                    neighbor_cell_i = cell_neighbors[cell_i]

                    for j_cell in neighbor_cell_i
                        for j_res in cell_particles[j_cell]
                            if j_res < i_res
                                continue
                            end
                            if cg_chains[cg_chain_id[j_res]].moltype != MOL_RNA
                                continue
                            end
                            if cg_chain_id[j_res] == i_chain
                                continue
                            end
                            if cg_bead_name[j_res] == "RP" && ! cgRNA_use_phosphate_go
                                continue
                            end
                            coor_j = cg_bead_coor[:, j_res]
                            adist, nhb  = compute_RNA_native_contact(cg_residues[i_res].atoms,
                                                                     cg_residues[j_res].atoms,
                                                                     aa_atom_name,
                                                                     aa_coor)
                            if adist > RNA_GO_ATOMIC_CUTOFF
                                continue
                            end
                            native_dist = compute_distance(coor_i, coor_j)
                            if cg_bead_name[i_res] == "RB" && cg_bead_name[j_res] == "RB"
                                if nhb == 2
                                    # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                    # push!(top_cg_RNA_base_pair, tmp_top_cnt)
                                    # push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_2HB)
                                    push!(top_cg_RNA_base_pair[i_res], [j_res, native_dist, RNA_BPAIR_EPSILON_2HB])
                                elseif nhb >= 3
                                    # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                    # push!(top_cg_RNA_base_pair, tmp_top_cnt)
                                    # push!(param_cg_RNA_e_base_pair, RNA_BPAIR_EPSILON_3HB)
                                    push!(top_cg_RNA_base_pair[i_res], [j_res, native_dist, RNA_BPAIR_EPSILON_3HB])
                                else
                                    # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                    # push!(top_cg_RNA_other_contact, tmp_top_cnt)
                                    # push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER["BB"])
                                    push!(top_cg_RNA_other_contact[i_res], [j_res, native_dist, RNA_PAIR_EPSILON_OTHER["BB"]])
                                end
                            else
                                contact_type = cg_bead_name[i_res][end] * cg_bead_name[j_res][end]
                                # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                                # push!(top_cg_RNA_other_contact, tmp_top_cnt)
                                # push!(param_cg_RNA_e_other_contact, RNA_PAIR_EPSILON_OTHER[contact_type])
                                push!(top_cg_RNA_other_contact[i_res], [j_res, native_dist, RNA_PAIR_EPSILON_OTHER[contact_type]])
                            end
                        end
                    end
                end
            end
            if verbose
                print("\n              ... inter-molecular contacts: DONE! \n")
            end
        end

        if verbose
            println("------------------------------------------------------------")
            @printf("          > Total number of RNA contacts:     %12d  \n",
                    sum(length.(top_cg_RNA_base_stack)) +
                    sum(length.(top_cg_RNA_base_pair)) +
                    sum(length.(top_cg_RNA_other_contact)) )
        end

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

    if ff_pro_rna == FF_pro_RNA_Go && num_chain_RNA > 0 && num_chain_pro > 0 && !do_test_local_only
        i_step += 1
        if verbose
            println("============================================================")
            println("> Step $(i_step): Generating protein-RNA native contacts.")

            @printf("%11s Calculating protein-RNA contacts... \n", " ")
            @printf("              ... progress: %32s", " ")
        end
        for i_chain in 1:aa_num_chain
            # -----------------
            # show progress bar
            # -----------------
            if verbose
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            end
            # ------------------

            chain_pro = cg_chains[i_chain]

            if chain_pro.moltype != MOL_PROTEIN
                continue
            end


            Threads.@threads for i_res in chain_pro.first : chain_pro.last
                cell_i = cell_index_cg_bead[i_res]
                neighbor_cell_i = cell_neighbors[cell_i]
                coor_i = cg_bead_coor[:, i_res]

                for j_cell in neighbor_cell_i
                    for j_res in cell_particles[j_cell]
                        if cg_chains[cg_chain_id[j_res]].moltype != MOL_RNA
                            continue
                        end
                        if cg_bead_name[j_res] == "RP" && ! cgRNA_use_phosphate_go
                            continue
                        end
                        if !is_protein_RNA_native_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                            continue
                        end
                        coor_j = cg_bead_coor[:, j_res]
                        native_dist = compute_distance(coor_i, coor_j)
                        if cg_bead_name[j_res] == "RS"
                            # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            # push!(top_cg_pro_RNA_contact, tmp_top_cnt)
                            # push!(param_cg_pro_RNA_e_contact, PRO_RNA_GO_EPSILON_S)
                            push!(top_cg_pro_RNA_contact[i_res], [j_res, native_dist, PRO_RNA_GO_EPSILON_S])
                        elseif cg_bead_name[j_res] == "RB"
                            # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            # push!(top_cg_pro_RNA_contact, tmp_top_cnt)
                            # push!(param_cg_pro_RNA_e_contact, PRO_RNA_GO_EPSILON_B)
                            push!(top_cg_pro_RNA_contact[i_res], [j_res, native_dist, PRO_RNA_GO_EPSILON_B])
                        elseif cg_bead_name[j_res] == "RP" && ! cgRNA_use_phosphate_go
                            # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                            # push!(top_cg_pro_RNA_contact, tmp_top_cnt)
                            # push!(param_cg_pro_RNA_e_contact, PRO_RNA_GO_EPSILON_P)
                            push!(top_cg_pro_RNA_contact[i_res], [j_res, native_dist, PRO_RNA_GO_EPSILON_P])
                        end
                    end
                end
            end
        end

        if verbose
            println("\n------------------------------------------------------------")
            @printf("          > Total number of protein-RNA contacts: %8d  \n",
                    sum( length.(top_cg_pro_RNA_contact) ) )
        end
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
        if verbose
            println("============================================================")
            println("> Step $(i_step): Generating PWMcos parameters.")
        end

        # ------------------------------------------------
        #        Step 7.1: determine protein-DNA contacts
        # ------------------------------------------------
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).1: determine contacts between protein and DNA.")
        end

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
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).2: read in position frequency matrix (PFM).")
        end

        pwmcos_pwm, pwmcos_chain_a, pwmcos_chain_b = read_modified_pfm(pfm_filename)
        num_pwmcos_terms = length(pwmcos_chain_a)


        # ------------------------------------------------
        #        Step 7.2: Read in PFM and convert to PWM
        # ------------------------------------------------
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).3: decomposing PWM.")
        end

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

        if verbose
            println(">           ... DONE!")
        end
    end


    # ===================================================================
    # PWMcos-ns parameters: protein-DNA sequence-non-specific interaction
    # ===================================================================
    #    ______        ____  __
    #   |  _ \ \      / /  \/  | ___ ___  ___       _ __  ___
    #   | |_) \ \ /\ / /| |\/| |/ __/ _ \/ __|_____| '_ \/ __|
    #   |  __/ \ V  V / | |  | | (_| (_) \__ \_____| | | \__ \
    #   |_|     \_/\_/  |_|  |_|\___\___/|___/     |_| |_|___/
    # 
    # ============================================================
    if ff_pro_dna == FF_PWMcos_ns

        # pwmcos_ns_native_contacts = []

        if num_chain_pro == 0
            error("Cannot generate PWMcos parameters without protein...")
        end

        i_step += 1
        if verbose
            println("============================================================")
            println("> Step $(i_step): Generating PWMcos-ns parameters.")
        end

        # ------------------------------------------------
        #        Step 7.1: determine protein-DNA contacts
        # ------------------------------------------------
        if verbose
            println("------------------------------------------------------------")
            println(">      $(i_step).1: determine contacts between CA and DP.")
        end

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

                    for j_res in chain_DNA.first + 2 : chain_DNA.last - 2
                        if cg_bead_name[j_res] != "DP"
                            continue
                        end
                        if !is_PWMcos_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                            continue
                        end

                        coor_dna_j       = cg_bead_coor[:, j_res]
                        coor_dna_S       = cg_bead_coor[:, j_res + 1]

                        vec0   = coor_pro_i - coor_dna_j
                        vec1   = coor_dna_S - coor_dna_j
                        vec3   = coor_pro_N - coor_pro_C
                        r0     = norm(vec0)
                        theta1 = compute_vec_angle(vec0, vec1)
                        theta3 = compute_vec_angle(vec0, vec3)

                        # push!(pwmcos_ns_native_contacts, (i_res - chain_pro.first + 1,
                        #                                   r0,
                        #                                   theta1,
                        #                                   theta3))

                        if do_debug
                            println("PWMcos | pro ===> ", i_res - chain_pro.first + 1,
                                    " DNA ===> ", j_res, " : ", cg_resid_index[j_res],
                                    " r0 = ", r0,
                                    " theta1 = ", theta1,
                                    " theta3 = ", theta3)
                        end

                        tmp_top_pwmcosns = CGTopPWMcos(i_res - chain_pro.first + 1,
                                                       r0, theta1, 0.0, theta3,
                                                       0.0, 0.0, 0.0, 0.0)
                        push!(top_cg_pro_DNA_pwmcosns, tmp_top_pwmcosns)

                    end
                end
            end
        end

        if verbose
            println(">           ... DONE!")
        end
    end


    # ============================================================
    # Protein-DNA structure-based Go-like interaction
    # ============================================================
    #                        ____  _   _    _       ____
    #  _ __  _ __ ___       |  _ \| \ | |  / \     / ___| ___
    # | '_ \| '__/ _ \ _____| | | |  \| | / _ \   | |  _ / _ \
    # | |_) | | | (_) |_____| |_| | |\  |/ ___ \  | |_| | (_) |
    # | .__/|_|  \___/      |____/|_| \_/_/   \_\  \____|\___/
    # |_|
    # ============================================================
    if ff_pro_dna == FF_pro_DNA_Go

        if num_chain_pro == 0
            error("Cannot generate protein-DNA parameters without protein...")
        end
        if num_chain_DNA == 0
            error("Cannot generate protein-DNA parameters without DNA...")
        end

        i_step += 1
        if verbose
            println("============================================================")
            println("> Step $(i_step): Generating protein-DNA Go-like interactions.")
        end

        for i_chain in 1:aa_num_chain
            chain_pro = cg_chains[i_chain]

            if chain_pro.moltype != MOL_PROTEIN
                continue
            end

            Threads.@threads for i_res in chain_pro.first : chain_pro.last
                cell_i = cell_index_cg_bead[i_res]
                neighbor_cell_i = cell_neighbors[cell_i]
                coor_i = cg_bead_coor[:, i_res]

                for j_cell in neighbor_cell_i
                    for j_res in cell_particles[j_cell]

                        if cg_chains[cg_chain_id[j_res]].moltype != MOL_DNA
                            continue
                        end
                        if !is_protein_DNA_Go_contact(cg_residues[i_res].atoms, cg_residues[j_res].atoms, aa_atom_name, aa_coor)
                            continue
                        end

                        coor_j = cg_bead_coor[:, j_res]
                        native_dist = compute_distance(coor_i, coor_j)
                        # tmp_top_cnt = CGTopContact(i_res, j_res, native_dist)
                        # push!(top_cg_pro_DNA_contact, tmp_top_cnt)
                        push!(top_cg_pro_DNA_contact[i_res], [j_res, native_dist])

                    end
                end
            end
        end

        if verbose
            println("\n------------------------------------------------------------")
            @printf("          > Total number of protein-DNA contacts: %8d  \n",
                    sum( length.(top_cg_pro_DNA_contact) ) )
        end
    end



    # =============================
    # Make a new topology structure
    # =============================

    top_default_params             = GenTopDefault(0, 0, false, 0.0, 0.0) 
    top_default_atomtype           = Vector{GenTopAtomType}(undef, 0)
    top_default_CGDNA_bp           = Vector{GenTopCGDNABasepairType}(undef, 0)
    top_default_CGDNA_bs           = Vector{GenTopCGDNABasestackType}(undef, 0)
    top_default_CGDNA_cs           = Vector{GenTopCGDNABasecrossType}(undef, 0)
    top_default_CGDNA_exv          = Vector{GenTopCGDNAExvType}(undef, 0)
    top_default_CGPro_flx_angle    = Vector{GenTopCGProAICGFlexAngleType}(undef, 0)
    top_default_CGPro_flx_dihedral = Vector{GenTopCGProAICGFlexDihedralType}(undef, 0)

    global_index_2_local_index     = Vector{Int}(undef, 0)
    global_index_2_local_molid     = Vector{Int}(undef, 0)
    top_atoms                      = Vector{GenTopAtom}(undef, 0)
    top_bonds                      = Vector{GenTopBond}(undef, 0)
    top_angles                     = Vector{GenTopAngle}(undef, 0)
    top_dihedrals                  = Vector{GenTopDihedral}(undef, 0)
    top_pairs                      = Vector{GenTopPair}(undef, 0)
    top_exclusions                 = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos                     = Vector{GenTopPWMcos}(undef, 0)
    top_pwmcosns                   = Vector{GenTopPWMcos}(undef, 0)
    top_idr_hps                    = Vector{GenTopRegion}(undef, 0)
    top_idr_kh                     = Vector{GenTopRegion}(undef, 0)
    top_mol_list                   = Vector{GenTopMolList}(undef, 0)

    # ---------
    # [ atoms ]
    # ---------
    for i_bead in 1 : cg_num_particles
        a_type = cg_bead_type[i_bead]
        r_indx = cg_resid_index[i_bead]
        r_name = cg_resid_name[i_bead]
        a_name = cg_bead_name[i_bead]
        f_type = AICG_ATOM_FUNC_NR
        charge = cg_bead_charge[i_bead]
        mass   = cg_bead_mass[i_bead]
        c_id   = cg_chain_id[i_bead]
        s_name = cg_seg_name[i_bead]
        new_atom = GenTopAtom(i_bead, a_type, r_indx, r_name, a_name, f_type, charge, mass, c_id, s_name)
        push!(top_atoms, new_atom)
        push!(global_index_2_local_index, i_bead)
        push!(global_index_2_local_molid, c_id)
    end

    # ---------
    # [ bonds ]
    # ---------
    # AICG2+ bonds
    if ff_pro == FF_pro_AICG2p
        for bond in top_cg_pro_bonds
            new_bond = GenTopBond(bond.i, bond.j, AICG_BOND_FUNC_TYPE, bond.r0, AICG_BOND_K)
            push!(top_bonds, new_bond)
        end
    # Clementi Go bonds
    elseif ff_pro == FF_pro_Clementi_Go
        for bond in top_cg_pro_bonds
            new_bond = GenTopBond(bond.i, bond.j, CCGO_BOND_FUNC_TYPE, bond.r0, CCGO_BOND_K)
            push!(top_bonds, new_bond)
        end
    end

    # 3SPN.2C bonds
    if ff_dna == FF_DNA_3SPN2C && gen_3spn_itp > 0
        for bond in top_cg_DNA_bonds
            new_bond = GenTopBond(bond.i, bond.j, DNA3SPN_BOND_FUNC4_TYPE, bond.r0, DNA3SPN_BOND_K_2)
            push!(top_bonds, new_bond)
        end
    end

    # Structure-based RNA bonds
    if ff_rna == FF_RNA_HT
        for ( i_bond, bond ) in enumerate( top_cg_RNA_bonds )
            new_bond = GenTopBond(bond.i, bond.j, RNA_BOND_FUNC_TYPE, bond.r0, param_cg_RNA_k_bonds[i_bond])
            push!(top_bonds, new_bond)
        end
    end

    # ----------
    # [ angles ]
    # ----------
    # AICG2+ angles
    if ff_pro == FF_pro_AICG2p
        # AICG2+ 1-3
        if length(top_cg_pro_aicg13) > 0
            for ( i_13, a13 ) in enumerate( top_cg_pro_aicg13 )
                if  in(a13.i, AICG2p_flexible_local) || in(a13.j, AICG2p_flexible_local) || in(a13.k, AICG2p_flexible_local) ||
                    in(a13.i, HPS_IDR_region) || in(a13.j, HPS_IDR_region) || in(a13.k, HPS_IDR_region) ||
                    in(a13.i, KH_IDR_region) || in(a13.j, KH_IDR_region) || in(a13.k, KH_IDR_region)
                    continue
                end
                new_angle = GenTopAngle(a13.i, a13.j, a13.k, AICG_ANG_G_FUNC_TYPE, a13.a0, param_cg_pro_e_13[i_13], AICG_13_SIGMA)
                push!(top_angles, new_angle)
            end
        end
        # AICG2+ flexible
        if length(top_cg_pro_angles) > 0
            for ang in top_cg_pro_angles
                if  ( in(ang.i, HPS_IDR_region) && in(ang.j, HPS_IDR_region) && in(ang.k, HPS_IDR_region) ) || 
                    ( in(ang.i, KH_IDR_region) && in(ang.j, KH_IDR_region) && in(ang.k, KH_IDR_region) )
                    continue
                end
                new_angle = GenTopAngle(ang.i, ang.j, ang.k, AICG_ANG_F_FUNC_TYPE, 0.0, 0.0, 0.0)
                push!(top_angles, new_angle)
            end
        end
    # Clementi Go angles
    elseif ff_pro == FF_pro_Clementi_Go
        for ang in top_cg_pro_angles
            if  in(ang.i, HPS_IDR_region) || in(ang.j, HPS_IDR_region) || in(ang.k, HPS_IDR_region) || 
                in(ang.i, KH_IDR_region) || in(ang.j, KH_IDR_region) || in(ang.k, KH_IDR_region)
                continue
            end
            new_angle = GenTopAngle(ang.i, ang.j, ang.k, CCGO_ANG_FUNC_TYPE, ang.a0, CCGO_ANGL_K, 0.0)
            push!(top_angles, new_angle)
        end
    end

    # 3SPN.2C angles
    if ff_dna == FF_DNA_3SPN2C && gen_3spn_itp > 0
        for ( i_ang, ang ) in enumerate( top_cg_DNA_angles )
            new_angle = GenTopAngle(ang.i, ang.j, ang.k, DNA3SPN_ANG_FUNC_TYPE, ang.a0, param_cg_DNA_k_angles[i_ang], 0.0)
            push!(top_angles, new_angle)
        end
    end

    # RNA structure-based angles
    if ff_rna == FF_RNA_HT
        for ( i_ang, ang ) in enumerate( top_cg_RNA_angles )
            new_angle = GenTopAngle(ang.i, ang.j, ang.k, RNA_ANG_FUNC_TYPE, ang.a0, param_cg_RNA_k_angles[i_ang], 0.0)
            push!(top_angles, new_angle)
        end
    end

    # -------------
    # [ dihedrals ]
    # -------------
    function is_dihedral_dangerous(dih::CGTopDihedral)
        coor1 = cg_bead_coor[:, dih.i]
        coor2 = cg_bead_coor[:, dih.j]
        coor3 = cg_bead_coor[:, dih.k]
        coor4 = cg_bead_coor[:, dih.l]
        ang1 = compute_angle(coor1, coor2, coor3)
        ang2 = compute_angle(coor2, coor3, coor4)
        if ang1 > DIHEDRAL_SAFE_CUTOFF || ang2 > DIHEDRAL_SAFE_CUTOFF
            return true
        end
        return false
    end
    # AICG2+ dihedrals
    if ff_pro == FF_pro_AICG2p
        # AICG2+ Gaussian dihedrals
        for ( i_dih, dih ) in enumerate( top_cg_pro_aicg14 )
            if  in(dih.i, AICG2p_flexible_local) || in(dih.j, AICG2p_flexible_local) ||
                in(dih.k, AICG2p_flexible_local) || in(dih.l, AICG2p_flexible_local)
                continue
            elseif  in(dih.i, HPS_IDR_region) || in(dih.j, HPS_IDR_region) ||
                in(dih.k, HPS_IDR_region) || in(dih.l, HPS_IDR_region) || 
                in(dih.i, KH_IDR_region) || in(dih.j, KH_IDR_region) ||
                in(dih.k, KH_IDR_region) || in(dih.l, KH_IDR_region) 
                continue
            end
            # if is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_GAUS_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = AICG_DIH_G_FUNC_TYPE
            end
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          dih.t0, param_cg_pro_e_14[i_dih], AICG_14_SIGMA, 0)
            push!(top_dihedrals, new_dihedral)
        end
        # AICG2+ flexible dihedrals
        for dih in top_cg_pro_dihedrals
            if  ( in(dih.i, HPS_IDR_region) && in(dih.j, HPS_IDR_region) &&
                in(dih.k, HPS_IDR_region) && in(dih.l, HPS_IDR_region) ) ||
                ( in(dih.i, KH_IDR_region) && in(dih.j, KH_IDR_region) &&
                in(dih.k, KH_IDR_region) && in(dih.l, KH_IDR_region) )
                continue
            end
            dih_func_type = DIHEDRAL_TABU_MOD_TYPE[use_safe_dihedral]
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type, 0.0, 0.0, 0.0, 0)
            push!(top_dihedrals, new_dihedral)
        end
    # Clementi Go dihedral
    elseif ff_pro == FF_pro_Clementi_Go
        for dih in top_cg_pro_dihedrals
            if  in(dih.i, HPS_IDR_region) || in(dih.j, HPS_IDR_region) ||
                in(dih.k, HPS_IDR_region) || in(dih.l, HPS_IDR_region) ||
                in(dih.i, KH_IDR_region) || in(dih.j, KH_IDR_region) ||
                in(dih.k, KH_IDR_region) || in(dih.l, KH_IDR_region)
                continue
            end
            # if is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_PERI_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = CCGO_DIH_P_FUNC_TYPE
            end
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          dih.t0 - 180.0, CCGO_DIHE_K_1, 0.0, 1)
            push!(top_dihedrals, new_dihedral)
        end
        for dih in top_cg_pro_dihedrals
            if  in(dih.i, HPS_IDR_region) || in(dih.j, HPS_IDR_region) ||
                in(dih.k, HPS_IDR_region) || in(dih.l, HPS_IDR_region) ||
                in(dih.i, KH_IDR_region) || in(dih.j, KH_IDR_region) ||
                in(dih.k, KH_IDR_region) || in(dih.l, KH_IDR_region)
                continue
            end
            # if is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_PERI_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = CCGO_DIH_P_FUNC_TYPE
            end
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          3 * dih.t0 - 180.0, CCGO_DIHE_K_3, 0.0, 3)
            push!(top_dihedrals, new_dihedral)
        end
    end

    # 3SPN.2C dihedrals
    if ff_dna == FF_DNA_3SPN2C && gen_3spn_itp > 0
        # 3SPN.2C Gaussian dihedrals
        for dih in top_cg_DNA_dih_Gaussian
            # if use_safe_dihedral > 0 && is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_GAUS_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = DNA3SPN_DIH_G_FUNC_TYPE
            end
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          dih.t0, DNA3SPN_DIH_G_K, DNA3SPN_DIH_G_SIGMA, 0)
            push!(top_dihedrals, new_dihedral)
        end

        # 3SPN.2C Periodic dihedrals
        for dih in top_cg_DNA_dih_periodic
            # if use_safe_dihedral > 0 && is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_PERI_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = DNA3SPN_DIH_P_FUNC_TYPE
            end
            n_dih_tmp = DNA3SPN_DIH_P_FUNC_PERI
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          n_dih_tmp * dih.t0 - 180.0, DNA3SPN_DIH_P_K, 0.0, n_dih_tmp)
            push!(top_dihedrals, new_dihedral)
        end
    end

    # RNA structure-based Periodic dihedrals
    if ff_rna == FF_RNA_HT
        for ( i_dih, dih ) in enumerate( top_cg_RNA_dihedrals )
            # if is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_PERI_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = RNA_DIH_FUNC_TYPE
            end
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          dih.t0 - 180.0, param_cg_RNA_k_dihedrals[i_dih], 0.0, 1)
            push!(top_dihedrals, new_dihedral)
        end
        for ( i_dih, dih ) in enumerate( top_cg_RNA_dihedrals )
            # if is_dihedral_dangerous(dih)
            if true
                dih_func_type = DIHEDRAL_PERI_MOD_TYPE[use_safe_dihedral]
            else
                dih_func_type = RNA_DIH_FUNC_TYPE
            end
            new_dihedral = GenTopDihedral(dih.i, dih.j, dih.k, dih.l, dih_func_type,
                                          3 * dih.t0 - 180.0, param_cg_RNA_k_dihedrals[i_dih] / 2, 0.0, 3)
            push!(top_dihedrals, new_dihedral)
        end
    end


    # ---------
    # [ pairs ]
    # ---------
    # protein Go-type native contacts
    if ff_pro == FF_pro_AICG2p
        for i_res in 1:cg_num_particles
            for c in top_cg_pro_go_contact[i_res]
                if  in(i_res, AICG2p_flexible_nonlocal) || in(c[1], AICG2p_flexible_nonlocal) ||
                    in(i_res, HPS_IDR_region) || in(c[1], HPS_IDR_region) ||
                    in(i_res, KH_IDR_region) || in(c[1], KH_IDR_region)
                    continue
                end
                new_pair = GenTopPair(i_res, c[1], AICG_CONTACT_FUNC_TYPE, c[2], c[3])
                push!(top_pairs, new_pair)
            end
        end
    # Clementi Go native contacts
    elseif ff_pro == FF_pro_Clementi_Go
        for i_res in 1:cg_num_particles
            for c in top_cg_pro_go_contact[i_res]
                if  in(i_res, HPS_IDR_region) || in(c[1], HPS_IDR_region) ||
                    in(i_res, KH_IDR_region) || in(c[1], KH_IDR_region)
                    continue
                end
                new_pair = GenTopPair(i_res, c[1], CCGO_CONTACT_FUNC_TYPE, c[2], CCGO_NATIVE_EPSILON * ccgo_contact_scale)
                push!(top_pairs, new_pair)
            end
        end
    end

    # RNA HT-type native contacts
    if ff_rna == FF_RNA_HT
        for i_res in 1:cg_num_particles
            for c in top_cg_RNA_base_stack[i_res]
                new_pair = GenTopPair(i_res, c[1], RNA_CONTACT_FUNC_TYPE, c[2], c[3])
                push!(top_pairs, new_pair)
            end
            for c in top_cg_RNA_base_pair[i_res]
                new_pair = GenTopPair(i_res, c[1], RNA_CONTACT_FUNC_TYPE, c[2], c[3])
                push!(top_pairs, new_pair)
            end
            for c in top_cg_RNA_other_contact[i_res]
                new_pair = GenTopPair(i_res, c[1], RNA_CONTACT_FUNC_TYPE, c[2], c[3])
                push!(top_pairs, new_pair)
            end
        end
    end

    # protein-RNA native contacts
    if ff_pro_rna == FF_pro_RNA_Go
        for i_res in 1:cg_num_particles
            for c in top_cg_pro_RNA_contact[i_res]
                new_pair = GenTopPair(i_res, c[1], RNA_CONTACT_FUNC_TYPE, c[2], c[3])
                push!(top_pairs, new_pair)
            end
        end
    end


    # protein-DNA native contacts
    if ff_pro_dna == FF_pro_DNA_Go
        for i_res in 1:cg_num_particles
            for c in top_cg_pro_DNA_contact[i_res]
                new_pair = GenTopPair(i_res, c[1], CCGO_CONTACT_FUNC_TYPE, c[2], CCGO_NATIVE_EPSILON * ccgo_contact_scale)
                push!(top_pairs, new_pair)
            end
        end
    end


    # --------------
    # [ exclusions ]
    # --------------
    # contact pairs
    # for c in top_pairs
    #     i_exc = c.i
    #     j_exc = c.j
    #     new_exc = GenTopExclusion(i_exc, j_exc)
    #     push!(top_exclusions, new_exc)
    # end


    # ----------
    # [ pwmcos ]
    # ----------
    for p in top_cg_pro_DNA_pwmcos
        new_pwmcos = GenTopPWMcos(p.i, PWMCOS_FUNC_TYPE, p.r0, p.t1, p.t2, p.t3,
                                  p.eA, p.eC, p.eG, p.eT, pwmcos_gamma, pwmcos_epsil)
        push!(top_pwmcos, new_pwmcos)
    end

    for p in top_cg_pro_DNA_pwmcosns
        new_pwmcos = GenTopPWMcos(p.i, PWMCOS_FUNC_TYPE, p.r0, p.t1, p.t2, p.t3,
                                  p.eA, p.eC, p.eG, p.eT, pwmcos_gamma, pwmcos_epsil)
        push!(top_pwmcosns, new_pwmcos)
    end


    # ---------------------
    # [ cg_IDR_HPS_region ]
    # ---------------------
    if has_toml_mod
        if haskey(ff_detail_config, "IDR")
            if haskey(ff_detail_config["IDR"], "HPS_region")
                index_string = ff_detail_config["IDR"]["HPS_region"]
                hps_words = split(index_string, r"\s*,\s*", keepempty=false)
                for w in hps_words
                    idxwords = split(w, r"\s*to\s*", keepempty=false)
                    i = parse(Int, idxwords[1])
                    if length(idxwords) > 1
                        j = parse(Int, idxwords[2])
                    else
                        j = i
                    end
                    new_idr = GenTopRegion(i, j)
                    push!(top_idr_hps, new_idr)
                end
            end
        end
    end

    # ---------------------
    # [ cg_IDR_KH_region ]
    # ---------------------
    if has_toml_mod
        if haskey(ff_detail_config, "IDR")
            if haskey(ff_detail_config["IDR"], "KH_region")
                index_string = ff_detail_config["IDR"]["KH_region"]
                kh_words = split(index_string, r"\s*,\s*", keepempty=false)
                for w in kh_words
                    idxwords = split(w, r"\s*to\s*", keepempty=false)
                    i = parse(Int, idxwords[1])
                    if length(idxwords) > 1
                        j = parse(Int, idxwords[2])
                    else
                        j = i
                    end
                    new_idr = GenTopRegion(i, j)
                    push!(top_idr_kh, new_idr)
                end
            end
        end
    end


    mol_name = pdb_name[1:end-4]
    mytop = GenTopology(mol_name, cg_num_particles,
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
                        top_pwmcosns,
                        top_idr_hps,
                        top_idr_kh,
                        top_mol_list)
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
            @printf(log_file, " - Number of protein contacts:     %12d  \n", sum( length.(top_cg_pro_go_contact) ))
        end
        if num_chain_RNA > 0
            @printf(log_file, " - Number of RNA contacts:         %12d  \n", sum( length.(top_cg_RNA_base_stack) )
                    + sum( length.(top_cg_RNA_base_pair) ) + sum( length.(top_cg_RNA_other_contact) ) )
        end
        if num_chain_RNA > 0 && num_chain_pro > 0
            @printf(log_file, " - Number of protein-RNA contacts: %12d  \n", sum(length.(top_cg_pro_RNA_contact)) )
        end
        println(log_file, "================================================================================")

        close(log_file)
    end



    return ( mytop, myconf )

end
