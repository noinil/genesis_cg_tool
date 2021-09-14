#!/usr/bin/env julia

using Printf
using ArgParse

include("../../../src/lib/gcj.jl")

function read_DNA_standard_base(fname, base_type)
    if base_type == "purine"
        # purine: 9 atoms
        list_key_atoms = [" N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N1 ", " C2 ", " N3 ", " C4 "]
    elseif base_type == "pyrimidine"
        # pyrimidine: 6 atoms
        list_key_atoms = [" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "]
    end
    pdb_atom_names = []
    pdb_coors = []
    for line in eachline(fname)
        if startswith(line, "ATOM")
            atom_name = line[13:16]
            if !(atom_name in list_key_atoms)
                continue
            end
            push!(pdb_atom_names, strip(atom_name))
            x = parse(Float64, line[31:38])
            y = parse(Float64, line[39:46])
            z = parse(Float64, line[47:54])
            push!(pdb_coors, [x, y, z])
        end
    end
    coor_array = zeros(Float64, (3, length(pdb_coors)))
    for j in 1:length(pdb_coors)
        coor_array[:, j] = pdb_coors[j][:]
    end
    return [pdb_atom_names, coor_array]
end

function analyze_NA_structure(args)

    lib_path = @__DIR__

    # ----
    # args
    # ----
    pdb_name = get(args, "PDB", "")
    out_pref = get(args, "output", "_DNA_analysis_")
    is_mmCIF = get(args, "mmCIF", false)
    do_debug = get(args, "debug", false)

    # ======================
    # read in standard bases
    # ======================
    std_base_A = read_DNA_standard_base(lib_path * "/lib/aa_A.pdb", "purine")
    std_base_C = read_DNA_standard_base(lib_path * "/lib/aa_C.pdb", "pyrimidine")
    std_base_G = read_DNA_standard_base(lib_path * "/lib/aa_G.pdb", "purine")
    std_base_T = read_DNA_standard_base(lib_path * "/lib/aa_T.pdb", "pyrimidine")
    map_base_atom_coors = Dict("DA" => std_base_A[2],
                               "DC" => std_base_C[2],
                               "DG" => std_base_G[2],
                               "DT" => std_base_T[2])
    map_base_atom_names = Dict("DA" => std_base_A[1],
                               "DC" => std_base_C[1],
                               "DG" => std_base_G[1],
                               "DT" => std_base_T[1])

    # ====================
    # prepare output files
    # ====================
    of_base_fname  = out_pref * "_base_reference_frames.dat"
    of_bp_fname    = out_pref * "_basepair_reference_frames.dat"
    of_parm_fname  = out_pref * "_structure_parameters.dat"
    if do_debug
        of_base_frames = open(of_base_fname, "w")
        of_bp_frames   = open(of_bp_fname, "w")
    end
    of_parameters  = open(of_parm_fname, "w")


    ###########################################################################
    #                               Analyze PDB                               #
    ###########################################################################
    if is_mmCIF
        cif_data = read_mmCIF(pdb_name)
        aa_molecule = mmCIF_to_AAMolecule(cif_data)
    else
        aa_molecule = read_PDB(pdb_name)
    end

    # get number of bp
    num_residues_chain_A = length(aa_molecule.chains[1].residues)
    num_residues_chain_B = length(aa_molecule.chains[2].residues)
    if num_residues_chain_A != num_residues_chain_B
        println("Inconsistent number of residues in chains 1 and 2.")
        exit()
    else
        num_bp = num_residues_chain_A
    end

    # ---------------
    # data structures
    # ---------------
    bp_frame_orign = zeros(Float64, (3, num_bp))
    bp_frame_basis = zeros(Float64, (3, 3, num_bp))
    bp_names = ["" for i in 1:num_bp]

    bp_params = zeros(Float64, (6, num_bp)) # base-pair quantities
    bs_params = zeros(Float64, (6, num_bp)) # base-step quantities

    # =============================
    # get base and base-pair frames
    # =============================
    for ibp in 1:num_bp
        i_base1 = aa_molecule.chains[1].residues[ibp]
        i_base2 = aa_molecule.chains[2].residues[num_bp - ibp + 1]
        resname_base1 = aa_molecule.residues[i_base1].name
        resname_base2 = aa_molecule.residues[i_base2].name
        bp_names[ibp] = resname_base1[end] * "-" * resname_base2[end]
        # -------------------
        # get coors of base 1
        # -------------------
        tmp_dict = map_base_atom_names[resname_base1]
        num_atom_base1 = length(tmp_dict)
        pdb_coor_base1 = zeros(Float64, (3, num_atom_base1))
        for iatm in aa_molecule.residues[i_base1].atoms
            atom_name = aa_molecule.atom_names[iatm]
            if atom_name in tmp_dict
                ii = findfirst(x->x==atom_name, tmp_dict)
                pdb_coor_base1[:, ii] = aa_molecule.atom_coors[:, iatm]
            end
        end
        # -------------------
        # get coors of base 2
        # -------------------
        tmp_dict = map_base_atom_names[resname_base2]
        num_atom_base2 = length(tmp_dict)
        pdb_coor_base2 = zeros(Float64, (3, num_atom_base2))
        for iatm in aa_molecule.residues[i_base2].atoms
            atom_name = aa_molecule.atom_names[iatm]
            if atom_name in tmp_dict
                ii = findfirst(x->x==atom_name, tmp_dict)
                pdb_coor_base2[:, ii] = aa_molecule.atom_coors[:, iatm]
            end
        end

        # ---------------
        # superpose bases
        # ---------------
        base1_fit = compute_superimposition_transformation(map_base_atom_coors[resname_base1], pdb_coor_base1)
        base2_fit = compute_superimposition_transformation(map_base_atom_coors[resname_base2], pdb_coor_base2)

        # -------------------------
        # get base reference frames
        # -------------------------
        base1_orign = base1_fit.translation
        base1_basis = base1_fit.rotation
        base2_orign = base2_fit.translation
        base2_basis = base2_fit.rotation
        if base1_basis[:, 3]' * base2_basis[:, 3] < 0
            base2_basis *= Diagonal([1, -1, -1]);
        end

        if do_debug
            @printf(of_base_frames, "=============================================\n")
            @printf(of_base_frames, "... %5d %s ...\n", ibp, bp_names[ibp])
            @printf(of_base_frames, "---------------------------------------------\n")
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base1_orign[1], base1_orign[2], base1_orign[3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base1_basis[1,1], base1_basis[1,2], base1_basis[1,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base1_basis[2,1], base1_basis[2,2], base1_basis[2,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base1_basis[3,1], base1_basis[3,2], base1_basis[3,3])
            @printf(of_base_frames, "---------------------------------------------\n")
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base2_orign[1], base2_orign[2], base2_orign[3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base2_basis[1,1], base2_basis[1,2], base2_basis[1,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base2_basis[2,1], base2_basis[2,2], base2_basis[2,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", base2_basis[3,1], base2_basis[3,2], base2_basis[3,3])
        end

        # ===================================
        # calculate base-pair reference frame
        # ===================================
        # base-pair frame origin!!!
        bp_frame_orign[:, ibp] = ( base1_orign .+ base2_orign ) .* 0.5

        # base-pair frame basis!!!
        if true                 # give the correct result as latest 3DNA
            base1_z = base1_basis[:, 3]
            base2_z = base2_basis[:, 3]
            bp_hinge = normalize(cross(base2_z, base1_z))
            bp_BucTwist_angle = acosd(base2_z' * base1_z)
            base1_basis_new = rotation_matrix_around_axis(bp_hinge, -bp_BucTwist_angle/2) * base1_basis
            base2_basis_new = rotation_matrix_around_axis(bp_hinge, +bp_BucTwist_angle/2) * base2_basis
        else                    # OLD STATEGY: align y-axis
            # base1_y = base1_basis[:, 2]
            # base2_y = base2_basis[:, 2]
            # bp_hinge = normalize(cross(base2_y, base1_y))
            # bp_BucOpen_angle = acosd(base2_y' * base1_y)
            # base1_basis_new = rotation_matrix_around_axis(bp_hinge, -bp_BucOpen_angle/2) * base1_basis
            # base2_basis_new = rotation_matrix_around_axis(bp_hinge, +bp_BucOpen_angle/2) * base2_basis
        end

        bp_x = normalize(base1_basis_new[:, 1] + base2_basis_new[:, 1])
        bp_y = normalize(base1_basis_new[:, 2] + base2_basis_new[:, 2])
        bp_z = normalize(base1_basis_new[:, 3] + base2_basis_new[:, 3])

        bp_frame_basis[:, :, ibp] = hcat(bp_x, bp_y, bp_z)

        # ----------------------------
        # analyze base pair quantities
        # ----------------------------
        if true                 # gives correct results as 3DNA...
            # translation
            (tmp_shear, tmp_stretch, tmp_stagger) = (base1_orign - base2_orign)' * bp_frame_basis[:, :, ibp]
            # opening
            opening_sign = cross(base2_basis_new[:, 2], base1_basis_new[:, 2])' * bp_z > 0 ? 1 : -1
            tmp_opening = acosd(base2_basis_new[:, 1]' * base1_basis_new[:, 1]) * opening_sign
            # phase
            phase_sign = cross(bp_hinge, bp_y)' * bp_z > 0 ? 1 : -1
            tmp_phase  = acosd(bp_y' * bp_hinge) * phase_sign
            # twist
            tmp_twist  = bp_BucTwist_angle * cosd(tmp_phase)
            # buckle
            tmp_buckle = bp_BucTwist_angle * sind(tmp_phase)
        else                    # OLD STRATEGY...
            # twist_sign = cross(base2_basis_new[:, 1], base1_basis_new[:, 1])'
            # * bp_y > 0 ? 1 : -1
            # tmp_twist = acosd(base2_basis_new[:, 1]' * base1_basis_new[:, 1])
            # * twist_sign
            # phase_sign = cross(bp_hinge, bp_x)' * bp_y > 0 ? 1 : -1
            # tmp_phase  = acosd(bp_x' * bp_hinge) * phase_sign
            # tmp_buckle  = bp_BucOpen_angle * cosd(tmp_phase)
            # tmp_opening = bp_BucOpen_angle * sind(tmp_phase)
        end
        # fill to array of params
        bp_params[:, ibp] = [tmp_shear, tmp_stretch, tmp_stagger, tmp_buckle, tmp_twist, tmp_opening]
    end

    # ====================
    # get base-step params
    # ====================
    for ibp in 2:num_bp
        bp1_orign = bp_frame_orign[:, ibp - 1]
        bp1_basis = bp_frame_basis[:, :, ibp - 1]
        bp2_orign = bp_frame_orign[:, ibp]
        bp2_basis = bp_frame_basis[:, :, ibp]

        # bs (base-step) origin is the average of bp1 and bp2
        bs_orign = ( bp1_orign + bp2_orign ) .* 0.5

        # bp orientation is a bit complicated...
        bp1_z = bp1_basis[:, 3]
        bp2_z = bp2_basis[:, 3]
        bs_hinge = normalize(cross(bp1_z, bp2_z))
        bs_RollTilt_angle = acosd(bp1_z' * bp2_z)
        bp1_basis_new = rotation_matrix_around_axis(bs_hinge, +bs_RollTilt_angle/2) * bp1_basis
        bp2_basis_new = rotation_matrix_around_axis(bs_hinge, -bs_RollTilt_angle/2) * bp2_basis
        bs_x = normalize(bp1_basis_new[:, 1] + bp2_basis_new[:, 1])
        bs_y = normalize(bp1_basis_new[:, 2] + bp2_basis_new[:, 2])
        bs_z = normalize(bp1_basis_new[:, 3] + bp2_basis_new[:, 3])

        bs_frame_basis = hcat(bs_x, bs_y, bs_z)

        # ----------------------------
        # analyze base-step quantities
        # ----------------------------
        # translation
        (tmp_shift, tmp_slide, tmp_rise) = (bp2_orign - bp1_orign)' * bs_frame_basis
        # twist
        twist_sign = cross(bp1_basis_new[:, 2], bp2_basis_new[:, 2])' * bs_z > 0 ? 1 : -1
        tmp_twist  = acosd(bp1_basis_new[:, 1]' * bp2_basis_new[:, 1]) * twist_sign
        # phase
        phase_sign = cross(bs_hinge, bs_y)' * bs_z > 0 ? 1 : -1
        tmp_phase  = acosd(bs_y' * bs_hinge) * phase_sign
        # roll
        tmp_roll   = bs_RollTilt_angle * cosd(tmp_phase)
        # tilt
        tmp_tilt   = bs_RollTilt_angle * sind(tmp_phase)
        # fill to array of params
        bs_params[:, ibp] = [tmp_shift, tmp_slide, tmp_rise, tmp_tilt, tmp_roll, tmp_twist]
    end


    # ---------------------------
    # output the reference frames
    # ---------------------------
    if do_debug
        for ibp in 1:num_bp
            @printf(of_bp_frames, "... %5d %s ...\n", ibp, bp_names[ibp])
            @printf(of_bp_frames, "%10.4f%10.4f%10.4f\n", bp_frame_orign[1,ibp], bp_frame_orign[2,ibp], bp_frame_orign[3,ibp])
            @printf(of_bp_frames, "%10.4f%10.4f%10.4f\n", bp_frame_basis[1,1,ibp], bp_frame_basis[1,2,ibp], bp_frame_basis[1,3,ibp])
            @printf(of_bp_frames, "%10.4f%10.4f%10.4f\n", bp_frame_basis[2,1,ibp], bp_frame_basis[2,2,ibp], bp_frame_basis[2,3,ibp])
            @printf(of_bp_frames, "%10.4f%10.4f%10.4f\n", bp_frame_basis[3,1,ibp], bp_frame_basis[3,2,ibp], bp_frame_basis[3,3,ibp])
        end
    end


    # -------------------------------
    # output the structure quantities
    # -------------------------------
    @printf(of_parameters, "%3s %10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
            "#  ",
            "Shear", "Stretch", "Stagger", "Buckle", "Prop-Tw", "Opening",
            "Shift", "Slide", "Rise", "Tilt", "Roll", "Twist")
    for ibp in 1:num_bp
        bp_shear, bp_stretch, bp_stagger, bp_buckle, bp_twist, bp_opening = bp_params[:, ibp]
        bs_shift, bs_slide, bs_rise, bs_tilt, bs_roll, bs_twist = bs_params[:, ibp]
        @printf(of_parameters, "%3s %10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
                bp_names[ibp],
                bp_shear, bp_stretch, bp_stagger, bp_buckle, bp_twist, bp_opening,
                bs_shift, bs_slide, bs_rise, bs_tilt, bs_roll, bs_twist)
    end


    # ==============
    # ending work...
    # ==============
    if do_debug
        close(of_base_frames)
        close(of_bp_frames)
    end
    close(of_parameters)

end

if abspath(PROGRAM_FILE) == @__FILE__
    args = ArgParseSettings()

    @add_arg_table args begin

        "PDB"
        help     = "PDB file of DNA structures."
        required = true
        arg_type = String

        "--mmCIF"
        help     = "PDB file of DNA structures."
        action   = :store_true

        "--output", "-o"
        help     = "Output name prefix."
        arg_type = String
        default  = "_DNA_analysis_"

        "--debug"
        help     = "DEBUG."
        action   = :store_true
    end

    analyze_NA_structure(parse_args(args))
end
