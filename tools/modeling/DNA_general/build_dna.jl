#!/usr/bin/env julia

using Printf
using ArgParse

include("../../../src/lib/biomath.jl")

function read_DNA_standard_base(fname)
    pdb_atom_names = []
    pdb_coors = []
    for line in eachline(fname)
        if startswith(line, "ATOM")
            push!(pdb_atom_names, line[13:16])
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

function generate_NA_structure(args)

    lib_path = @__DIR__
    # ==============
    # read in params
    # ==============
    par_fname = get(args, "param", "")
    pdb_fname = get(args, "output", "_DNA_constructed_.pdb")
    do_debug  = get(args, "debug", false)

    bp_names = []
    bp_parms = []
    bs_parms = []
    if length(par_fname) == 0
        println("\e[0;31;40m ERROR: Please specify the DNA structural parameter file. \e[0m")
        exit()
    end
    for line in eachline(par_fname)
        if line[1] in "ACGT"
            words = split(line)
            local_bp_parms = [parse(Float64, words[i]) for i in 2:7]
            local_bs_parms = [parse(Float64, words[i]) for i in 8:13]
            push!(bp_names, words[1])
            push!(bp_parms, local_bp_parms)
            push!(bs_parms, local_bs_parms)
        end
    end

    # ------------------
    # read in structures
    # ------------------
    std_base_A = read_DNA_standard_base(lib_path * "/lib/Atomic_A.pdb")
    std_base_C = read_DNA_standard_base(lib_path * "/lib/Atomic_C.pdb")
    std_base_G = read_DNA_standard_base(lib_path * "/lib/Atomic_G.pdb")
    std_base_T = read_DNA_standard_base(lib_path * "/lib/Atomic_T.pdb")
    map_base_atom_coors = Dict('A' => std_base_A[2],
                               'C' => std_base_C[2],
                               'G' => std_base_G[2],
                               'T' => std_base_T[2])
    map_base_atom_names = Dict('A' => std_base_A[1],
                               'C' => std_base_C[1],
                               'G' => std_base_G[1],
                               'T' => std_base_T[1])


    ###########################################################################
    #                          Build base-pair frames                         #
    ###########################################################################
    if do_debug
        of_bp_frames = open("DEBUG_basepair_reference_frames.dat", "w")
        of_base_frames = open("DEBUG_base_reference_frames.dat", "w")
    end

    num_bp = length(bp_parms)
    if do_debug
        @printf(of_bp_frames, "%5d base-pairs\n", num_bp)
    end

    bp_frame_orign = zeros(Float64, (3, num_bp))
    bp_frame_basis = zeros(Float64, (3, 3, num_bp))

    # ----------------------
    # set the frame for BP 1
    # ----------------------
    bp_frame_basis[:, :, 1] = diagm([1., 1., 1.])

    # --------------------
    # propagate to all bps
    # --------------------
    for ibp in 2:num_bp
        # parameters
        tmp_shift = bs_parms[ibp][1]
        tmp_slide = bs_parms[ibp][2]
        tmp_rise  = bs_parms[ibp][3]
        tmp_tilt  = bs_parms[ibp][4]
        tmp_roll  = bs_parms[ibp][5]
        tmp_twist = bs_parms[ibp][6]
        # angles
        tmp_RollTilt = sqrt(tmp_tilt * tmp_tilt + tmp_roll * tmp_roll)
        tmp_phi_sign = tmp_tilt > 0 ? 1 : -1
        tmp_phi      = acosd(tmp_roll / tmp_RollTilt) * tmp_phi_sign

        # -----------------------------------------------
        # calculate frames with respect to the current bp
        # -----------------------------------------------
        # basis of middle frame
        local_middle_basis = rotation_matrix_around_z(tmp_twist/2 - tmp_phi) * rotation_matrix_around_y(tmp_RollTilt/2) * rotation_matrix_around_z(tmp_phi)
        # basis of next bp frame
        local_bp2_basis    = rotation_matrix_around_z(tmp_twist/2 - tmp_phi) * rotation_matrix_around_y(tmp_RollTilt)   * rotation_matrix_around_z(tmp_twist/2 + tmp_phi)
        # origin of next bp
        local_bp2_orign    = local_middle_basis * [tmp_shift, tmp_slide, tmp_rise]

        # ----------------------------------
        # remap frames to global coordinates
        # ----------------------------------
        global_bp1_orign = bp_frame_orign[:, ibp - 1]
        global_bp1_basis = bp_frame_basis[:, :, ibp - 1]
        global_bp2_orign = global_bp1_orign + global_bp1_basis * local_bp2_orign
        global_bp2_basis = global_bp1_basis * local_bp2_basis

        bp_frame_orign[:, ibp] = global_bp2_orign
        bp_frame_basis[:, :, ibp] = global_bp2_basis
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


    ###########################################################################
    #                           Build DNA structures                          #
    ###########################################################################
    # --------------------------------
    # count number of atoms in strands
    # --------------------------------
    # strand A and B
    num_atom_strand_A = 0
    num_atom_strand_B = 0
    for ibp in 1:num_bp
        num_atom_strand_A += length(map_base_atom_names[bp_names[ibp][1]])
        num_atom_strand_B += length(map_base_atom_names[bp_names[ibp][3]])
    end
    coors_strand_A = zeros(Float64, (3, num_atom_strand_A))
    coors_strand_B = zeros(Float64, (3, num_atom_strand_B))
    aname_strand_A = ["" for j in 1:num_atom_strand_A] # atom names
    aname_strand_B = ["" for j in 1:num_atom_strand_B] # atom names
    rname_strand_A = ["" for j in 1:num_atom_strand_A] # resid names
    rname_strand_B = ["" for j in 1:num_atom_strand_B] # resid names
    resid_strand_A = [0 for j in 1:num_atom_strand_A]  # resid index
    resid_strand_B = [0 for j in 1:num_atom_strand_B]  # resid index

    i_start_strand_A = 1
    i_end_strand_B   = num_atom_strand_B
    for ibp in 1:num_bp
        # =======================
        # Build local base frames
        # =======================
        # parameters
        tmp_shear    = bp_parms[ibp][1]
        tmp_stretch  = bp_parms[ibp][2]
        tmp_stagger  = bp_parms[ibp][3]
        tmp_buckle   = bp_parms[ibp][4]
        tmp_twist    = bp_parms[ibp][5]
        tmp_opening  = bp_parms[ibp][6]
        # angles
        if true                 # give the correct output as 3DNA v2.4
            tmp_BucTwist = sqrt(tmp_buckle * tmp_buckle + tmp_twist * tmp_twist)
            tmp_phi_sign = tmp_buckle > 0 ? 1 : -1
            tmp_phi      = acosd(tmp_twist / tmp_BucTwist) * tmp_phi_sign
        else                    # OLD STRATEGY
            # tmp_BuckOpen = sqrt(tmp_buckle * tmp_buckle + tmp_opening * tmp_opening)
            # tmp_phi_sign = tmp_opening > 0 ? 1 : -1
            # tmp_phi      = acosd(tmp_buckle / tmp_BuckOpen) * tmp_phi_sign
        end

        # -----------------------------------------------
        # calculate frames with respect to the current bp
        # -----------------------------------------------
        if true                 # give the same structure as 3DNA v2.4
            # basis of base 1
            local_base1_basis = rotation_matrix_around_z(-tmp_phi) * rotation_matrix_around_y(+tmp_BucTwist / 2) * rotation_matrix_around_z(+tmp_phi + tmp_opening / 2)
            # basis of base 2
            local_base2_basis = rotation_matrix_around_z(-tmp_phi) * rotation_matrix_around_y(-tmp_BucTwist / 2) * rotation_matrix_around_z(+tmp_phi - tmp_opening / 2)
        else                    # OLD STRATEGY
            # local_base1_basis = rotation_matrix_around_y(-tmp_phi) * rotation_matrix_around_x( tmp_BuckOpen / 2) * rotation_matrix_around_y(tmp_phi + tmp_twist / 2)
            # local_base2_basis = rotation_matrix_around_y(-tmp_phi) * rotation_matrix_around_x(-tmp_BuckOpen / 2) * rotation_matrix_around_y(tmp_phi - tmp_twist / 2)
        end
        # rotate base 2 around x-axis by 180
        local_base2_basis = local_base2_basis * diagm([1, -1, -1])

        # origin of next bp
        global_bp_orign = bp_frame_orign[:, ibp]
        global_bp_basis = bp_frame_basis[:, :, ibp]
        local_base1_orign = global_bp_basis * [tmp_shear, tmp_stretch, tmp_stagger] ./ 2

        # ----------------------------------
        # remap frames to global coordinates
        # ----------------------------------
        global_base1_orign = global_bp_orign + local_base1_orign
        global_base1_basis = global_bp_basis * local_base1_basis
        global_base2_orign = global_bp_orign - local_base1_orign
        global_base2_basis = global_bp_basis * local_base2_basis

        if do_debug
            @printf(of_base_frames, "=============================================\n")
            @printf(of_base_frames, "... %5d %s ...\n", ibp, bp_names[ibp])
            @printf(of_base_frames, "---------------------------------------------\n")
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base1_orign[1], global_base1_orign[2], global_base1_orign[3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base1_basis[1,1], global_base1_basis[1,2], global_base1_basis[1,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base1_basis[2,1], global_base1_basis[2,2], global_base1_basis[2,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base1_basis[3,1], global_base1_basis[3,2], global_base1_basis[3,3])
            @printf(of_base_frames, "---------------------------------------------\n")
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base2_orign[1], global_base2_orign[2], global_base2_orign[3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base2_basis[1,1], global_base2_basis[1,2], global_base2_basis[1,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base2_basis[2,1], global_base2_basis[2,2], global_base2_basis[2,3])
            @printf(of_base_frames, "%10.4f%10.4f%10.4f\n", global_base2_basis[3,1], global_base2_basis[3,2], global_base2_basis[3,3])
        end

        # =======================
        # fill in the coordinates
        # =======================
        b_name_1 = bp_names[ibp][1]
        b_name_2 = bp_names[ibp][3]
        coor_std_base1 = map_base_atom_coors[b_name_1]
        coor_std_base2 = map_base_atom_coors[b_name_2]

        coor_new_base1 = global_base1_basis * coor_std_base1 .+ global_base1_orign
        coor_new_base2 = global_base2_basis * coor_std_base2 .+ global_base2_orign

        num_atom_base1 = length(map_base_atom_names[b_name_1])
        num_atom_base2 = length(map_base_atom_names[b_name_2])

        coors_strand_A[:, i_start_strand_A:i_start_strand_A + num_atom_base1 - 1] = coor_new_base1[:, :]
        coors_strand_B[:, i_end_strand_B - num_atom_base2 + 1:i_end_strand_B]     = coor_new_base2[:, :]
        aname_strand_A[i_start_strand_A:i_start_strand_A + num_atom_base1 - 1] = map_base_atom_names[b_name_1]
        aname_strand_B[i_end_strand_B - num_atom_base2 + 1:i_end_strand_B]     = map_base_atom_names[b_name_2]
        rname_strand_A[i_start_strand_A:i_start_strand_A + num_atom_base1 - 1] .= " D" * b_name_1
        rname_strand_B[i_end_strand_B - num_atom_base2 + 1:i_end_strand_B]     .= " D" * b_name_2
        resid_strand_A[i_start_strand_A:i_start_strand_A + num_atom_base1 - 1] .= ibp
        resid_strand_B[i_end_strand_B - num_atom_base2 + 1:i_end_strand_B]     .= 2 * num_bp - ibp + 1

        i_start_strand_A += num_atom_base1
        i_end_strand_B   -= num_atom_base2
    end

    # =============
    # output to PDB
    # =============
    new_PDB_file = open(pdb_fname, "w")
    # --------
    # strand A
    # --------
    for iatm in 1:num_atom_strand_A
        atm_name = aname_strand_A[iatm]
        x = coors_strand_A[1, iatm]
        y = coors_strand_A[2, iatm]
        z = coors_strand_A[3, iatm]
        res_name = rpad(rname_strand_A[iatm], 4)
        res_id   = resid_strand_A[iatm]
        @printf(new_PDB_file, "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s \n",
                iatm, atm_name, res_name,
                "A",
                res_id,
                x, y, z,
                1.0,
                1.0,
                "DNA1",
                strip(atm_name)[1])
    end
    println(new_PDB_file, "TER")

    # --------
    # strand B
    # --------
    for iatm in 1:num_atom_strand_B
        atm_name = aname_strand_B[iatm]
        x = coors_strand_B[1, iatm]
        y = coors_strand_B[2, iatm]
        z = coors_strand_B[3, iatm]
        res_name = rpad(rname_strand_B[iatm], 4)
        res_id   = resid_strand_B[iatm]
        @printf(new_PDB_file, "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s \n",
                iatm + num_atom_strand_A, atm_name, res_name,
                "B",
                res_id,
                x, y, z,
                1.0,
                1.0,
                "DNA2",
                strip(atm_name)[1])
    end
    println(new_PDB_file, "END")

    # ==============
    # ending work...
    # ==============
    close(new_PDB_file)

    if do_debug
        close(of_bp_frames)
        close(of_base_frames)
    end

end

if abspath(PROGRAM_FILE) == @__FILE__

    args = ArgParseSettings()

    @add_arg_table args begin
        "--param", "-p"
        help     = "Parameter file of base-pair and base-step local structures."
        arg_type = String
        default  = ""

        "--output", "-o"
        help     = "PDB file as output."
        arg_type = String
        default  = "_DNA_constructed_.pdb"

        "--debug"
        help = "DEBUG."
        action = :store_true
    end

    generate_NA_structure(parse_args(args))
end
