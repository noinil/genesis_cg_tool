###############################################################################
#                                 _           _                               #
#                              __| |  ___  __| |                              #
#                             / _` | / __|/ _` |                              #
#                            | (_| || (__| (_| |                              #
#                             \__,_| \___|\__,_|                              #
#                                                                             #
###############################################################################

using Printf

struct DCD_trajectory
    traj_type::String
    traj_frames::Int
    traj_first_step::Int
    traj_output_interval::Int
    traj_steps::Int
    boundary_type::Int
    traj_format::Int
    md_doc::Vector{String}
    num_atoms::Int
    boundary_box_size::Array{Float64, 2}
    conformations::Vector{Conformation}
end


function read_dcd(dcd_filename::String, args::Dict{String, <:Any}=Dict{String, Any}())

    verbose  = get(args, "verbose", false)
    i_frame_begin = get(args, "begin", 1)
    i_frame_end   = get(args, "end", -1)
    i_frame_step  = get(args, "step", 1)

    dcd_file = open(dcd_filename, "r")

    if verbose
        println("============================================================")
        println("> Open DCD file: ", dcd_filename)
    end

    # ========================
    # Read in head information
    # ========================

    # DCD first block
    block_size_0 = read(dcd_file, Int32)
    if block_size_0 != 84
        error("ERROR: wrong DCD format!")
    end

    # "CORD" or "VELD"
    tmp_char_array = Array{Char}(undef, 0)
    for i = 1:4
        c = read(dcd_file, Char)
        push!(tmp_char_array, c)
    end
    file_type = String(tmp_char_array)

    # simulation info
    tmp_int_array = Array{Int32}(undef, 0)
    for i = 1:20
        m = read(dcd_file, Int32)
        push!(tmp_int_array, m)
    end
    n_frames      = tmp_int_array[1]
    n_ts_first    = tmp_int_array[2]
    n_ts_interval = tmp_int_array[3]
    n_ts_all      = tmp_int_array[4]
    bc_flag       = tmp_int_array[11]
    traj_format   = tmp_int_array[20]

    bc_type       = bc_flag == 0 ? "No boundary" : "Periodic boundary"

    block_size_1 = read(dcd_file, Int32)
    if block_size_1 != 84
        error("ERROR: wrong DCD format!")
    end

    # ======================
    # Read in MD system info
    # ======================

    block_size_0 = read(dcd_file, Int32)

    dcd_doc_lines = Vector{String}(undef, 0)
    n_doc_line = read(dcd_file, Int32)
    for i in 1:n_doc_line
        tmp_char_array = Array{Char}(undef, 0)
        for j = 1:80
            c = read(dcd_file, Char)
            push!(tmp_char_array, c)
        end
        doc_line = String(tmp_char_array)
        push!(dcd_doc_lines, doc_line)
    end

    block_size_1 = read(dcd_file, Int32)

    # =======================
    # Read in particle number
    # =======================

    block_size_0 = read(dcd_file, Int32)

    n_particles  = read(dcd_file, Int32)

    block_size_1 = read(dcd_file, Int32)

    # ========================
    # Brief output of DCD info
    # ========================

    if verbose
        println("> DCD information:")
        println("  > File type           : $(file_type)")
        println("  --------------------------------------------------")
        println("  > Number of particles : $(n_particles)")
        println("  > Boundary codition   : $(bc_type)")
        println("  --------------------------------------------------")
        println("  > Number of frames    : $(n_frames)")
        println("    > First MD step       : $(n_ts_first)")
        println("    > Total MD steps      : $(n_ts_all)")
        println("    > Output interval     : $(n_ts_interval)")
        println("  > Trajectory format     : $(traj_format)")

        println("============================================================")
        println("> DCD information:")
    end

    if verbose
        for line in dcd_doc_lines
            println("  > |", line, "|")
        end
    end

    # ======================
    # Read in Coordinates!!!
    # ======================
    data_block_size = n_particles * 4 # Single precision (4-byte), Float32

    data_frame_size = (data_block_size + 8) * 3
    if bc_flag == 1
        data_frame_size += 8 * 6 + 8
    end

    if i_frame_end < 1
        i_frame_end = n_frames
    end
    i_frame_read_indices = [i_frame_begin:i_frame_step:i_frame_end...]

    boundary_conditions = zeros(Float64, 6, length(i_frame_read_indices))
    conformations = Vector{Conformation}(undef, 0)
    is_broken_trajectory = false
    i_read_frame = 0
    for t in 1 : n_frames
        if !(t in i_frame_read_indices)
            skip(dcd_file, data_frame_size)
            continue
        end
        if verbose
            println(" ~~~~> Reading Frame: ", t)
        end
        i_read_frame += 1
        # ----------------
        # Read in box info
        # ----------------
        if bc_flag == 1
            block_size_0 = read(dcd_file, Int32)
            bc_size = zeros(Float64, 6)
            for i in 1 : 6
                l = read(dcd_file, Float64)
                bc_size[i] = l
            end
            block_size_1 = read(dcd_file, Int32)
            # store box information
            boundary_conditions[:, i_read_frame] = bc_size
        end

        # ------------------
        # Read in coors/vels
        # ------------------
        coors = zeros(Float64, (3, n_particles))
        for dim = 1 : 3
            block_size_0 = read(dcd_file, Int32)
            for i  in 1 : n_particles
                x = read(dcd_file, Float32)
                coors[dim, i] = x
            end
            block_size_1 = read(dcd_file, Int32)
            if block_size_0 != data_block_size || block_size_1 != data_block_size
                println("[1;31m WARNING: wrong block size in coordinate reading![0m")
                println("         Incomplete trajectory?")
                is_broken_trajectory = true
                break
            end
        end
        if is_broken_trajectory
            break
        end
        push!(conformations, Conformation(n_particles, coors))
        if eof(dcd_file)
            break
        end
    end

    if verbose
        println("------------------------------------------------------------")
        println(">[1;32m FINISH[0m reading the DCD file. Have fun!")
        println("============================================================")
        println("         ")
    end

    close(dcd_file)

    new_trajectory = DCD_trajectory(file_type,
                                    i_read_frame,
                                    n_ts_first,
                                    n_ts_interval,
                                    n_ts_all,
                                    bc_flag,
                                    traj_format,
                                    dcd_doc_lines,
                                    n_particles,
                                    boundary_conditions,
                                    conformations)

    return new_trajectory

end


function write_dcd(dcd_trajectory::DCD_trajectory, dcd_filename::AbstractString)

    dcd_file = open(dcd_filename, "w")

    n_atom = dcd_trajectory.num_atoms

    # ======================
    # Write head information
    # ======================
    #
    # block size
    write(dcd_file, Int32(84))

    # "CORD" or "VELD"
    for i = 1:4
        write(dcd_file, dcd_trajectory.traj_type[i])
    end

    # simulation info
    tmp_int_array = zeros(Int32, 20)
    tmp_int_array[1]  = dcd_trajectory.traj_frames
    tmp_int_array[2]  = dcd_trajectory.traj_first_step
    tmp_int_array[3]  = dcd_trajectory.traj_output_interval
    tmp_int_array[4]  = dcd_trajectory.traj_steps
    tmp_int_array[11] = dcd_trajectory.boundary_type
    tmp_int_array[20] = dcd_trajectory.traj_format
    for i = 1:20
        write(dcd_file, tmp_int_array[i])
    end

    # block size
    write(dcd_file, Int32(84))

    # ===========================
    # Write MD information string
    # ===========================
    #
    # block size calculation
    n_doc_line = length(dcd_trajectory.md_doc)
    block_size = 4 + 80 * n_doc_line

    # block size
    write(dcd_file, Int32(block_size))

    # write dcd documentation
    write(dcd_file, Int32(n_doc_line))
    for i in 1:n_doc_line - 1
        for j = 1:80
            write(dcd_file, dcd_trajectory.md_doc[i][j])
        end
    end
    new_doc_line = rpad("REMARKS ** MODIFIED BY JULIA CG TOOL **", 80)
    for j = 1:80
        write(dcd_file, new_doc_line[j])
    end

    # block size
    write(dcd_file, Int32(block_size))

    # =====================
    # Write number of atoms
    # =====================
    #
    # block size
    write(dcd_file, Int32(4))

    # write num_atoms
    write(dcd_file, Int32(n_atom))

    # block size
    write(dcd_file, Int32(4))

    # ===================
    # Write trajectory!!!
    # ===================
    #
    coor_block_size = 4 * n_atom

    # loop over frames
    for i_frame in 1:length( dcd_trajectory.conformations )
        # -------------------------
        # Write boundary conditions
        # -------------------------
        if dcd_trajectory.boundary_type == 1
            # block size
            write(dcd_file, Int32(48))

            for j_dim in 1:6
                write(dcd_file, Float64(dcd_trajectory.boundary_box_size[j_dim, i_frame]))
            end

            # block size
            write(dcd_file, Int32(48))
        end

        # ------------
        # coordinate x
        # ------------
        #
        # block size
        write(dcd_file, Int32(coor_block_size))

        # write X
        for j in 1:n_atom
            write(dcd_file, Float32(dcd_trajectory.conformations[i_frame].coors[1, j]))
        end

        # block size
        write(dcd_file, Int32(coor_block_size))

        # ------------
        # coordinate y
        # ------------
        #
        # block size
        write(dcd_file, Int32(coor_block_size))

        # write Y
        for j in 1:n_atom
            write(dcd_file, Float32(dcd_trajectory.conformations[i_frame].coors[2, j]))
        end

        # block size
        write(dcd_file, Int32(coor_block_size))

        # ------------
        # coordinate z
        # ------------
        #
        # block size
        write(dcd_file, Int32(coor_block_size))

        # write Z
        for j in 1:n_atom
            write(dcd_file, Float32(dcd_trajectory.conformations[i_frame].coors[3, j]))
        end

        # block size
        write(dcd_file, Int32(coor_block_size))

    end

    close(dcd_file)

end
