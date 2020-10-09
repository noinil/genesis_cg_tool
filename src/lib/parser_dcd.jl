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
    md_doc::Vector{String}
    num_atoms::Int
    boundary_box_size::Array{Float64, 2}
    conformations::Vector{Conformation}
end

                   
function read_dcd(dcd_filename::String, args::Dict{String, <:Any}=Dict{String, Any}())

    verbose  = get(args, "verbose", false)

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
        if length(strip(doc_line)) == 0
            continue
        end
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

        println("============================================================")
        println("> DCD information:")
    end

    if verbose
        for line in dcd_doc_lines
            println("  > ", replace(line, "REMARKS " => ""))
        end
    end

    # ======================
    # Read in Coordinates!!!
    # ======================

    boundary_conditions = zeros(Float64, 6, n_frames)
    conformations = Vector{Conformation}(undef, 0)

    data_block_size = n_particles * 4 # Single precision (4-byte), Float32
    is_broken_trajectory = false
    for t in 1 : n_frames
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
            boundary_conditions[:, t] = bc_size
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
                                    n_frames,
                                    n_ts_first,
                                    n_ts_interval,
                                    n_ts_all,
                                    bc_flag,
                                    dcd_doc_lines,
                                    n_particles,
                                    boundary_conditions,
                                    conformations)

    return new_trajectory

end


function write_dcd(dcd_trajectory::DCD_trajectory, dcd_filename::AbstractString)

    dcd_file = open(dcd_filename, "w")

    # ======================
    # Write head information
    # ======================
    # block size
    write(dcd_file, Int32(84))

    # "CORD" or "VELD"
    for i = 1:4
        write(dcd_file, dcd_trajectory.traj_type[i])
    end

    # simulation info
    tmp_int_array = Array{Int32}(undef, 0)
    tmp_int_array[1]  = dcd_trajectory.traj_frames
    tmp_int_array[2]  = dcd_trajectory.traj_first_step
    tmp_int_array[3]  = dcd_trajectory.traj_output_interval
    tmp_int_array[4]  = dcd_trajectory.traj_steps
    tmp_int_array[11] = dcd_trajectory.boundary_type
    for i = 1:20
        write(dcd_file, tmp_int_array[i])
    end

    # block size
    write(dcd_file, Int32(84))

    # ===========================
    # Write MD information string
    # ===========================

    n_doc_line = length(dcd_trajectory.md_doc)
    block_size = 4 + 80 * n_doc_line

    # block size
    write(dcd_file, Int32(block_size))

    # write dcd documentation
    for i in 1:n_doc_line
        for j = 1:80
            write(dcd_file, dcd_trajectory.md_doc[i][j])
        end
    end

    # block size
    write(dcd_file, Int32(block_size))

end
