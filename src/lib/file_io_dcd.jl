###############################################################################
#                                 _           _                               #
#                              __| |  ___  __| |                              #
#                             / _` | / __|/ _` |                              #
#                            | (_| || (__| (_| |                              #
#                             \__,_| \___|\__,_|                              #
#                                                                             #
###############################################################################

using Printf
                   
function read_dcd(dcd_filename::String)

    dcd_file = open(dcd_filename, "r")

    println("============================================================")
    println("> Open DCD file: ", dcd_filename)

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
        for i = 1:80
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

    for line in dcd_doc_lines
        println("  > ", replace(line, "REMARKS " => ""))
    end

    # ======================
    # Read in Coordinates!!!
    # ======================

    conformations = Vector{Conformation}(undef, 0)

    data_block_size = n_particles * 4 # Single precision, Float32
    for t in 1 : n_frames
        # ----------------
        # Read in box info
        # ----------------
        if bc_flag == 1
            bc_size = zeros(Float64, 6)
            for i  in 1 : n_particles
                l = read(dcd_file, Float64)
                bc_size[i] = l
            end
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
                error("ERROR: wrong block size in coordinate reading!")
            end
        end
        push!(conformations, Conformation(n_particles, coors))

        if eof(dcd_file)
            break
        end
    end

    println("------------------------------------------------------------")
    println(">[1;32m FINISH[0m reading the DCD file. Have fun!")
    println("============================================================")
    println("         ")

    close(dcd_file)

    return conformations

end
