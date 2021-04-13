#!/usr/bin/env julia

using Printf
using ArgParse

function read_genesis_log(genesis_log_filename::AbstractString)

    # ----------------------------
    # read in lines from data file
    # ----------------------------
    #
    info_line  = ""
    data_lines = []
    for line in eachline(genesis_log_filename)
        if startswith(line, "INFO:")
            if length(info_line) > 0
                push!(data_lines, line)
            else
                keywords = split(line)
                if !isdigit(keywords[1][1])
                    info_line = line
                end
            end
        end
    end

    num_steps = length(data_lines)

    # ----------------
    # Extract keywords
    # ----------------
    keywords  = split(info_line)[2:end]

    num_keywords = length(keywords)
    keyword_dict = Dict(keywords .=> [i for i = 1 : num_keywords])

    # println(keyword_dict)

    # ------------------------
    # Extract time series data
    # ------------------------
    data_matrix = zeros(Float64, (num_steps, num_keywords))

    for ( i, line ) in enumerate( data_lines )
        words = split(line)[2:end]
        new_vec = [parse(Float64, w) for w in words]
        data_matrix[i, :] = new_vec
    end

    return (keyword_dict, data_matrix)
end

function read_genesis_remd_log(genesis_remd_log_filename::AbstractString)

    par_dim_1 = []
    i_step    = []
    exchange_data = []
    param_data    = []
    replica_data  = []

    # ----------------------------
    # read in lines from data file
    # ----------------------------
    #
    rep2par_map_lines = []
    par2rep_map_lines = []
    exchange_info_lines = []

    is_setup_block = false
    is_log_block   = false
    is_exchange_block = false
    for line in eachline(genesis_remd_log_filename)
        if startswith(line, "Setup_Remd>")
            is_setup_block = true
        end
        if startswith(line, "[STEP5]")
            is_setup_block = false
            is_log_block = true
        end
        if is_setup_block
            words = split(line)
            if length(words) > 0 && words[1] == "Dim"
                par_dim_1 = [parse(Float64, w) for w in words[4:end]]
            end
        end
        if is_log_block
            words = split(line)
            if length(words) < 1
                continue
            end
            if words[1] == "RepIDtoParmID:"
                push!(rep2par_map_lines, line)
            elseif words[1] == "ParmIDtoRepID:"
                push!(par2rep_map_lines, line)
            elseif words[1] == "REMD>" && words[2] == "Step:"
                push!(i_step, parse(Int, words[3]))
            elseif words[1] == "Replica" && words[2] == "ExchangeTrial"
                is_exchange_block = true
                exchange_info_lines = []
                continue
            elseif words[1] == "Parameter"
                is_exchange_block = false
            end
            if is_exchange_block
                push!(exchange_info_lines, line)
            end
        end
    end

    # -----------------
    # basic information
    # -----------------
    num_frames = length(rep2par_map_lines)
    num_params = length(par_dim_1)

    # ================
    # acceptance ratio
    # ================
    exchange_from  = Vector{Int}(undef, 0)
    exchange_to    = Vector{Int}(undef, 0)
    exchange_ratio = Vector{Float64}(undef, 0)
    for line in exchange_info_lines
        words   = split(line)
        e_from  = parse(Int, words[2])
        e_to    = parse(Int, words[4])
        e_ratio = parse(Float64, words[6]) / parse(Float64, words[8])
        if e_from < e_to
            push!(exchange_from, e_from)
            push!(exchange_to, e_to)
            push!(exchange_ratio, e_ratio)
        end
    end
    sorted_indx = sortperm(exchange_from)
    push!(exchange_data, exchange_from[sorted_indx])
    push!(exchange_data, exchange_to[sorted_indx])
    push!(exchange_data, exchange_ratio[sorted_indx])

    # ===============================
    # specific param on different rep
    # ===============================
    param_data = zeros(Int, (num_frames, 1 + num_params))
    for ( i, line ) in enumerate( par2rep_map_lines )
        words = split(line)[2:end]
        new_vec = [parse(Int, w) for w in words]
        param_data[i, 1] = i_step[i]
        param_data[i, 2:end] = new_vec[:]
    end

    # =====================================
    # specific replica with different param
    # =====================================
    replica_data = zeros(Int, (num_frames, 1 + num_params))
    for ( i, line ) in enumerate( rep2par_map_lines )
        words = split(line)[2:end]
        new_vec = [parse(Int, w) for w in words]
        replica_data[i, 1] = i_step[i]
        replica_data[i, 2:end] = new_vec[:]
    end

    return (par_dim_1, exchange_data, param_data, replica_data)

end

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--remd", "-R"
        help     = "Analyze REMD results."
        action   = :store_true

        "--remd-main-log"
        help     = "File name of the master log of REMD."
        arg_type = String
        default  = ""
    end

    return parse_args(s)
end
