#!/usr/bin/env julia

using Printf

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
