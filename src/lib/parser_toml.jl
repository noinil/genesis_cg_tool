###############################################################################
#                           _____ ___  __  __ _                               #
#                          |_   _/ _ \|  \/  | |                              #
#                            | || | | | |\/| | |                              #
#                            | || |_| | |  | | |___                           #
#                            |_| \___/|_|  |_|_____|                          #
#                                                                             #
###############################################################################

function read_TOML_normal_line(toml_line::AbstractString)
    # remove comments
    if occursin('#', toml_line)
        sep = findfirst("#", toml_line)
        toml_line = toml_line[1 : sep[1] - 1]
    end

    # split by first "="
    sep = findfirst("=", toml_line)
    if isnothing( sep )
        error("ERROR: wrong format in toml: $toml_line")
    end

    t_key_str = strip( toml_line[1 : sep[1] - 1])
    t_val_str = strip( toml_line[sep[end] + 1 : end] )

    t_key = strip( t_key_str, ['\'', '\"'])

    # t_val type: String, Int, Float, Bool
    if t_val_str == "true"
        t_val = true
    elseif t_val_str == "false"
        t_val = false
    elseif t_val_str[1] == '\"' || t_val_str[1] == '\''
        t_val = strip(t_val_str, ['\'', '\"'])
    elseif occursin(r"^[+-]?[1-9][0-9_]*$", t_val_str)
        t_int_str = replace(t_val_str, "_" => "")
        t_val = parse(Int, t_int_str)
    elseif occursin(r"^[+-]?(([0-9]*\.[0-9]+)|([0-9]+\.))([eE][+-]?[0-9]+)?$", t_val_str)
        t_val = parse(Float64, t_val_str)
    elseif t_val_str[1] == '['
        try
            t_expression = Meta.parse(t_val_str)
            t_val = eval(t_expression)
        catch
            error("BUG: cannot understand complex array format in TOML: $toml_line")
        end
    else
        error("ERROR: data type not supported by TOML: $toml_line")
    end

    return (t_key, t_val)
end
                         
"A rather simple subroutine to parse TOML data."
function read_TOML(toml_filename::String)

    # TODO: date
    # TODO: dotted key
    # TODO: multi-line string, array, table...
    # TODO: array

    root_dict = Dict()

    tmp_dict = root_dict
    for line in eachline(toml_filename)
        toml_line = strip(line)

        if length(toml_line) == 0
            continue
        end

        if toml_line[1] == '['
            if toml_line[2] == '['
                sep = findfirst("]]", toml_line)
                dict_name = strip(toml_line[3 : sep[1] - 1])
                if haskey(root_dict, dict_name)
                    push!(root_dict[dict_name], Dict())
                    tmp_dict = root_dict[dict_name][end]
                else
                    root_dict[dict_name] = []
                    push!(root_dict[dict_name], Dict())
                    tmp_dict = root_dict[dict_name][end]
                end
            else
                sep = findfirst("]", toml_line)
                dict_name = strip(toml_line[2 : sep[1] - 1])
                root_dict[dict_name] = Dict()
                tmp_dict = root_dict[dict_name]
            end
        elseif occursin('=', toml_line)
            t_key, t_val = read_TOML_normal_line(toml_line)
            tmp_dict[t_key] = t_val
        end
    end

    return root_dict
end
