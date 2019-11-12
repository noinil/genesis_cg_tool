###############################################################################
#                           _____ ___  __  __ _                               #
#                          |_   _/ _ \|  \/  | |                              #
#                            | || | | | |\/| | |                              #
#                            | || |_| | |  | | |___                           #
#                            |_| \___/|_|  |_|_____|                          #
#                                                                             #
###############################################################################

function parse_TOML_normal_line(toml_line::String)
    # remove comments
    if occursin('#', toml_line)
        sep = findfirst("#", toml_line)
        toml_line = toml_line[1 : sep[1] - 1]
    end

    # split by first "="
    sep = findfirst("=", toml_line)
    if isnothing( sep )
        error("ERROR: wrong format in toml.")
    end

    t_key_str = strip( toml_line[1 : sep[1] - 1])
    t_val_str = strip( toml_line[sep[end] + 1 : end] )

    t_key = strip( t_key_str, ['\'', '\"'])

    # t_val type: String, Int, Float, Bool
    if t_val_str == "ture"
        t_val = true
    elseif t_val_str == "false"
        t_val = false
    elseif t_val_str[1] == '\"' || t_val_str[1] == '\''
        t_val = strip(t_val_str, ['\'', '\"'])
    elseif occursin(r"^[+-]?[1-9][0-9_]*$", t_val_str)
        t_int_str = replace(t_val_str, " " => "")
        t_val = parse(Int, t_int_str)
    elseif occursin(r"^[+-]?(([0-9]*\.[0-9]+)|([0-9]+\.))([eE][+-]?[0-9]+)?$", t_val_str)
        t_val = parse(Float64, t_int_str)
    elseif t_val_str[1] == '['
        try
            t_expression = Meta.parse(t_val_str)
            t_val = eval(t_expression)
        catch
            error("BUG: cannot understand complex array format in TOML.")
        end
    else
        error("ERROR: data type not supported by TOML!")
    end

    return (t_key, t_val)
end
                         
"A rather simple subroutine to parse TOML data."
function read_TOML(toml_filename::String)

    # TODO: date
    # TODO: multi-line string
    # TODO: array

    for line in eachline(toml_filename)
        toml_line = strip(line)

        if toml_line[1] == '['
            # do parse block 
        elseif occursin('=', toml_line)
            # do parse normal line
        end
    end

end
