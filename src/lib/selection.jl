###############################################################################
#                                  selections                                 #
###############################################################################

function parse_selection(s::AbstractString)
    v = Vector{Int}(undef, 0)
    sel_pieces = split(s, r"\s*,\s*", keepempty=false)
    for sel in sel_pieces
        if occursin("to", sel)
            bounds = split(sel, r"\s*to\s*", limit=2)
            first = parse(Int, bounds[1])
            last  = parse(Int, bounds[2])
            for i in first : last
                push!(v, i)
            end
        else
            i = parse(Int, sel)
            push!(v, i)
        end
    end
    return v
end
