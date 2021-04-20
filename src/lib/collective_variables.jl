###############################################################################
#                             Collective Variables                            #
###############################################################################
# compute_nativeness(top, conformation, args)
# compute_center_of_mass(atom_indices, top, conformation, args)
###############################################################################

# =============
# CG nativeness
# =============

"""
    compute_nativeness(t, c, args)

# Arguments
- `t`: GenTopology
- `c`: Conformation
- `args`: other arguments

## Contact type
- 1: cannonical contact type: when r < r0 * cutoff contact is formed
- 2: (1/N) * sum (1 / (1 + exp(beta * (r - r0 * cutoff))))
    (Robert B. Best et al. PNAS 2013)
"""
function compute_nativeness(t::GenTopology, c::Conformation, args::Dict{String, <:Any}=Dict{String, Any}())
    verbose = get(args, "verbose", false)
    q_type  = get(args, "type", 1)
    cutoff  = get(args, "cutoff", 1.2)
    beta    = get(args, "beta", 5.0)

    region  = get(args, "region", [])

    # count number of native contacts
    num_native_contacts = 0
    for p in t.top_pairs
        if p.function_type == 1 || p.function_type == 2
            if length(region) > 0 && (!in(p.i, region) || !in(p.j, region))
                continue
            end
            num_native_contacts += 1
        end
    end

    # loop over all the native contacts
    num_correct_contact = 0
    num_contact_type_2  = 0.0
    for p in t.top_pairs
        if p.function_type == 1 || p.function_type == 2
            r0 = p.r0
            i  = p.i
            j  = p.j
            if length(region) > 0 && (!in(i, region) || !in(j, region))
                continue
            end
            r  = compute_distance(c.coors[:, i], c.coors[:, j])
            if q_type == 1      # simple type
                if r <= r0 * cutoff
                    num_correct_contact += 1
                end
            elseif q_type == 2  # complex type
                num_contact_type_2 += 1 / (1 + exp( beta * (r - cutoff * r0)))
            end
        end
    end

    if q_type == 1      # simple type
        return num_correct_contact / num_native_contacts
    elseif q_type == 2  # complex type
        return num_contact_type_2 / num_native_contacts
    end
end


# ==============
# Center of mass
# ==============

"""
    compute_center_of_mass(idx, t, c, args)

# Arguments
- `idx`: indices of particles
- `t`: GenTopology
- `c`: Conformation
- `args`: other arguments

"""
function compute_center_of_mass(atom_indices::Vector{Int}, t::GenTopology, c::Conformation, args::Dict{String, <:Any}=Dict{String, Any}())

    verbose = get(args, "verbose", false)

    total_mass = 0
    tmp_coor   = zeros(Float64, 3)

    for i in atom_indices
        a_mass      = RES_MASS_DICT[t.top_atoms[i].residue_name]
        a_coor      = c.coors[:, i]
        total_mass += a_mass
        tmp_coor   += a_coor * a_mass
    end
    com = tmp_coor / total_mass

    return com
end
