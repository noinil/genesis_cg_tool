###############################################################################
#                   Functions compting geometric quantities                   #
###############################################################################

using LinearAlgebra

# ===================
# Geometric Functions
# ===================

# --------
# Distance
# --------

function compute_distance(coor1::Vector{<:Real}, coor2::Vector{<:Real})
    d = coor1 - coor2
    return norm(d)
end

# -----
# Angle
# -----

function compute_angle(coor1::Vector{<:Real}, coor2::Vector{<:Real}, coor3::Vector{<:Real})
    v1 = coor1 - coor2
    v2 = coor3 - coor2
    n1 = norm(v1)
    n2 = norm(v2)
    return acos( dot(v1, v2) / n1 / n2) / pi * 180.0
end

function compute_vec_angle(vec1::Vector{<:Real}, vec2::Vector{<:Real})
    n1 = norm(vec1)
    n2 = norm(vec2)
    return acos( dot(vec1, vec2) / n1 / n2) / pi * 180.0
end

# --------
# Dihedral
# --------

function compute_dihedral(coor1::Vector{<:Real}, coor2::Vector{<:Real}, coor3::Vector{<:Real}, coor4::Vector{<:Real})
    v12   = coor2 - coor1
    v23   = coor3 - coor2
    v34   = coor4 - coor3
    c123  = cross(v12, v23)
    c234  = cross(v23, v34)
    nc123 = norm(c123)
    nc234 = norm(c234)
    dih   = acos( dot(c123, c234) / nc123 / nc234)
    c1234 = cross(c123, c234)
    judge = dot(c1234, v23)
    dih   = judge < 0 ? - dih : dih
    return dih / pi * 180.0
end



# ===================
# Physical properties
# ===================

# --------------
# Center of mass
# --------------

function compute_center_of_mass(atom_indices::Vector{Int}, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
    total_mass      = 0
    tmp_coor        = zeros(Float64, 3)
    for i in atom_indices
        a_mass      = ATOM_MASS_DICT[atom_names[i][1]]
        a_coor      = atom_coors[:, i]
        total_mass += a_mass
        tmp_coor   += a_coor * a_mass
    end
    com = tmp_coor / total_mass
    return com
end

function centeroid(coors::Array{<:Real, 2})
    num_coor = size(coors, 2)
    coor_centroid = zeros(Float64, 3)
    for i_bead in 1 : num_coor
        coor_centroid .+= coors[:, i_bead]
    end
    coor_centroid ./= num_coor
    return coor_centroid
end

function radius_of_gyration(coors::Array{<:Real, 2})
    num_coor = size(coors, 2)
    coor_centroid = zeros(Float64, 3)
    for i_bead in 1 : num_coor
        coor_centroid .+= coors[:, i_bead]
    end
    coor_centroid ./= num_coor

    dist_sq_sum = 0
    for i_bead in 1 : num_coor
        v = coors[:, i_bead] - coor_centroid
        dist_sq_sum += v' * v
    end
    rg = sqrt(dist_sq_sum / num_coor)
end

