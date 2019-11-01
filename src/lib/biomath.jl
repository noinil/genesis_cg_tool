module BioMath

using LinearAlgebra

# ===================
# Geometric Functions
# ===================

# --------
# Distance
# --------

function compute_distance(coor1, coor2)
    d = coor1 - coor2
    return norm(d)
end

# -----
# Angle
# -----

function compute_angle(coor1, coor2, coor3)
    v1 = coor1 - coor2
    v2 = coor3 - coor2
    n1 = norm(v1)
    n2 = norm(v2)
    return acos( dot(v1, v2) / n1 / n2) / pi * 180.0
end

function compute_vec_angle(vec1, vec2)
    n1 = norm(vec1)
    n2 = norm(vec2)
    return acos( dot(vec1, vec2) / n1 / n2) / pi * 180.0
end

# --------
# Dihedral
# --------

function compute_dihedral(coor1, coor2, coor3, coor4)
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



end
