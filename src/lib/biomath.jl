###############################################################################
#                   Functions compting geometric quantities                   #
###############################################################################

using LinearAlgebra

# ======================
# Mathematical Structure
# ======================

struct Quaternion
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

struct GeoTransformation
    rotation::Array{<:Real, 2}
    translation::Array{<:Real, 2}
end


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
    return acosd(clamp( dot(v1, v2) / n1 / n2, -1.0, 1.0 ))
end

function compute_vec_angle(vec1::Vector{<:Real}, vec2::Vector{<:Real})
    n1 = norm(vec1)
    n2 = norm(vec2)
    return acosd(clamp( dot(vec1, vec2) / n1 / n2, -1.0, 1.0 ))
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
    dih   = acosd(clamp( dot(c123, c234) / nc123 / nc234, -1.0, 1.0 ))
    c1234 = cross(c123, c234)
    judge = dot(c1234, v23)
    dih   = judge < 0 ? - dih : dih
    return dih
end

# --------
# Centroid
# --------

function centroid(coors::Array{<:Real, 2})
    num_coor = size(coors, 2)
    coor_centroid = zeros(Float64, 3)
    for i_bead in 1 : num_coor
        coor_centroid .+= coors[:, i_bead]
    end
    coor_centroid ./= num_coor
    return coor_centroid
end

# ------------------
# Radius of gyration
# ------------------

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

# ---------------
# Superimposition
# ---------------

"""
    compute_superimposition_transformation(coors_group_1, coors_group_2)

Find out the transformation (rotation + translation) to superimpose Group 1 onto Group 2.

# Arguments
- `coors_group_1`: Group of particles to be moved;
- `coors_group_2`: Group of particles used as target.
"""
function compute_superimposition_transformation(coors_group_1::Array{<:Real, 2}, coors_group_2::Array{<:Real, 2})
    coor_size = size(coors_group_1)[2]

    if coor_size != size(coors_group_2)[2]
        error("Can not perform superimposition for conformations with different size.")
    end

    # Step 1: scaling group 2 to group 3 to math group 1
    measure_group_1 = 0
    measure_group_2 = 0
    for i in 1:coor_size - 1
        measure_group_1 += norm(coors_group_1[:, i] - coors_group_1[:, i + 1])
        measure_group_2 += norm(coors_group_2[:, i] - coors_group_2[:, i + 1])
    end
    measure_scale = measure_group_1 / measure_group_2
    coors_group_3 = coors_group_2 .* measure_scale

    # Step 2: compute centroids
    #
    coor_centroid_1 = sum(coors_group_1, dims=2) .* (1 / coor_size)
    coor_centroid_3 = sum(coors_group_3, dims=2) .* (1 / coor_size)

    # Step 3: shift coordinates to centroid
    #
    coors_shift_1 = coors_group_1 .- coor_centroid_1
    coors_shift_3 = coors_group_3 .- coor_centroid_3

    # SVD
    #
    s = svd(coors_shift_1 * coors_shift_3')

    # rotation
    #
    d = det(s.V * s.U') < 0.0 ? -1.0 : 1.0
    m = diagm([1, 1, d])
    rotation_matrix = s.V * m * s.U'

    # translation
    #
    translation_matrix = ( coor_centroid_3 ./ measure_scale ) - ( rotation_matrix * coor_centroid_1 )

    # final RMSD fit
    #
    fit = GeoTransformation(rotation_matrix, translation_matrix)

    return fit
end

function apply_transformation(t::GeoTransformation, coors_group_old::Array{<:Real, 2})
    coors_group_new = t.rotation * coors_group_old .+ t.translation

    return coors_group_new
end

# ------------
# Compute RMSD
# ------------

function compute_rmsd(coors_group_1::Array{<:Real, 2}, coors_group_2::Array{<:Real, 2})
    coor_size = size(coors_group_1)[2]

    if coor_size != size(coors_group_2)[2]
        error("Can not perform superimposition for conformations with different size.")
    end

    # -----------------------
    # perform superimposition
    # -----------------------
    fit = compute_superimposition_transformation(coors_group_1, coors_group_2)
    coors_group_3 = apply_transformation(fit, coors_group_1)

    d = sum((coors_group_2 - coors_group_3).^2)
    rmsd = sqrt(d / coor_size)

    return rmsd
end

# --------------------------
# Generate a random rotation
# --------------------------
function generate_random_rotation()
    # step 1: random axis
    rand_vec = rand(Float64, 3) .- 0.5
    rand_vec_norm = sqrt(rand_vec' * rand_vec)
    axis = rand_vec ./ rand_vec_norm
    (ax, ay, az) = axis

    outer_product_matrix = axis * axis'
    cross_product_matrix = [0 -az ay; az 0 -ax; -ay ax 0]

    # step 2: random θ
    theta = rand() * 2 * pi
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    R = cos_theta * I + (1 - cos_theta) * outer_product_matrix + sin_theta * cross_product_matrix

    # return the rotation matrix
    return R
end

# ------------------------------
# Rotate by theta around an axis
# ------------------------------
function rotation_matrix_around_axis(axis::Vector{<:Real}, theta)
    cost = cosd(theta)
    sint = sind(theta)
    _1_m_cost = 1 - cost
    ux = axis[1]
    uy = axis[2]
    uz = axis[3]
    rotmat = reshape([cost + ux*ux*_1_m_cost,
                      uy*ux*_1_m_cost + uz*sint,
                      uz*ux*_1_m_cost - uy*sint,
                      ux*uy*_1_m_cost - uz*sint,
                      cost + uy*uy*_1_m_cost,
                      uz*uy*_1_m_cost + ux*sint,
                      ux*uz*_1_m_cost + uy*sint,
                      uy*uz*_1_m_cost - ux*sint,
                      cost + uz*uz*_1_m_cost], (3, 3))
    return rotmat
end

function rotation_matrix_around_x(theta)
    cost = cosd(theta)
    sint = sind(theta)
    rotmat = [1 0 0; 0 cost -sint; 0 sint cost]
    return rotmat
end

function rotation_matrix_around_y(theta)
    cost = cosd(theta)
    sint = sind(theta)
    rotmat = [cost 0 sint; 0 1 0; -sint 0 cost]
    return rotmat
end

function rotation_matrix_around_z(theta)
    cost = cosd(theta)
    sint = sind(theta)
    rotmat = [cost -sint 0; sint cost 0; 0 0 1]
    return rotmat
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
