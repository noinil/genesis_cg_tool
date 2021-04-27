###############################################################################
#                                 Conformation                                #
###############################################################################

struct Conformation
    num_particle::Int
    coors::Array{Float64, 2}
end

function centroid(c::Conformation)
    coor_centroid = zeros(Float64, 3)
    for i_bead in 1 : c.num_particle
        coor_centroid .+= c.coors[:, i_bead]
    end
    coor_centroid ./= c.num_particle
    return coor_centroid
end

function radius_of_gyration(c::Conformation)
    coor_centroid = zeros(Float64, 3)
    for i_bead in 1 : c.num_particle
        coor_centroid .+= c.coors[:, i_bead]
    end
    coor_centroid ./= c.num_particle

    dist_sq_sum = 0
    for i_bead in 1 : c.num_particle
        v = c.coors[:, i_bead] - coor_centroid
        dist_sq_sum += v' * v
    end
    rg = sqrt(dist_sq_sum / c.num_particle)
end


function radius_of_circumshpere(c::Conformation)
    coor_centroid = zeros(Float64, 3)
    for i_bead in 1 : c.num_particle
        coor_centroid .+= c.coors[:, i_bead]
    end
    coor_centroid ./= c.num_particle

    tmp_dist = 0
    for i_bead in 1 : c.num_particle
        v = c.coors[:, i_bead] - coor_centroid
        v_norm_sqr = sqrt( v' * v )
        tmp_dist = v_norm_sqr > tmp_dist ? v_norm_sqr : tmp_dist
    end
    rc = tmp_dist
end
