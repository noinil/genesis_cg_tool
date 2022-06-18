#!/usr/bin/env julia

using LinearAlgebra
using Printf

function icosahedron_construction(r_P)

    # key coordinate
    ξ = (sqrt(5) + 1) / 2
    # circum-radius
    r_c = sqrt(ξ * ξ + 1)
    # scaling factor
    s_r = r_P / r_c

    # --------
    # vertices
    # --------
    nodes = zeros(Float64, (3, 12))
    nodes[:,  1] = [ 0,  1,  ξ]
    nodes[:,  2] = [ 0, -1,  ξ]
    nodes[:,  3] = [ ξ,  0,  1]
    nodes[:,  4] = [-ξ,  0,  1]
    nodes[:,  5] = [ 1,  ξ,  0]
    nodes[:,  6] = [-1,  ξ,  0]
    nodes[:,  7] = [-1, -ξ,  0]
    nodes[:,  8] = [ 1, -ξ,  0]
    nodes[:,  9] = [ ξ,  0, -1]
    nodes[:, 10] = [-ξ,  0, -1]
    nodes[:, 11] = [ 0,  1, -ξ]
    nodes[:, 12] = [ 0, -1, -ξ]
    nodes .*= s_r

    # -----
    # edges
    # -----
    edge_idx_list = []
    i_count = 1
    for i in 1:12
        for j in i + 1:12
            coor_1 = nodes[:, i]
            coor_2 = nodes[:, j]
            distance = norm(coor_1 - coor_2)
            # 1.05 * edge_length = 1.05 * 2 * s_r
            if distance < 2.1 * s_r
                push!(edge_idx_list, [i, j])
            end
        end
    end

    # -----
    # faces
    # -----
    face_idx_list = []
    i_count = 1
    for i_edge in 1:length(edge_idx_list)
        i1 = edge_idx_list[i_edge][1]
        i2 = edge_idx_list[i_edge][2]
        for i3 in i2 + 1:12
            if [i1, i3] in edge_idx_list && [i2, i3] in edge_idx_list
                push!(face_idx_list, [i1, i2, i3])
            end
        end
    end

    # return
    return (nodes, edge_idx_list, face_idx_list)

end

function determine_subdivision_level(r_P)

    inv_ξ = 2 / (sqrt(5) + 1)
    n_s = Int(ceil(2 * atan(inv_ξ) * r_P / 25))
    d_0 = 2 * r_P * atan(inv_ξ) / n_s
    sigma = d_0 / sqrt(3)

    @printf("Polyhedron radius: %8.3f Å  |  # subdivision level: %5d  |  σ ~ %8.3f Å \n", r_P, n_s, sigma)

    return n_s

end

function subdivided_polyhedron_construction(f_polyhedron_construction, r_P)

    # construct polyhedron
    (polyhedron_vertices_coors, edge_list, face_list) = f_polyhedron_construction(r_P)
    # determine number of subdivision
    n_subdivision = determine_subdivision_level(r_P)

    # calculate number of vertices
    n_new_points_per_edge = (n_subdivision - 1) * length(edge_list)
    n_new_points_per_face = Int((n_subdivision - 1) * (n_subdivision - 2) / 2 * length(face_list))
    n_vertices_total = size(polyhedron_vertices_coors)[2] + n_new_points_per_edge + n_new_points_per_face

    # new coordinates
    new_polyhedron_coors = zeros(Float64, (3, n_vertices_total))
    new_polyhedron_coors[:, 1:size(polyhedron_vertices_coors)[2]] = polyhedron_vertices_coors

    # ===========
    # subdivision
    # ===========
    i_vertex = size(polyhedron_vertices_coors)[2] + 1
    # -------------
    # step 1: edges
    # -------------
    for i_edge in 1:length(edge_list)
        v1 = edge_list[i_edge][1]
        v2 = edge_list[i_edge][2]
        coor_v1 = polyhedron_vertices_coors[:, v1]
        coor_v2 = polyhedron_vertices_coors[:, v2]
        e1 = (coor_v2 - coor_v1) ./ n_subdivision
        for i in 1:n_subdivision - 1
            coor_new = coor_v1 .+ (i .* e1)
            coor_new .*= r_P / norm(coor_new)
            new_polyhedron_coors[:, i_vertex] = coor_new
            i_vertex += 1
        end
    end

    # -------------
    # step 2: faces
    # -------------
    for i_face in 1:length(face_list)
        v1 = face_list[i_face][1]
        v2 = face_list[i_face][2]
        v3 = face_list[i_face][3]
        coor_v1 = polyhedron_vertices_coors[:, v1]
        coor_v2 = polyhedron_vertices_coors[:, v2]
        coor_v3 = polyhedron_vertices_coors[:, v3]
        e1 = (coor_v2 - coor_v1) ./ n_subdivision
        e2 = (coor_v3 - coor_v2) ./ n_subdivision
        for i in 2:n_subdivision - 1
            for j in 1:i - 1
                coor_new = coor_v1 .+ (i .* e1) .+ (j .* e2)
                coor_new .*= r_P / norm(coor_new)
                new_polyhedron_coors[:, i_vertex] = coor_new
                i_vertex += 1
            end
        end
    end

    return new_polyhedron_coors

end


if abspath(PROGRAM_FILE) == @__FILE__
    # unit: Angstrom
    r_P = 120

    # icosahedron_construction(r_P)
    # determine_subdivision_level(r_P)
    my_polyhedron_coors = subdivided_polyhedron_construction(icosahedron_construction, r_P)

    # test
    io = open("test.pdb", "w")
    n_atom = size(my_polyhedron_coors)[2]
    for i_atom in 1:n_atom
        coor = my_polyhedron_coors[:, i_atom]
        @printf(io, "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s%2d \n",
                i_atom, "AU", "ICO", "A", i_atom, coor[1], coor[2], coor[3],
                1.0, 1.0, "NANOP", "AU", 0)
    end
    close(io)
end
