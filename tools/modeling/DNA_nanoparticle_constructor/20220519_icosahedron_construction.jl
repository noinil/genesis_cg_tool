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
    # nodes[:,  1] = [ 0,  1,  ξ]
    # nodes[:,  2] = [ 0, -1,  ξ]
    # nodes[:,  3] = [ ξ,  0,  1]
    # nodes[:,  4] = [-ξ,  0,  1]
    # nodes[:,  5] = [ 1,  ξ,  0]
    # nodes[:,  6] = [-1,  ξ,  0]
    # nodes[:,  7] = [ 1, -ξ,  0]
    # nodes[:,  8] = [-1, -ξ,  0]
    # nodes[:,  9] = [ ξ,  0, -1]
    # nodes[:, 10] = [-ξ,  0, -1]
    # nodes[:, 11] = [ 0,  1, -ξ]
    # nodes[:, 12] = [ 0, -1, -ξ]
    nodes[:,  1] = [ 0,  1,  ξ] # 1
    nodes[:,  2] = [ 0, -1, -ξ] # 12
    nodes[:,  3] = [ 0, -1,  ξ] # 2
    nodes[:,  4] = [ 0,  1, -ξ] # 11
    nodes[:,  5] = [ ξ,  0,  1] # 3
    nodes[:,  6] = [-ξ,  0, -1] # 10
    nodes[:,  7] = [-ξ,  0,  1] # 4
    nodes[:,  8] = [ ξ,  0, -1] # 9
    nodes[:,  9] = [ 1,  ξ,  0] # 5
    nodes[:, 10] = [-1, -ξ,  0] # 8
    nodes[:, 11] = [-1,  ξ,  0] # 6
    nodes[:, 12] = [ 1, -ξ,  0] # 7
    nodes .*= s_r

    # -----
    # edges
    # -----
    edge_idx_list = []
    i_count = 1
    for i in 1:11
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
    # estimated upper limit of d_0 is 15Å
    n_s = Int(ceil(2 * atan(inv_ξ) * r_P / 15))
    d_0 = 2 * r_P * atan(inv_ξ) / n_s
    sigma = d_0 / sqrt(3) * 2

    @printf("Polyhedron radius: %8.3f Å  |  # subdivision level: %5d  |  σ ~ %8.3f Å \n", r_P, n_s, sigma)

    return n_s

end

function subdivided_polyhedron_construction(f_polyhedron_construction, r_P, DNA_density, linker_length)

    # construct polyhedron
    (polyhedron_vertices_coors, edge_list, face_list) = f_polyhedron_construction(r_P)
    n_poly_vert = size(polyhedron_vertices_coors)[2]
    n_poly_edge = length(edge_list)
    n_poly_face = length(face_list)

    # determine number of subdivision
    n_subdivision = determine_subdivision_level(r_P)

    # calculate number of vertices
    n_new_points_per_edge = n_subdivision - 1
    n_new_points_all_edge = n_new_points_per_edge * n_poly_edge
    n_new_points_per_face = Int((n_subdivision - 1) * (n_subdivision - 2) / 2)
    n_new_points_all_face = n_new_points_per_face * n_poly_face
    n_vertices_total = n_poly_vert + n_new_points_all_edge + n_new_points_all_face

    # ===============
    # DATA structures
    # ===============
    # new coordinates
    new_polyhedron_coors = zeros(Float64, (3, n_vertices_total))
    # flags for DNA-planting particles
    new_polyhedron_flags = zeros(Int, n_vertices_total)

    # ===========
    # subdivision
    # ===========
    # vertices coordinates (first n_poly_vert)
    new_polyhedron_coors[:, 1:n_poly_vert] = polyhedron_vertices_coors
    i_node = n_poly_vert + 1
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
            new_polyhedron_coors[:, i_node] = coor_new
            i_node += 1
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
                new_polyhedron_coors[:, i_node] = coor_new
                i_node += 1
            end
        end
    end


    # =================================
    # Determine vertices to connect DNA
    # =================================
    A_nanop = pi * r_P * r_P
    N_DNA   = Int(ceil(A_nanop / DNA_density))
    if N_DNA <= n_poly_vert
        new_polyhedron_flags[1:N_DNA] = 1
    else
        new_polyhedron_flags[1:n_poly_vert] = 1
        # non-poly-vertex DNA=connecting particles
        n_DNA_node_on_face_edge = N_DNA - n_poly_vert
        n_DNA_node_on_face = Int(ceil(n_DNA_node_on_face_edge * n_new_points_all_face / (n_new_points_all_face + n_new_points_all_edge)))
        n_DNA_node_on_edge = n_DNA_node_on_face_edge - n_DNA_node_on_face
        # ----------------------------
        # put extra DNA nodes on edges
        # ----------------------------
        n_interval_tmp = div(n_new_points_all_edge, n_DNA_node_on_edge)
        new_polyhedron_flags[n_poly_vert+1:n_interval_tmp:n_poly_vert+n_interval_tmp*n_DNA_node_on_edge] .= 1
        # ----------------------------
        # put extra DNA nodes on faces
        # ----------------------------
        n_interval_tmp = div(n_new_points_all_face, n_DNA_node_on_face)
        n_shift_tmp    = n_poly_vert + n_poly_edge
        new_polyhedron_flags[n_shift_tmp+1:n_interval_tmp:n_shift_tmp+n_interval_tmp*n_DNA_node_on_face] .= 1
    end

    return (new_polyhedron_coors, new_polyhedron_flags)

end


if abspath(PROGRAM_FILE) == @__FILE__
    # unit: Angstrom
    r_P = 120

    # unit: nm^-2
    DNA_density = 0.2

    # unit: 1
    len_linker = 5

    # icosahedron_construction(r_P)
    # determine_subdivision_level(r_P)
    my_polyhedron_coors = subdivided_polyhedron_construction(icosahedron_construction, r_P, DNA_density, len_linker)

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
