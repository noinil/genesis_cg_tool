#!/usr/bin/env julia

include("../../../src/lib/gcj.jl")
using ArgParse
using LinearAlgebra
using Printf

function icosahedron_construction(r_P)

    println(" ------------------------------------------------------------ ")
    println(" Constructing icosahedron...")

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

    println(" > DONE! ")

    # return
    return (nodes, edge_idx_list, face_idx_list)

end

function determine_subdivision_level(r_P)

    println(" ------------------------------------------------------------ ")
    println(" Calculating subdivision level...")
    inv_ξ = 2 / (sqrt(5) + 1)
    # estimated upper limit of d_0 is 15Å
    n_s = Int(ceil(2 * atan(inv_ξ) * r_P / 15))
    d_0 = 2 * r_P * atan(inv_ξ) / n_s
    sigma = d_0 / sqrt(3) * 2

    @printf(" > Polyhedron radius : %8.3f Å  |  # subdivision level: %8d  \n", r_P, n_s)
    @printf(" >      Estimated d0 : %8.3f Å  |                   σ : %8.3f Å \n", d_0, sigma)
    println(" > DONE! ")

    return n_s

end

function subdivided_polyhedron_construction(f_polyhedron_construction, r_P, DNA_density)

    # construct polyhedron
    (polyhedron_vertices_coors, edge_list, face_list) = f_polyhedron_construction(r_P)
    n_poly_vert = size(polyhedron_vertices_coors)[2]
    n_poly_edge = length(edge_list)
    n_poly_face = length(face_list)

    # determine number of subdivision
    n_subdivision = determine_subdivision_level(r_P)

    println(" ------------------------------------------------------------ ")
    # calculate number of vertices
    n_new_points_per_edge = n_subdivision - 1
    n_new_points_all_edge = n_new_points_per_edge * n_poly_edge
    n_new_points_per_face = Int((n_subdivision - 1) * (n_subdivision - 2) / 2)
    n_new_points_all_face = n_new_points_per_face * n_poly_face
    n_nodes_total = n_poly_vert + n_new_points_all_edge + n_new_points_all_face

    @printf(" > Adding %d edge nodes and %d face nodes. \n", n_new_points_all_edge, n_new_points_all_face)
    @printf(" > Totally %d beads in the polyhedron. \n", n_nodes_total)

    # ===============
    # DATA structures
    # ===============
    # new coordinates
    new_polyhedron_coors = zeros(Float64, (3, n_nodes_total))
    # flags for DNA-planting particles
    new_polyhedron_flags = zeros(Int, n_nodes_total)

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
    A_nanop = 4 * pi * r_P * r_P
    N_DNA   = Int(ceil(A_nanop * DNA_density))
    println(" ------------------------------------------------------------ ")
    @printf(" > Marking %d DNA-connecting nodes in nanoparticle\n", N_DNA)
    if N_DNA <= n_poly_vert
        new_polyhedron_flags[1:N_DNA] .= 1
        @printf(" >> %d vertices marked \n", N_DNA)
    else
        new_polyhedron_flags[1:n_poly_vert] .= 1
        @printf(" >> %d vertices marked \n", n_poly_vert)
        # non-poly-vertex DNA=connecting particles
        n_DNA_node_on_face_edge = N_DNA - n_poly_vert
        n_DNA_node_on_face = Int(ceil(n_DNA_node_on_face_edge * n_new_points_all_face / (n_new_points_all_face + n_new_points_all_edge)))
        n_DNA_node_on_edge = n_DNA_node_on_face_edge - n_DNA_node_on_face
        @printf(" >> %d edge-nodes marked \n", n_DNA_node_on_edge)
        @printf(" >> %d face-nodes marked \n", n_DNA_node_on_face)
        # ----------------------------
        # put extra DNA nodes on edges
        # ----------------------------
        n_DNA_per_edge, n_DNA_extra_edge = divrem(n_DNA_node_on_edge, n_poly_edge)
        n_tmp_edge_node_count = n_poly_vert
        for k_edge in 1:n_poly_edge
            if k_edge <= n_DNA_extra_edge
                n_tmp = n_DNA_per_edge + 1
            else
                n_tmp = n_DNA_per_edge
            end
            if n_tmp >= 1
                n_interval_tmp = div(n_new_points_per_edge, n_tmp + 1)
                new_polyhedron_flags[n_tmp_edge_node_count + n_interval_tmp:n_interval_tmp:n_tmp_edge_node_count + n_interval_tmp*n_tmp] .= 1
                n_tmp_edge_node_count += n_new_points_per_edge
            end
        end
        # ----------------------------
        # put extra DNA nodes on faces
        # ----------------------------
        n_DNA_per_face, n_DNA_extra_face = divrem(n_DNA_node_on_face, n_poly_face)
        n_tmp_face_node_count = n_poly_vert + n_new_points_all_edge
        for k_face in 1:n_poly_face
            if k_face <= n_DNA_extra_face
                n_tmp = n_DNA_per_face + 1
            else
                n_tmp = n_DNA_per_face
            end
            if n_tmp >= 1
                n_interval_tmp = div(n_new_points_per_face, n_tmp + 1)
                new_polyhedron_flags[n_tmp_face_node_count + n_interval_tmp:n_interval_tmp:n_tmp_face_node_count + n_interval_tmp*n_tmp] .= 1
                n_tmp_face_node_count += n_new_points_per_face
            end
        end
    end

    println(" > DONE! ")

    return (new_polyhedron_coors, new_polyhedron_flags)

end


function generate_nanoparticle_with_DNA(nanoP_coors, nanoP_flags, linker_length, DNA_name, out_name, do_debug)

    println(" ------------------------------------------------------------ ")
    println(" Constructing nanoparticle + DNA...")

    if linker_length < 3
        println(" Please use a longer linker (at least 3 beads)!")
        exit()
    end

    # =============================
    # Make a new topology structure
    # =============================

    top_default_params             = GenTopDefault(0, 0, false, 0.0, 0.0)
    top_default_atomtype           = Vector{GenTopAtomType}(undef, 0)
    top_default_CGDNA_bp           = Vector{GenTopCGDNABasepairType}(undef, 0)
    top_default_CGDNA_bs           = Vector{GenTopCGDNABasestackType}(undef, 0)
    top_default_CGDNA_cs           = Vector{GenTopCGDNABasecrossType}(undef, 0)
    top_default_CGDNA_exv          = Vector{GenTopCGDNAExvType}(undef, 0)
    top_default_CGPro_flx_angle    = Vector{GenTopCGProAICGFlexAngleType}(undef, 0)
    top_default_CGPro_flx_dihedral = Vector{GenTopCGProAICGFlexDihedralType}(undef, 0)

    global_index_2_local_index     = Vector{Int}(undef, 0)
    global_index_2_local_molid     = Vector{Int}(undef, 0)
    top_atoms                      = Vector{GenTopAtom}(undef, 0)
    top_bonds                      = Vector{GenTopBond}(undef, 0)
    top_angles                     = Vector{GenTopAngle}(undef, 0)
    top_dihedrals                  = Vector{GenTopDihedral}(undef, 0)
    top_pairs                      = Vector{GenTopPair}(undef, 0)
    top_exclusions                 = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos                     = Vector{GenTopPWMcos}(undef, 0)
    top_pwmcosns                   = Vector{GenTopPWMcos}(undef, 0)
    top_idr_hps                    = Vector{GenTopRegion}(undef, 0)
    top_idr_kh                     = Vector{GenTopRegion}(undef, 0)
    top_mol_list                   = Vector{GenTopMolList}(undef, 0)

    # ---------------------
    # read DNA top and coor
    # ---------------------
    tmp_top_name = @sprintf("%s.top", DNA_name)
    tmp_crd_name = @sprintf("%s.gro", DNA_name)
    println(" > Reading topology and coordinates from: ", tmp_top_name, " and ", tmp_crd_name)
    meta_DNA_top = read_grotop(tmp_top_name)
    meta_DNA_cnf = read_grocrd(tmp_crd_name)
    n_atom_DNA = meta_DNA_top.num_atom

    # ====================================
    # prepare coordinates of all particles
    # ====================================
    N_DNA = count(>(0), nanoP_flags)
    println(" > Number of DNAs to add: ", N_DNA)
    n_atom_nanoP = length(nanoP_flags)
    n_atom_total = 1 + n_atom_nanoP + N_DNA * (linker_length + n_atom_DNA)
    println(" > Total number of CG beads:", n_atom_total, " : ")
    @printf("   >> %12d for nanoparticle; \n", 1 + n_atom_nanoP)
    @printf("   >> %12d for linkers; \n", N_DNA * linker_length)
    @printf("   >> %12d for DNAs; \n", N_DNA * n_atom_DNA)

    # -----------------
    # nanoP coordinates
    # -----------------
    atom_coors = zeros(Float64, (3, n_atom_total))
    atom_coors[:, 1] = [0, 0, 0] # the central bead of nanoparticle!
    atom_coors[:, 2:1 + n_atom_nanoP] = nanoP_coors[:, :]

    # -----------------------
    # adding linkers and DNAs
    # -----------------------
    println(" > Adding coordinates of linkers and DNAs")
    i_atom = 1 + n_atom_nanoP
    for (i_nanoP, f_nanoP) in enumerate(nanoP_flags)
        if f_nanoP > 0
            anchor_coor = nanoP_coors[:, i_nanoP]
            v0 = anchor_coor / norm(anchor_coor) # direction vector
            if v0[3] > 0.5
                v_tmp = [v0[1], v0[3], -v0[2]]
            else
                v_tmp = [v0[2], -v0[1], v0[3]]
            end
            vx = cross(v_tmp, v0) # normal vector 1
            vx0 = vx / norm(vx)
            vy = cross(v0, vx0)    # normal vector 2
            vy0 = vy / norm(vy)
            frame_local = [vx0 vy0 v0]

            # adding linker coordinates...
            for j in 1:linker_length
                i_atom += 1
                anchor_coor += v0 * 3.33
                atom_coors[:, i_atom] = anchor_coor
            end
            anchor_coor .+= v0 * 3.33

            # adding DNA coordinates...
            dna_coors_0 = meta_DNA_cnf.coors .- meta_DNA_cnf.coors[:, 1]
            dna_coors_1 = (frame_local * dna_coors_0) .+ anchor_coor
            atom_coors[:, i_atom + 1:i_atom + n_atom_DNA] = dna_coors_1
            i_atom += n_atom_DNA
        end
    end


    # ========================
    # prepare topology for all
    # ========================
    println(" > Preparing topology of the whole structure...")
    # ---------
    # [ atoms ]
    # ---------
    println("   >> adding atoms")
    for i_bead in 1 : 1 + n_atom_nanoP
        # new_atom = GenTopAtom(i_bead, a_type, r_indx, r_name, a_name, f_type, charge, mass, c_id, s_name)
        new_atom = GenTopAtom(i_bead, "NP", 1, "NP", "NP", AICG_ATOM_FUNC_NR, 0, 200, 1, "NANOP")
        push!(top_atoms, new_atom)
        push!(global_index_2_local_index, i_bead)
        push!(global_index_2_local_molid, 1)
    end
    n_tmp_res = 1
    n_tmp_bead = 1 + n_atom_nanoP
    cid = 1
    for i_DNA in 1:N_DNA
        cid += 1
        for j_bead in 1 : linker_length
            n_tmp_res += 1
            n_tmp_bead += 1
            # new_atom = GenTopAtom(i_bead, a_type, r_indx, r_name, a_name, f_type, charge, mass, c_id, s_name)
            new_atom = GenTopAtom(n_tmp_bead, "LNK", n_tmp_res, "LNK", "LNK", AICG_ATOM_FUNC_NR, 0, 100, 1, "THIOL")
            push!(top_atoms, new_atom)
            push!(global_index_2_local_index, j_bead)
            push!(global_index_2_local_molid, 1)
        end
        n_dna_res = 0
        for k_bead in 1 : n_atom_DNA
            n_tmp_bead += 1
            # new_atom = GenTopAtom(i_bead, a_type, r_indx, r_name, a_name, f_type, charge, mass, c_id, s_name)
            my_atom = meta_DNA_top.top_atoms[k_bead]
            if my_atom.chain_id == 1
                new_chain_id = 1
            else
                new_chain_id = cid
            end
            n_dna_res = my_atom.residue_index
            new_atom = GenTopAtom(n_tmp_bead, my_atom.atom_type,
                                  n_tmp_res + n_dna_res,
                                  my_atom.residue_name,
                                  my_atom.atom_name,
                                  AICG_ATOM_FUNC_NR,
                                  my_atom.charge,
                                  my_atom.mass,
                                  new_chain_id, my_atom.seg_name)
            push!(top_atoms, new_atom)
            push!(global_index_2_local_index, k_bead)
            push!(global_index_2_local_molid, new_chain_id)
        end
        n_tmp_res += n_dna_res
    end
    # ---------
    # [ bonds ]
    # ---------
    println("   >> adding bonds")
    # nanoparticle
    r_P = norm(atom_coors[:, 1] - atom_coors[:, 2])
    # @printf(" Radius of nanoparticle: %12.3f \n", r_P)
    # 1: central bead to all the other beads
    for i_bead in 2:1 + n_atom_nanoP
        new_bond = GenTopBond(1, i_bead, AICG_BOND_FUNC_TYPE, r_P, 80 * JOU2CAL)
        push!(top_bonds, new_bond)
    end
    # 2: between every two neighboring beads
    # 2.1: determine the minimum distance
    coor1 = atom_coors[:, 2]
    d_min = 1e10
    for j_nanop in 3:1+n_atom_nanoP
        coor2 = atom_coors[:, j_nanop]
        d_tmp = norm(coor1 - coor2)
        d_min = d_tmp < d_min ? d_tmp : d_min
    end
    @printf("     >>> Distance between every two beads in nanoparticle = %12.3f Å \n", d_min)
    # 2.2: add nanop-nanop bonds
    for j_bead in 2:1+n_atom_nanoP
        coor1 = atom_coors[:, j_bead]
        for k_bead in j_bead + 1:1+n_atom_nanoP
            coor2 = atom_coors[:, k_bead]
            d_tmp = norm(coor1 - coor2)
            if d_tmp < 1.5 * d_min
                new_bond = GenTopBond(j_bead, k_bead, AICG_BOND_FUNC_TYPE, d_tmp, 80 * JOU2CAL)
                push!(top_bonds, new_bond)
            end
        end
    end
    # 3: B-L; L-L; L-L
    n_DNA_connecting_node = 0
    for (i_nanoP, f_nanoP) in enumerate(nanoP_flags)
        if f_nanoP > 0
            n_DNA_connecting_node += 1
            j_shift = 1 + n_atom_nanoP + (n_DNA_connecting_node - 1) * (linker_length + n_atom_DNA)
            # add B-L
            i_bead = 1 + i_nanoP
            j_bead = j_shift + 1
            d_tmp = compute_distance(atom_coors[:, i_bead], atom_coors[:, j_bead])
            new_bond = GenTopBond(i_bead, j_bead, DNA3SPN_BOND_FUNC4_TYPE, d_tmp, 0.6 * JOU2CAL)
            push!(top_bonds, new_bond)
            # add L-L and L-S
            for k in 1:linker_length
                i_bead = j_shift + k
                j_bead = i_bead  + 1
                d_tmp = compute_distance(atom_coors[:, i_bead], atom_coors[:, j_bead])
                new_bond = GenTopBond(i_bead, j_bead, DNA3SPN_BOND_FUNC4_TYPE, d_tmp, 0.6 * JOU2CAL)
                push!(top_bonds, new_bond)
            end
        end
    end
    # 4: DNA bonds
    for i_dna in 1:N_DNA
        i_shift = 1 + n_atom_nanoP + (i_dna - 1) * (linker_length + n_atom_DNA) + linker_length
        n_bond_DNA = length(meta_DNA_top.top_bonds)
        for i_bond in 1:n_bond_DNA
            b = meta_DNA_top.top_bonds[i_bond]
            new_bond = GenTopBond(i_shift + b.i, i_shift + b.j, b.function_type, b.r0, b.coef)
            push!(top_bonds, new_bond)
        end
    end

    # ----------
    # [ angles ]
    # ----------
    println("   >> adding angles")
    # 1: B-B-L; B-L-L; L-L-L; L-L-S; L-S-Base
    n_DNA_connecting_node = 0
    for (i_nanoP, f_nanoP) in enumerate(nanoP_flags)
        if f_nanoP > 0
            n_DNA_connecting_node += 1
            j_shift = 1 + n_atom_nanoP + (n_DNA_connecting_node - 1) * (linker_length + n_atom_DNA)
            # add B-B-L
            i_bead = 1 + i_nanoP
            j_bead = j_shift + 1
            coor1 = atom_coors[:, i_bead]
            for h_bead in 2:1+n_atom_nanoP
                if h_bead == i_bead
                    continue
                end
                coor2 = atom_coors[:, h_bead]
                d_tmp = norm(coor1 - coor2)
                if d_tmp < 1.5 * d_min
                    new_angle = GenTopAngle(h_bead, i_bead, j_bead, DNA3SPN_ANG_FUNC_TYPE, 105.19, 200.0 * JOU2CAL, 0.0)
                    push!(top_angles, new_angle)
                end
            end
            # add B-L-L
            new_angle = GenTopAngle(i_bead, j_bead, j_bead + 1, DNA3SPN_ANG_FUNC_TYPE, 178.0, 200.0 * JOU2CAL, 0.0)
            push!(top_angles, new_angle)
            # add L-L-L
            for k in 1:linker_length - 2
                j_bead = j_shift + k
                new_angle = GenTopAngle(j_bead, j_bead + 1, j_bead + 2, DNA3SPN_ANG_FUNC_TYPE, 176.0, 200.0 * JOU2CAL, 0.0)
                push!(top_angles, new_angle)
            end
            # add L-L-S
            k = linker_length - 1
            j_bead = j_shift + k
            new_angle = GenTopAngle(j_bead, j_bead + 1, j_bead + 2, DNA3SPN_ANG_FUNC_TYPE, 178.0, 200.0 * JOU2CAL, 0.0)
            push!(top_angles, new_angle)
            # add L-S-Base
            k = linker_length
            j_bead = j_shift + k
            new_angle = GenTopAngle(j_bead, j_bead + 1, j_bead + 2, DNA3SPN_ANG_FUNC_TYPE, 103.28, 200.0 * JOU2CAL, 0.0)
            push!(top_angles, new_angle)
            # add L-S-P
            new_angle = GenTopAngle(j_bead, j_bead + 1, j_bead + 3, DNA3SPN_ANG_FUNC_TYPE, 134.03, 200.0 * JOU2CAL, 0.0)
            push!(top_angles, new_angle)
        end
    end
    # 2: DNA
    for i_dna in 1:N_DNA
        i_shift = 1 + n_atom_nanoP + (i_dna - 1) * (linker_length + n_atom_DNA) + linker_length
        n_angle_DNA = length(meta_DNA_top.top_angles)
        for i_angle in 1:n_angle_DNA
            a = meta_DNA_top.top_angles[i_angle]
            new_angle = GenTopAngle(i_shift + a.i, i_shift + a.j, i_shift + a.k, DNA3SPN_ANG_FUNC_TYPE, a.a0, a.coef, 0.0)
            push!(top_angles, new_angle)
        end
    end

    # -------------
    # [ dihedrals ]
    # -------------
    println("   >> adding dihedrals")
    # DNA only...
    for i_dna in 1:N_DNA
        i_shift = 1 + n_atom_nanoP + (i_dna - 1) * (linker_length + n_atom_DNA) + linker_length
        n_dihedral_DNA = length(meta_DNA_top.top_dihedrals)
        for i_dihedral in 1:n_dihedral_DNA
            d = meta_DNA_top.top_dihedrals[i_dihedral]
            new_dihedral = GenTopDihedral(i_shift + d.i, i_shift + d.j, i_shift + d.k, i_shift + d.l, d.function_type,
                                          d.d0, d.coef, d.w, d.n)
            push!(top_dihedrals, new_dihedral)
        end
    end


    # ===================
    # Assemble everything
    # ===================
    println(" ------------------------------------------------------------ ")
    println(" > Output files...")
    if length(out_name) > 0
        mol_name = out_name
    else
        mol_name = @sprintf("NANOPARTICLE_%s", DNA_name)
    end
    mytop = GenTopology(mol_name, n_atom_total,
                        top_default_params,
                        top_default_atomtype,
                        top_default_CGDNA_bp,
                        top_default_CGDNA_bs,
                        top_default_CGDNA_cs,
                        top_default_CGDNA_exv,
                        top_default_CGPro_flx_angle,
                        top_default_CGPro_flx_dihedral,
                        global_index_2_local_index,
                        global_index_2_local_molid,
                        top_atoms,
                        top_bonds,
                        top_angles,
                        top_dihedrals,
                        top_pairs,
                        top_exclusions,
                        top_pwmcos,
                        top_pwmcosns,
                        top_idr_hps,
                        top_idr_kh,
                        top_mol_list)
    myconf = Conformation(n_atom_total, atom_coors)

    # -----------------
    # output some files
    # -----------------
    write_grotop(mytop, mol_name)
    write_grocrd(mytop, myconf, mol_name)
    write_mmCIF(mytop, myconf, mol_name)

    pdb_out_args = Dict("cgconnect"=>true)
    write_pdb(mytop, myconf, mol_name, pdb_out_args)

    # -------------
    # return values
    # -------------

end


if abspath(PROGRAM_FILE) == @__FILE__


    struct_args = ArgParseSettings()

    @add_arg_table struct_args begin
        "--r_P", "-r"
        help     = "Radius of nanoparticle (Å)"
        arg_type = Float64
        default  = 100.0

        "--DNA_density", "-d"
        help     = "Area density of DNA (nm^-2)"
        arg_type = Float64
        default  = 0.2

        "--linker_length", "-l"
        help     = "Length of linker (>=3)"
        arg_type = Int
        default  = 5

        "--DNA_file_name", "-D"
        help     = "File name for DNA topology and coordinates (basename of .top, .gro)"
        arg_type = String
        default  = ""

        "--out_file_name", "-o"
        help     = "File name for your new pet"
        arg_type = String
        default  = ""

        "--log", "-L"
        help = "Debug mode"
        action = :store_true

        "--debug"
        help = "Debug mode"
        action = :store_true
    end

    main_args = parse_args(struct_args)

    r_P = main_args["r_P"]
    DNA_density = main_args["DNA_density"] * 0.01 # unit: Å^-2
    len_linker = main_args["linker_length"]
    DNA_name = get(main_args, "DNA_file_name", "bdna_cg")
    out_name = get(main_args, "out_file_name", "bdna_cg")
    do_log = main_args["log"]
    do_debug = main_args["debug"]

    # =============
    # main function
    # =============
    @printf(" ============================================================ \n")
    @printf(" DNA-nanoparticle generated with GENESIS-cg-tool  \n\n")
    @printf(" Command line parameters:  \n")
    @printf(" Radius of nanoparticle (Å)  > %12.3f \n", r_P)
    @printf(" Area density of DNA (nm^-2) > %12.3f \n", DNA_density * 100)
    @printf(" Length of linker            > %12d   \n", len_linker)
    @printf(" ------------------------------------------------------------ \n")

    (my_polyhedron_coors, my_polyhedron_flags) = subdivided_polyhedron_construction(icosahedron_construction, r_P, DNA_density)
    generate_nanoparticle_with_DNA(my_polyhedron_coors, my_polyhedron_flags, len_linker, DNA_name, out_name, do_debug)

    @printf(" ============================================================ \n")

end
