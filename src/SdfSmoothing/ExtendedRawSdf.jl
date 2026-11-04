# Raw SDF of the ρₜ iso-contour alone, extended beyond the design domain D.
#
# `evalDistances` counts the faces of D as a part of the surface. Where the body touches a face of
# D, the raw SDF is therefore zero all over the face, with kinks along the edges of the contact. The
# Gaussian RBF field is not monotonic between the nodes, so it rings around this zero. The clip
# min(φ̃ - c, d_D) removes only the positive lobes, and the negative ones stay as pits in the face.
#
# Here the SDF is split into two parts:
# - the smoothed part is the SDF of the iso-contour alone, which crosses the faces of D smoothly;
# - the other part is d_D (DomainDistance.jl), which is not smoothed.
# The body is {min(φ̃ - c, d_D) ≥ 0}: φ̃ smooths only the iso-contour and d_D draws the faces of D
# exactly. The clip removes everything beyond D, so the values there only have to continue those on
# the faces of D. Each node outside D takes the value at its nearest point of the boundary of D,
# i.e. the SDF is constant along the normals of D.

"""
    extended_raw_sdf(mesh, grid, points, ρₙ, ρₜ)

Raw SDF of the ρₜ iso-contour alone, extended beyond the design domain D, the region of `mesh` (see
the top of the file). It is meant for `RBFs_smoothing(...; isocontour_only = true)`.

1. d_D with the cap 8 h. Nodes beyond the cap are not extended: a jump of the interpolated data
   8 nodes away changes the field at D by about 0.44⁸ ≈ 1e-3 of the jump.
2. Nodes outside D within the cap move to their nearest point of the boundary of D.
3. Each point is used once, however many nodes share it.
4. Distances from the points to the iso-contour (`evalDistances` without the faces of D).
5. Signs of the nodes (`Sign_Detection`).
6. Signs of the moved nodes: the density at their point decides.
7. A summary and the signed distances.

# Arguments
- `mesh::Mesh`: Hexahedral mesh of the design domain (the SIMP mesh)
- `grid::Grid`: Grid structure
- `points::Matrix`: Grid nodes (3 × `grid.ngp`, `generateGridPoints`)
- `ρₙ::Vector{Float64}`: Nodal densities
- `ρₜ::Float64`: Threshold density of the iso-contour

# Returns
- Vector{Float64}: Signed distances at the grid nodes, positive inside. Nodes farther than 2.5 h
  from the iso-contour keep the placeholders ±1e10 of `evalDistances` (see `fill_far_distances`)
"""
function extended_raw_sdf(
    mesh::Mesh,
    grid::Grid,
    points::Matrix,
    ρₙ::Vector{Float64},
    ρₜ::Float64,
)
    ngp = size(points, 2)

    # 1. Distance of the design domain up to the reach of the extension:
    d_D = DomainDistance(mesh, 8 * grid.cell_size)

    # 2. Nodes outside D move to their nearest boundary point; `boundary_elements` holds the element
    #    of that point (0 = the node does not move):
    moved_points = copy(points)
    boundary_elements = zeros(Int, ngp)
    @threads for v = 1:ngp
        x = Point(points[1, v], points[2, v], points[3, v])
        if !inside_domain(d_D, x)
            _, point, triangle = nearest_boundary_point(d_D, x)
            if triangle > 0
                moved_points[:, v] = point
                boundary_elements[v] = d_D.triangle_elements[triangle]
            end
        end
    end

    # 3. Each point once, identified by its coordinates rounded to 1e-9 h. Nodes sharing a point
    #    must not differ: the edge of the band of `evalDistances` is sensitive to round-off. Done
    #    serially to be deterministic; the nodes that did not move go first to keep their values.
    moved_nodes = findall(>(0), boundary_elements)
    resolution = 1e-9 * grid.cell_size
    point_numbers = Dict{NTuple{3,Int},Int}()
    point_of_node = zeros(Int, ngp)
    first_nodes = Int[]     # the first node of each point
    for v in [findall(==(0), boundary_elements); moved_nodes]
        key = ntuple(a -> round(Int, moved_points[a, v] / resolution), 3)
        point_of_node[v] = get!(point_numbers, key) do
            push!(first_nodes, v)
            length(first_nodes)
        end
    end
    unique_points = moved_points[:, first_nodes]

    # 4. Distances to the iso-contour, the faces of D are not a part of the surface:
    unique_distances, _ = evalDistances(mesh, grid, unique_points, ρₙ, ρₜ; domain_faces = false)
    dist = unique_distances[point_of_node]

    # 5. `Sign_Detection` finds the candidate elements by the index of a grid node, so it needs the
    #    original nodes, not the moved points:
    signs = Sign_Detection(mesh, grid, points, ρₙ, ρₜ)

    # 6. Moved nodes: the density at their point, in the element of its triangle, once per point:
    point_elements = zeros(Int, length(first_nodes))
    point_elements[point_of_node[moved_nodes]] = boundary_elements[moved_nodes]
    point_signs = zeros(length(first_nodes))
    not_found = zeros(Bool, length(first_nodes))    # points without local coordinates
    @threads for p in findall(>(0), point_elements)
        element_nodes = mesh.IEN[:, point_elements[p]]
        found, ξ = find_local_coordinates(mesh.X[:, element_nodes], unique_points[:, p])
        ρ = mesh.sfce(ξ) ⋅ ρₙ[element_nodes]
        point_signs[p] = found && ρ >= ρₜ ? 1.0 : -1.0  # not found: outside, as in `Sign_Detection`
        not_found[p] = !found
    end
    signs[moved_nodes] = point_signs[point_of_node[moved_nodes]]

    # 7. Summary and the result:
    @printf("Raw SDF of the iso-contour beyond the domain: %d of %d nodes moved to %d points\n",
        length(moved_nodes), ngp, count(>(0), point_elements))
    @printf("  points without local coordinates (taken as outside): %d\n", count(not_found))
    return dist .* signs
end
