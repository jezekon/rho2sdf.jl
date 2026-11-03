"""
    export_sdf_results(sdf, grid, taskName)

Save the raw SDF and the grid to JLD2 files for the mesher (Implicit2TetMesh.jl,
`generate_tetrahedral_mesh(grid_file, sdf_file)`), which interpolates the values trilinearly.

# Arguments
- `sdf::Vector{Float64}`: SDF at the grid nodes, positive inside, without placeholders (`fill_far_distances`)
- `grid::Grid`: Grid on which the SDF is defined
- `taskName::String`: Base name for the output files

# Returns
- Tuple(String, String): Names of the SDF file and of the grid file
"""
function export_sdf_results(sdf::Vector{Float64}, grid::Grid, taskName::String)
    B = round(grid.cell_size, digits=4)
    dims = Tuple(grid.N .+ 1)

    # The variable names fine_sdf and fine_grid are those the mesher reads.
    # SDF values as a 3D array:
    fine_sdf = reshape(sdf, dims)

    # Node coordinates [x, y, z] as a 3D array (same points as generateGridPoints):
    fine_grid = Array{Vector{Float64},3}(undef, dims)
    for k = 0:grid.N[3]
        for j = 0:grid.N[2]
            for i = 0:grid.N[1]
                fine_grid[i+1, j+1, k+1] = grid.AABB_min .+ grid.cell_size .* [i, j, k]
            end
        end
    end

    sdf_file = "Z_$(taskName)_SDF_B-$(B).jld2"
    grid_file = "Z_$(taskName)_Grid_B-$(B).jld2"
    @save sdf_file fine_sdf
    @save grid_file fine_grid

    println("Saved raw SDF to $sdf_file and grid to $grid_file")
    return (sdf_file, grid_file)
end
