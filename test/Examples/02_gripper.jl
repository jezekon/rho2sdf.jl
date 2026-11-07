# =============================================================================
# GRIPPER - SDF from the SIMP densities
# =============================================================================
#
# Description:
#   Smooth geometry of the optimized gripper: nodal densities, the SDF of the
#   iso-contour ρ_t = 0.5 on a Cartesian grid, Gaussian RBF interpolation, and the
#   clip by the design domain, min(φ̃ - c, d_D), with the shift c keeping the
#   volume Σ ρ_e v_e. The script checks the result against a reference run.
#
# Input: data/gripper_simp.mat
#   SIMP result of EasySIMP.jl, test/Examples/04_gripper_complex.jl
#   (half model, symmetry plane x = 0, 16,490 hexahedra with a median edge of
#   4 mm, volume fraction 0.3, filter radius 1.5 × element size).
#   Bounding box: x 0..75, y -75..75, z -90..145 mm.
#
# Settings:
#   - Threshold density ρ_t = 0.5
#   - Grid step B = 3.0 mm: h = 3.0128, 155,040 grid nodes. B = 3.0 came from a
#     study of the grid step: the lowest deformation energy and less waviness of
#     the surface than B = 2.0, with three times fewer grid nodes.
#   - RBF interpolation of the SDF of the iso-contour alone, clipped by D
#
# Outputs (in ./results/02_gripper):
#   - gripper_nodal_densities.vtu                  nodal densities
#   - gripper_SDF_CellSize-3.0128_IsoContour.vti   raw SDF of the iso-contour
#   - Z_gripper_RBF_B-3.0128_Interpolation.jld2    RBF field for Implicit2TetMesh.jl
#   - gripper_B-3.0128_RBF_Interpolation.vti       the field sampled at h/2 (ParaView)
#
# Reference run (Julia 1.11.5, 4 threads, 1 BLAS thread):
#   clip shift c = -0.3669154072106465, volume 317,846.6 mm³ (target 317,846.6 mm³)
#
# =============================================================================

using Rho2sdf
using MAT
using LinearAlgebra
using Test

# One BLAS thread as in the reference run
LinearAlgebra.BLAS.set_num_threads(1)

# -----------------------------------------------------------------------------
# 1. SIMP RESULT
# -----------------------------------------------------------------------------
println("Reading the SIMP result...")
data = matread(joinpath(@__DIR__, "data", "gripper_simp.mat"))
(X, IEN, rho) = MeshInformations(data)
println("  ✓ $(length(IEN)) elements, $(length(X)) nodes")

# -----------------------------------------------------------------------------
# 2. OPTIONS
# -----------------------------------------------------------------------------
options = Rho2sdfOptions(
    threshold_density = 0.5,        # iso-contour ρ_t
    grid_step = 3.0,                # grid step B [mm], no prompt
    export_nodal_densities = true,  # nodal densities to ParaView
    export_raw_sdf = true,          # raw SDF of the iso-contour to ParaView
    rbf_smoothing = :interpolation, # Gaussian RBF interpolation of the SDF
    domain_clip = true,             # clip by the design domain, volume kept
    rbf_isocontour_only = true,     # smooth only the iso-contour, the clip draws the faces of D
)

# -----------------------------------------------------------------------------
# 3. RUN
# -----------------------------------------------------------------------------
results_dir = "./results/02_gripper"
mkpath(results_dir)
(sdf_grid, sdf_dists, rbf) = cd(results_dir) do
    rho2sdf("gripper", X, IEN, rho; options = options)
end

println("\n" * "="^80)
println("SDF COMPLETED")
println("="^80)
println("Grid: h = $(sdf_grid.cell_size), $(sdf_grid.ngp) nodes")
println("Target volume: $(rbf.target_volume)")
println("Clip shift: $(rbf.clip_shift)")
println("Clipped volume: $(rbf.clip_volume)")
println("\nResults saved to: $results_dir")
println("="^80)

# -----------------------------------------------------------------------------
# 4. COMPARISON WITH THE REFERENCE RUN
# -----------------------------------------------------------------------------
@testset "Gripper vs. reference run" begin
    @test sdf_grid.ngp == 155040
    @test isapprox(sdf_grid.cell_size, 3.0128205128205128; rtol = 1e-12)
    @test isapprox(rbf.target_volume, 317846.5817959004; rtol = 1e-9)
    @test isapprox(rbf.clip_shift, -0.3669154072106465; rtol = 1e-6)
    @test isapprox(rbf.clip_volume, rbf.target_volume; rtol = 1e-6)
end
