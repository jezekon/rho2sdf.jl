# =============================================================================
# CANTILEVER BEAM - SDF from the SIMP densities
# =============================================================================
#
# Description:
#   Smooth geometry of the optimized cantilever: nodal densities, the SDF of the
#   iso-contour ρ_t = 0.5 on a Cartesian grid, Gaussian RBF interpolation, and the
#   clip by the design domain, min(φ̃ - c, d_D), with the shift c keeping the
#   volume Σ ρ_e v_e. The script checks the result against a reference run.
#
# Input: data/cantilever_beam_simp.mat
#   SIMP result of EasySIMP.jl, test/Examples/01_basic_cantilever.jl
#   (60 × 20 × 4 mm, 1 mm hexahedra, 4,800 elements, volume fraction 0.4,
#   filter radius 1.5 × element size).
#
# Settings:
#   - Threshold density ρ_t = 0.5
#   - Grid step B = 1.0 mm (the element size): h = 1.0, 21,708 grid nodes
#   - RBF interpolation of the SDF of the iso-contour alone, clipped by D
#
# Outputs (in ./results/01_cantilever_beam):
#   - cantilever_nodal_densities.vtu               nodal densities
#   - cantilever_SDF_CellSize-1.0_IsoContour.vti   raw SDF of the iso-contour
#   - Z_cantilever_RBF_B-1.0_Interpolation.jld2    RBF field for Implicit2TetMesh.jl
#   - cantilever_B-1.0_RBF_Interpolation.vti       the field sampled at h/2 (ParaView)
#
# Reference run (Julia 1.11.5, 4 threads, 1 BLAS thread):
#   clip shift c = -0.026284727513354193, volume 1920.0 mm³ (target 1920.0 mm³)
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
data = matread(joinpath(@__DIR__, "data", "cantilever_beam_simp.mat"))
(X, IEN, rho) = MeshInformations(data)
println("  ✓ $(length(IEN)) elements, $(length(X)) nodes")

# -----------------------------------------------------------------------------
# 2. OPTIONS
# -----------------------------------------------------------------------------
options = Rho2sdfOptions(
    threshold_density = 0.5,        # iso-contour ρ_t
    grid_step = 1.0,                # grid step B [mm], no prompt
    export_nodal_densities = true,  # nodal densities to ParaView
    export_raw_sdf = true,          # raw SDF of the iso-contour to ParaView
    rbf_smoothing = :interpolation, # Gaussian RBF interpolation of the SDF
    domain_clip = true,             # clip by the design domain, volume kept
    rbf_isocontour_only = true,     # smooth only the iso-contour, the clip draws the faces of D
)

# -----------------------------------------------------------------------------
# 3. RUN
# -----------------------------------------------------------------------------
results_dir = "./results/01_cantilever_beam"
mkpath(results_dir)
(sdf_grid, sdf_dists, rbf) = cd(results_dir) do
    rho2sdf("cantilever", X, IEN, rho; options = options)
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
@testset "Cantilever vs. reference run" begin
    @test sdf_grid.ngp == 21708
    @test sdf_grid.cell_size ≈ 1.0
    @test isapprox(rbf.target_volume, 1919.999999460622; rtol = 1e-9)
    @test isapprox(rbf.clip_shift, -0.026284727513354193; rtol = 1e-6)
    @test isapprox(rbf.clip_volume, rbf.target_volume; rtol = 1e-6)
end
