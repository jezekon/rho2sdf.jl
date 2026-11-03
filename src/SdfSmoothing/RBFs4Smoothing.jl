# RBF smoothing of the raw SDF.
#
# The smoothed SDF is a sum of Gaussians centred at the grid nodes x_i:
#     φ̃(x) = Σ_i w_i exp(-(|x - x_i| / h)²),   h = grid.cell_size,
# where kernel values ≤ `threshold` are dropped (the kernel reaches √(-ln threshold)·h ≈ 2.6 h).
# The body is {φ̃ ≥ shift}, positive inside as in the rest of rho2sdf.
#
# - Interpolation: the weights solve K w = φ, so φ̃ equals the raw SDF φ at the grid nodes.
# - Approximation: w = φ / S, S the sum of the kernel over the grid. φ̃ is smoother but does not
#   pass through φ; it keeps the scale of a distance.
#
# The field is not sampled back on a grid: its definition (weights, grid, kernel, shift) is
# exported and the mesher (Implicit2TetMesh.jl) evaluates it wherever it needs.
#
# Positions in this file are in grid units: node (i, j, k) (0-based) lies at (i, j, k), so the
# kernel is exp(-(Δi² + Δj² + Δk²)). It factorizes into exp(-Δi²)·exp(-Δj²)·exp(-Δk²).

"""
    kernel_range(position, reach, n)

Indices (0-based) of the grid nodes along one axis whose kernel can reach `position`
(|index - position| ≤ `reach`), limited to the grid with `n` nodes.
"""
function kernel_range(position::Float64, reach::Float64, n::Int)
    return max(0, ceil(Int, position - reach)):min(n - 1, floor(Int, position + reach))
end

"""
    kernel_offsets(threshold)

Offsets (Δi, Δj, Δk) between two grid nodes for which the kernel exceeds `threshold`, and the
kernel values. Ordered with Δk outermost and Δi innermost, i.e. by the linear node index.
"""
function kernel_offsets(threshold::Float64)
    max_offset = floor(Int, sqrt(-log(threshold)))
    offsets = NTuple{3,Int}[]
    kernel_values = Float64[]
    for Δk = -max_offset:max_offset
        for Δj = -max_offset:max_offset
            for Δi = -max_offset:max_offset
                kernel_value = exp(-(Δi^2 + Δj^2 + Δk^2))
                if kernel_value > threshold
                    push!(offsets, (Δi, Δj, Δk))
                    push!(kernel_values, kernel_value)
                end
            end
        end
    end
    return offsets, kernel_values
end

"""
    rbf_kernel_matrix(dims, threshold)

Sparse kernel matrix K_ij = exp(-|x_i - x_j|²/h²) of the grid nodes, values ≤ `threshold`
dropped. All nodes have the same neighbours (`kernel_offsets`), so each column is filled
directly from the offsets.
"""
function rbf_kernel_matrix(dims::NTuple{3,Int}, threshold::Float64)
    nx, ny, nz = dims
    n = nx * ny * nz
    offsets, kernel_values = kernel_offsets(threshold)

    # Compressed sparse column format: the rows and values of column `col` are
    # rowval[colptr[col]:colptr[col+1]-1] and nzval[colptr[col]:colptr[col+1]-1]
    colptr = Vector{Int}(undef, n + 1)
    rowval = Int[]
    nzval = Float64[]
    sizehint!(rowval, n * length(offsets))
    sizehint!(nzval, n * length(offsets))

    colptr[1] = 1
    for k = 0:nz-1
        for j = 0:ny-1
            for i = 0:nx-1
                col = 1 + i + nx * (j + ny * k)
                # Neighbours inside the grid (the offsets are ordered, so the rows are sorted):
                for (offset, kernel_value) in zip(offsets, kernel_values)
                    ni, nj, nk = i + offset[1], j + offset[2], k + offset[3]
                    if 0 <= ni < nx && 0 <= nj < ny && 0 <= nk < nz
                        push!(rowval, 1 + ni + nx * (nj + ny * nk))
                        push!(nzval, kernel_value)
                    end
                end
                colptr[col+1] = length(rowval) + 1
            end
        end
    end
    return SparseMatrixCSC(n, n, colptr, rowval, nzval)
end

"""
    rbf_value(weights, threshold, position)

RBF sum Σ w_i exp(-|x - x_i|²/h²) at `position` (grid units, see the top of the file).
"""
function rbf_value(weights::Array{Float64,3}, threshold::Float64, position::NTuple{3,Float64})
    nx, ny, nz = size(weights)
    x, y, z = position
    reach = sqrt(-log(threshold))   # kernel reach in cell sizes

    i_range = kernel_range(x, reach, nx)
    j_range = kernel_range(y, reach, ny)
    k_range = kernel_range(z, reach, nz)

    value = 0.0
    @inbounds for k in k_range
        gz = exp(-(k - z)^2)
        for j in j_range
            gyz = gz * exp(-(j - y)^2)
            gyz > threshold || continue     # all i would be below the threshold
            for i in i_range
                g = gyz * exp(-(i - x)^2)
                if g > threshold
                    value += weights[i+1, j+1, k+1] * g
                end
            end
        end
    end
    return value
end

"""
    rbf_sample(weights, threshold, refine)

RBF sum (`rbf_value`) at the nodes of the grid refined `refine` times (spacing h/refine).
Returns an array of size `refine .* (size(weights) .- 1) .+ 1`. Faster than calling `rbf_value`
at every node, because the 1D kernel values exp(-Δ²) are shared along the grid lines.
"""
function rbf_sample(weights::Array{Float64,3}, threshold::Float64, refine::Int)
    nx, ny, nz = size(weights)
    mx, my, mz = refine * (nx - 1) + 1, refine * (ny - 1) + 1, refine * (nz - 1) + 1
    reach = sqrt(-log(threshold))   # kernel reach in cell sizes
    values = zeros(mx, my, mz)

    @threads :static for sk = 1:mz
        # 1D kernel values along x, y, z (only the entries in the kernel range are used)
        gx = zeros(nx)
        gy = zeros(ny)
        gz = zeros(nz)

        z = (sk - 1) / refine
        k_range = kernel_range(z, reach, nz)
        for k in k_range
            gz[k+1] = exp(-(k - z)^2)
        end
        for sj = 1:my
            y = (sj - 1) / refine
            j_range = kernel_range(y, reach, ny)
            for j in j_range
                gy[j+1] = exp(-(j - y)^2)
            end
            for si = 1:mx
                x = (si - 1) / refine
                i_range = kernel_range(x, reach, nx)
                for i in i_range
                    gx[i+1] = exp(-(i - x)^2)
                end

                value = 0.0
                @inbounds for k in k_range
                    for j in j_range
                        gyz = gy[j+1] * gz[k+1]
                        gyz > threshold || continue
                        for i in i_range
                            g = gx[i+1] * gyz
                            if g > threshold
                                value += weights[i+1, j+1, k+1] * g
                            end
                        end
                    end
                end
                values[si, sj, sk] = value
            end
        end
    end
    return values
end

"""
    RBFs_smoothing(mesh::Mesh, dist::Vector, grid::Grid, is_interpolation::Bool;
                   threshold = 1e-3, volume_refine = 8, cg_reltol = 1e-12)

Smooth the raw SDF by Gaussian RBFs (see the top of the file) and find the shift that keeps the
volume of the SIMP result.

1. Placeholders of the far nodes are replaced (`fill_far_distances`).
2. Weights: interpolation solves K w = φ by CG (relative tolerance `cg_reltol`), approximation
   takes w = φ / S.
3. The shift is found by bisection so that the volume of {φ̃ ≥ shift} (`RBFBandVolume`) equals
   mesh.V_frac·mesh.V_domain.

# Arguments
- `mesh::Mesh`: Original mesh with known volume
- `dist::Vector`: Raw SDF at the grid nodes (positive inside, placeholders ±1e10)
- `grid::Grid`: Grid structure
- `is_interpolation::Bool`: Interpolation (true) or approximation (false)

# Returns
- NamedTuple: the RBF field (weights, origin, h, sigma, threshold, shift, interpolation) and
  checks (target and achieved volume, deviation from the raw SDF near the surface)
"""
function RBFs_smoothing(
    mesh::Mesh,
    dist::Vector,
    grid::Grid,
    is_interpolation::Bool;
    threshold::Float64 = 1e-3,
    volume_refine::Int = 8,
    cg_reltol::Float64 = 1e-12,
)
    name = is_interpolation ? "Interpolation" : "Approximation"
    dims = Tuple(grid.N .+ 1)
    h = grid.cell_size

    # 1. Raw SDF without placeholders:
    sdf = fill_far_distances(Float64.(dist))

    # 2. RBF weights:
    println("Computing RBF weights ($name)...")
    if is_interpolation
        K = rbf_kernel_matrix(dims, threshold)
        w = cg(K, sdf; reltol = cg_reltol, maxiter = 10_000)
    else
        _, kernel_values = kernel_offsets(threshold)
        w = sdf ./ sum(kernel_values)
    end
    weights = reshape(w, dims)

    # Check: deviation of φ̃ from the raw SDF at the nodes near the surface
    rbf_at_nodes = vec(rbf_sample(weights, threshold, 1))
    near_surface = abs.(sdf) .<= 2.5 * h
    nodal_error = maximum(abs.(rbf_at_nodes[near_surface] .- sdf[near_surface]))

    # 3. Shift that keeps the volume, found by bisection:
    println("Computing the shift of the RBF field to meet the volume condition...")
    target_volume = mesh.V_frac * mesh.V_domain
    volume_of = RBFBandVolume(weights, threshold, h, reshape(sdf, dims); refine = volume_refine)

    # The volume decreases as the shift grows
    shift_low, shift_high = -h, h
    if !(volume_of(shift_low).volume >= target_volume >= volume_of(shift_high).volume)
        error("RBFs_smoothing: target volume $target_volume outside the bracket of the shift")
    end

    shift = 0.0
    achieved = (volume = 0.0, volume_fine = 0.0, volume_coarse = 0.0)
    for _ = 1:60
        shift = (shift_low + shift_high) / 2
        achieved = volume_of(shift)
        if abs(achieved.volume - target_volume) <= 1e-9 * target_volume
            break
        end
        if achieved.volume > target_volume
            shift_low = shift
        else
            shift_high = shift
        end
    end

    # The volume is piecewise constant in the shift, so 1e-9 is not always reached
    relative_error = abs(achieved.volume - target_volume) / target_volume
    if relative_error > 1e-6
        @warn @sprintf("RBF volume shift: volume %.6f deviates from the target %.6f by %.1e (relative)",
            achieved.volume, target_volume, relative_error)
    end

    @printf("RBF %s: %d nodes, h = %.4f\n", lowercase(name), prod(dims), h)
    @printf("  max deviation from the raw SDF near the surface: %.2e\n", nodal_error)
    @printf("  shift = %.6f, volume %.4f (target %.4f; sub-grid h/%d: %.4f, h/%d: %.4f)\n",
        shift, achieved.volume, target_volume,
        volume_refine, achieved.volume_fine, volume_refine ÷ 2, achieved.volume_coarse)

    return (
        weights = weights,
        origin = Float64.(grid.AABB_min),
        h = h,
        sigma = h,
        threshold = threshold,
        shift = shift,
        interpolation = is_interpolation,
        target_volume = target_volume,
        volume = achieved.volume,
        volume_sampled = (achieved.volume_fine, achieved.volume_coarse),
        volume_refine = volume_refine,
        nodal_error = nodal_error,
    )
end

"""
    export_rbf_results(rbf, grid, taskName)

Save the RBF field of `RBFs_smoothing` for the mesher (`load_rbf_source` of Implicit2TetMesh.jl)
to `Z_<taskName>_RBF_B-<h>_<Interpolation|Approximation>.jld2`, and φ̃ - shift sampled at h/2 to
a VTI for viewing in ParaView.

# Returns
- String: Name of the JLD2 file
"""
function export_rbf_results(rbf, grid::Grid, taskName::String)
    name = rbf.interpolation ? "Interpolation" : "Approximation"
    B = round(rbf.h, digits = 4)

    # The variable names are those the mesher reads
    file = "Z_$(taskName)_RBF_B-$(B)_$(name).jld2"
    jldsave(file;
        weights = rbf.weights,
        origin = rbf.origin,
        h = rbf.h,
        sigma = rbf.sigma,
        threshold = rbf.threshold,
        shift = rbf.shift,
        inside = "positive",
        smoothing = lowercase(name),
        target_volume = rbf.target_volume,
        volume = rbf.volume,
    )

    # Picture of the field for ParaView (not used by the mesher)
    field = rbf_sample(rbf.weights, rbf.threshold, 2) .- rbf.shift
    exportSdfToVTI("$(taskName)_B-$(B)_RBF_$(name).vti", grid, vec(field), "distance", 2)

    println("Saved RBF definition to $file")
    return file
end
