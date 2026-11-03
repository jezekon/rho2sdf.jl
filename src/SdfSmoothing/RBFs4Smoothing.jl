# RBF smoothing of the raw SDF.
#
# The smoothed SDF is a sum of Gaussians centred at the grid nodes x_i:
#     φ̃(x) = Σ_i w_i exp(-(|x - x_i| / h)²),   h = grid.cell_size,
# where kernel values ≤ `threshold` are dropped (the kernel reaches √(-ln threshold)·h ≈ 2.6 h).
# The body is {φ̃ ≥ shift}, positive inside as in the rest of rho2sdf.
#
# φ̃ can reach slightly beyond the design domain D (the SIMP mesh). With the domain clip (default)
# the body is {φ ≥ 0}, φ = min(φ̃ - clip_shift, d_D), d_D the signed distance of D
# (DomainDistance.jl), and clip_shift keeps the volume of this clipped body.
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
    find_shift(volume_of, target_volume, h, body)

Shift in [-h, h] for which the volume `volume_of(shift)` (`RBFBandVolume`) equals the target,
found by bisection (the volume decreases as the shift grows). `body` names the body in the
messages.

# Returns
- Tuple: the shift and `volume_of(shift)`
"""
function find_shift(volume_of::RBFBandVolume, target_volume::Float64, h::Float64, body::String)
    shift_low, shift_high = -h, h
    if !(volume_of(shift_low).volume >= target_volume >= volume_of(shift_high).volume)
        error("RBFs_smoothing ($body): target volume $target_volume outside the bracket of the shift")
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
        @warn @sprintf(
            "RBF volume shift (%s): volume %.6f deviates from the target %.6f by %.1e (relative)",
            body, achieved.volume, target_volume, relative_error)
    end
    return shift, achieved
end

"""
    RBFs_smoothing(mesh::Mesh, dist::Vector, grid::Grid, is_interpolation::Bool;
                   threshold = 1e-3, volume_refine = 8, cg_reltol = 1e-12, domain_clip = true)

Smooth the raw SDF by Gaussian RBFs (see the top of the file) and find the shift that keeps the
volume of the SIMP result.

1. Placeholders of the far nodes are replaced (`fill_far_distances`).
2. Weights: interpolation solves K w = φ by CG (relative tolerance `cg_reltol`), approximation
   takes w = φ / S.
3. The shift is found by bisection (`find_shift`) so that the volume of {φ̃ ≥ shift}
   (`RBFBandVolume`) equals mesh.V_frac·mesh.V_domain.
4. Domain clip (`domain_clip = true`): the body is clipped by the design domain, the region of
   `mesh`: φ = min(φ̃ - clip_shift, d_D) (`DomainDistance`). clip_shift is found as in step 3 so
   that the volume of {φ ≥ 0} equals the same target. The weights and the shift of step 3 do not
   change; the volume of {φ̃ ≥ shift} outside the domain is reported.

# Arguments
- `mesh::Mesh`: Original mesh with known volume, also the design domain of the clip
- `dist::Vector`: Raw SDF at the grid nodes (positive inside, placeholders ±1e10)
- `grid::Grid`: Grid structure
- `is_interpolation::Bool`: Interpolation (true) or approximation (false)
- `domain_clip::Bool`: Clip the body by the design domain (step 4)

# Returns
- NamedTuple: the RBF field (weights, origin, h, sigma, threshold, shift, interpolation), checks
  (target and achieved volume, deviation from the raw SDF near the surface) and the clip
  (domain_clip, clip_shift, clip_volume, volume_outside, domain, domain_time; `nothing` without
  the clip)
"""
function RBFs_smoothing(
    mesh::Mesh,
    dist::Vector,
    grid::Grid,
    is_interpolation::Bool;
    threshold::Float64 = 1e-3,
    volume_refine::Int = 8,
    cg_reltol::Float64 = 1e-12,
    domain_clip::Bool = true,
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
    raw_sdf = reshape(sdf, dims)
    band_width = 2.5        # half-width of the band of RBFBandVolume in cell sizes
    volume_of = RBFBandVolume(weights, threshold, h, raw_sdf; refine = volume_refine, band_width)
    shift, achieved = find_shift(volume_of, target_volume, h, "unclipped body")

    # 4. Clip by the design domain, φ = min(φ̃ - clip_shift, d_D), and the shift of the clipped
    #    body:
    if domain_clip
        println("Clipping the RBF field by the design domain...")
        # |d_D| is needed exactly up to the band half-width, beyond it only its sign
        domain_time = @elapsed d_D = DomainDistance(mesh, (band_width + 1) * h)
        origin = Point(grid.AABB_min...)
        clipped_volume_of = RBFBandVolume(volume_of, weights, threshold, raw_sdf, d_D, origin;
            band_width)
        clip_shift, clip_achieved = find_shift(clipped_volume_of, target_volume, h, "clipped body")

        # Part of {φ̃ ≥ shift} outside the domain, removed by the clip:
        volume_outside = achieved.volume - clipped_volume_of(shift).volume
    else
        d_D = domain_time = clip_shift = clip_achieved = volume_outside = nothing
    end

    @printf("RBF %s: %d nodes, h = %.4f\n", lowercase(name), prod(dims), h)
    @printf("  max deviation from the raw SDF near the surface: %.2e\n", nodal_error)
    @printf("  shift = %.6f, volume %.4f (target %.4f; sub-grid h/%d: %.4f, h/%d: %.4f)\n",
        shift, achieved.volume, target_volume,
        volume_refine, achieved.volume_fine, volume_refine ÷ 2, achieved.volume_coarse)
    if domain_clip
        @printf("  domain clip (%d boundary triangles, set-up %.2f s):\n",
            length(d_D.triangles), domain_time)
        @printf("  clip_shift = %.6f, volume %.4f (sub-grid h/%d: %.4f, h/%d: %.4f)\n",
            clip_shift, clip_achieved.volume,
            volume_refine, clip_achieved.volume_fine, volume_refine ÷ 2, clip_achieved.volume_coarse)
        @printf("  volume of {φ̃ ≥ shift} outside the domain: %.4f (%.3f %% of the target)\n",
            volume_outside, 100 * volume_outside / target_volume)
    end

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
        domain_clip = domain_clip,
        clip_shift = clip_shift,
        clip_volume = domain_clip ? clip_achieved.volume : nothing,
        clip_volume_sampled =
            domain_clip ? (clip_achieved.volume_fine, clip_achieved.volume_coarse) : nothing,
        volume_outside = volume_outside,
        domain = d_D,
        domain_time = domain_time,
    )
end

"""
    export_rbf_results(rbf, grid, taskName)

Save the RBF field of `RBFs_smoothing` for the mesher (`load_rbf_source` of Implicit2TetMesh.jl)
to `Z_<taskName>_RBF_B-<h>_<Interpolation|Approximation>.jld2`, and the field to be meshed,
sampled at h/2, to a VTI for viewing in ParaView: φ̃ - shift, or min(φ̃ - clip_shift, d_D) with the
domain clip.

The keys of the unclipped field keep their meaning with the clip: `shift` is the shift of φ̃
alone, so a reader that does not know the clip meshes the unclipped body, with its own volume.
The clip adds `clip_shift`, `clip_volume`, `volume_outside` and the design domain:
`domain_nodes` (3 × nodes) and `domain_elements` (8 × elements, node order of `Mesh`), from which
d_D is built as in `DomainDistance`. `domain_clip` tells whether the clip was used.

# Returns
- String: Name of the JLD2 file
"""
function export_rbf_results(rbf, grid::Grid, taskName::String)
    name = rbf.interpolation ? "Interpolation" : "Approximation"
    B = round(rbf.h, digits = 4)

    # The variable names are those the mesher reads
    file = "Z_$(taskName)_RBF_B-$(B)_$(name).jld2"
    fields = (
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
        domain_clip = rbf.domain_clip,
    )
    if rbf.domain_clip
        clip_fields = (
            clip_shift = rbf.clip_shift,
            clip_volume = rbf.clip_volume,
            volume_outside = rbf.volume_outside,
            domain_nodes = rbf.domain.nodes,
            domain_elements = rbf.domain.elements,
        )
        fields = merge(fields, clip_fields)
    end
    jldsave(file; fields...)

    # Picture of the field to be meshed, for ParaView (not used by the mesher)
    field = rbf_sample(rbf.weights, rbf.threshold, 2)
    if rbf.domain_clip
        domain_field = domain_sample(rbf.domain, Point(rbf.origin...), rbf.h / 2, size(field))
        field = min.(field .- rbf.clip_shift, domain_field)
    else
        field = field .- rbf.shift
    end
    exportSdfToVTI("$(taskName)_B-$(B)_RBF_$(name).vti", grid, vec(field), "distance", 2)

    println("Saved RBF definition to $file")
    return file
end
