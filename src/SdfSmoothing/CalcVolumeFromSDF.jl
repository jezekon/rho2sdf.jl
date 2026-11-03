"""
    RBFBandVolume(weights, threshold, h, raw_sdf; refine = 8, band_width = 2.5)

Volume of the body {φ̃ ≥ shift} of the RBF field φ̃ (`RBFs_smoothing`) for a given shift.
Usage: `volume_of = RBFBandVolume(...)`, then `volume_of(shift)`.

Evaluating φ̃ is expensive, so it is done once here and every `volume_of(shift)` only compares
the stored values with the shift:
1. The grid cells are sorted by the raw SDF at their corners. Cells deeper than `band_width`·h
   inside count fully, cells deeper than `band_width`·h outside are skipped (the smoothing and
   the shift move the surface by less than h). The remaining cells form the band around the surface.
2. In each band cell, φ̃ is sampled on a sub-grid with spacing h/`refine`.
3. `volume_of(shift)` adds up the inside parts of the sub-cells (`band_volume`).

The body clipped by the design domain has a second constructor (below).

# Arguments
- `weights::Array{Float64,3}`: RBF weights at the grid nodes
- `threshold::Float64`: Kernel cut-off (`RBFs_smoothing`)
- `h::Float64`: Grid cell size
- `raw_sdf::Array{Float64,3}`: Raw SDF at the grid nodes, without placeholders
- `refine::Int`: Number of sub-cells along a cell edge (even)
- `band_width::Float64`: Half-width of the band in cell sizes
"""
struct RBFBandVolume
    h::Float64                                # grid cell size
    refine::Int                               # sub-cells along a cell edge
    n_inside_cells::Int                       # cells fully inside the body
    band_cells::Vector{NTuple{3,Int}}         # grid cells (i, j, k) of the band
    band_samples::Vector{Array{Float64,3}}    # φ̃ on the sub-grid of each band cell
    domain_samples::Vector{Array{Float64,3}}  # d_D on the same sub-grid, NOT_CLIPPED if not clipped
end

# Domain samples of a band cell that the clip does not change
const NOT_CLIPPED = Array{Float64,3}(undef, 0, 0, 0)

"""
    cell_samples(field, cell, refine)

`field` on the sub-grid with spacing h/`refine` of the grid cell `cell` = (i, j, k). `field` takes
the position in grid units (first grid node at 0).
"""
function cell_samples(field::Function, cell::NTuple{3,Int}, refine::Int)
    i, j, k = cell
    samples = Array{Float64,3}(undef, refine + 1, refine + 1, refine + 1)
    for sk = 0:refine
        for sj = 0:refine
            for si = 0:refine
                position = (i - 1 + si / refine, j - 1 + sj / refine, k - 1 + sk / refine)
                samples[si+1, sj+1, sk+1] = field(position)
            end
        end
    end
    return samples
end

function RBFBandVolume(
    weights::Array{Float64,3},
    threshold::Float64,
    h::Float64,
    raw_sdf::Array{Float64,3};
    refine::Int = 8,
    band_width::Float64 = 2.5,
)
    # The coarse sub-grid of the Richardson extrapolation takes every second sample:
    iseven(refine) || error("RBFBandVolume: refine must be even")

    # 1. Sort the grid cells by the raw SDF at their 8 corners:
    nx, ny, nz = size(raw_sdf)
    n_inside_cells = 0
    band_cells = NTuple{3,Int}[]
    for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:nx-1
                corner_min, corner_max = extrema(@view raw_sdf[i:i+1, j:j+1, k:k+1])
                if corner_min > band_width * h
                    n_inside_cells += 1             # deep inside
                elseif corner_max >= -band_width * h
                    push!(band_cells, (i, j, k))    # near the surface
                end                                 # else deep outside
            end
        end
    end

    # 2. φ̃ on the sub-grid of each band cell:
    rbf_field = position -> rbf_value(weights, threshold, position)
    band_samples = Vector{Array{Float64,3}}(undef, length(band_cells))
    @threads :static for b in eachindex(band_cells)
        band_samples[b] = cell_samples(rbf_field, band_cells[b], refine)
    end

    domain_samples = fill(NOT_CLIPPED, length(band_cells))
    return RBFBandVolume(h, refine, n_inside_cells, band_cells, band_samples, domain_samples)
end

"""
    RBFBandVolume(unclipped, weights, threshold, raw_sdf, d_D, origin; band_width = 2.5)

Volume of the body clipped by the design domain, {min(φ̃ - shift, d_D) ≥ 0}, for a given shift.
`unclipped` is the `RBFBandVolume` of the same field without the clip; its samples of φ̃ are reused.

1. d_D is evaluated at the grid nodes. The grid cells are sorted by the raw SDF and d_D at their
   corners: cells deeper than `band_width`·h inside both count fully, cells deeper than
   `band_width`·h outside either are skipped. The remaining cells form the band.
2. In each band cell, φ̃ is taken from `unclipped` (sampled here if it is not there) and d_D is
   sampled on the same sub-grid. Cells whose corners all lie deeper than `band_width`·h inside the
   domain need no d_D: within them d_D > (`band_width` - √3/2)·h, more than φ̃ - shift near the
   surface of the body, so the clip changes nothing there.
3. `volume_of(shift)` adds up the inside parts of the sub-cells (`band_volume`).

# Arguments
- `unclipped::RBFBandVolume`: Volume of the unclipped body (constructor above)
- `weights::Array{Float64,3}`, `threshold::Float64`, `raw_sdf::Array{Float64,3}`: As above
- `d_D::DomainDistance`: Signed distance of the design domain (its cap must exceed `band_width`·h)
- `origin::Point`: Position of the first grid node
- `band_width::Float64`: Half-width of the band in cell sizes
"""
function RBFBandVolume(
    unclipped::RBFBandVolume,
    weights::Array{Float64,3},
    threshold::Float64,
    raw_sdf::Array{Float64,3},
    d_D::DomainDistance,
    origin::Point;
    band_width::Float64 = 2.5,
)
    h, refine = unclipped.h, unclipped.refine
    if d_D.cap <= band_width * h
        error("RBFBandVolume: the cap of d_D must exceed the band half-width $(band_width)·h")
    end

    # 1. d_D at the grid nodes, then the grid cells sorted by the raw SDF and d_D at their corners:
    nx, ny, nz = size(raw_sdf)
    node_distances = domain_sample(d_D, origin, h, (nx, ny, nz))

    n_inside_cells = 0
    band_cells = NTuple{3,Int}[]
    for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:nx-1
                sdf_min, sdf_max = extrema(@view raw_sdf[i:i+1, j:j+1, k:k+1])
                domain_min, domain_max = extrema(@view node_distances[i:i+1, j:j+1, k:k+1])
                if sdf_min > band_width * h && domain_min > band_width * h
                    n_inside_cells += 1             # deep inside the body and the domain
                elseif sdf_max >= -band_width * h && domain_max >= -band_width * h
                    push!(band_cells, (i, j, k))    # near the surface of the body or the domain
                end                                 # else deep outside one of them
            end
        end
    end

    # 2. φ̃ and d_D on the sub-grid of each band cell:
    unclipped_index = Dict(zip(unclipped.band_cells, eachindex(unclipped.band_cells)))
    rbf_field = position -> rbf_value(weights, threshold, position)
    domain_field = position -> d_D(origin + h * Point(position))
    band_samples = Vector{Array{Float64,3}}(undef, length(band_cells))
    domain_samples = Vector{Array{Float64,3}}(undef, length(band_cells))
    @threads :static for b in eachindex(band_cells)
        cell = band_cells[b]
        if haskey(unclipped_index, cell)
            band_samples[b] = unclipped.band_samples[unclipped_index[cell]]
        else
            band_samples[b] = cell_samples(rbf_field, cell, refine)
        end

        i, j, k = cell
        if minimum(@view node_distances[i:i+1, j:j+1, k:k+1]) <= band_width * h
            domain_samples[b] = cell_samples(domain_field, cell, refine)
        else
            domain_samples[b] = NOT_CLIPPED     # deep inside the domain
        end
    end

    return RBFBandVolume(h, refine, n_inside_cells, band_cells, band_samples, domain_samples)
end

"""
    sample_value(samples, distances, i, j, k, shift)

Sample (i, j, k) of a band cell, to be compared with the shift: φ̃, or in a clipped cell
min(φ̃, d_D + shift), which is ≥ shift exactly where min(φ̃ - shift, d_D) ≥ 0.
"""
@inline function sample_value(
    samples::Array{Float64,3},
    distances::Array{Float64,3},
    i::Int,
    j::Int,
    k::Int,
    shift::Float64,
)
    if isempty(distances)
        return samples[i, j, k]
    end
    return min(samples[i, j, k], distances[i, j, k] + shift)
end

"""
    band_volume(V::RBFBandVolume, shift, step)

Volume of {φ̃ ≥ shift}: the cells inside plus the inside parts of the band cells. The band cells
are split into sub-cells whose corners are every `step`-th sample (spacing `step`·h/refine).
A sub-cell with all corners ≥ shift counts fully. If the surface passes through it, the inside
part is integrated by 3×3×3 Gauss quadrature of the trilinear interpolation of the corners.
With the domain clip, the corners of the clipped cells are min(φ̃, d_D + shift) (`sample_value`).
"""
function band_volume(V::RBFBandVolume, shift::Float64, step::Int)
    # Gauss points and weights on [0, 1]:
    gp, gw = FastGaussQuadrature.gausslegendre(3)
    ξ = (gp .+ 1) ./ 2
    w = gw ./ 2

    n_sub = V.refine ÷ step         # sub-cells along a cell edge
    sub_volume = (V.h / n_sub)^3

    cell_volumes = zeros(length(V.band_samples))
    @threads :static for b in eachindex(V.band_samples)
        samples = V.band_samples[b]
        distances = V.domain_samples[b]
        volume = 0.0
        for sk = 0:n_sub-1
            for sj = 0:n_sub-1
                for si = 0:n_sub-1
                    # Corner values of the sub-cell, cxyz = corner at x, y, z ∈ {0, 1}:
                    i, j, k = 1 + si * step, 1 + sj * step, 1 + sk * step
                    c000 = sample_value(samples, distances, i, j, k, shift)
                    c100 = sample_value(samples, distances, i + step, j, k, shift)
                    c010 = sample_value(samples, distances, i, j + step, k, shift)
                    c110 = sample_value(samples, distances, i + step, j + step, k, shift)
                    c001 = sample_value(samples, distances, i, j, k + step, shift)
                    c101 = sample_value(samples, distances, i + step, j, k + step, shift)
                    c011 = sample_value(samples, distances, i, j + step, k + step, shift)
                    c111 = sample_value(samples, distances, i + step, j + step, k + step, shift)

                    corner_min = min(c000, c100, c010, c110, c001, c101, c011, c111)
                    corner_max = max(c000, c100, c010, c110, c001, c101, c011, c111)

                    if corner_min >= shift
                        # Sub-cell fully inside
                        volume += sub_volume
                    elseif corner_max >= shift
                        # Surface passes through the sub-cell: inside fraction by Gauss quadrature
                        fraction = 0.0
                        for kq = 1:3, jq = 1:3, iq = 1:3
                            x, y, z = ξ[iq], ξ[jq], ξ[kq]

                            # Trilinear interpolation of the corner values:
                            c00 = c000 * (1 - x) + c100 * x
                            c10 = c010 * (1 - x) + c110 * x
                            c01 = c001 * (1 - x) + c101 * x
                            c11 = c011 * (1 - x) + c111 * x
                            c0 = c00 * (1 - y) + c10 * y
                            c1 = c01 * (1 - y) + c11 * y
                            value = c0 * (1 - z) + c1 * z

                            if value >= shift
                                fraction += w[iq] * w[jq] * w[kq]
                            end
                        end
                        volume += fraction * sub_volume
                    end
                end
            end
        end
        cell_volumes[b] = volume
    end

    return V.n_inside_cells * V.h^3 + sum(cell_volumes)
end

"""
    (V::RBFBandVolume)(shift) -> (volume, volume_fine, volume_coarse)

Volume of {φ̃ ≥ shift}, or of {min(φ̃ - shift, d_D) ≥ 0} for the clipped body. `volume_fine` and
`volume_coarse` are computed on the sub-grids h/refine and 2h/refine. Their error falls with the
square of the spacing, so the Richardson extrapolation (4·fine - coarse)/3 removes its main part.
"""
function (V::RBFBandVolume)(shift::Float64)
    volume_fine = band_volume(V, shift, 1)
    volume_coarse = band_volume(V, shift, 2)
    volume = (4 * volume_fine - volume_coarse) / 3
    return (volume = volume, volume_fine = volume_fine, volume_coarse = volume_coarse)
end
