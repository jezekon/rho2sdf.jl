## Overview
This package converts density-based topology optimization results (from SIMP method) into smooth geometry defined by Signed Distance Functions (SDFs). The transformation precisely preserves the volume of the original optimized structure while eliminating jagged boundaries and intermediate density regions, making the resulting geometry suitable for manufacturing and further analysis.

<p align="center">
  <img src="doc/Cantilever_beam-RAW.png" width="45%" alt="Raw topology optimization result" />
  <img src="doc/Cantilever_beam-smooth_SDF_Approx.png" width="45%" alt="Smoothed SDF approximation" />
</p>

## Function Input
Main function for transforming TO results into an implicit geometric representation using SDF:
```julia
rho2sdf(taskName, X, IEN, rho; options=Rho2sdfOptions())
```

### Parameters:
- `taskName::String`: Base name for output files
- `X::Vector{Vector{Float64}}`: Mesh node coordinates
- `IEN::Vector{Vector{Int64}}`: Element connectivity (indices of element nodes)
- `rho::Vector{Float64}`: Element densities
- `options::Rho2sdfOptions`: Configuration options (optional)

### Return Value:
- `Tuple`: (sdf_grid, sdf_dists, rbf)
  - `sdf_grid`: grid on which the raw SDF is computed
  - `sdf_dists`: raw SDF at the grid nodes (positive inside)
  - `rbf`: RBF field (weights, kernel, volume shift and diagnostics), `nothing` for `rbf_smoothing=:none`

### Output Files:
Depending on `rbf_smoothing`, one of the following inputs for the mesher ([Implicit2TetMesh.jl](https://github.com/jezekon/Implicit2TetMesh.jl)) is written:

| `rbf_smoothing` | File | Content | Mesher |
|---|---|---|---|
| `:none` | `Z_<taskName>_SDF_B-<h>.jld2`, `Z_<taskName>_Grid_B-<h>.jld2` | raw SDF on the grid (`fine_sdf`, `fine_grid`) | `generate_tetrahedral_mesh(grid_file, sdf_file, ...)` |
| `:interpolation`, `:approximation` | `Z_<taskName>_RBF_B-<h>_<Interpolation\|Approximation>.jld2` | RBF field: weights, origin, `h`, `sigma`, `threshold`, `shift`; with the domain clip also `clip_shift` and the design domain (below) | `load_rbf_source(file)`, then `generate_tetrahedral_mesh(source, ...; dx=source.h)` with `cut_points=:bisection` |

The RBF field is never sampled back on a grid: the mesher evaluates the smooth field wherever it needs it. The file `<taskName>_B-<h>_RBF_<...>.vti` (the field to be meshed, sampled at h/2; clipped with the domain clip) is only for viewing in Paraview.

### Domain clip
The smoothed field φ̃ can reach slightly beyond the design domain D (the input hex mesh) where the body touches the domain boundary. With `domain_clip=true` (default) the body is

    φ(x) = min(φ̃(x) - clip_shift, d_D(x)) ≥ 0,

where d_D is the signed distance of D (positive inside) and `clip_shift` is found so that the volume of this clipped body equals the volume of the SIMP result. d_D is the distance to the boundary faces of the hex mesh (faces of one element only, each split into 4 triangles around its centre), positive inside the elements; it is exact for planar faces.

The JLD2 keeps every key of the unclipped field with its meaning: `shift` is the volume-keeping shift of φ̃ alone, so a reader that does not know the clip meshes the unclipped body. The clip adds:

| Key | Content |
|---|---|
| `domain_clip` | `true` if the clip was used (also written as `false` without it) |
| `clip_shift` | shift of the clipped field φ |
| `clip_volume` | volume of {φ ≥ 0} (equals `target_volume`) |
| `volume_outside` | volume of {φ̃ ≥ shift} outside D, which the clip removes |
| `domain_nodes`, `domain_elements` | the design domain: hex mesh nodes (3 × nodes) and elements (8 × elements, node order of `Mesh`) |

To mesh the clipped body, the mesher has to evaluate min(φ̃ - `clip_shift`, d_D) with d_D built from `domain_nodes` and `domain_elements` as described above. Implicit2TetMesh.jl does this: `load_rbf_source(file)` returns its `ClippedRBFSDF` for a file with `domain_clip = true` (branch `dev`), with the same d_D (a port of `DomainDistance`).

## Rho2sdfOptions

The `Rho2sdfOptions` struct allows for customization of the SDF generation process:

```julia
Rho2sdfOptions(;
    threshold_density=0.5,          # value for isocontour (0, 1)
    sdf_grid_setup=:manual,         # manual/automatic grid setup
    export_nodal_densities=true,    # export nodal field to Paraview
    export_raw_sdf=true,            # export non-smoothed SDF to Paraview
    rbf_smoothing=:interpolation,   # none/interpolation/approximation of SDF using RBFs
    domain_clip=true                # clip the RBF field by the design domain
)
```

### Options:

#### threshold_density::Union{Float64, Nothing}
- Threshold density value in range [0,1] used for isosurface generation
- If `nothing`, it will be automatically calculated from volume fraction

#### sdf_grid_setup::Symbol
- SDF grid step configuration. Valid values:
  - `:automatic`: Uses original mesh properties for grid step
  - `:manual` (default): Uses interactive grid setup allowing user configuration

#### export_nodal_densities::Bool
- Controls whether nodal densities are exported to VTU format
- Default: `false`

#### export_raw_sdf::Bool
- Controls whether raw SDF values are exported to VTI format
- Default: `false`

#### rbf_smoothing::Symbol
- Smoothing of the raw SDF with Gaussian RBFs (σ = grid step) centred at the grid nodes. Valid values:
  - `:none`: No smoothing, the raw SDF on the grid is the output
  - `:interpolation` (default): RBF interpolation, the field passes through the raw SDF at the grid nodes
  - `:approximation`: RBF approximation (normalized Gaussian quasi-interpolation), smoother than the interpolation, does not pass through the nodal values
- For both RBF variants, the field is shifted so that the volume enclosed by its zero level equals the volume of the SIMP result

#### domain_clip::Bool
- Clips the RBF field by the design domain, min(φ̃ - clip_shift, d_D), with `clip_shift` keeping the volume (see [Domain clip](#domain-clip))
- Used by `:interpolation` and `:approximation` only; `:none` exports the raw SDF unchanged
- Default: `true`

## Example Usage

```julia
# Basic usage with default options
result = rho2sdf("chapadlo", X, IEN, rho)

# Custom options
options = Rho2sdfOptions(
    threshold_density=0.5,
    sdf_grid_setup=:manual,
    export_nodal_densities=true,
    export_raw_sdf=true,
    rbf_smoothing=:approximation,
    domain_clip=false               # the old behaviour: no clip
)
result = rho2sdf("chapadlo", X, IEN, rho, options=options)
```

## TODO List

- [ ] Improvements to the calculation of nodal densities
- [ ] Use IPopt library for points where local coordinates were not found
- [ ] Extend implementation to include tetrahedral elements
