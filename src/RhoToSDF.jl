"""
    Rho2sdfOptions

Configuration options for density to SDF conversion.

# Fields
- `threshold_density::Union{Float64, Nothing}`: Threshold density in range [0,1]. If `nothing`, calculated automatically.
- `sdf_grid_setup::Symbol`: Grid setup method: `:manual` (the grid step is asked for in the terminal)
  or `:automatic` (the median edge of the input mesh). Not used when `grid_step` is set.
- `grid_step::Union{Float64, Nothing}`: Grid step B; the cell size is max extent / floor(max extent / B),
  as with `:manual`, but without the prompt. Default `nothing`: `sdf_grid_setup` decides.
- `export_nodal_densities::Bool`: Whether to export nodal densities to VTU.
- `export_raw_sdf::Bool`: Whether to export raw SDF to VTI.
- `rbf_smoothing::Symbol`: Smoothing of the raw SDF:
    * `:none`: no smoothing, the raw SDF on the grid goes to the mesher (`export_sdf_results`).
    * `:interpolation`: Gaussian RBF interpolation of the raw SDF (`RBFs_smoothing`).
    * `:approximation`: Gaussian RBF approximation of the raw SDF (`RBFs_smoothing`).
  With RBFs, the field itself (weights, shift) goes to the mesher (`export_rbf_results`).
- `domain_clip::Bool`: Clip the RBF field by the design domain (the input mesh),
  min(φ̃ - clip_shift, d_D), with clip_shift keeping the volume (`RBFs_smoothing`). Default `true`.
  Used by `:interpolation` and `:approximation` only; `:none` exports the raw SDF unchanged.
- `rbf_sharpening::Bool`: One sharpening step of the approximation weights: the field follows the
  raw SDF more closely at thin members and corners, flat parts stay as they are (`RBFs_smoothing`).
  Default `false`. Used by `:approximation` only.
- `rbf_isocontour_only::Bool`: Smooth only the SDF of the iso-contour, continued beyond the design
  domain (`extended_raw_sdf`); the faces of the domain are drawn by the clip, not smoothed, so the
  field leaves no pits in them. `false` smooths the whole raw SDF, which includes the faces of the
  domain (files with the suffix `_WholeSDF`). Default `true`. Used by `:interpolation` and
  `:approximation` only, and it needs `domain_clip`; `:none` exports the raw SDF unchanged.
"""
struct Rho2sdfOptions
    threshold_density::Union{Float64, Nothing}
    sdf_grid_setup::Symbol
    grid_step::Union{Float64, Nothing}
    export_nodal_densities::Bool
    export_raw_sdf::Bool
    rbf_smoothing::Symbol
    domain_clip::Bool
    rbf_sharpening::Bool
    rbf_isocontour_only::Bool

    # Constructor with default values
    function Rho2sdfOptions(;
        threshold_density=nothing,
        sdf_grid_setup=:manual,
        grid_step=nothing,
        export_nodal_densities=false,
        export_raw_sdf=false,
        rbf_smoothing=:interpolation,
        domain_clip=true,
        rbf_sharpening=false,
        rbf_isocontour_only=true
    )
        # Validate threshold density if provided
        if threshold_density !== nothing
            if !(0.0 <= threshold_density <= 1.0)
                @warn "Threshold density $threshold_density is outside the valid range [0.0, 1.0]. Will use automatic calculation instead."
                threshold_density = nothing
            elseif threshold_density == 0.0 || threshold_density == 1.0
                @warn "Using extreme threshold density value: $threshold_density"
            end
        end
        
        # Validate sdf_grid_setup and set default if invalid
        if !(sdf_grid_setup in [:manual, :automatic])
            @warn "Invalid sdf_grid_setup: $sdf_grid_setup. Must be either :manual or :automatic. Using default :manual instead."
            sdf_grid_setup = :manual
        end
        
        # Validate rbf_smoothing and set default if invalid
        if !(rbf_smoothing in [:none, :interpolation, :approximation])
            @warn "Invalid rbf_smoothing: $rbf_smoothing. Must be :none, :interpolation or :approximation. Using default :interpolation instead."
            rbf_smoothing = :interpolation
        end

        # A fixed grid step must be positive
        if grid_step !== nothing && !(grid_step > 0)
            error("Rho2sdfOptions: grid_step must be positive, got $grid_step")
        end

        # Sharpening changes the approximation weights only
        if rbf_sharpening && rbf_smoothing != :approximation
            @warn "rbf_sharpening is used by rbf_smoothing = :approximation only, not by :$rbf_smoothing. Ignoring it."
            rbf_sharpening = false
        end

        # Without the clip nothing would draw the faces of the design domain
        if rbf_isocontour_only && !domain_clip && rbf_smoothing != :none
            @warn "rbf_isocontour_only needs domain_clip = true. Smoothing the whole raw SDF instead."
            rbf_isocontour_only = false
        end

        new(threshold_density, sdf_grid_setup, grid_step === nothing ? nothing : Float64(grid_step), export_nodal_densities, export_raw_sdf, rbf_smoothing, domain_clip, rbf_sharpening, rbf_isocontour_only)
    end
end

"""
    rho2sdf(taskName, X, IEN, rho; options=Rho2sdfOptions())

Convert element densities to Signed Distance Function representation.

# Arguments
- `taskName::String`: Base name for output files
- `X::Vector{Vector{Float64}}`: Mesh node coordinates
- `IEN::Vector{Vector{Int64}}`: Element connectivity
- `rho::Vector{Float64}}`: Element densities
- `options::Rho2sdfOptions`: Configuration options (optional)

# Returns
- `Tuple`: (sdf_grid, sdf_dists, rbf). `sdf_dists` is the raw SDF at the grid nodes, `rbf` the
  RBF field of `RBFs_smoothing` (`nothing` for `rbf_smoothing = :none`). `sdf_dists`, and the file
  of `export_raw_sdf`, hold what went to the RBFs: with `rbf_isocontour_only`, the SDF of the
  iso-contour alone (`extended_raw_sdf`, file with the suffix `_IsoContour`).

# Example
```julia
# Basic usage with default options
result = rho2sdf("chapadlo", X, IEN, rho)

# Custom options
options = Rho2sdfOptions(
    threshold_density=0.5,
    grid_step=1.0,
    export_nodal_densities=true,
    rbf_smoothing=:approximation,
    domain_clip=true,
    rbf_isocontour_only=true
)
result = rho2sdf("chapadlo", X, IEN, rho, options=options)
"""

function rho2sdf(taskName::String, X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64};
    options::Rho2sdfOptions=Rho2sdfOptions())
    # Generate FEM mesh structure
    mesh = Mesh(X, IEN, rho, hex8_shape)
    
    # Setup SDF grid: a fixed grid step, or the selected method
    sdf_grid = if options.grid_step !== nothing
        noninteractive_sdf_grid_setup(mesh, options.grid_step)
    elseif options.sdf_grid_setup == :manual
        interactive_sdf_grid_setup(mesh)
    else
        # Use automatic setup with default cell size factor
        noninteractive_sdf_grid_setup(mesh)
    end

    # Map elemental densities to the nodes using least squares
    ρₙ = DenseInNodes(mesh, rho)
    
    # Determine threshold density (ρₜ)
    ρₜ = if options.threshold_density === nothing
        # Automatically calculate threshold based on volume
        find_threshold_for_volume(mesh, ρₙ)
    else
        options.threshold_density
    end
    
    # Export nodal densities if requested
    if options.export_nodal_densities
        VTK_CODE = 12  # Hexahedron code per VTK specification
        exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)
    end
    
    # Generate grid points for the SDF calculation
    points = generateGridPoints(sdf_grid)
    
    # Calculate SDF from densities. By default the RBFs smooth only the SDF of the iso-contour and
    # the clip draws the faces of the design domain
    isocontour_only = options.rbf_smoothing != :none && options.rbf_isocontour_only
    if isocontour_only
        sdf_dists = extended_raw_sdf(mesh, sdf_grid, points, ρₙ, ρₜ)
    else
        (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
        signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
        sdf_dists = dists .* signs  # Combine distances with signs to create signed distance field
    end

    # Export raw SDF to VTK if requested; the SDF of the iso-contour goes to a file of its own
    if options.export_raw_sdf
        B = round(sdf_grid.cell_size, digits=4)
        suffix = isocontour_only ? "_IsoContour" : ""
        exportSdfToVTI(taskName * "_SDF_CellSize-" * string(B) * suffix * ".vti", sdf_grid, sdf_dists, "distance")
    end
    
    # Raw SDF only: exported on the grid for the mesher
    if options.rbf_smoothing == :none
        export_sdf_results(fill_far_distances(sdf_dists), sdf_grid, taskName)
        return (sdf_grid, sdf_dists, nothing)
    end

    # RBF smoothing: the RBF field itself is exported and evaluated in the mesher
    rbf = RBFs_smoothing(mesh, sdf_dists, sdf_grid, options.rbf_smoothing == :interpolation;
        domain_clip = options.domain_clip, sharpening = options.rbf_sharpening, isocontour_only)
    export_rbf_results(rbf, sdf_grid, taskName)

    # Return all relevant results
    return (sdf_grid, sdf_dists, rbf)

end


