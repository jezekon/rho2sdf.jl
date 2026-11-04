# using Base: print_without_params
using Test
using Rho2sdf
using Rho2sdf.TerminalUtils
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using Rho2sdf.SdfSmoothing: DomainDistance, Point, rbf_value, rbf_kernel_matrix, kernel_offsets,
  nearest_boundary_point, inside_domain
using MAT
using JLD2
using LinearAlgebra
using BenchmarkTools

@testset "Rho2sdf.jl" begin

  # @time @testset "PrimitiveGeometriesTest" begin include("PrimitiveGeometriesTest/runtests.jl") end
  # @time @testset "MeshGridTest" begin include("MeshGridTest/runtests.jl") end
  # @time @testset "SignedDistancesTest" begin include("SignedDistancesTest/runtests.jl") end
  # end
  # exit()
  #     
  # # Data from Matlab:
  # taskName = "chapadlo"

  RUN_PLANE = false
  RUN_BLOCK = false
  RUN_SPHERE = false
  RUN_BEAM = false
  RUN_CHAPADLO = false
  RUN_CHAPADLO_cele = false
  RUN_RBF_SMOOTHING = true

  if (RUN_RBF_SMOOTHING)
    @testset "RBF smoothing" begin
      # Design domain: the block [-1, 1] × [-0.5, 0.5]² of 8 × 4 × 4 elements. Body: the bar
      # [-1, 0] × [-0.25, 0.25]², which touches the domain face x = -1.
      X, IEN, _, _ = selectPrimitiveGeometry("block", [8, 4, 4])
      centres = [sum(X[element]) / 8 for element in IEN]
      rho = [c[1] < 0 && abs(c[2]) < 0.25 && abs(c[3]) < 0.25 ? 1.0 : 0.0 for c in centres]
      mesh = Mesh(X, IEN, rho, hex8_shape)

      # Distance from x to the axis-aligned box or rectangle [low, high] (zero inside)
      box_gap(x, low, high) = norm(max.(low .- x, 0.0, x .- high))

      @testset "Signed distance of the domain" begin
        # d_D against the exact signed distance of the block at a lattice of points inside and
        # outside it, including points on element faces (x = 0, y = 0, z = 0) and domain faces
        # (z = ±0.5)
        cap = 0.4
        d_D = DomainDistance(mesh, cap)
        low, high = [-1.0, -0.5, -0.5], [1.0, 0.5, 0.5]
        largest_error = 0.0
        for x in range(-1.3, 1.3, 23)
          for y in range(-0.8, 0.8, 15)
            for z in (-0.7, -0.5, -0.31, 0.0, 0.17, 0.5, 0.66)
              p = [x, y, z]
              if box_gap(p, low, high) > 0
                exact = -box_gap(p, low, high)
              else
                exact = minimum(min.(p .- low, high .- p))
              end
              largest_error = max(largest_error, abs(d_D(Point(p...)) - clamp(exact, -cap, cap)))
            end
          end
        end
        @test largest_error <= 1e-12
      end

      # Analytic raw SDF of the bar: distance to its faces without the end x = -1 on the domain
      # boundary (`evalDistances` counts that one too), positive in the bar. It jumps from + to -
      # across x = -1, so the RBF field reaches beyond the domain face.
      iso_faces = (([0.0, -0.25, -0.25], [0.0, 0.25, 0.25]),
                   ([-1.0, -0.25, -0.25], [0.0, -0.25, 0.25]),
                   ([-1.0, 0.25, -0.25], [0.0, 0.25, 0.25]),
                   ([-1.0, -0.25, -0.25], [0.0, 0.25, -0.25]),
                   ([-1.0, -0.25, 0.25], [0.0, 0.25, 0.25]))
      function bar_sdf(x)
        distance = minimum(box_gap(x, low, high) for (low, high) in iso_faces)
        in_bar = box_gap(x, [-1.0, -0.25, -0.25], [0.0, 0.25, 0.25]) == 0
        return in_bar ? distance : -distance
      end

      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, 16, 3)
      points = generateGridPoints(sdf_grid)
      raw_sdf = [bar_sdf(points[:, n]) for n in 1:sdf_grid.ngp]

      @testset "Clip keeps the volume" begin
        # Without the clip, rho2sdf ca529c2 gave (shift, volume, sum of the weights). The last
        # bits differ between machines, so they are compared with tolerances; the volume only to
        # 1e-5, as the bisection stops anywhere within the tolerance of the volume condition.
        reference = Dict(true => (0.008666676631240672, 0.25000012558674123, -721.699707129952),
                         false => (-0.014575590606290718, 0.24999950463007617, -616.4109287386418))

        for is_interpolation in (true, false)
          clipped = RBFs_smoothing(mesh, raw_sdf, sdf_grid, is_interpolation)
          unclipped = RBFs_smoothing(mesh, raw_sdf, sdf_grid, is_interpolation; domain_clip = false)

          # Without the clip: the result of ca529c2. With it: the same weights and shift of φ̃
          reference_shift, reference_volume, reference_weights = reference[is_interpolation]
          @test isapprox(unclipped.shift, reference_shift; atol = 1e-10)
          @test isapprox(unclipped.volume, reference_volume; rtol = 1e-5)
          @test isapprox(sum(unclipped.weights), reference_weights; rtol = 1e-12)
          @test clipped.weights == unclipped.weights
          @test (clipped.shift, clipped.volume) == (unclipped.shift, unclipped.volume)

          # The clipped body keeps the volume. The volume is piecewise constant in the shift, with
          # steps of about 6e-6 on this coarse grid, and the bisection cannot get below them.
          @test abs(clipped.clip_volume / clipped.target_volume - 1) <= 1e-5
          @test clipped.volume_outside > 0

          # Just beyond the face x = -1, φ̃ - shift > 0 at some points (the overshoot) while the
          # clipped field φ = min(φ̃ - clip_shift, d_D) is ≤ 0 at all of them
          beyond_face = [Point(-1.02, y, z) for y in range(-0.2, 0.2, 9) for z in range(-0.2, 0.2, 9)]
          grid_position(x) = Tuple((x .- clipped.origin) ./ clipped.h)
          rbf_field = [rbf_value(clipped.weights, clipped.threshold, grid_position(x))
                       for x in beyond_face]
          clipped_field = [min(f - clipped.clip_shift, clipped.domain(x))
                           for (f, x) in zip(rbf_field, beyond_face)]
          @test any(rbf_field .- clipped.shift .> 0)
          @test all(clipped_field .<= 0)
        end
      end

      @testset "Sharpening" begin
        # φ̃ follows the raw SDF more closely, so the shift is smaller; the volume is kept and
        # interpolation does not use it
        approximation = RBFs_smoothing(mesh, raw_sdf, sdf_grid, false)
        sharpened = RBFs_smoothing(mesh, raw_sdf, sdf_grid, false; sharpening = true)
        @test !approximation.sharpening && sharpened.sharpening
        @test sharpened.nodal_error < approximation.nodal_error
        @test abs(sharpened.shift) < abs(approximation.shift)
        @test abs(sharpened.clip_volume / sharpened.target_volume - 1) <= 1e-5
        @test RBFs_smoothing(mesh, raw_sdf, sdf_grid, true; sharpening = true).weights ==
              RBFs_smoothing(mesh, raw_sdf, sdf_grid, true).weights

        # One step with the kernel matrix K, the formula at the top of RBFs4Smoothing.jl
        K = rbf_kernel_matrix(Tuple(sdf_grid.N .+ 1), sharpened.threshold)
        S = sum(kernel_offsets(sharpened.threshold)[2])
        @test vec(sharpened.weights) ≈ vec(approximation.weights) .+
              (raw_sdf .- K * vec(approximation.weights)) ./ S

        # With and without sharpening the results go to different files (JLD2 and VTI)
        mktempdir() do dir
          cd(dir) do
            approximation_file = export_rbf_results(approximation, sdf_grid, "bar")
            sharpened_file = export_rbf_results(sharpened, sdf_grid, "bar")
            @test approximation_file != sharpened_file
            @test length(readdir()) == 4
            @test !load(approximation_file, "sharpening") && load(sharpened_file, "sharpening")
          end
        end

        # Off by default; outside the approximation it is ignored with a warning
        @test !Rho2sdfOptions().rbf_sharpening
        @test Rho2sdfOptions(rbf_smoothing = :approximation, rbf_sharpening = true).rbf_sharpening
        ignored = @test_logs (:warn, r"rbf_sharpening") Rho2sdfOptions(rbf_sharpening = true)
        @test !ignored.rbf_sharpening
      end

      @testset "Nearest boundary point and position in the domain" begin
        # Against the exact nearest point of the block. It is unique outside the block; inside it,
        # d_D above checks the distance.
        cap = 0.4
        d_D = DomainDistance(mesh, cap)
        low, high = [-1.0, -0.5, -0.5], [1.0, 0.5, 0.5]
        largest_error = 0.0
        position_right = true
        element_right = true
        for x in range(-1.3, 1.3, 23)
          for y in range(-0.8, 0.8, 15)
            for z in (-0.7, -0.5, -0.31, 0.0, 0.17, 0.5, 0.66)
              p = [x, y, z]
              gap = box_gap(p, low, high)
              position_right &= inside_domain(d_D, Point(p...)) == (gap == 0)
              if 0 < gap < cap
                distance, point, triangle = nearest_boundary_point(d_D, Point(p...))
                largest_error = max(largest_error, norm(point - clamp.(p, low, high)),
                                    abs(distance - gap))

                # The point lies on the element of its triangle
                corners = X[IEN[d_D.triangle_elements[triangle]]]
                element_low = reduce((a, b) -> min.(a, b), corners)
                element_high = reduce((a, b) -> max.(a, b), corners)
                element_right &= box_gap(point, element_low, element_high) <= 1e-12
              end
            end
          end
        end
        @test largest_error <= 1e-12
        @test position_right
        @test element_right

        # No triangle within the cap, inside or outside the bins
        for p in (Point(-1.5, 0.0, 0.0), Point(10.0, 0.0, 0.0))
          @test nearest_boundary_point(d_D, p) == (cap, p, 0)
          @test !inside_domain(d_D, p)
        end
      end

      # The raw SDF of the iso-contour below is built from the nodal densities, as rho2sdf builds
      # it, not from the analytic `bar_sdf`
      ρₙ = DenseInNodes(mesh, rho)
      ρₜ = 0.5
      dims = Tuple(sdf_grid.N .+ 1)
      whole_distances, _ = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      iso_distances, _ = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ; domain_faces = false)
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      whole_sdf = whole_distances .* signs
      iso_sdf = extended_raw_sdf(mesh, sdf_grid, points, ρₙ, ρₜ)

      # Node numbers (1-based) of the domain face x = -1 and of the block along y, z: the grid has
      # h = 0.125 and 3 cells of margin. The section of the bar is |y|, |z| < 0.25.
      face_i = 4
      block_jk = 4:12
      section_jk = 7:9

      @testset "Distances at any points, without the faces of the domain" begin
        # Fewer points than grid nodes give the same values
        some = 1:3:sdf_grid.ngp
        @test evalDistances(mesh, sdf_grid, points[:, some], ρₙ, ρₜ)[1] == whole_distances[some]

        # The faces of the domain are a part of the surface by default, and only there
        on_face(distances) = reshape(distances, dims)[face_i, section_jk, section_jk]
        @test all(on_face(whole_distances) .== 0)
        @test all(on_face(iso_distances) .> 0)
      end

      @testset "Raw SDF of the iso-contour, extended beyond the domain" begin
        d_D = DomainDistance(mesh, 8 * sdf_grid.cell_size)
        in_domain = [inside_domain(d_D, Point(points[:, v]...)) for v in 1:sdf_grid.ngp]

        # In the domain: the SDF of the iso-contour
        @test iso_sdf[in_domain] == (iso_distances .* signs)[in_domain]

        # Beyond the face x = -1: the values on the face, constant along x, and the sign says
        # whether the node is in front of the bar
        grid_sdf = reshape(iso_sdf, dims)
        @test all(grid_sdf[i, block_jk, block_jk] == grid_sdf[face_i, block_jk, block_jk]
                  for i in 1:face_i-1)
        @test all(grid_sdf[1, section_jk, section_jk] .> 0)

        # Nodes that share their nearest boundary point have the same value, and each of them the
        # value of evalDistances at the point (the edge of its band is sensitive to round-off, so
        # the placeholders are left out)
        values_at_point = Dict{NTuple{3,Int},Vector{Float64}}()
        nearest_points = copy(points)
        for v in findall(.!in_domain)
          _, point, triangle = nearest_boundary_point(d_D, Point(points[:, v]...))
          if triangle > 0
            nearest_points[:, v] = point
            key = ntuple(a -> round(Int, point[a] / (1e-9 * sdf_grid.cell_size)), 3)
            push!(get!(values_at_point, key, Float64[]), iso_sdf[v])
          end
        end
        @test !isempty(values_at_point)
        @test all(allequal, values(values_at_point))
        at_nearest_points, _ = evalDistances(mesh, sdf_grid, nearest_points, ρₙ, ρₜ;
                                             domain_faces = false)
        computed = (abs.(iso_sdf) .< 1e9) .& (at_nearest_points .< 1e9)
        @test maximum(abs.(abs.(iso_sdf[computed]) .- at_nearest_points[computed])) <= 1e-9
      end

      @testset "Only the clipped body of the iso-contour" begin
        # Interpolation, approximation and sharpening; the tolerance of the volume as in "Clip keeps
        # the volume".
        for (is_interpolation, sharpening) in ((true, false), (false, false), (false, true))
          rbf = RBFs_smoothing(mesh, iso_sdf, sdf_grid, is_interpolation; sharpening,
                               isocontour_only = true)
          @test rbf.isocontour_only
          @test abs(rbf.clip_volume / rbf.target_volume - 1) <= 1e-5
          @test rbf.shift === rbf.volume === rbf.volume_sampled === rbf.volume_outside === nothing

          # The face x = -1 lies inside the field over the section of the bar, so the clip draws it
          section = range(-0.125, 0.125, 9)
          on_face = [Point(-1.0, y, z) for y in section for z in section]
          grid_position(x) = Tuple((x .- rbf.origin) ./ rbf.h)
          @test all(rbf_value(rbf.weights, rbf.threshold, grid_position(x)) - rbf.clip_shift > 0
                    for x in on_face)
        end

        @test_throws ErrorException RBFs_smoothing(mesh, iso_sdf, sdf_grid, true;
                                                   domain_clip = false, isocontour_only = true)
      end

      @testset "No pits in a solid face of the domain" begin
        # Body: the half of the block x < 0, so the face x = -1 is solid over several elements. The
        # whole SDF is zero all over it and φ̃ rings around this zero, which leaves pits where
        # φ̃ - clip_shift < 0 (the bar above is too thin: there the shift outweighs the ringing).
        # The SDF of the iso-contour has no zero there, so φ̃ - clip_shift stays well above zero.
        half_rho = [c[1] < 0 ? 1.0 : 0.0 for c in centres]
        half_mesh = Mesh(X, IEN, half_rho, hex8_shape)
        half_ρₙ = DenseInNodes(half_mesh, half_rho)
        half_whole_sdf = evalDistances(half_mesh, sdf_grid, points, half_ρₙ, ρₜ)[1] .*
                         Sign_Detection(half_mesh, sdf_grid, points, half_ρₙ, ρₜ)
        half_iso_sdf = extended_raw_sdf(half_mesh, sdf_grid, points, half_ρₙ, ρₜ)

        section = range(-0.4, 0.4, 41)
        on_face = [Point(-1.0, y, z) for y in section for z in section]
        function shifted_field(rbf)
          grid_position(x) = Tuple((x .- rbf.origin) ./ rbf.h)
          return [rbf_value(rbf.weights, rbf.threshold, grid_position(x)) - rbf.clip_shift
                  for x in on_face]
        end

        @test any(<(0), shifted_field(RBFs_smoothing(half_mesh, half_whole_sdf, sdf_grid, true)))
        for (is_interpolation, sharpening) in ((true, false), (false, false), (false, true))
          rbf = RBFs_smoothing(half_mesh, half_iso_sdf, sdf_grid, is_interpolation; sharpening,
                               isocontour_only = true)
          @test all(>(0.1), shifted_field(rbf))
        end
      end

      @testset "Export of the clipped body and the names of the files" begin
        mesher_keys = ("weights", "origin", "h", "sigma", "threshold", "inside", "domain_clip",
                       "clip_shift", "domain_nodes", "domain_elements")
        unclipped_keys = ("shift", "volume", "volume_outside")
        mktempdir() do dir
          cd(dir) do
            # The iso-contour alone: no keys of the unclipped body, so a reader that does not know
            # the clip fails. The files keep the plain name.
            rbf = RBFs_smoothing(mesh, iso_sdf, sdf_grid, true; isocontour_only = true)
            file = export_rbf_results(rbf, sdf_grid, "bar")
            stored = load(file)
            @test file == "Z_bar_RBF_B-0.125_Interpolation.jld2"
            @test isfile("bar_B-0.125_RBF_Interpolation.vti")
            @test stored["isocontour_only"]
            @test all(haskey(stored, key) for key in mesher_keys)
            @test !any(haskey(stored, key) for key in unclipped_keys)

            # The whole SDF: all keys, and files that the run above does not overwrite
            rbf = RBFs_smoothing(mesh, whole_sdf, sdf_grid, true)
            file = export_rbf_results(rbf, sdf_grid, "bar")
            stored = load(file)
            @test file == "Z_bar_RBF_B-0.125_Interpolation_WholeSDF.jld2"
            @test isfile("bar_B-0.125_RBF_Interpolation_WholeSDF.vti")
            @test !stored["isocontour_only"]
            @test all(haskey(stored, key) for key in mesher_keys)
            @test all(haskey(stored, key) for key in unclipped_keys)

            sharpened = RBFs_smoothing(mesh, whole_sdf, sdf_grid, false; sharpening = true)
            @test export_rbf_results(sharpened, sdf_grid, "bar") ==
                  "Z_bar_RBF_B-0.125_Approximation_Sharpened_WholeSDF.jld2"
            @test isfile("bar_B-0.125_RBF_Approximation_Sharpened_WholeSDF.vti")

            # Without the clip: the whole SDF, the keys of the unclipped body only
            unclipped = RBFs_smoothing(mesh, whole_sdf, sdf_grid, false; domain_clip = false)
            stored = load(export_rbf_results(unclipped, sdf_grid, "bar"))
            @test haskey(stored, "shift") && !haskey(stored, "clip_shift")
          end
        end
      end

      @testset "Option rbf_isocontour_only" begin
        @test Rho2sdfOptions().rbf_isocontour_only
        @test !Rho2sdfOptions(rbf_isocontour_only = false).rbf_isocontour_only

        # Nothing would draw the faces of the domain without the clip
        no_clip = @test_logs (:warn, r"rbf_isocontour_only") Rho2sdfOptions(domain_clip = false)
        @test !no_clip.rbf_isocontour_only

        # :none exports the raw SDF unchanged, so the option does not matter there
        @test_logs Rho2sdfOptions(rbf_smoothing = :none, domain_clip = false)
      end

      @testset "rho2sdf" begin
        # The branches of rho2sdf(): the raw SDF alone, the iso-contour alone (the default) and the
        # whole SDF. Each run goes to its own directory, to see the files it writes.
        function run_rho2sdf(; settings...)
          return mktempdir() do dir
            cd(dir) do
              options = Rho2sdfOptions(; threshold_density = ρₜ, sdf_grid_setup = :automatic,
                                       export_raw_sdf = true, settings...)
              return (rho2sdf("bar", X, IEN, rho; options)..., readdir())
            end
          end
        end

        # The grid of the runs, for the comparison with the functions called directly
        grid = noninteractive_sdf_grid_setup(mesh)
        grid_points = generateGridPoints(grid)
        expected_whole = evalDistances(mesh, grid, grid_points, ρₙ, ρₜ)[1] .*
                         Sign_Detection(mesh, grid, grid_points, ρₙ, ρₜ)

        _, raw, rbf, files = run_rho2sdf(rbf_smoothing = :none)
        @test rbf === nothing
        @test raw == expected_whole
        @test !any(contains("RBF"), files)
        @test any(contains("_SDF_CellSize-"), files) && !any(endswith("_IsoContour.vti"), files)

        _, raw, rbf, files = run_rho2sdf()
        @test rbf.isocontour_only && rbf.shift === nothing
        @test raw == extended_raw_sdf(mesh, grid, grid_points, ρₙ, ρₜ)
        @test abs(rbf.clip_volume / rbf.target_volume - 1) <= 1e-5
        @test any(endswith("_Interpolation.jld2"), files)
        @test any(endswith("_IsoContour.vti"), files)

        _, raw, rbf, files = run_rho2sdf(rbf_isocontour_only = false)
        @test !rbf.isocontour_only && rbf.shift !== nothing
        @test raw == expected_whole
        @test rbf.weights == RBFs_smoothing(mesh, expected_whole, grid, true).weights
        @test any(endswith("_Interpolation_WholeSDF.jld2"), files)
        @test any(contains("_SDF_CellSize-"), files) && !any(endswith("_IsoContour.vti"), files)
      end
    end
  end

  if (RUN_PLANE)
    @testset "Plane" begin

      taskName = "plane"
      N = 15  # Number of cells along the longest side
      ρₜ = 0.5 # Threshold density (isosurface level)

      X = [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
      ]
      IEN = [[1, 2, 3, 4, 5, 6, 7, 8]]
      # ρₙ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1]
      # ρₙ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 0.5]
      # ρₙ = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0] # dve roviny 4 uzlové
      # ρₙ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1]
      # ρₙ = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1]
      ρₙ = [1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1]
      # ρₙ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] # 1, 7

      ## Generate FEM mesh structure:
      mesh = Mesh(X, IEN, rho, hex8_shape)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      ## Export to VTK:
      exportStructuredPointsToVTK(taskName * "_SDF.vtk", sdf_grid, sdf_dists, "distance")

      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists

    end
  end

  if (RUN_BLOCK)
    @testset "Block" begin
      ## Inputs:
      taskName = "block"

      N = 10  # Number of cells along the longest side
      ρₜ = 0.5 # Threshold density (isosurface level)

      (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [2, 1, 1])
      # ρₙ = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]

      mesh = Mesh(X, IEN, rho, hex8_shape)
      # ρₙ = DenseInNodes(mesh, rho) # LSQ

      # Modif ρₙ:
      ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]
      
      # ρₜ= find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      ## Export to VTK:
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      rbf = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp) # interpolation == true, aproximation == false
      export_rbf_results(rbf, sdf_grid, taskName)

      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists
      # @save "Z_$(taskName)_rho.jld2" ρₙ

    end
  end


  if (RUN_SPHERE)
    @testset "Sphere" begin
      ## Inputs:
      taskName = "sphere"

      N = 10  # Number of cells along the longest side
      # ρₜ = 0.5 # Threshold density (isosurface level)
      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      # data = matread("test/" * taskName * ".mat")
      (X, IEN, rho) = MeshInformations(data)

      ## Generate FEM mesh structure:
      mesh = Mesh(X, IEN, rho, hex8_shape)

      ## Map elemental densities to the nodes:
      ρₙ = DenseInNodes(mesh, rho) # LSQ
      #ρₙ = elementToNodalValues(mesh, rho) # average

      ρₜ= @time find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export to VTK:
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      rbf = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp) # interpolation == true, aproximation == false
      export_rbf_results(rbf, sdf_grid, taskName)
     
      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists
      # @save "Z_$(taskName)_rho.jld2" ρₙ

    end
  end

  if (RUN_BEAM)
    @testset "Chapadlo" begin
      ## Inputs:
      taskName = "cantilever_beam_vfrac_04"
      
      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshInformations(data)

      # Custom options
      options = Rho2sdfOptions(
          threshold_density=0.5,
          sdf_grid_setup=:automatic,
          export_nodal_densities=false,
          export_raw_sdf=false,
          rbf_smoothing=:interpolation
      )

      result = rho2sdf("beam", X, IEN, rho, options=options)

    end
  end

  if (RUN_CHAPADLO)
    @testset "Chapadlo" begin
      ## Inputs:
      taskName = "chapadlo"
      N = 60  # Number of cells along the longest side
      # ρₜ = 0.5 # Threshold density (isosurface level)

      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshInformations(data)
      #Z,idx_Z = findall(x->X[3,i] > 50 for i in [1:size(X,2)])

      ## Generate FEM mesh structure:
      mesh = Mesh(X, IEN, rho, hex8_shape)

      ## Map elemental densities to the nodes:
      ρₙ = DenseInNodes(mesh, rho) # LSQ
      #ρₙ = elementToNodalValues(mesh, rho) # average

      ρₜ= @time find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      X_min[2] = 0
      X_min[3] = 50
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export to VTK:
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      rbf = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp) # interpolation == true, aproximation == false
      export_rbf_results(rbf, sdf_grid, taskName)
    
      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists
      # @save "Z_$(taskName)_rho.jld2" ρₙ

    end
  end

  if (RUN_CHAPADLO_cele)
    @testset "Chapadlo" begin
      ## Inputs:
      taskName = "chapadlo"
      # ρₜ = 0.5 # Threshold density (isosurface level)

      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshInformations(data)

      ## Generate FEM mesh structure:
      mesh = Mesh(X, IEN, rho, hex8_shape)

      ## Grid:
      # sdf_grid = interactive_sdf_grid_setup(mesh)
      sdf_grid = noninteractive_sdf_grid_setup(mesh)
      points = generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## Map elemental densities to the nodes:
      ρₙ = DenseInNodes(mesh, rho) # LSQ
      
      ρₜ= find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      ## Export to VTK:
      B = round(sdf_grid.cell_size, digits=4)
      exportSdfToVTI(taskName * "_SDF_CellSize-" * string(B) * ".vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      rbf = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp) # interpolation == true, aproximation == false
      export_rbf_results(rbf, sdf_grid, taskName)

      # @save "Z_$(taskName)_cele_xp_B-$(B).jld2" xp
      # @save "Z_$(taskName)_cele_SDF_B-$(B).jld2" sdf_dists
      # @save "Z_$(taskName)_cele_Grid_B-$(B).jld2" sdf_grid
      # @save "Z_$(taskName)_cele_Points_B-$(B).jld2" points
      # @save "Z_$(taskName)_cele_Mesh_B-$(B).jld2" mesh
      #
    end
  end
end
