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
using Rho2sdf.SdfSmoothing: DomainDistance, Point, rbf_value
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
  RUN_DOMAIN_CLIP = true

  if (RUN_DOMAIN_CLIP)
    @testset "Domain clip" begin
      # Design domain: the block [-1, 1] × [-0.5, 0.5]² of 8 × 4 × 4 elements (built here: the
      # primitive "block" fills a Dict from several threads). Body: the bar
      # [-1, 0] × [-0.25, 0.25]², which touches the domain face x = -1.
      node(i, j, k) = 1 + i + 9 * (j + 5 * k)
      X = [[-1.0 + 0.25 * i, -0.5 + 0.25 * j, -0.5 + 0.25 * k] for k = 0:4 for j = 0:4 for i = 0:8]
      IEN = [[node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k),
              node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1),
              node(i, j + 1, k + 1)]
             for k = 0:3 for j = 0:3 for i = 0:7]
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

      @testset "Clip keeps the volume" begin
        # Raw SDF as rho2sdf builds it: distance to the iso-contour, i.e. the faces of the bar
        # except its end x = -1, which lies on the domain boundary; positive in the bar. It jumps
        # from + to - across x = -1, so the RBF field reaches beyond the domain face.
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

          # The clipped body keeps the volume
          @test abs(clipped.clip_volume / clipped.target_volume - 1) <= 1e-6
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
     
      # fig = visualize_stable_isosurface(fine_LSF)
      # display(fig)

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
      IEN = [subvector .- 1 for subvector in IEN]  # This will work

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
