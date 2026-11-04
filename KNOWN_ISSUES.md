# Known issues

Open issues found in the code review of rho2sdf.jl (review ids R…) and during the work on the domain clip (D…): none at the moment, apart from the limitations below.

Limitations of the domain clip (`domain_clip`, `src/SdfSmoothing/DomainDistance.jl`):
- d_D approximates a warped (bilinear) element face by 4 triangles around the face centre; it is
  exact for planar faces only.
- Implicit2TetMesh.jl evaluates d_D with its own copy of `DomainDistance` (`ClippedRBFSDF`, branch
  `dev`): a change of d_D here has to be made there too, or the meshed body loses the volume.
- With the raw SDF of the iso-contour alone (`rbf_isocontour_only`,
  `src/SdfSmoothing/ExtendedRawSdf.jl`), the unclipped field is not a body: φ̃ continues beyond the
  domain up to the edge of the grid. The body is given by `clip_shift` and the domain only, `shift`,
  `volume` and `volume_outside` are not in the JLD2 file.
- The extension beyond the domain (the value at the nearest point of its boundary) jumps behind the
  concave edges of the domain, where the nearest point changes. The clip removes these values, but
  the field keeps a part of the jump through the weights of the nodes around it. To verify: whether
  this shapes the surface near such edges, with a field fitted without the values beyond the domain
  and on a mesh finer than h. Example on the gripper (h = 2 mm): a bump of 0.3 mm over 2 mm at
  (60.3, -10.5, 30), next to a node with +1.69 and a neighbour with -1.22; 0.45 % of the surface has
  a higher bump over the same 2 mm, so it is not an outlier.
- The edge of the band of 2.5 h around the iso-contour is sensitive to round-off: nodes between
  2.5 h and 4.3 h from the iso-contour get the distance or the placeholder, depending on the
  rounding of the indices of the mini AABB grid (`evalDistances`).

Limitations of the nodal densities (`src/MeshGrid/NodalDensities.jl`):
- `NodalDensityLeastSquares` does not clamp the fit value, so a nodal density can leave [0, 1], although
  the element densities from SIMP lie in it. `KeepRange` exists but nothing calls it. simp2sdf.jl has
  the same file and clamps the fit value, this is the one deliberate difference between them.
- The effect of the clamp on the final geometry is not known. To verify: compare the geometry with
  and without the clamp, here also through `find_threshold_for_volume`, which takes the threshold from
  the nodal densities. Example on the gripper of simp2sdf.jl (`test/TestData/Gripper.vtu`, 20,542
  nodes): the clamp would change one node by more than 1e-12 (1.0047 → 1.0), the other 543 nodes above
  1 differ by round-off only. There are no values below 0 on this mesh.

Resolved:
- R2 (weights of `FilterForNodalDensity` described as distance-weighted, while they are
  1 − L/(1.2 L_max)): the code is unchanged, its docstring now gives the weights.
- R5 (only `ret == :FAILURE` of SLSQP was checked): `compute_coords_on_iso` warns for every
  return code other than `:SUCCESS`, `:STOPVAL_REACHED`, `:FTOL_REACHED` and `:XTOL_REACHED`,
  with the code in the message. Sphere and beam tests: 198× `:FAILURE` (warned before too),
  58× `:ROUNDOFF_LIMITED` (new). The effect of these foot points on the field is still not
  quantified.
- R6 (`Sign_Detection` found the candidate elements of every grid node by a scan over the AABBs
  of all elements): the candidate elements are collected once, each element testing only the grid
  nodes of its mini AABB grid (`calculateMiniAABB_grid`, one cell of margin), in the same ascending
  order. Gripper (76,050 grid nodes): 18.7 s → 0.4 s, the signs are unchanged.
- R8 (`cantilever_beam.jl` subtracted 1 from `IEN` after `MeshInformations` had added 1):
  `MeshInformations` now finds out whether `IEN` numbers the nodes from 0 or from 1, the example
  `.mat` files use both. The subtraction in `cantilever_beam.jl` and in the beam test is removed.
- D1 (element volumes summed with `atomic_add!` inside `@threads`, so the last bits of `V_domain`,
  `V_frac` and the target volume changed from run to run): `calculate_mesh_volume` and
  `calculate_isocontour_volume` store the element volumes and sum them after the loop.
- D2 (`TestGeometryBlock` and `TestGeometryCube` filled the `Dict` `node_map` inside `@threads`
  and crashed on several threads): the node numbers are computed by a function instead.

Resolved by the removal of the grid-sampled RBF path:
- R7 (Float32 kernel matrix and weights, default `cg` tolerance): the kernel matrix is now assembled
  in Float64 from the lattice offsets (`rbf_kernel_matrix`), CG runs with reltol 1e-12.
- The shift bisection of `RBFs_smoothing` stopped after 60 steps without a warning when the
  relative tolerance 1e-9 was not reached. It now warns when the volume deviates by more than 1e-6.
