# Known issues

Open issues found in the code review of rho2sdf.jl (review ids R…) and during the work on the domain clip (D…): none at the moment.

Limitations of the domain clip (`domain_clip`, `src/SdfSmoothing/DomainDistance.jl`):
- d_D approximates a warped (bilinear) element face by 4 triangles around the face centre; it is
  exact for planar faces only.
- Implicit2TetMesh.jl evaluates d_D with its own copy of `DomainDistance` (`ClippedRBFSDF`, branch
  `dev`): a change of d_D here has to be made there too, or the meshed body loses the volume.

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
