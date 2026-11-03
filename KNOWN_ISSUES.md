# Known issues

Open issues found in the code review of rho2sdf.jl (review ids R…) and during the work on the domain clip (D…).

| ID | Where | Issue | Effect |
|---|---|---|---|
| R2 | `FilterForNodalDensity` (`src/MeshGrid/NodalDensities.jl`) | Nodes with 2–3 elements get the weights 1 − L/(1.2 L_max), where L is the distance from the node to the element centre. These are not inverse-distance weights. | Same result on uniform meshes, slightly different on irregular boundary edges. Align either the code or its description. |
| R5 | `compute_coords_on_iso` (`src/SignedDistances/ComputeCoordsOnIso.jl`) | Only `ret == :FAILURE` of SLSQP is checked. `:ROUNDOFF_LIMITED` and `:MAXEVAL_REACHED` pass without a warning. | A non-converged foot point gives a wrong distance with no warning. The effect on the field has not been quantified (it is local). |
| R6 | `Sign_Detection` (`src/SignedDistances/SignDetection.jl`) | For every grid node, the candidate elements are found by a linear scan over the AABBs of all elements. | Correct but O(grid nodes × elements): 457,560 × 16,490 checks on the gripper. A k-d tree or a bin grid would remove it. |
| R8 | `cantilever_beam.jl` | Subtracts 1 from `IEN` after `MeshInformations` has already added 1. | Correct only for `.mat` inputs with 0-based `IEN`. The example `.mat` files use different index bases. |

Limitations of the domain clip (`domain_clip`, `src/SdfSmoothing/DomainDistance.jl`):
- d_D approximates a warped (bilinear) element face by 4 triangles around the face centre; it is
  exact for planar faces only.
- Implicit2TetMesh.jl evaluates d_D with its own copy of `DomainDistance` (`ClippedRBFSDF`, branch
  `dev`): a change of d_D here has to be made there too, or the meshed body loses the volume.

Resolved:
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
