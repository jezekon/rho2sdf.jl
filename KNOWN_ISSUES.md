# Known issues

Open issues found in the code review of rho2sdf.jl (review ids R…). None of them is fixed yet.

| ID | Where | Issue | Effect |
|---|---|---|---|
| R2 | `FilterForNodalDensity` (`src/MeshGrid/NodalDensities.jl`) | Nodes with 2–3 elements get the weights 1 − L/(1.2 L_max), where L is the distance from the node to the element centre. These are not inverse-distance weights. | Same result on uniform meshes, slightly different on irregular boundary edges. Align either the code or its description. |
| R5 | `compute_coords_on_iso` (`src/SignedDistances/ComputeCoordsOnIso.jl`) | Only `ret == :FAILURE` of SLSQP is checked. `:ROUNDOFF_LIMITED` and `:MAXEVAL_REACHED` pass without a warning. | A non-converged foot point gives a wrong distance with no warning. The effect on the field has not been quantified (it is local). |
| R6 | `Sign_Detection` (`src/SignedDistances/SignDetection.jl`) | For every grid node, the candidate elements are found by a linear scan over the AABBs of all elements. | Correct but O(grid nodes × elements): 457,560 × 16,490 checks on the gripper. A k-d tree or a bin grid would remove it. |
| R8 | `cantilever_beam.jl` | Subtracts 1 from `IEN` after `MeshInformations` has already added 1. | Correct only for `.mat` inputs with 0-based `IEN`. The example `.mat` files use different index bases. |

Resolved by the removal of the grid-sampled RBF path:
- R7 (Float32 kernel matrix and weights, default `cg` tolerance): the kernel matrix is now assembled
  in Float64 from the lattice offsets (`rbf_kernel_matrix`), CG runs with reltol 1e-12.
- The shift bisection of `RBFs_smoothing` stopped after 60 steps without a warning when the
  relative tolerance 1e-9 was not reached. It now warns when the volume deviates by more than 1e-6.
