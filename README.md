# Full Potential Solver Notes

This document summarizes how the Chapel `fullPotentialSolver` code is organized and what numerical method it is actually implementing today.

## Overview

The solver is a 2D cell-centered finite-volume/full-potential code on an unstructured mesh read from CGNS. The primary unknown is the velocity potential `phi` stored at cell centers. The code reconstructs velocity from `grad(phi)`, computes density from an isentropic full-potential relation, forms mass fluxes across faces, and drives the steady residual to zero with an approximate Newton method whose linear system is solved with PETSc GMRES.

At a high level, the steady algorithm is:

1. Read mesh and build connectivity, including ghost cells on boundaries.
2. Precompute geometry, least-squares gradient weights, wake bookkeeping, and Jacobian sparsity.
3. Initialize `phi`, reconstruct velocity, density, face fluxes, and residual.
4. Repeatedly:
   - assemble an analytical-but-approximate Jacobian,
   - solve `J * dphi = -R` with GMRES,
   - update `phi`,
   - recompute circulation from the trailing-edge jump,
   - reevaluate the residual and aerodynamic outputs.

## Source Map

- `src/main.chpl`: top-level driver.
- `src/input.chpl`: runtime/config parameters and freestream initialization.
- `src/mesh.chpl`: CGNS mesh reader, edge/element/node connectivity, ghost-cell construction.
- `src/leastSquaresGradient.chpl`: precomputed least-squares gradient operators.
- `src/spatialDiscretization.chpl`: geometry metrics, wake/Kutta setup, residual assembly, wall/farfield treatment, aerodynamic coefficients.
- `src/temporalDiscretization.chpl`: nonlinear iteration driver, Jacobian assembly, PETSc linear solve.
- `src/linearAlgebra.chpl`: helper norms and PETSc wrappers.
- `src/gmres.chpl`: native Chapel GMRES implementation with Jacobi/ILU(0), currently auxiliary.

## Build and Run

From the `fullPotentialSolver` directory:

```bash
make
./bin/main --help
./bin/main --GRID_FILENAME=pre/naca0012_129x129.cgns
```

In practice this code expects the external Chapel/C/HPC dependencies referenced in [`Makefile`](/home/user/test/fullPotentialSolver/Makefile) to be available, including HDF5, CGNS, PETSc, MPI, MKL, and the local `champs` support modules.

## Governing Model

The discretization matches a steady compressible full-potential formulation written as a flux balance:

`div( rho(grad(phi)) * grad(phi) ) = 0`

with density recovered from an isentropic Bernoulli relation. In code, the cell-centered density is computed as:

`rho = [1 + (gamma - 1)/2 * M_inf^2 * (1 - u^2 - v^2)]^(1/(gamma-1))`

where:

- `u = dphi/dx`
- `v = dphi/dy`
- the freestream is nondimensionalized so `|V_inf| = 1`

This makes the problem nonlinear because `rho` depends on `grad(phi)`.

## Mesh and Control Volumes

`MeshData` in [`src/mesh.chpl`](/home/user/test/fullPotentialSolver/src/mesh.chpl) reads:

- node coordinates,
- element-to-node connectivity,
- wall boundary edges,
- farfield boundary edges.

From that it builds:

- unique edges (`edge2node_`),
- element-to-edge connectivity (`elem2edge_`),
- edge-to-element connectivity (`edge2elem_`),
- element-neighbor structure (`esuel_`),
- node-neighbor support tables (`esup*`, `psup*`).

Boundary faces are converted into ghost-cell neighbors by replacing the missing side of a boundary edge with a new ghost element index. The ghost cell centroid is mirrored across the physical boundary face. This lets the flux and gradient loops treat interior and boundary faces through the same edge-based connectivity.

The control volume is cell-centered:

- unknown `phi` lives at element centroids,
- residual is one scalar per physical element,
- face fluxes are accumulated with sign according to edge orientation.

## Spatial Discretization

### 1. Cell geometry and face metrics

`spatialDiscretization.initializeMetrics()` computes:

- cell centroids,
- cell volumes/areas,
- face centroids,
- face lengths,
- oriented face normals,
- interpolation weights from element centroids to face centroids.

It also precomputes for each face:

- `t_IJ`: unit vector from cell `I` to cell `J`,
- `invL_IJ = 1 / |x_J - x_I|`,
- correction coefficients `n / (n . t_IJ)`.

Those coefficients are later reused both in the flux formula and in the Jacobian.

### 2. Gradient reconstruction

The production path uses `LeastSquaresGradientQR`, not the older normal-equation version.

For each cell, the gradient is reconstructed from neighboring cell values using a QR-based least-squares method in the style of Blazek:

`grad(phi)_i = sum_j w_ij * (phi_j - phi_i)`

Important implementation details:

- the weights depend only on mesh geometry and are precomputed once,
- final face-based weights are stored for both directions of each face,
- `sumWx_` and `sumWy_` are also stored to simplify Jacobian assembly.

### 3. Wake treatment and Kutta bookkeeping

`initializeKuttaCells()` identifies the trailing-edge node and constructs a wake line extending downstream. Cells are tagged as:

- `+1` above the wake,
- `-1` below the wake,
- `0` elsewhere.

For faces that cross the wake, the code enforces continuity of the physical potential by applying a circulation jump when evaluating neighbor values:

- upper looking to lower: `phi_lower + Gamma`
- lower looking to upper: `phi_upper - Gamma`

This correction is used consistently in:

- gradient reconstruction,
- face-potential differences,
- Jacobian sensitivity to circulation.

### 4. Boundary conditions

Ghost cells are used for both `phi` and velocity:

- Wall:
  - `phi_ghost = phi_interior` to impose zero normal potential gradient.
  - velocity is mirrored: `V_g = V_i - 2 (V_i . n) n`.
- Farfield:
  - `phi_ghost` is set from either freestream potential or the analytical cylinder solution.
  - velocity ghost state is chosen so that the face-average velocity matches the prescribed farfield value.

### 5. Face velocity and flux

After reconstructing cell velocities, `computeFaceProperties()` forms a corrected face velocity.

The basic idea is:

- interpolate cell-centered reconstructed velocity to the face,
- compare that interpolated directional derivative with the direct potential difference across the two cells,
- correct the interpolated velocity so its normal component is consistent with the direct potential jump.

In code, this is the deferred-correction style formula:

`V_face = V_avg - (V_avg . t_IJ - dphi/dl) * n / (n . t_IJ)`

where:

- `V_avg` is the distance-weighted average of cell velocities,
- `dphi/dl = (phi_J - phi_I) / |x_J - x_I|`,
- the wake jump is applied to `phi_J - phi_I` when the face crosses the wake.

Then the isentropic face density is computed from the corrected face velocity.

### 6. Artificial density / transonic stabilization

`artificialDensity()` adds a Jameson-like upwind blend in transonic regions by replacing the isentropic face density with

`rho_face = rho_isen - mu * (rho_isen - rho_upwind)`

where:

- `rho_upwind` comes from the upwind cell based on `V_face . n`,
- `mu = MU_C * max(0, M^2 - M_C^2)`.

So:

- subsonic faces keep the isentropic density,
- supersonic/transonic faces blend toward an upwind density for robustness.

This is the main shock-capturing/stabilizing mechanism in the current code.

### 7. Residual

For each face, the continuity flux is

`F_f = rho_face * (V_face . n_f) * |f|`

and the cell residual is the signed sum of its face fluxes.

`computeResiduals()` also evaluates a separate Kutta residual:

`R_Gamma = Gamma - Gamma_computed`

where `Gamma_computed` is obtained from the extrapolated upper/lower trailing-edge potentials.

Important: that Kutta residual is monitored, but there is no separate circulation unknown in the linear system right now.

## Nonlinear Solve Strategy

Despite the module name `temporalDiscretization`, the implemented path in [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl) is a steady nonlinear solve, not a physical time integrator.

The iteration in `solve()` is:

1. Evaluate the current residual `R(phi)`.
2. Assemble an approximate Jacobian `J`.
3. Solve `J * delta_phi = -R` with PETSc GMRES.
4. Update `phi <- phi + omega * delta_phi`.
5. Recompute circulation from the TE potential jump.
6. Re-run the residual assembly and check convergence.

The stopping logic uses:

- relative residual tolerance `CONV_TOL`,
- absolute residual tolerance `CONV_ATOL`,
- max iterations `IT_MAX`.

The code currently applies a fixed relaxation `OMEGA`; there is scaffolding for line search and other continuation features, but the steady path shown in `solve()` uses the direct update.

## Linear Solve

The production linear solve path is PETSc GMRES:

- matrix type: sequential AIJ,
- preconditioner selected by `GMRES_PRECON`,
- restart, tolerances, and max iterations taken from `input.chpl`.

The solver builds a sparse matrix with one row per physical cell. There is no extra row/column for circulation in the current assembled system.

There is also a native Chapel GMRES implementation in [`src/gmres.chpl`](/home/user/test/fullPotentialSolver/src/gmres.chpl), including Jacobi and ILU(0), but the current steady solve uses the PETSc wrapper in [`src/linearAlgebra.chpl`](/home/user/test/fullPotentialSolver/src/linearAlgebra.chpl).

## Jacobian Construction

### What is differentiated

`computeJacobian()` differentiates the discrete residual with respect to neighboring cell potentials using the precomputed least-squares weights and face metrics.

For each cell and each of its faces, the code differentiates the face contribution with respect to:

- the current cell potential,
- the neighboring cell potential,
- indirectly, the trailing-edge potentials through circulation sensitivity.

### Interior-face structure

For an interior face, the face contribution is split into:

1. a gradient/reconstruction part coming from the reconstructed face velocity,
2. a direct potential-difference term coming from `(phi_J - phi_I) / L_IJ`.

The gradient part uses:

- the QR least-squares weights for the current cell and neighbor,
- the sums of weights `sumWx_`, `sumWy_`,
- the effective corrected normal `m = n - k t` with `k = 1/(n.t)`.

The direct term contributes:

- a negative amount to the diagonal,
- a positive amount to the off-diagonal.

### Wall-face derivative

At a wall face, the Jacobian uses the mirrored ghost relation `phi_g = phi_i`, so only the interior-cell derivative appears. The effective normal is projected tangentially to reflect the slip-wall condition.

### Farfield-face derivative

At a farfield face, the ghost potential is fixed by the boundary condition, so the only retained derivative is the direct interior-cell term. The farfield state is treated as constant.

### Circulation sensitivity

The code computes `d(grad phi)/dGamma` for wake-affected cells in `computeGradientSensitivity()`. Then, in `computeJacobian()`, it forms `dR/dGamma` and converts it into sensitivity with respect to the upper and lower trailing-edge cell potentials using

- `Gamma = phi_upperTE - phi_lowerTE`
- `dGamma/dphi_upperTE = +1`
- `dGamma/dphi_lowerTE = -1`

So the wake/circulation effect enters the matrix through extra contributions added to those two physical columns.

### Important limitation: this is not a fully exact Newton Jacobian

The current Jacobian is analytical for the discretized flux form it differentiates, but it is still approximate in an important sense:

- `rho_face` is treated as frozen when assembling `J`,
- the derivative of density with respect to velocity/potential is not included,
- the derivative of the artificial-density switching/blending with respect to `phi` is not included,
- circulation is not solved as a separate unknown; it is updated explicitly after the linear solve.

So the method is best described as a sparse analytical approximate Newton or defect-correction Jacobian, not a fully consistent exact Jacobian of the whole nonlinear residual.

### AD exact reduced phi-only Jacobian option

There is now also a reduced exact path selected with `JACOBIAN_TYPE = "ad_reduced_exact"`.

A detailed implementation note for this path lives in [`AD_IMPLEMENTATION.md`](/home/user/test/fullPotentialSolver/AD_IMPLEMENTATION.md#L1). That document explains:

- how Enzyme is inserted into the Chapel build,
- which scalar kernels are differentiated,
- how the reduced Jacobian `dR/dphi - dR/dGamma * (dK/dphi)/(dK/dGamma)` is assembled,
- how circulation is re-synchronized after each phi update,
- and how the row-check / finite-difference verification hooks work.

There is also a derivation note for a future hand-coded exact Jacobian in [`EXACT_ANALYTICAL_JACOBIAN.md`](/home/user/test/fullPotentialSolver/EXACT_ANALYTICAL_JACOBIAN.md#L1). That note writes out the full chain rule needed to match the AD Jacobian analytically.

The codebase also now includes row-local exact analytical derivative routines used only for verification. They are compared against the reduced AD Jacobian and finite differences through the existing `CHECK_AD_ROW` diagnostics, which makes it possible to validate the derivation before attempting a full analytical matrix assembly.

This mode keeps `Gamma` constrained by the Kutta equation instead of storing it as an explicit linear-solve unknown. Writing the coupled exact system as

- `R(phi, Gamma) = 0`,
- `K(phi, Gamma) = 0`,

and assuming `dK/dGamma != 0`, the exact Jacobian of the reduced phi-only nonlinear map is the Schur-complement form

- `J_reduced = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`.

In [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl), Enzyme is used to differentiate both:

- each cell residual row `R_i(phi, Gamma)`,
- the scalar Kutta residual `K(phi, Gamma)`.

The code then assembles the reduced row

- `dR_i/dphi - (dR_i/dGamma)/(dK/dGamma) * dK/dphi`

directly into the phi-only PETSc matrix. After each Newton update of `phi`, `Gamma` is re-synchronized by solving the scalar Kutta equation `K(phi, Gamma) = 0` with a short Newton iteration before recomputing the residual.

So `ad_reduced_exact` is mathematically equivalent to eliminating `Gamma` from the coupled exact system, while keeping the linear system size at `N_phi` instead of `N_phi + 1`.

### Hand-coded exact reduced phi-only Jacobian option

There is now a second reduced exact path selected with `JACOBIAN_TYPE = "analytical_reduced_exact"`.

This mode assembles the same reduced phi-only Jacobian as `ad_reduced_exact`,

- `J_reduced = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`,

but computes all row derivatives from hand-coded analytical chain-rule routines instead of Enzyme. The implementation reuses the same row stencil and reduced Schur-complement structure, so it is directly comparable to the AD path.

The current implementation assembles those derivatives face by face on the local row stencil, in the same overall spirit as the original hand-coded `analytical` Jacobian, but with the full exact chain rule carried through density, transonic switching, and Kutta reduction.

The expensive row-local coefficient work in this path is now parallelized with Chapel `forall` loops. PETSc insertion remains sequential, but the analytical row assembly itself is computed concurrently across rows before being flushed into the matrix.

In the current codebase, this mode has been validated row-by-row against AD and finite differences through the `CHECK_AD_ROW` diagnostics. On the tested rows, the hand-coded reduced-exact row matched the AD row to machine precision.

## What `run()` Actually Does

The full residual evaluation sequence in `spatialDiscretization.run()` is:

1. update ghost-cell `phi`,
2. reconstruct velocity from `phi`,
3. compute density and the transonic switch `mu`,
4. update ghost-cell velocities,
5. compute face properties,
6. apply artificial-density stabilization,
7. compute continuity fluxes,
8. sum the cell residuals and Kutta residual.

This is the best single routine to read if you want the discrete operator in execution order.

## Outputs

The solver writes:

- cell-centered solution fields,
- wall-face fields,
- wake-face fields,
- convergence history,
- aerodynamic coefficients `Cl`, `Cd`, `Cm`,
- circulation history.

Output is written to CGNS through `writeCGNS.chpl`.

## Configuration Notes

The main runtime controls are defined in [`src/input.chpl`](/home/user/test/fullPotentialSolver/src/input.chpl) and sample values live in [`input.txt`](/home/user/test/fullPotentialSolver/input.txt).

Important active controls include:

- mesh file and element type,
- Mach number and angle of attack,
- shock-stabilization parameters `MU_C`, `MACH_C`,
- Jacobian mode through `JACOBIAN_TYPE` (`analytical`, `ad_reduced_exact`, or `analytical_reduced_exact`),
- GMRES tolerances and preconditioner,
- relaxation factor `OMEGA`,
- nonlinear backtracking controls `LINE_SEARCH`, `MAX_LINE_SEARCH`, and `SUFFICIENT_DECREASE`,
- convergence tolerances.

When `LINE_SEARCH=true`, `SUFFICIENT_DECREASE` controls the acceptance threshold for a trial step. Values above `1.0` permit residual growth, while values below `1.0` require an actual decrease.

There is also configuration scaffolding for:

- line search,
- adaptive upwinding,
- adaptive beta,
- adaptive GMRES forcing,
- dual-time stepping / unsteady flow,
- native GMRES selection,
- numerical Jacobian selection.

However, based on the current steady driver, many of those options are not fully wired into the active solve path yet. The current production path is the steady PETSc-GMRES solve with the analytical approximate Jacobian described above.

## Practical Interpretation

If you want a short mental model of the current code, it is:

- a cell-centered unstructured full-potential solver,
- using QR least-squares gradients,
- ghost-cell wall and farfield boundaries,
- a wake potential jump to model circulation,
- Jameson-style density upwinding for transonic robustness,
- either an analytical sparse approximate Jacobian or the exact phi-only Enzyme-reduced Jacobian,
- and a PETSc GMRES linear solve inside a steady nonlinear iteration.
