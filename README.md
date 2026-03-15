# Full Potential Solver

This directory contains a 2D compressible full-potential solver written in Chapel for cell-centered unstructured meshes read from CGNS.

The solver is built around a nonlinear residual in the velocity potential `phi`, a wake/circulation correction `Gamma`, and a Newton-Krylov solve driven by PETSc GMRES.

## Current Solver Structure

The active production path is:

1. Read the mesh and build ghost cells on wall and farfield boundaries.
2. Precompute geometry, face metrics, least-squares weights, wake bookkeeping, and Jacobian sparsity data.
3. Reconstruct velocity from `grad(phi)`.
4. Recover density from the isentropic full-potential relation.
5. Build corrected face velocities and continuity fluxes.
6. Add the transonic artificial-density stabilization.
7. Assemble the residual and Kutta relation.
8. Solve the steady nonlinear problem with Newton updates and PETSc GMRES.

## Source Map

- `src/main.chpl`: top-level driver.
- `src/input.chpl`: runtime parameters and freestream initialization.
- `src/mesh.chpl`: CGNS mesh reader and connectivity/ghost-cell construction.
- `src/leastSquaresGradient.chpl`: QR least-squares gradient operators.
- `src/spatialDiscretization.chpl`: geometry, wake tagging, boundary treatment, residual assembly, wall/wake output, aerodynamic coefficients.
- `src/temporalDiscretization.chpl`: nonlinear solve driver, PETSc setup, line search, omega scheduling, shared Jacobian helpers.
- `src/temporalDiscretizationAnalyticalApprox.chpl`: original analytical approximate Jacobian.
- `src/temporalDiscretizationAD.chpl`: Enzyme-based reduced exact Jacobian.
- `src/temporalDiscretizationAnalyticalExact.chpl`: hand-coded reduced exact Jacobian.
- `src/linearAlgebra.chpl`: norms and PETSc helper routines.
- `src/gmres.chpl`: native Chapel GMRES implementation.

## Spatial Discretization

### Unknown and control volumes

- The unknown `phi` is stored at cell centers.
- One residual equation is assembled for each physical element.
- Boundary faces are paired with ghost cells so the face loops can stay edge-based.

### Gradient reconstruction

The production gradient path uses QR least-squares reconstruction. For each cell,

`grad(phi)_i = sum_j W_ij * (phi_j^corr - phi_i)`

where `phi_j^corr` includes the wake jump when the edge crosses the wake.

The geometric QR weights are precomputed once and then reused in the residual and all Jacobian modes.

### Wake and circulation

The wake is introduced by tagging cells above and below the trailing-edge cut. Across a wake-crossing face, the solver corrects the neighbor potential by `+Gamma` or `-Gamma` so the physical potential remains continuous.

This affects:

- the least-squares gradient reconstruction,
- the direct face potential jump,
- the Kutta relation,
- the Jacobian sensitivity with respect to circulation.

### Face velocity

The corrected face velocity follows the deferred-correction form

`V_face = V_avg - (V_avg . t_IJ - dphi/dl) * corrCoeff`

with:

- `V_avg`: weighted average of the two cell velocities,
- `dphi/dl`: direct potential jump divided by cell-center distance,
- `corrCoeff = n / (n . t_IJ)`.

This is the core face formula used in both the primal residual and the exact Jacobians.

### Density and transonic stabilization

The isentropic density is computed from the full-potential Bernoulli relation:

`rho = [1 + (gamma-1)/2 * M_inf^2 * (1 - u^2 - v^2)]^(1/(gamma-1))`

Then `artificialDensity()` applies a Jameson-like transonic blend:

`rho_face = rho_isen - mu * (rho_isen - rho_upwind)`

with

`mu = MU_C * max(0, M^2 - M_C^2)`

using a hard upwind selection based on the sign of `V_face . n`.

### Residual

Each face contributes

`F_f = rho_face * (V_face . n_f) * A_f`

and each cell residual is the signed sum of its incident face fluxes.

The solver also evaluates a scalar Kutta residual that constrains circulation through the extrapolated upper/lower trailing-edge potentials.

## Nonlinear Solve

The steady solve lives in `temporalDiscretization.solve()`. The active algorithm is:

1. Evaluate the residual.
2. Assemble the selected Jacobian.
3. Solve `J * delta_phi = -R` with PETSc GMRES.
4. Apply a relaxed Newton update.
5. Re-synchronize `Gamma` from the Kutta condition.
6. Recompute the primal residual.
7. Accept or backtrack the step with the line search.

The solve uses:

- relative tolerance `CONV_TOL`,
- absolute tolerance `CONV_ATOL`,
- nonlinear backtracking via `LINE_SEARCH`, `MAX_LINE_SEARCH`, `SUFFICIENT_DECREASE`,
- optional omega scheduling for difficult early iterations.

## Jacobian Modes

### `JACOBIAN_TYPE="analytical"`

This is the original sparse analytical approximate Jacobian. It differentiates the gradient reconstruction and direct face-potential terms but freezes the nonlinear density response and the exact reduced Kutta coupling. It is cheap to assemble and often robust, but it is not the exact Jacobian of the implemented residual.

### `JACOBIAN_TYPE="ad_reduced_exact"`

This is the Enzyme reference Jacobian. It differentiates:

- each scalar residual row `R_i(phi, Gamma)`,
- the scalar Kutta residual `K(phi, Gamma)`,

and then assembles the reduced phi-only Jacobian

`J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`.

It is the most direct exact differentiation of the implemented discrete operator, but it is slower because Enzyme is invoked row by row.

Detailed note: [`AD_IMPLEMENTATION.md`](/home/user/test/fullPotentialSolver/AD_IMPLEMENTATION.md#L1)

### `JACOBIAN_TYPE="analytical_reduced_exact"`

This is the production exact Jacobian. It assembles the same reduced matrix as `ad_reduced_exact`, but with hand-coded chain-rule derivatives.

The implementation is optimized with:

- sparse row stencils,
- precomputed cell velocity-sensitivity templates,
- active-index loops inside each face contribution,
- parallel row assembly with Chapel `forall`.

This mode has been checked row-by-row against AD and finite differences and is much faster than the AD assembly.

Detailed note: [`EXACT_ANALYTICAL_JACOBIAN.md`](/home/user/test/fullPotentialSolver/EXACT_ANALYTICAL_JACOBIAN.md#L1)

## Exact Reduced Jacobian in One Equation

The nonlinear system is

- `R(phi, Gamma) = 0`
- `K(phi, Gamma) = 0`

where `K` is the Kutta relation.

Eliminating `Gamma` gives the phi-only exact Jacobian:

`J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`

Both exact modes assemble this same reduced matrix. The difference is only how the derivatives are obtained:

- AD with Enzyme,
- or hand-coded analytical chain rule.

## Practical Notes on Robustness

For difficult transonic cases on fine meshes, the early Newton steps can be too aggressive even with the exact Jacobian. The solver therefore supports:

- backtracking line search,
- fixed `OMEGA`,
- residual-triggered omega switching,
- CFD-style iteration-based `OMEGA_RAMP`.

The current fine-grid transonic workflow generally benefits from a small initial `omega` and a gradual ramp toward `1.0`.

## Build and Run

From `fullPotentialSolver/`:

```bash
make
./bin/main -f input.txt
```

Useful runtime overrides:

```bash
./bin/main -f input.txt --JACOBIAN_TYPE=analytical
./bin/main -f input.txt --JACOBIAN_TYPE=ad_reduced_exact
./bin/main -f input.txt --JACOBIAN_TYPE=analytical_reduced_exact
./bin/main -f input.txt --PROFILE_ITERATION_TIMINGS=true
```

## Documentation Map

- [`README.md`](/home/user/test/fullPotentialSolver/README.md#L1): solver overview and current implementation status.
- [`AD_IMPLEMENTATION.md`](/home/user/test/fullPotentialSolver/AD_IMPLEMENTATION.md#L1): how Enzyme is wired in and how the AD reduced Jacobian is assembled.
- [`EXACT_ANALYTICAL_JACOBIAN.md`](/home/user/test/fullPotentialSolver/EXACT_ANALYTICAL_JACOBIAN.md#L1): detailed derivation and implementation of the analytical reduced exact Jacobian.
- [`AD_JACOBIAN_PLAN.md`](/home/user/test/fullPotentialSolver/AD_JACOBIAN_PLAN.md#L1): archived design notes from the earlier AD development phase.
