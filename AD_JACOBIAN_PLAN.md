# AD Full Jacobian Plan

This note records what "full Jacobian with automatic differentiation" should mean for the current `fullPotentialSolver` codebase, and what has to change to make it real.

## Why the current Jacobian is not full

The existing Jacobian in [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L198) is a hand-derived approximate Jacobian. It only inserts entries for:

- the diagonal `dR_i / dphi_i`,
- immediate face neighbors `dR_i / dphi_j`,
- two trailing-edge columns used to represent circulation sensitivity.

That is already useful, but it is not the full derivative of the implemented residual operator in [`src/spatialDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L952).

The true residual row `R_i(phi)` depends on more than the one-ring face neighbors because:

1. `grad(phi)` is reconstructed by QR least squares.
   - `grad(phi)_i` depends on all one-ring neighbors of `i`.
   - `grad(phi)_j` on a face `(i,j)` depends on all one-ring neighbors of `j`.
   - therefore the flux on face `(i,j)` depends on a two-ring stencil around `i`.

2. `rho_face` depends on `uFace` and `vFace`.
   - so the chain rule introduces extra terms even if the flux formula itself looks local.

3. `artificialDensity()` uses:
   - upwind switching based on `V_face . n`,
   - `mu = MU_C * max(0, M^2 - M_C^2)`,
   - upwind cell density.
   This adds further state dependence that is currently frozen in the analytical Jacobian.

4. wake/circulation corrections enter both:
   - the reconstructed gradient,
   - the direct `phi` jump across wake-crossing faces.

So a true assembled Jacobian will be larger-stencil and denser than the current matrix.

## What "full Jacobian" should differentiate

The target operator is the actual discrete residual evaluation:

1. update ghost-cell `phi`
2. reconstruct cell gradients / velocities
3. compute cell density and `mu`
4. update ghost-cell velocities
5. compute face velocity
6. compute face density
7. apply artificial density blending
8. compute face fluxes
9. sum cell residuals

If we say "full Jacobian", we should mean:

`J_ij = d R_i / d phi_j`

for the residual that comes out of that entire chain, not a partially frozen linearization.

## Practical AD strategies

There are three realistic strategies.

### 1. Full assembled Jacobian by AD over local row kernels

This is the best match to "full Jacobian".

Idea:

- Extract a pure scalar residual kernel `R_i(local_phi_patch)` for one cell.
- The local patch includes:
  - cell `i`,
  - face neighbors of `i`,
  - neighbors of those neighbors,
  - TE cells if wake coupling reaches the row.
- Use AD repeatedly with one seed per local DOF to assemble the full row.

Why this fits the solver:

- `R_i` is scalar, so forward-mode AD is straightforward.
- The local stencil is still sparse, so we avoid differentiating the whole global state at once.
- The resulting row automatically includes:
  - LS-gradient couplings,
  - density chain rule,
  - artificial-density dependence,
  - wake-jump dependence.

Main cost:

- one AD call per nonzero in the row stencil,
- row stencils are wider than current preallocation,
- some refactoring is required to isolate a side-effect-free row residual function.

### 2. Matrix-free AD Jacobian-vector product

This is cheaper and often the best engineering choice, but it is not an explicitly assembled full Jacobian.

Idea:

- expose a vector residual `R(phi)`,
- use AD to compute `J * v`,
- switch GMRES to matrix-free mode.

Pros:

- no sparse assembly,
- no row-by-row seeding,
- naturally includes the full chain rule.

Cons:

- not an assembled Jacobian,
- harder to inspect/debug sparsity,
- preconditioning becomes the main challenge.

### 3. Global assembled Jacobian by AD over the entire residual vector

This is the most literal interpretation, but it is the least attractive here.

It would require:

- differentiating the full residual vector routine,
- seeding every global DOF,
- handling a very large amount of intermediate activity.

For this codebase, this is likely too intrusive and too expensive compared with row-local assembly.

## Recommended path

For this solver, the recommended path is:

`assembled row-local AD Jacobian`

because it preserves the existing PETSc matrix + GMRES workflow while finally differentiating the real discrete operator.

## Required refactor before AD

Before Enzyme can be used cleanly, the residual logic needs a more AD-friendly shape.

### A. Separate geometry from state

The following are already geometry-only and should remain frozen inputs:

- centroids,
- volumes,
- face normals and areas,
- `t_IJ`, `invL_IJ`, correction coefficients,
- LS-gradient weights,
- wake connectivity metadata.

That part is good news: the code already precomputes most of it.

### B. Create a pure row residual kernel

Today the residual is evaluated through global mutable arrays inside `spatialDiscretization`.

For AD we want a function conceptually like:

```chapel
proc residualRowFromPatch(row: int, patchPhi: [] real(64)): real(64)
```

that:

- gathers local `phi` values from a patch,
- reconstructs only the quantities needed for row `i`,
- returns only `R_i`.

The more this function avoids mutating shared arrays, the easier it will be to differentiate robustly.

### C. Build an explicit row stencil map

The current matrix preallocation assumes a narrow row pattern. A full AD Jacobian needs a real stencil builder.

For each cell `i`, build:

- `rowStencil_[i] = sorted list of global phi indices that can affect R_i`

At minimum this should include:

- `i`,
- all face neighbors of `i`,
- all neighbors of those face neighbors,
- TE cells used by wake sensitivity if needed.

This stencil map is needed for:

- matrix preallocation,
- local patch extraction,
- AD seeding order,
- debugging the resulting sparsity.

## Important non-smooth pieces

AD will differentiate the executed branch, not a smoothed mathematical idealization.

That matters for:

- upwind selection in `artificialDensity()` from `V_face . n`,
- `max(0, M^2 - M_C^2)`,
- wake crossing conditionals,
- any future shock sensors.

This is still valid for Newton-like methods away from switching surfaces, but it means:

- the Jacobian may change abruptly when the upwind direction changes,
- convergence may still benefit from relaxation or continuation.

## Enzyme-specific implementation shape

At implementation time, the simplest likely pattern is:

1. add a small AD wrapper module with Enzyme entry points,
2. create a scalar row-residual kernel,
3. for each row and each stencil DOF:
   - seed one tangent,
   - call forward-mode AD,
   - store the returned derivative into PETSc.

Conceptually:

```chapel
J[row, col] = d residualRow(row, patchPhi) / d patchPhi[localCol]
```

mapped back to the global column index.

## Expected solver impact

Once a true AD Jacobian is used, we should expect:

- wider row sparsity than the current matrix,
- more expensive assembly,
- fewer Newton iterations if the Jacobian is materially more consistent,
- better agreement between linear model and residual update,
- a clearer path to comparing "approximate analytical" versus "full AD" Jacobians.

## Suggested implementation order

1. Build row stencil extraction utilities.
2. Write a scalar row residual kernel that reproduces the current residual row.
3. Add an AD proof of concept for one row and compare against finite differences.
4. Assemble the reduced exact PETSc matrix from AD row derivatives plus the scalar Kutta Schur complement.
5. Add a runtime switch:
   - `JACOBIAN_TYPE=analytical`
   - `JACOBIAN_TYPE=ad_reduced_exact`
6. Compare:
   - row sparsity,
   - assembly cost,
   - GMRES iterations,
   - nonlinear convergence rate.

## Bottom line

Using AD here is worthwhile mainly because the current hand Jacobian is structurally incomplete relative to the implemented residual. The most practical way to get a true phi-only exact Jacobian is not to differentiate the entire global solver in one shot, but to:

- extract a scalar row residual,
- differentiate it over an explicit local stencil,
- combine it with the scalar Kutta derivatives through a Schur complement,
- assemble those derivatives into PETSc.
