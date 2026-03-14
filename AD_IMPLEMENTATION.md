# Automatic Differentiation Implementation

This note explains how automatic differentiation is implemented in `fullPotentialSolver`, what is differentiated, how the reduced exact Jacobian is derived, and how that Jacobian is used inside the nonlinear solve.

The active AD mode is:

- `JACOBIAN_TYPE = "ad_reduced_exact"`

in [`src/input.chpl`](/home/user/test/fullPotentialSolver/src/input.chpl#L71).

This is a phi-only exact Jacobian obtained by differentiating the discrete residual with Enzyme and eliminating circulation through the Kutta relation.

## 1. Build-time wiring

Enzyme is inserted into the Chapel build in [`Makefile`](/home/user/test/fullPotentialSolver/Makefile#L64).

The build sequence is:

1. Chapel emits LLVM bitcode with `--driver-compilation-phase`.
2. LLVM `opt` loads the Enzyme plugin and runs `-passes=enzyme`.
3. The transformed bitcode is assembled back to `.bc`.
4. Chapel finishes linking with `--driver-makebinary-phase`.

The relevant lines are:

- plugin discovery at [`Makefile:60`](/home/user/test/fullPotentialSolver/Makefile#L60)
- Chapel split compile at [`Makefile:66`](/home/user/test/fullPotentialSolver/Makefile#L66)
- Enzyme pass at [`Makefile:67`](/home/user/test/fullPotentialSolver/Makefile#L67)
- final Chapel link at [`Makefile:70`](/home/user/test/fullPotentialSolver/Makefile#L70)

So Enzyme operates on Chapel-generated LLVM IR rather than on source-level Chapel AST.

## 2. What is being differentiated

The implementation does not differentiate the whole solver driver or the full `spatialDiscretization.run()` routine directly. Instead, it differentiates two scalar functions in [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl):

- [`residualRowForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L207): one cell residual row `R_i(phi, Gamma)`
- [`kuttaResidualForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L289): the scalar Kutta residual `K(phi, Gamma)`

This is the central design choice.

The primal residual is still evaluated by [`spatialDiscretization.run()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L971), but the Jacobian comes from differentiating scalar kernels that reproduce the same discrete operator.

## 3. Why the code uses scalar AD kernels

Differentiating the whole primal residual pipeline directly would be difficult because it:

- mutates many arrays,
- mixes geometry/state updates,
- computes an entire residual vector when only one row derivative is needed,
- would be expensive to differentiate globally.

The scalar-kernel approach avoids that. Each AD function:

- takes immutable geometry/discretization metadata through `sd`,
- takes the active state as `phi` plus `Gamma`,
- returns a single scalar residual.

This makes Enzyme differentiation much easier and keeps the Jacobian assembly naturally row-oriented.

## 4. External Enzyme interface

The Enzyme ABI is declared at [`src/temporalDiscretization.chpl:28`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L28):

```chapel
extern {
    int enzyme_dup;
    int enzyme_const;
    void __enzyme_autodiff(void*, ...);
}
```

The important conventions are:

- `enzyme_const`: argument is treated as constant
- `enzyme_dup`: argument is active and accompanied by a derivative storage argument

So an Enzyme call of the form

```chapel
__enzyme_autodiff(f,
                  enzyme_const, a,
                  enzyme_dup, x, dx,
                  enzyme_dup, y, dy);
```

means:

- differentiate `f(a, x, y)`
- `a` is constant
- `x` and `y` are active inputs
- write `df/dx` into `dx`
- write `df/dy` into `dy`

## 5. Discrete primal operator being matched

The primal residual evaluation order in [`spatialDiscretization.run()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L971) is:

1. update ghost `phi`
2. reconstruct velocity from `phi`
3. compute density and transonic switch `mu`
4. update ghost velocities
5. compute corrected face properties
6. apply artificial-density stabilization
7. compute face fluxes
8. sum cell residuals and Kutta residual

The main corresponding primal routines are:

- [`computeVelocityFromPhiLeastSquaresQR()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L740)
- [`computeDensityFromVelocity()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L744)
- [`computeFaceProperties()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L811)
- [`artificialDensity()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L874)
- [`computeResiduals()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L943)

The AD kernels in `temporalDiscretization.chpl` are written to reproduce this same operator, but in pure scalar form.

## 6. Helper routines used by the AD kernels

Several small routines are there purely to make the AD kernels self-contained and pointer-based:

- [`phiFromPtr()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L34)
- [`gammaFromPtr()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L38)
- [`ghostPhiForFace()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L50)
- [`computeVelocityForCellFromPtr()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L70)
- [`computeDensityAndMuFromVelocity()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L113)
- [`computeGhostVelocityForFace()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L129)
- [`computeRhoMuForElemOnFace()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L164)

These routines do two things:

- they avoid dependence on mutable global solver arrays during differentiation,
- they reconstruct exactly the local quantities needed for one residual evaluation.

## 7. The row residual being differentiated

For a cell `i`, the differentiated row function is:

```text
R_i(phi, Gamma) = sum over faces f in cell i of sign(i,f) * rho_face * (V_face · n_f) * A_f
```

implemented in [`residualRowForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L207).

### 7.1 Gradient reconstruction

Inside [`computeVelocityForCellFromPtr()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L70), the velocity is reconstructed from the QR least-squares weights:

```text
grad(phi)_i = sum over faces/neighbor relations of w_if * (phi_j - phi_i)
```

If the neighbor crosses the wake, the code adds or subtracts `Gamma` from the neighbor potential before forming the difference. This is how circulation enters both the direct wake jump and the reconstructed gradients.

### 7.2 Face state reconstruction

For each face in [`residualRowForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L207):

1. Reconstruct both side velocities `V_1`, `V_2`
2. Form the weighted average velocity
3. Reconstruct the direct potential jump with wake correction
4. Apply deferred correction

The corrected face velocity is:

```text
V_face = V_avg - delta * corrCoeff
delta  = V_avg · t - dphi/dl
```

This is the same formula used in the primal face-property computation.

### 7.3 Density and transonic stabilization

The face density is first computed from the isentropic Bernoulli relation:

```text
rho_isen = (1 + c * (1 - u_face^2 - v_face^2))^(1/(gamma-1))
```

Then the transonic stabilizing blend is applied:

```text
rho_face = rho_isen - mu_upwind * (rho_isen - rho_upwind)
```

where:

- `rho_upwind` is a smooth or hard upwind blend of the neighboring cell densities
- `mu_upwind` is the corresponding blend of the cell switch `mu`

This matches the logic in [`artificialDensity()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L874).

The AD kernel uses the same hard transonic activation `max(0, M^2-M_c^2)` and the same hard upwind selector based on `v·n` as the primal discretization, so the differentiated operator stays consistent with the primal one.

## 8. The Kutta residual being differentiated

The Kutta residual is defined in [`kuttaResidualForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L289):

```text
K(phi, Gamma) = Gamma - Gamma_computed(phi)
```

with

```text
Gamma_computed(phi) = phi_upper^ext - phi_lower^ext
```

where the upper and lower trailing-edge values are extrapolated from the cell-centered potential using the local reconstructed velocity.

So the full nonlinear discrete system is:

```text
R(phi, Gamma) = 0
K(phi, Gamma) = 0
```

## 9. Reduced exact phi-only Jacobian

The solver keeps only `phi` in the linear system. `Gamma` is eliminated through the Kutta equation.

Starting from:

```text
K(phi, Gamma(phi)) = 0
```

differentiation gives:

```text
dK/dphi + dK/dGamma * dGamma/dphi = 0
```

therefore:

```text
dGamma/dphi = -(dK/dphi) / (dK/dGamma)
```

Now differentiate `R(phi, Gamma(phi))`:

```text
dR_reduced/dphi = dR/dphi + dR/dGamma * dGamma/dphi
```

which becomes:

```text
J_reduced = dR/dphi - (dR/dGamma) * (dK/dphi) / (dK/dGamma)
```

This is the exact phi-only Jacobian of the reduced nonlinear map.

## 10. How `computeADReducedExactJacobian()` assembles it

The implementation is in [`computeADReducedExactJacobian()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L615).

The assembly sequence is:

1. Copy the current global `phi` field into a flat array `phiGlobal`
2. Copy the current circulation into `gammaGlobal`
3. Differentiate the scalar Kutta residual once
4. For each cell row, differentiate the scalar row residual once
5. Form the reduced row using the Schur-complement formula
6. Insert those entries into PETSc

### 10.1 Kutta differentiation

The first Enzyme call produces:

- `kuttaDphi = dK/dphi`
- `kuttaDgamma = dK/dGamma`

This is done once because the Kutta residual is scalar and the same for every reduced row.

### 10.2 Per-row differentiation

For each row `i`, Enzyme produces:

- `rowDphi = dR_i/dphi`
- `rowDgamma = dR_i/dGamma`

The reduced row entry for column `j` is then computed as:

```text
J_ij = rowDphi[j] - (rowDgamma / kuttaDgamma) * kuttaDphi[j]
```

That line is the concrete implementation of the reduced exact Jacobian.

## 11. Why the matrix uses an explicit stencil

Although Enzyme returns derivatives with respect to the full `phiGlobal` array, most entries are structurally zero or negligibly small.

To avoid storing a dense row, the code builds a sparse support using:

- [`buildRowStencil()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L1048)
- [`buildKuttaStencil()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L1081)

The row stencil includes:

- the row cell,
- its one-ring neighbors,
- neighbors of those neighbors,
- wake-related TE influence cells when needed.

The Kutta stencil is the union of the upper-TE and lower-TE row stencils.

The PETSc matrix is preallocated for that wider support in [`initializeJacobian()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L326).

## 12. How `Gamma` is kept consistent after a phi update

Even though the matrix solve is phi-only, the primal operator still needs a consistent current circulation.

That is handled in [`enforceConsistentGammaForCurrentPhi()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L583).

This routine performs a small scalar Newton solve on:

```text
K(phi, Gamma) = 0
```

using Enzyme again to compute `dK/dGamma`.

So after every update of `phi`, the solver recomputes the circulation that makes the Kutta relation hold for the current potential field.

## 13. Where the AD Jacobian enters the nonlinear solve

The nonlinear loop is in [`solve()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L1276).

For `JACOBIAN_TYPE="ad_reduced_exact"` the iteration is:

1. compute the reduced exact Jacobian
2. assemble the right-hand side `-R(phi)`
3. solve the linear system with PETSc GMRES
4. update `phi`
5. enforce consistent `Gamma(phi)`
6. rerun the primal residual through `spatialDiscretization.run()`
7. accept or backtrack the step with line search if needed

So AD is used only for Jacobian assembly and for scalar Kutta resynchronization. The residual itself is still always recomputed by the primal solver path.

## 14. AD row-check and validation tools

The code includes built-in row validation controls near the top of [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L18):

- `CHECK_AD_ROW`
- `AD_ROW`
- `AD_ROW_CHECK_ONLY`
- `AD_FD_EPS`
- `AD_ROW_PRINT_TOL`
- `AD_ROW_MAX_PRINT`

The validation routines are:

- [`runADRowCheck()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L1112)
- `runADKuttaCheck()` in the same section

`runADRowCheck()` compares, for one selected row:

- the AD derivative
- a centered finite-difference derivative
- the matrix entry currently assembled in PETSc

This is the best debugging entry point if the AD Jacobian ever looks suspicious.

## 15. Why this approach works well for this solver

This AD design is a good fit because:

- the residual is naturally face-local and row-assembling,
- the true stencil is wider than the old analytical approximation,
- the Kutta relation can be eliminated exactly with a scalar Schur complement,
- PETSc still receives an assembled sparse matrix,
- the primal solver structure does not need to be rewritten around matrix-free products.

In short, the code gets an exact phi-only Jacobian of the implemented discrete operator without turning the whole solver into a giant differentiated program.

## 16. Practical mental model

The easiest way to think about the implementation is:

1. the primal code defines the nonlinear map,
2. `residualRowForAD()` and `kuttaResidualForAD()` restate that map in scalar form,
3. Enzyme provides `dR_i/dphi`, `dR_i/dGamma`, `dK/dphi`, and `dK/dGamma`,
4. the solver assembles

```text
J_reduced = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)
```

5. PETSc solves the resulting phi-only Newton step.

That is the whole AD implementation in one line.
