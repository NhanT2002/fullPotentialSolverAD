# Automatic Differentiation Implementation

This note describes the Enzyme-based reduced exact Jacobian used by

- `JACOBIAN_TYPE="ad_reduced_exact"`

and how it relates to the hand-coded exact analytical Jacobian.

## 1. Role of the AD Jacobian

In the current solver, the AD Jacobian serves two purposes:

1. it is a fully discrete exact reference for the reduced phi-only Jacobian,
2. it remains available as a runtime Jacobian mode for verification and debugging.

The analytical reduced exact Jacobian is the faster production implementation, but the AD path is still the cleanest reference because it differentiates the actual discrete operator.

## 2. Build-Time Wiring

Enzyme is inserted into the Chapel build through the split compile/link flow in [`Makefile`](/home/user/test/fullPotentialSolver/Makefile#L1):

1. Chapel emits LLVM bitcode.
2. `opt` loads the Enzyme plugin and runs the Enzyme pass.
3. The transformed IR is assembled back to bitcode.
4. Chapel finishes the final link step.

So Enzyme acts on Chapel-generated LLVM IR, not on source-level Chapel code.

## 3. What Is Differentiated

The implementation differentiates two scalar functions from [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L1):

- `residualRowForAD(sd, row, phi, Gamma, n)`
- `kuttaResidualForAD(sd, phi, Gamma, n)`

These functions reproduce the same discrete operator as the primal solver, but in scalar form.

That means AD is not differentiating the whole nonlinear driver or `spatialDiscretization.run()` directly. Instead, it differentiates the row-local residual operator and the scalar Kutta operator that define the Newton system.

## 4. Why Scalar Kernels Are Used

Differentiating the entire primal residual assembly would be awkward and expensive because it:

- mutates many arrays,
- computes all rows even when only one row derivative is needed,
- mixes geometry/state updates with the actual residual definition.

The scalar-kernel design keeps the differentiated function small and local:

- `sd` carries geometry and precomputed coefficients as constant data,
- `phi` and `Gamma` are the active variables,
- each call returns one scalar residual.

That structure maps naturally to row-by-row sparse matrix assembly.

## 5. Enzyme Interface

The Chapel-side Enzyme ABI is declared in `temporalDiscretization.chpl` with:

```chapel
extern {
    int enzyme_dup;
    int enzyme_const;
    void __enzyme_autodiff(void*, ...);
}
```

The conventions are:

- `enzyme_const`: argument is treated as constant
- `enzyme_dup`: argument is active and paired with derivative storage

So a call of the form

```chapel
__enzyme_autodiff(f,
                  enzyme_const, a,
                  enzyme_dup, x, dx);
```

means “differentiate `f(a, x)` with respect to `x` and store the result in `dx`”.

## 6. Discrete Operator Reproduced by the AD Kernels

The primal path in `spatialDiscretization.run()` does the following:

1. reconstruct cell velocities from `phi`,
2. compute density and transonic switch `mu`,
3. update ghost-cell velocities,
4. reconstruct corrected face velocities,
5. apply artificial-density stabilization,
6. compute face fluxes,
7. sum cell residuals and the Kutta residual.

The AD kernels reproduce this same operator with local helper routines:

- ghost-state accessors,
- least-squares gradient reconstruction,
- Bernoulli density recovery,
- wall/farfield ghost velocity treatment,
- wake jump logic,
- face density blending.

So the AD Jacobian is exact for the implemented discrete residual on the active branch of the hard switches.

## 7. Coupled Nonlinear System

The discrete steady system is:

- `R(phi, Gamma) = 0`
- `K(phi, Gamma) = 0`

where:

- `R_i` is the flux-balance residual for cell `i`,
- `K` is the scalar Kutta residual.

The AD kernels provide:

- `dR_i/dphi`
- `dR_i/dGamma`
- `dK/dphi`
- `dK/dGamma`

## 8. Reduced Phi-Only Jacobian

The linear system solved by the code keeps only `phi` as an unknown. `Gamma` is eliminated through the Kutta relation.

Differentiating `K(phi, Gamma(phi)) = 0` gives:

`dGamma/dphi = -(dK/dphi)/(dK/dGamma)`

Substituting that into the residual derivative gives the reduced exact Jacobian:

`J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`

This is the matrix assembled in the AD path.

## 9. AD Assembly Algorithm

`computeADReducedExactJacobian()` follows this sequence:

1. Copy the current `phi` state into a flat active array.
2. Copy the current circulation into an active scalar.
3. Differentiate `kuttaResidualForAD()` once to obtain:
   - `kuttaDphi`
   - `kuttaDgamma`
4. For each residual row, differentiate `residualRowForAD()` once to obtain:
   - `rowDphi`
   - `rowDgamma`
5. Assemble the reduced row:

   `J_ij = rowDphi[j] - (rowDgamma / kuttaDgamma) * kuttaDphi[j]`

6. Insert those entries into PETSc on the expected sparse stencil.

## 10. Why the Matrix Is Sparse

Although AD returns derivatives with respect to the full `phi` array, the discrete row only depends on:

- the row cell,
- first-ring neighbors,
- second-ring neighbors through least-squares gradients,
- trailing-edge support entering through the Kutta reduction.

So the solver builds an explicit row stencil and only inserts those columns into PETSc.

This is the main reason the AD path is feasible at all: the differentiated row is mathematically wider than the original approximate Jacobian, but it is still sparse enough to assemble as a matrix.

## 11. Gamma Re-Synchronization

Even in the reduced phi-only solve, the primal residual still needs a consistent current `Gamma`.

After every Newton update of `phi`, the solver calls `enforceConsistentGammaForCurrentPhi()`, which solves the scalar Kutta equation

`K(phi, Gamma) = 0`

with a short Newton iteration. This keeps the primal residual consistent with the reduced Jacobian formulation.

## 12. Why AD Is Slower Than the Analytical Exact Jacobian

The AD path is slower because it does generic differentiation row by row:

- Enzyme is called once per row residual and once for the Kutta residual,
- the differentiated row lives in dense derivative arrays before sparse extraction,
- the implementation cannot exploit as much problem-specific sparsity as the hand-coded exact path.

By contrast, the analytical exact Jacobian uses:

- precomputed sparse velocity-sensitivity templates,
- active-index loops,
- hand-coded chain-rule expressions,
- row-parallel assembly with less generic derivative overhead.

So AD is best understood as the exact reference implementation, while the analytical exact Jacobian is the optimized production implementation.

## 13. Validation Hooks

The code includes row-check diagnostics controlled by:

- `CHECK_AD_ROW`
- `AD_ROW`
- `AD_ROW_CHECK_ONLY`

These compare:

- AD derivative,
- finite-difference derivative,
- and, when available, the assembled matrix row.

That machinery is what made it possible to derive and trust the hand-coded exact Jacobian.

## 14. Relationship to the Analytical Exact Jacobian

Both exact modes target the same reduced matrix:

`J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`

The difference is only how the derivatives are produced:

- `ad_reduced_exact`: Enzyme computes them automatically,
- `analytical_reduced_exact`: they are written out by hand and assembled with sparse optimized loops.

The exact analytical derivation and implementation details are documented in [`EXACT_ANALYTICAL_JACOBIAN.md`](/home/user/test/fullPotentialSolver/EXACT_ANALYTICAL_JACOBIAN.md#L1).
