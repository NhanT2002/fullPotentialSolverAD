# Exact Analytical Jacobian Derivation

This note documents the hand-coded exact analytical Jacobian used by

- `JACOBIAN_TYPE="analytical_reduced_exact"`

and explains how it is derived from the discrete residual implemented by the solver.

The most important point is that this Jacobian is not a continuous PDE Jacobian. It is the exact Jacobian of the **implemented discrete operator** on the active branch of the hard switches.

## 1. Nonlinear System

The steady discrete system is written as

- `R(phi, Gamma) = 0`
- `K(phi, Gamma) = 0`

where:

- `phi` is the cell-centered potential,
- `Gamma` is the wake circulation jump,
- `R_i` is the residual row for cell `i`,
- `K` is the scalar Kutta residual.

The solver ultimately keeps only `phi` in the linear system, so the exact analytical Jacobian is assembled in reduced form.

## 2. Reduced Phi-Only Exact Jacobian

Starting from the Kutta relation

`K(phi, Gamma(phi)) = 0`

we differentiate to get

`dGamma/dphi = -(dK/dphi)/(dK/dGamma)`

Then for the residual:

`dR_red/dphi = dR/dphi + dR/dGamma * dGamma/dphi`

which gives the reduced exact Jacobian

`J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`

This is exactly the same reduced matrix assembled by the AD implementation.

So the job of the analytical derivation is to obtain:

- `dR_i/dphi_j`
- `dR_i/dGamma`
- `dK/dphi_j`
- `dK/dGamma`

and then assemble the Schur-complement correction above.

## 3. Discrete Residual Row

For a cell `i`,

`R_i = res_scale * sum_{f in i} sign(i,f) * F_f`

where the face flux is

`F_f = A_f * rho_f * q_f`

with:

- `A_f`: face area,
- `q_f = V_f . n_f`,
- `V_f = (u_f, v_f)`,
- `rho_f`: stabilized face density.

The derivation therefore reduces to differentiating each face contribution.

## 4. Cell Velocity Reconstruction

For each cell `e`, the least-squares QR gradient operator gives

`u_e = sum_l Gx_{e,l} * phi_l + Gx_{e,Gamma} * Gamma`

`v_e = sum_l Gy_{e,l} * phi_l + Gy_{e,Gamma} * Gamma`

where:

- `Gx_{e,l}`, `Gy_{e,l}` are geometric coefficients from the QR reconstruction,
- `Gx_{e,Gamma}`, `Gy_{e,Gamma}` come from wake-crossing gradient sensitivity.

That means the exact cell-velocity derivatives are linear:

`du_e/dphi_j = Gx_{e,j}`

`dv_e/dphi_j = Gy_{e,j}`

`du_e/dGamma = Gx_{e,Gamma}`

`dv_e/dGamma = Gy_{e,Gamma}`

In the implementation, these coefficients are not rebuilt from scratch during every row assembly. They are cached in sparse template form and then injected into the current row stencil.

## 5. Face Potential Jump

For a face between states `1` and `2`,

`DeltaPhi_f = phi_2 - phi_1 + s_f * Gamma`

where `s_f` is:

- `+1` for one wake-crossing orientation,
- `-1` for the opposite wake-crossing orientation,
- `0` otherwise.

So the direct derivatives are

`d(DeltaPhi_f)/dphi_j = delta_{2j} - delta_{1j}`

`d(DeltaPhi_f)/dGamma = s_f`

## 6. Face Velocity Reconstruction

The solver uses the corrected face velocity

`u_avg = w1 * u_1 + w2 * u_2`

`v_avg = w1 * v_1 + w2 * v_2`

`dphi/dl = invL_f * DeltaPhi_f`

`delta_f = u_avg * t_x + v_avg * t_y - dphi/dl`

`u_f = u_avg - delta_f * c_x`

`v_f = v_avg - delta_f * c_y`

where:

- `(t_x, t_y)` is the unit cell-to-cell direction,
- `(c_x, c_y)` is the deferred-correction coefficient.

For any active scalar `s` in `{phi_j, Gamma}`:

`du_avg/ds = w1 * du_1/ds + w2 * du_2/ds`

`dv_avg/ds = w1 * dv_1/ds + w2 * dv_2/ds`

`d(delta_f)/ds = t_x * du_avg/ds + t_y * dv_avg/ds - invL_f * d(DeltaPhi_f)/ds`

`du_f/ds = du_avg/ds - c_x * d(delta_f)/ds`

`dv_f/ds = dv_avg/ds - c_y * d(delta_f)/ds`

These are the exact analytical face-velocity derivatives used by the hand-coded exact Jacobian.

## 7. Face Normal Velocity

The flux uses

`q_f = u_f * n_x + v_f * n_y`

so

`dq_f/ds = n_x * du_f/ds + n_y * dv_f/ds`

This is the first part of the face-flux derivative.

## 8. Isentropic Density Derivative

The isentropic face density is

`rho_isen,f = B_f^p`

with

`B_f = 1 + a * (1 - u_f^2 - v_f^2)`

`a = ((gamma - 1)/2) * M_inf^2`

`p = 1/(gamma - 1)`

Differentiate:

`drho_isen,f/ds = p * B_f^(p-1) * dB_f/ds`

`dB_f/ds = -2 a (u_f * du_f/ds + v_f * dv_f/ds)`

so

`drho_isen,f/ds = -2 a p B_f^(p-1) * (u_f * du_f/ds + v_f * dv_f/ds)`

## 9. Cell Density and Mu Derivatives

For each cell,

`rho_e = [1 + a * (1 - u_e^2 - v_e^2)]^p`

so

`drho_e/ds = -2 a p B_e^(p-1) * (u_e * du_e/ds + v_e * dv_e/ds)`

The transonic switch is

`M_e^2 = M_inf^2 * (u_e^2 + v_e^2) * rho_e^(1-gamma)`

and

`mu_e = MU_C * max(0, M_e^2 - M_C^2)`

On the active branch:

- if `M_e^2 <= M_C^2`, then `dmu_e/ds = 0`
- if `M_e^2 > M_C^2`, then `dmu_e/ds = MU_C * d(M_e^2)/ds`

with

`d(M_e^2)/ds = M_inf^2 * [`

`2 u_e rho_e^(1-gamma) du_e/ds`

`+ 2 v_e rho_e^(1-gamma) dv_e/ds`

`+ (u_e^2 + v_e^2)(1-gamma) rho_e^(-gamma) drho_e/ds ]`

This is where the exact Jacobian differs strongly from the old analytical approximate Jacobian: the nonlinear density/switch response is carried explicitly.

## 10. Upwind Branch

The face stabilization uses a hard upwind branch selected by the sign of `q_f`.

On the active branch:

- the upwind cell index is treated as fixed,
- the derivative does not include a derivative of the branch selector itself.

So if the active upwind cell is `u`,

`rho_up,f = rho_u`

`mu_up,f = mu_u`

and

`drho_up,f/ds = drho_u/ds`

`dmu_up,f/ds = dmu_u/ds`

This is exactly the same piecewise derivative that AD computes on the executed branch.

## 11. Stabilized Face Density Derivative

The stabilized density is

`rho_f = rho_isen,f - mu_up,f * (rho_isen,f - rho_up,f)`

or equivalently

`rho_f = (1 - mu_up,f) rho_isen,f + mu_up,f rho_up,f`

Differentiate:

`drho_f/ds =`

`(1 - mu_up,f) drho_isen,f/ds`

`+ mu_up,f drho_up,f/ds`

`+ (rho_up,f - rho_isen,f) dmu_up,f/ds`

If `mu_up,f <= 0`, this reduces to

`drho_f/ds = drho_isen,f/ds`

## 12. Exact Face-Flux Derivative

The face flux is

`F_f = A_f * rho_f * q_f`

therefore

`dF_f/ds = A_f * [ q_f * drho_f/ds + rho_f * dq_f/ds ]`

This is the exact analytical face contribution used in the row assembly.

## 13. Exact Residual-Row Derivatives

For any scalar `s`,

`dR_i/ds = res_scale * sum_{f in i} sign(i,f) * dF_f/ds`

Choosing `s = phi_j` gives

`dR_i/dphi_j`

Choosing `s = Gamma` gives

`dR_i/dGamma`

Those are the coupled exact row derivatives before the Kutta reduction.

## 14. Exact Kutta Derivatives

The Kutta residual is

`K = res_scale * [ Gamma - phi_u^ext + phi_l^ext ]`

with extrapolated trailing-edge values

`phi_u^ext = phi_u + u_u * dSux + v_u * dSuy`

`phi_l^ext = phi_l + u_l * dSlx + v_l * dSly`

Therefore

`dK/dphi_j = res_scale * [`

`-delta_{u,j}`

`- dSux * du_u/dphi_j`

`- dSuy * dv_u/dphi_j`

`+ delta_{l,j}`

`+ dSlx * du_l/dphi_j`

`+ dSly * dv_l/dphi_j ]`

and

`dK/dGamma = res_scale * [`

`1`

`- dSux * du_u/dGamma`

`- dSuy * dv_u/dGamma`

`+ dSlx * du_l/dGamma`

`+ dSly * dv_l/dGamma ]`

So `dK/dGamma` is not just `1`; the TE extrapolation also depends on `Gamma` through the wake-sensitive gradients.

## 15. Final Reduced Exact Jacobian

Once the coupled exact derivatives are known, the final phi-only Jacobian is assembled as

`J_red(i,j) = dR_i/dphi_j - (dR_i/dGamma) * (dK/dphi_j) / (dK/dGamma)`

This is the exact reduced analytical Jacobian implemented by

- `JACOBIAN_TYPE="analytical_reduced_exact"`

and it is the same matrix targeted by the AD path.

## 16. Implemented Sparse Assembly Strategy

The current implementation does not build the analytical Jacobian with dense per-row algebra. It uses a sparse optimized strategy:

1. Build a row stencil containing the row cell, nearby support cells, and Kutta-support cells.
2. Reuse cached sparse cell velocity-sensitivity templates.
3. For each face in the row, form a small active index set of only the local stencil entries touched by that face.
4. Evaluate the exact chain rule only on those active local indices.
5. Accumulate:
   - local `dR_i/dphi`
   - local `dR_i/dGamma`
6. Apply the reduced Kutta correction on the merged row/Kutta support.
7. Insert the row into PETSc.

This is why the analytical exact Jacobian is both exact and much faster than the AD assembly.

## 17. Parallelism and Performance

The exact analytical assembly is parallelized across rows with Chapel `forall`. PETSc insertion is still performed after the row-local work is prepared, but the expensive derivative evaluation is row-parallel.

Performance improvements in the current implementation come from:

- sparse row stencils,
- cached velocity-sensitivity templates,
- active-index loops instead of dense row sweeps,
- avoiding repeated generic AD calls.

That is the main reason `analytical_reduced_exact` is substantially faster than `ad_reduced_exact` in practice.

## 18. Exactness and Branches

The analytical exact Jacobian is exact for the discrete operator on the current executed branch. Because the solver still contains hard switches:

- upwind branch from the sign of `q_f`,
- transonic activation through `max(0, M^2 - M_C^2)`,

the Jacobian is piecewise exact rather than globally smooth.

This is the same notion of exactness used by the AD implementation.

## 19. Validation Against AD

The row-check infrastructure was used to compare:

- hand-coded analytical exact derivatives,
- Enzyme derivatives,
- finite differences.

That validation was done before promoting the hand-coded path to a full runtime Jacobian mode.

So the intended interpretation is:

- AD is the exact reference implementation,
- analytical exact is the optimized implementation of the same reduced matrix.

## 20. Bottom Line

The analytical exact Jacobian is derived by differentiating the full discrete face-flux residual and Kutta relation through:

- least-squares velocity reconstruction,
- wake-corrected face potential jumps,
- deferred-correction face velocity,
- Bernoulli density,
- transonic switch `mu`,
- upwind density blending,
- trailing-edge extrapolation,

and then eliminating `Gamma` with the reduced formula

`J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)`.

That is the exact mathematical object now assembled in the fast hand-coded production Jacobian.
