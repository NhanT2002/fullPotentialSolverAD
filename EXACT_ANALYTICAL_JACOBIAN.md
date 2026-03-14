# Exact Analytical Jacobian Derivation

This note derives the exact analytical Jacobian of the current discrete full-potential operator, in the same reduced phi-only form used by the Enzyme implementation.

The goal is not to replace the working AD Jacobian immediately. The goal is to make the full chain rule explicit so we can:

1. implement a hand-coded exact Jacobian if desired,
2. compare that hand-coded Jacobian against the existing AD Jacobian,
3. use the AD Jacobian as the reference for debugging every analytical term.

## 1. Nonlinear system being linearized

The discrete steady system solved by the code is:

```text
R(phi, Gamma) = 0
K(phi, Gamma) = 0
```

where:

- `R_i(phi, Gamma)` is the cell residual from the face-flux balance
- `K(phi, Gamma)` is the Kutta residual

The current AD implementation eliminates `Gamma` and assembles the reduced exact phi-only Jacobian:

```text
J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)
```

So an exact analytical Jacobian should target the same reduced matrix. If the derivation is correct, it should match the AD Jacobian row by row on the active discrete branch.

## 2. Residual row definition

For a cell `i`,

```text
R_i = s_i * sum_{f in i} F_f
```

where `s_i` is the sign induced by the orientation of each face and

```text
F_f = A_f * rho_f * q_f
q_f = V_f · n_f
```

with:

- `A_f` face area
- `rho_f` face density after transonic blending
- `V_f = (u_f, v_f)` corrected face velocity
- `n_f = (n_x, n_y)` face normal

This is the same operator evaluated in [`residualRowForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L207) and in the primal pipeline through [`computeFaceProperties()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L811), [`artificialDensity()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L874), and [`computeResiduals()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L943).

## 3. Face-velocity reconstruction

For a face between states `1` and `2`:

```text
u_avg = w1 * u_1 + w2 * u_2
v_avg = w1 * v_1 + w2 * v_2
```

The potential jump used by the direct correction is:

```text
DeltaPhi_f = phi_2 - phi_1 + s_Gamma,f * Gamma
```

where:

- `s_Gamma,f = +1` for an upper-to-lower wake crossing,
- `s_Gamma,f = -1` for a lower-to-upper wake crossing,
- `s_Gamma,f = 0` otherwise.

Then:

```text
dphi/dl = invL_f * DeltaPhi_f
delta_f = u_avg * t_x + v_avg * t_y - dphi/dl
u_f = u_avg - delta_f * c_x
v_f = v_avg - delta_f * c_y
```

where:

- `(t_x, t_y)` is the unit vector along the cell-to-cell direction,
- `(c_x, c_y) = corrCoeff` is the deferred-correction coefficient.

It is often cleaner to work directly with:

```text
q_f = V_f · n
    = m_x * u_avg + m_y * v_avg + d_f * DeltaPhi_f
```

where:

```text
m = n - k t
k = 1 / (n · t)
d_f = k * invL_f
```

This is exactly the same simplification already described in the approximate Jacobian comments in [`computeJacobian()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L665).

## 4. Cell-velocity derivatives

The reconstructed cell velocity is a QR least-squares gradient:

```text
u_e = sum_{neighbors l of e} W^x_{e,l} * (phi_l^corr - phi_e)
v_e = sum_{neighbors l of e} W^y_{e,l} * (phi_l^corr - phi_e)
```

where `phi_l^corr` includes the wake jump when the edge `(e,l)` crosses the wake.

Therefore the exact derivatives are purely geometric coefficients:

```text
du_e/dphi_j = G^x_{e,j}
dv_e/dphi_j = G^y_{e,j}
du_e/dGamma = G^x_{e,Gamma}
dv_e/dGamma = G^y_{e,Gamma}
```

with:

- `G^x_{e,j}` and `G^y_{e,j}` assembled from the QR weights
- `G^x_{e,Gamma}` and `G^y_{e,Gamma}` given by the wake-crossing gradient sensitivities already computed in [`computeGradientSensitivity()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L441)

For a hand-coded exact Jacobian, the cleanest formulation is to explicitly view the gradient operator as a sparse linear map in `phi` and `Gamma`.

## 5. Derivative of face normal velocity

For any active variable `s` in `{phi_j, Gamma}`:

```text
dq_f/ds = m_x * du_avg/ds + m_y * dv_avg/ds + d_f * d(DeltaPhi_f)/ds
```

with:

```text
du_avg/ds = w1 * du_1/ds + w2 * du_2/ds
dv_avg/ds = w1 * dv_1/ds + w2 * dv_2/ds
```

and:

```text
d(DeltaPhi_f)/dphi_j = delta_{2j} - delta_{1j}
d(DeltaPhi_f)/dGamma = s_Gamma,f
```

So the exact face-normal-velocity derivative is already manageable once the cell-gradient derivatives are known.

## 6. Isentropic density derivative

The face isentropic density is:

```text
rho_isen,f = B_f^p
B_f = 1 + a * (1 - u_f^2 - v_f^2)
a = ((gamma - 1)/2) * M_inf^2
p = 1 / (gamma - 1)
```

Differentiate:

```text
drho_isen,f/ds = p * B_f^(p-1) * dB_f/ds
```

with:

```text
dB_f/ds = -2 a (u_f * du_f/ds + v_f * dv_f/ds)
```

so:

```text
drho_isen,f/ds =
-2 a p B_f^(p-1) * (u_f * du_f/ds + v_f * dv_f/ds)
```

This is one of the main terms missing from the old approximate analytical Jacobian.

## 7. Cell density derivative

Each cell density is:

```text
rho_e = (1 + a * (1 - u_e^2 - v_e^2))^p
```

so by the same chain rule:

```text
drho_e/ds =
-2 a p B_e^(p-1) * (u_e * du_e/ds + v_e * dv_e/ds)
```

These cell-density derivatives are needed because the face density blend depends on the upwind cell density.

## 8. Mach-squared and mu derivative

The code defines:

```text
M_e = M_inf * sqrt(u_e^2 + v_e^2) * rho_e^((1-gamma)/2)
```

from [`mach()`](/home/user/test/fullPotentialSolver/src/spatialDiscretization.chpl#L985), so:

```text
M_e^2 = M_inf^2 * (u_e^2 + v_e^2) * rho_e^(1-gamma)
```

Differentiate:

```text
d(M_e^2)/ds =
M_inf^2 * [ 2 u_e rho_e^(1-gamma) du_e/ds
          + 2 v_e rho_e^(1-gamma) dv_e/ds
          + (u_e^2 + v_e^2)(1-gamma) rho_e^(-gamma) drho_e/ds ]
```

The transonic switch is:

```text
mu_e = MU_C * max(0, M_e^2 - M_c^2)
```

Therefore, on a fixed active branch:

```text
dmu_e/ds =
0                                  if M_e^2 <= M_c^2
MU_C * d(M_e^2)/ds                 if M_e^2 >  M_c^2
```

This is piecewise exact, which matches what AD gives on the executed branch.

## 9. Upwind branch and face-density derivative

With the current hard upwind switch, the branch is:

```text
if q_f >= 0:
    upwind = elem1
else:
    upwind = elem2
```

On a fixed active branch, the selector itself is treated as constant. So:

```text
rho_up,f = rho_upwindCell
mu_up,f  = mu_upwindCell
```

and:

```text
drho_up,f/ds = drho_upwindCell/ds
dmu_up,f/ds  = dmu_upwindCell/ds
```

The blended face density is:

```text
rho_f = rho_isen,f - mu_up,f * (rho_isen,f - rho_up,f)
```

Equivalently:

```text
rho_f = (1 - mu_up,f) rho_isen,f + mu_up,f rho_up,f
```

Differentiate:

```text
drho_f/ds =
(1 - mu_up,f) drho_isen,f/ds
+ mu_up,f drho_up,f/ds
+ (rho_up,f - rho_isen,f) dmu_up,f/ds
```

Special case:

```text
if mu_up,f <= 0 then rho_f = rho_isen,f
```

so:

```text
drho_f/ds = drho_isen,f/ds
```

again on the active branch.

## 10. Exact face-flux derivative

The face flux is:

```text
F_f = A_f * rho_f * q_f
```

Therefore:

```text
dF_f/ds = A_f * [ q_f * drho_f/ds + rho_f * dq_f/ds ]
```

This is the core formula for an exact analytical Jacobian.

Everything below it is just expansion of:

- `dq_f/ds`
- `drho_f/ds`

through the least-squares reconstruction, Bernoulli density, transonic activation, and upwind density branch.

## 11. Exact residual-row derivative

For a cell row `i`:

```text
dR_i/ds = res_scale * sum_{f in i} sign(i,f) * dF_f/ds
```

where `s` is either:

- `phi_j`
- `Gamma`

This yields the exact coupled derivatives:

```text
dR_i/dphi_j
dR_i/dGamma
```

That is the analytical counterpart of what Enzyme currently returns as:

- `rowDphi`
- `rowDgamma`

inside [`computeADReducedExactJacobian()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L615).

## 12. Exact Kutta derivative

The Kutta residual is:

```text
K = res_scale * [ Gamma - phi_u^ext + phi_l^ext ]
```

with:

```text
phi_u^ext = phi_u + u_u * dSux + v_u * dSuy
phi_l^ext = phi_l + u_l * dSlx + v_l * dSly
```

Therefore:

```text
dK/dphi_j = res_scale * [
    -delta_{u,j}
    - dSux * du_u/dphi_j
    - dSuy * dv_u/dphi_j
    + delta_{l,j}
    + dSlx * du_l/dphi_j
    + dSly * dv_l/dphi_j
]
```

and:

```text
dK/dGamma = res_scale * [
    1
    - dSux * du_u/dGamma
    - dSuy * dv_u/dGamma
    + dSlx * du_l/dGamma
    + dSly * dv_l/dGamma
]
```

This is important: `dK/dGamma` is not just `res_scale`. The extrapolated trailing-edge velocities also depend on `Gamma` through the wake-corrected least-squares gradients.

That exact dependence is already captured automatically by AD in [`kuttaResidualForAD()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L289).

## 13. Reduced exact phi-only analytical Jacobian

Once the coupled exact derivatives are known, the phi-only Jacobian is:

```text
J_red(i,j) = dR_i/dphi_j - (dR_i/dGamma) * (dK/dphi_j) / (dK/dGamma)
```

This is the exact same reduced formula currently assembled by the AD implementation.

So if a hand-coded analytical Jacobian is implemented correctly, it should satisfy:

```text
J_analytical_exact(i,j) = J_AD(i,j)
```

up to roundoff, for all rows and all columns on the current active branch.

## 14. Why the old analytical Jacobian is not exact

The current non-AD analytical Jacobian in [`computeJacobian()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L665) already differentiates:

- the QR gradient reconstruction,
- the corrected face-velocity formula,
- the direct wake-jump term,
- an approximate circulation sensitivity route.

But it freezes or omits:

- `drho_f/ds`
- `drho_e/ds`
- `dmu_e/ds`
- the exact reduced Kutta Schur complement

So it is structurally useful, but not the exact Jacobian of the implemented nonlinear map.

## 15. How to implement the exact analytical Jacobian

The least risky implementation path is incremental.

### Step 1

Keep the current row-stencil logic and reduced phi-only matrix structure.

### Step 2

Introduce analytical helper routines for:

- `du_e/dphi_j`, `dv_e/dphi_j`
- `du_e/dGamma`, `dv_e/dGamma`
- `drho_e/ds`
- `dmu_e/ds`
- `dq_f/ds`
- `drho_f/ds`

### Step 3

Build a coupled exact row routine returning:

- `dR_i/dphi_j`
- `dR_i/dGamma`

### Step 4

Build exact Kutta derivatives:

- `dK/dphi_j`
- `dK/dGamma`

### Step 5

Assemble:

```text
J_red(i,j) = dR_i/dphi_j - (dR_i/dGamma) * (dK/dphi_j) / (dK/dGamma)
```

### Step 6

Compare every row against AD using the existing check framework in [`runADRowCheck()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L1112).

## 16. Comparison strategy against AD

Yes, the AD Jacobian is exactly what we should compare against.

The best validation ladder is:

1. Compare one analytical row against AD on a subsonic case.
2. Compare one wake-affected row against AD.
3. Compare the Kutta derivatives against AD.
4. Compare a full matrix row norm against AD.
5. Only then trust nonlinear convergence comparisons.

The reason this is a good strategy is that AD already differentiates the discrete code you are actually running. So the AD Jacobian is the most reliable reference for the analytical derivation.

There is now a row-local implementation of this comparison in [`src/temporalDiscretization.chpl`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl):

- [`residualRowDerivativeAnalytical()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L379)
- [`kuttaResidualDerivativeAnalytical()`](/home/user/test/fullPotentialSolver/src/temporalDiscretization.chpl#L492)

These routines hand-code the active-branch chain rule for one scalar seed and are used by the row-check diagnostics to compare:

- exact analytical row derivative,
- reduced AD row derivative,
- finite-difference row derivative.

So the project now has an implemented reference path for validating the analytical derivation row by row, even though the global matrix assembly is still driven by AD.

That is no longer just a validation-only path: the solver now also exposes a full global assembly mode,

- `JACOBIAN_TYPE = "analytical_reduced_exact"`,

which assembles the reduced phi-only matrix from hand-coded exact derivatives. The current implementation performs this assembly face by face on each row stencil, so it mirrors the structure of the original hand-coded `analytical` Jacobian much more closely than the earlier seed-based exact assembly prototype. This makes it possible to compare the global analytical reduced-exact solve directly against `ad_reduced_exact`.

## 17. Bottom line

Yes, the full exact analytical Jacobian can be derived for the current solver.

The exact reduced phi-only form is:

```text
J_red = dR/dphi - dR/dGamma * (dK/dphi) / (dK/dGamma)
```

and the missing work is to hand-code the full chain rule through:

- least-squares gradients,
- corrected face velocity,
- Bernoulli density,
- transonic switch `mu`,
- upwind density blending,
- and the Kutta extrapolation.

Once implemented, the correct reference for validation is the existing AD Jacobian.
