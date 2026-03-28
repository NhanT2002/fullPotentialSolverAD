module temporalDiscretizationAD {
use CTypes;
use temporalDiscretization;
use temporalDiscretizationADIBM;
use temporalDiscretizationAnalyticalExact;
use Set;

// AD-based reduced exact Jacobian plus derivative verification helpers.

proc temporalDiscretization.pickADRow(): int {
    if AD_ROW > 0 then
        return AD_ROW;

    var count = 0;
    for elem in this.spatialDisc_.fluid_dom {
        if this.spatialDisc_.wakeFaceIndexInfluenceOnElem_[elem] == -1 {
            count += 1;
            if count == this.spatialDisc_.fluid_dom.size / 2 then
                return elem;
        }
    }

    for elem in this.spatialDisc_.fluid_dom do
        return elem;

    return 1;
}

proc temporalDiscretization.enforceConsistentGammaForCurrentPhi(maxIts: int = 12,
                                                                tol: real(64) = 1.0e-14) {
    if this.spatialDisc_.isIBMFlow() {
        this.enforceIBMConsistentGammaForCurrentPhi(maxIts, tol);
        return;
    }

    const sd = this.spatialDisc_.borrow();
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    const n = this.spatialDisc_.nelemDomain_: c_int;
    var phiGlobal: [phiDom] real(64);
    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    var gamma = this.spatialDisc_.circulation_;
    for gammaIt in 1..maxIts {
        const residual = kuttaResidualForAD(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gamma), n);
        if abs(residual) < tol then break;

        var dgamma = 0.0;
        __enzyme_autodiff(c_ptrTo(kuttaResidualForAD): c_ptr(void),
                          enzyme_const,
                          sd,
                          enzyme_const,
                          c_ptrTo(phiGlobal[0]),
                          enzyme_dup,
                          c_ptrTo(gamma),
                          c_ptrTo(dgamma),
                          enzyme_const,
                          n);

        if abs(dgamma) < 1.0e-14 then break;
        gamma -= residual / dgamma;
    }

    this.spatialDisc_.circulation_ = gamma;
}

proc temporalDiscretization.computeADReducedExactJacobian() {
    if this.spatialDisc_.isIBMFlow() {
        this.computeIBMADReducedExactJacobian();
        return;
    }

    this.A_petsc.zeroEntries();

    const sd = this.spatialDisc_.borrow();
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    const n = this.spatialDisc_.nelemDomain_: c_int;
    var phiGlobal: [phiDom] real(64);
    var gammaGlobal = this.spatialDisc_.circulation_;
    var kuttaDphi: [phiDom] real(64) = 0.0;
    var kuttaDgamma = 0.0;
    const kuttaStencil = this.buildKuttaStencil();

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    __enzyme_autodiff(c_ptrTo(kuttaResidualForAD): c_ptr(void),
                      enzyme_const,
                      sd,
                      enzyme_dup,
                      c_ptrTo(phiGlobal[0]),
                      c_ptrTo(kuttaDphi[0]),
                      enzyme_dup,
                      c_ptrTo(gammaGlobal),
                      c_ptrTo(kuttaDgamma),
                      enzyme_const,
                      n);

    if abs(kuttaDgamma) < 1.0e-14 then
        halt("Reduced exact Jacobian failed: dR_gamma/dGamma is too small");

    for row in 1..this.spatialDisc_.nelemDomain_ {
        if AD_JACOBIAN_PROGRESS_FREQ > 0 &&
           (row == 1 || row % AD_JACOBIAN_PROGRESS_FREQ == 0 || row == this.spatialDisc_.nelemDomain_) {
            writeln("AD reduced exact Jacobian assembly row ", row, " / ", this.spatialDisc_.nelemDomain_);
        }

        var rowDphi: [phiDom] real(64) = 0.0;
        var rowDgamma = 0.0;
        __enzyme_autodiff(c_ptrTo(residualRowForAD): c_ptr(void),
                          enzyme_const,
                          sd,
                          enzyme_const,
                          row,
                          enzyme_dup,
                          c_ptrTo(phiGlobal[0]),
                          c_ptrTo(rowDphi[0]),
                          enzyme_dup,
                          c_ptrTo(gammaGlobal),
                          c_ptrTo(rowDgamma),
                          enzyme_const,
                          n);

        var stencilSet = new set(int);
        const rowStencil = this.buildRowStencil(row);
        for col in rowStencil do stencilSet.add(col);
        for col in kuttaStencil do stencilSet.add(col);

        for col in stencilSet {
            const reducedValue = rowDphi[col - 1] - (rowDgamma / kuttaDgamma) * kuttaDphi[col - 1];
            if abs(reducedValue) > AD_ROW_PRINT_TOL || col == row then
                this.A_petsc.add(row - 1, col - 1, reducedValue);
        }
    }

    this.A_petsc.assemblyComplete();
}

proc temporalDiscretization.runADRowCheck() {
    if this.spatialDisc_.isIBMFlow() {
        this.runIBMADRowCheck();
        return;
    }

    if AD_ROW < 0 {
        this.runADKuttaCheck();
        return;
    }

    const row = this.pickADRow();
    writeln("Running AD row Jacobian check for row ", row);

    const stencil = this.buildRowStencil(row);
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    var phiGlobal: [phiDom] real(64);
    var dphiGlobal: [phiDom] real(64) = 0.0;
    var dgammaGlobal = 0.0;
    var gammaGlobal = this.spatialDisc_.circulation_;
    const reducedExactMode = isReducedExactJacobianType(this.inputs_.JACOBIAN_TYPE_);
    var kuttaDphi: [phiDom] real(64) = 0.0;
    var kuttaDgamma = 0.0;
    var exactRowDgamma = 0.0;
    var exactKuttaDgamma = 0.0;

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    writeln("AD row check: evaluating primal residual");
    const sd = this.spatialDisc_.borrow();
    const baseResidual = residualRowForAD(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                          this.spatialDisc_.nelemDomain_: c_int);
    writeln("AD row check: primal residual done, calling Enzyme");
    __enzyme_autodiff(c_ptrTo(residualRowForAD): c_ptr(void),
                      enzyme_const,
                      sd,
                      enzyme_const,
                      row,
                      enzyme_dup,
                      c_ptrTo(phiGlobal[0]),
                      c_ptrTo(dphiGlobal[0]),
                      enzyme_dup,
                      c_ptrTo(gammaGlobal),
                      c_ptrTo(dgammaGlobal),
                      enzyme_const,
                      this.spatialDisc_.nelemDomain_: c_int);

    if reducedExactMode {
        __enzyme_autodiff(c_ptrTo(kuttaResidualForAD): c_ptr(void),
                          enzyme_const,
                          sd,
                          enzyme_dup,
                          c_ptrTo(phiGlobal[0]),
                          c_ptrTo(kuttaDphi[0]),
                          enzyme_dup,
                          c_ptrTo(gammaGlobal),
                          c_ptrTo(kuttaDgamma),
                          enzyme_const,
                          this.spatialDisc_.nelemDomain_: c_int);
        exactRowDgamma = residualRowDerivativeAnalytical(sd, row, c_ptrTo(phiGlobal[0]),
                                                         c_ptrTo(gammaGlobal), 0, true,
                                                         this.spatialDisc_.nelemDomain_: c_int);
        exactKuttaDgamma = kuttaResidualDerivativeAnalytical(sd, c_ptrTo(phiGlobal[0]),
                                                             c_ptrTo(gammaGlobal), 0, true,
                                                             this.spatialDisc_.nelemDomain_: c_int);
    }

    writeln("AD row base residual = ", baseResidual, " stencil size = ", stencil.size);

    var maxADFD = 0.0;
    var maxAnalyticalAD = 0.0;
    var maxExactAD = 0.0;
    var maxMatrixExact = 0.0;
    var printed = 0;

    for col in stencil {
        const exactCoupled = residualRowDerivativeAnalytical(sd, row, c_ptrTo(phiGlobal[0]),
                                                             c_ptrTo(gammaGlobal), col, false,
                                                             this.spatialDisc_.nelemDomain_: c_int);
        const exactKutta = if reducedExactMode then
                               kuttaResidualDerivativeAnalytical(sd, c_ptrTo(phiGlobal[0]),
                                                                 c_ptrTo(gammaGlobal), col, false,
                                                                 this.spatialDisc_.nelemDomain_: c_int)
                           else
                               0.0;
        const exact = if reducedExactMode then
                          exactCoupled - (exactRowDgamma / exactKuttaDgamma) * exactKutta
                      else
                          exactCoupled;

        const ad = if reducedExactMode then
                       dphiGlobal[col - 1] - (dgammaGlobal / kuttaDgamma) * kuttaDphi[col - 1]
                   else
                       dphiGlobal[col - 1];
        const analytical = this.A_petsc.get(row - 1, col - 1);

        const saved = phiGlobal[col - 1];
        var fd: real(64);
        if reducedExactMode {
            phiGlobal[col - 1] = saved + AD_FD_EPS;
            var gammaPlus = solveConsistentGammaForPhi(sd, c_ptrTo(phiGlobal[0]), gammaGlobal,
                                                       this.spatialDisc_.nelemDomain_: c_int);
            const rPlus = residualRowForAD(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaPlus),
                                           this.spatialDisc_.nelemDomain_: c_int);

            phiGlobal[col - 1] = saved - AD_FD_EPS;
            var gammaMinus = solveConsistentGammaForPhi(sd, c_ptrTo(phiGlobal[0]), gammaGlobal,
                                                        this.spatialDisc_.nelemDomain_: c_int);
            const rMinus = residualRowForAD(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaMinus),
                                            this.spatialDisc_.nelemDomain_: c_int);
            phiGlobal[col - 1] = saved;
            fd = (rPlus - rMinus) / (2.0 * AD_FD_EPS);
        } else {
            phiGlobal[col - 1] = saved + AD_FD_EPS;
            const rPlus = residualRowForAD(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                           this.spatialDisc_.nelemDomain_: c_int);
            phiGlobal[col - 1] = saved - AD_FD_EPS;
            const rMinus = residualRowForAD(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                            this.spatialDisc_.nelemDomain_: c_int);
            phiGlobal[col - 1] = saved;
            fd = (rPlus - rMinus) / (2.0 * AD_FD_EPS);
        }

        maxADFD = max(maxADFD, abs(ad - fd));
        maxAnalyticalAD = max(maxAnalyticalAD, abs(analytical - ad));
        maxExactAD = max(maxExactAD, abs(exact - ad));
        maxMatrixExact = max(maxMatrixExact, abs(analytical - exact));

        if (abs(ad) > AD_ROW_PRINT_TOL || abs(fd) > AD_ROW_PRINT_TOL ||
            abs(analytical) > AD_ROW_PRINT_TOL || abs(exact) > AD_ROW_PRINT_TOL) &&
           printed < AD_ROW_MAX_PRINT {
            writeln("  col ", col,
                    " matrix=", analytical,
                    " exact=", exact,
                    " ad=", ad,
                    " fd=", fd,
                    " |exact-ad|=", abs(exact - ad),
                    " |matrix-exact|=", abs(analytical - exact),
                    " |ad-fd|=", abs(ad - fd),
                    " |matrix-ad|=", abs(analytical - ad));
            printed += 1;
        }
    }

    var outsideCount = 0;
    var outsideMax = 0.0;
    for globalIdx in phiDom {
        const col = globalIdx + 1;
        var inStencil = false;
        for stencilCol in stencil do
            if stencilCol == col then
                inStencil = true;

        if !inStencil && abs(dphiGlobal[globalIdx]) > AD_ROW_PRINT_TOL {
            outsideCount += 1;
            outsideMax = max(outsideMax, abs(dphiGlobal[globalIdx]));
        }
    }

    writeln("AD row check summary: max |AD-FD| = ", maxADFD,
            " max |exact-AD| = ", maxExactAD,
            " max |matrix-exact| = ", maxMatrixExact,
            " max |matrix-AD| = ", maxAnalyticalAD,
            " d/dGamma = ", dgammaGlobal,
            " outside-stencil nonzeros = ", outsideCount,
            " outside-stencil max = ", outsideMax);
}

proc temporalDiscretization.runADKuttaCheck() {
    if this.spatialDisc_.isIBMFlow() {
        this.runIBMADKuttaCheck();
        return;
    }

    writeln("Running AD Kutta Jacobian check");

    const stencil = this.buildKuttaStencil();
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    var phiGlobal: [phiDom] real(64);
    var dphiGlobal: [phiDom] real(64) = 0.0;
    var gammaGlobal = this.spatialDisc_.circulation_;
    var dgammaGlobal = 0.0;

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    const sd = this.spatialDisc_.borrow();
    const baseResidual = kuttaResidualForAD(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                            this.spatialDisc_.nelemDomain_: c_int);
    __enzyme_autodiff(c_ptrTo(kuttaResidualForAD): c_ptr(void),
                      enzyme_const,
                      sd,
                      enzyme_dup,
                      c_ptrTo(phiGlobal[0]),
                      c_ptrTo(dphiGlobal[0]),
                      enzyme_dup,
                      c_ptrTo(gammaGlobal),
                      c_ptrTo(dgammaGlobal),
                      enzyme_const,
                      this.spatialDisc_.nelemDomain_: c_int);
    const exactDgamma = kuttaResidualDerivativeAnalytical(sd, c_ptrTo(phiGlobal[0]),
                                                          c_ptrTo(gammaGlobal), 0, true,
                                                          this.spatialDisc_.nelemDomain_: c_int);

    writeln("AD kutta base residual = ", baseResidual, " stencil size = ", stencil.size);

    var maxADFD = 0.0;
    var maxExactAD = 0.0;
    var printed = 0;
    for col in stencil {
        const exact = kuttaResidualDerivativeAnalytical(sd, c_ptrTo(phiGlobal[0]),
                                                        c_ptrTo(gammaGlobal), col, false,
                                                        this.spatialDisc_.nelemDomain_: c_int);
        const saved = phiGlobal[col - 1];
        phiGlobal[col - 1] = saved + AD_FD_EPS;
        const rPlus = kuttaResidualForAD(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                         this.spatialDisc_.nelemDomain_: c_int);
        phiGlobal[col - 1] = saved - AD_FD_EPS;
        const rMinus = kuttaResidualForAD(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                          this.spatialDisc_.nelemDomain_: c_int);
        phiGlobal[col - 1] = saved;

        const fd = (rPlus - rMinus) / (2.0 * AD_FD_EPS);
        const ad = dphiGlobal[col - 1];
        maxADFD = max(maxADFD, abs(ad - fd));
        maxExactAD = max(maxExactAD, abs(exact - ad));

        if (abs(ad) > AD_ROW_PRINT_TOL || abs(fd) > AD_ROW_PRINT_TOL || abs(exact) > AD_ROW_PRINT_TOL) &&
           printed < AD_ROW_MAX_PRINT {
            writeln("  kutta col ", col,
                    " exact=", exact,
                    " ad=", ad,
                    " fd=", fd,
                    " |exact-ad|=", abs(exact - ad),
                    " |ad-fd|=", abs(ad - fd));
            printed += 1;
        }
    }

    const savedGamma = gammaGlobal;
    gammaGlobal = savedGamma + AD_FD_EPS;
    const rPlus = kuttaResidualForAD(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                     this.spatialDisc_.nelemDomain_: c_int);
    gammaGlobal = savedGamma - AD_FD_EPS;
    const rMinus = kuttaResidualForAD(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                      this.spatialDisc_.nelemDomain_: c_int);
    gammaGlobal = savedGamma;
    const fdGamma = (rPlus - rMinus) / (2.0 * AD_FD_EPS);

    writeln("AD kutta check summary: max |AD-FD| = ", maxADFD,
            " max |exact-AD| = ", maxExactAD,
            " d/dGamma exact = ", exactDgamma,
            " d/dGamma ad = ", dgammaGlobal,
            " d/dGamma fd = ", fdGamma,
            " |exact-ad|_gamma = ", abs(exactDgamma - dgammaGlobal),
            " |ad-fd|_gamma = ", abs(dgammaGlobal - fdGamma));
}
}
