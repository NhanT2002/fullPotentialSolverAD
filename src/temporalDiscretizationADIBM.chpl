module temporalDiscretizationADIBM {
use CTypes;
use Math;
use Set;
use Sort;
use temporalDiscretization;
use spatialDiscretization;

inline proc temporalDiscretization.borrowIBMSpatialDisc(): borrowed spatialDiscretizationIBM {
    return try! (this.spatialDisc_: borrowed spatialDiscretizationIBM);
}

proc temporalDiscretization.pickIBMADRow(): int {
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

inline proc ibmGhostPhiForFace(sd: borrowed spatialDiscretizationIBM,
                               phiPtr: c_ptr(real(64)),
                               face: int,
                               interiorElem: int,
                               ghostElem: int): real(64) {
    if isWallFace(sd, face) {
        const ibmIndex = sd.getIBMEntryForWallFace(face);
        if ibmIndex <= 0 then
            halt("Missing IBM entry for wall face ", face, " while building the exact IBM Jacobian");

        var interpolatedPhi = 0.0;
        for stencilIndex in 1..sd.ibmInterpStencilSize_[ibmIndex] {
            const cellIndex = sd.ibmInterpStencilCellIndex_[ibmIndex, stencilIndex];
            const weight = sd.ibmInterpStencilWeight_[ibmIndex, stencilIndex];
            interpolatedPhi += weight * phiFromPtr(phiPtr, cellIndex);
        }
        return interpolatedPhi;
    }

    return ghostPhiForFace(sd, phiPtr, face, interiorElem, ghostElem);
}

proc computeVelocityForCellFromPtrIBM(sd: borrowed spatialDiscretizationIBM,
                                      phiPtr: c_ptr(real(64)),
                                      gammaPtr: c_ptr(real(64)),
                                      elem: int): (real(64), real(64)) {
    const phiI = phiFromPtr(phiPtr, elem);
    const gamma = gammaFromPtr(gammaPtr);
    var gx = 0.0;
    var gy = 0.0;
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];

    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        const elem1 = sd.mesh_.edge2elem_[1, face];
        const elem2 = sd.mesh_.edge2elem_[2, face];
        const neighbor = if elem1 == elem then elem2 else elem1;

        var phiJ: real(64);
        if neighbor <= sd.nelemDomain_ {
            phiJ = phiFromPtr(phiPtr, neighbor);

            const elemKuttaType = sd.kuttaCell_[elem];
            const neighborKuttaType = sd.kuttaCell_[neighbor];
            if elemKuttaType == 1 && neighborKuttaType == -1 then
                phiJ += gamma;
            else if elemKuttaType == -1 && neighborKuttaType == 1 then
                phiJ -= gamma;
        } else {
            phiJ = ibmGhostPhiForFace(sd, phiPtr, face, elem, neighbor);
        }

        const dphi = phiJ - phiI;
        if elem == elem1 {
            gx += sd.lsGradQR_!.wxFinal1_[face] * dphi;
            gy += sd.lsGradQR_!.wyFinal1_[face] * dphi;
        } else {
            gx += sd.lsGradQR_!.wxFinal2_[face] * dphi;
            gy += sd.lsGradQR_!.wyFinal2_[face] * dphi;
        }
    }

    return (gx, gy);
}

proc computeGhostVelocityForFaceIBM(sd: borrowed spatialDiscretizationIBM,
                                    phiPtr: c_ptr(real(64)),
                                    gammaPtr: c_ptr(real(64)),
                                    face: int,
                                    interiorElem: int): (real(64), real(64)) {
    if isWallFace(sd, face) {
        const ibmIndex = sd.getIBMEntryForWallFace(face);
        if ibmIndex <= 0 then
            halt("Missing IBM entry for wall face ", face, " while building the exact IBM Jacobian");

        const normalX = sd.ibmNormalX_[ibmIndex];
        const normalY = sd.ibmNormalY_[ibmIndex];
        const alpha = sd.ibmWallToGhostSignedDistance_[ibmIndex];
        const beta = sd.ibmWallToImageSignedDistance_[ibmIndex];
        var imageVelocityX = 0.0;
        var imageVelocityY = 0.0;

        for stencilIndex in 1..sd.ibmInterpStencilSize_[ibmIndex] {
            const cellIndex = sd.ibmInterpStencilCellIndex_[ibmIndex, stencilIndex];
            const weight = sd.ibmInterpStencilWeight_[ibmIndex, stencilIndex];
            const (uCell, vCell) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, cellIndex);
            imageVelocityX += weight * uCell;
            imageVelocityY += weight * vCell;
        }

        const normalVelImage = imageVelocityX * normalX + imageVelocityY * normalY;
        var tangentX = imageVelocityX - normalVelImage * normalX;
        var tangentY = imageVelocityY - normalVelImage * normalY;
        const tangentMag = sqrt(tangentX * tangentX + tangentY * tangentY);
        tangentX /= tangentMag;
        tangentY /= tangentMag;

        const tangentVelImage = imageVelocityX * tangentX + imageVelocityY * tangentY;
        const normalVelGhost = alpha / beta * normalVelImage;
        const tangentVelGhost = tangentVelImage;

        return (normalVelGhost * normalX + tangentVelGhost * tangentX,
                normalVelGhost * normalY + tangentVelGhost * tangentY);
    }

    const (uInt, vInt) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, interiorElem);
    const x = sd.faceCentroidX_[face];
    const y = sd.faceCentroidY_[face];
    var uFaceBC: real(64);
    var vFaceBC: real(64);

    if sd.farfieldIsCylinder_ {
        const R = sd.inputs_.CYLINDER_RADIUS_;
        const r2 = x * x + y * y;
        const R2_over_r2 = R * R / r2;
        const theta = atan2(y, x);
        const cosTheta = cos(theta);
        const sinTheta = sin(theta);

        const Vr = sd.inputs_.VEL_INF_ * (1.0 - R2_over_r2) * cosTheta;
        const Vtheta = -sd.inputs_.VEL_INF_ * (1.0 + R2_over_r2) * sinTheta;
        uFaceBC = Vr * cosTheta - Vtheta * sinTheta;
        vFaceBC = Vr * sinTheta + Vtheta * cosTheta;
    } else {
        uFaceBC = sd.inputs_.U_INF_;
        vFaceBC = sd.inputs_.V_INF_;
    }

    return (2.0 * uFaceBC - uInt, 2.0 * vFaceBC - vInt);
}

proc computeRhoMuForElemOnFaceIBM(sd: borrowed spatialDiscretizationIBM,
                                  phiPtr: c_ptr(real(64)),
                                  gammaPtr: c_ptr(real(64)),
                                  face: int,
                                  elem: int): (real(64), real(64)) {
    if elem <= sd.nelemDomain_ {
        const (u, v) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, elem);
        const (rho, mach, mu) = computeDensityAndMuFromVelocity(sd, u, v);
        return (rho, mu);
    }

    const elem1 = sd.mesh_.edge2elem_[1, face];
    const elem2 = sd.mesh_.edge2elem_[2, face];
    const interiorElem = if elem1 <= sd.nelemDomain_ then elem1 else elem2;

    if isWallFace(sd, face) {
        const (u, v) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, interiorElem);
        const (rho, mach, mu) = computeDensityAndMuFromVelocity(sd, u, v);
        return (rho, mu);
    }

    const (uGhost, vGhost) = computeGhostVelocityForFaceIBM(sd, phiPtr, gammaPtr, face, interiorElem);
    const (rhoGhost, machGhost, muGhost) = computeDensityAndMuFromVelocity(sd, uGhost, vGhost);
    return (rhoGhost, muGhost);
}

proc residualRowForADIBM(sd: borrowed spatialDiscretizationIBM,
                         elem: int,
                         phiPtr: c_ptr(real(64)),
                         gammaPtr: c_ptr(real(64)),
                         n: c_int): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    var res = 0.0;
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];

    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        const elem1 = sd.mesh_.edge2elem_[1, face];
        const elem2 = sd.mesh_.edge2elem_[2, face];
        const sign = if elem1 == elem then 1.0 else -1.0;

        var u1, v1, u2, v2: real(64);
        if elem1 <= sd.nelemDomain_ then
            (u1, v1) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, elem1);
        else
            (u1, v1) = computeGhostVelocityForFaceIBM(sd, phiPtr, gammaPtr, face, elem2);

        if elem2 <= sd.nelemDomain_ then
            (u2, v2) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, elem2);
        else
            (u2, v2) = computeGhostVelocityForFaceIBM(sd, phiPtr, gammaPtr, face, elem1);

        const uAvg = sd.weights1_[face] * u1 + sd.weights2_[face] * u2;
        const vAvg = sd.weights1_[face] * v1 + sd.weights2_[face] * v2;

        var phi1: real(64);
        var phi2: real(64);
        if elem1 <= sd.nelemDomain_ then
            phi1 = phiFromPtr(phiPtr, elem1);
        else
            phi1 = ibmGhostPhiForFace(sd, phiPtr, face, elem2, elem1);

        if elem2 <= sd.nelemDomain_ then
            phi2 = phiFromPtr(phiPtr, elem2);
        else
            phi2 = ibmGhostPhiForFace(sd, phiPtr, face, elem1, elem2);

        if elem1 <= sd.nelemDomain_ && elem2 <= sd.nelemDomain_ {
            const kuttaType1 = sd.kuttaCell_[elem1];
            const kuttaType2 = sd.kuttaCell_[elem2];
            if kuttaType1 == 1 && kuttaType2 == -1 then
                phi2 += gamma;
            else if kuttaType1 == -1 && kuttaType2 == 1 then
                phi2 -= gamma;
        }

        const dPhidl = (phi2 - phi1) * sd.invL_IJ_[face];
        const vDotT = uAvg * sd.t_IJ_x_[face] + vAvg * sd.t_IJ_y_[face];
        const delta = vDotT - dPhidl;
        const uFace = uAvg - delta * sd.corrCoeffX_[face];
        const vFace = vAvg - delta * sd.corrCoeffY_[face];
        var rhoFace = (1.0 + sd.gamma_minus_one_over_two_ * sd.inputs_.MACH_ * sd.inputs_.MACH_ *
                      (1.0 - uFace * uFace - vFace * vFace)) ** sd.one_over_gamma_minus_one_;

        const nx = sd.faceNormalX_[face];
        const ny = sd.faceNormalY_[face];
        const vDotN = uFace * nx + vFace * ny;

        const (rho1Blend, mu1Blend) = computeRhoMuForElemOnFaceIBM(sd, phiPtr, gammaPtr, face, elem1);
        const (rho2Blend, mu2Blend) = computeRhoMuForElemOnFaceIBM(sd, phiPtr, gammaPtr, face, elem2);
        const upwindWeight = if vDotN >= 0.0 then 1.0 else 0.0;
        const rhoUpwind = upwindWeight * rho1Blend + (1.0 - upwindWeight) * rho2Blend;
        const muUpwind = upwindWeight * mu1Blend + (1.0 - upwindWeight) * mu2Blend;
        if muUpwind > 0.0 then
            rhoFace = rhoFace - muUpwind * (rhoFace - rhoUpwind);

        res += sign * rhoFace * vDotN * sd.faceArea_[face];
    }

    return res * sd.res_scale_;
}

proc kuttaResidualForADIBM(sd: borrowed spatialDiscretizationIBM,
                           phiPtr: c_ptr(real(64)),
                           gammaPtr: c_ptr(real(64)),
                           n: c_int): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    const upperElem = sd.upperTEelem_;
    const lowerElem = sd.lowerTEelem_;
    const (uUpper, vUpper) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, upperElem);
    const (uLower, vLower) = computeVelocityForCellFromPtrIBM(sd, phiPtr, gammaPtr, lowerElem);

    const phiUpper = phiFromPtr(phiPtr, upperElem) +
                     (uUpper * sd.deltaSupperTEx_ + vUpper * sd.deltaSupperTEy_);
    const phiLower = phiFromPtr(phiPtr, lowerElem) +
                     (uLower * sd.deltaSlowerTEx_ + vLower * sd.deltaSlowerTEy_);
    const gammaComputed = phiUpper - phiLower;

    return (gamma - gammaComputed) * sd.res_scale_;
}

proc temporalDiscretization.buildIBMRowStencilUncached(row: int) {
    const sd = this.borrowIBMSpatialDisc();
    var stencilSet = new set(int);
    var cellsPending = new set(int);
    var cellsScanned = new set(int);

    const bodyStencil = this.buildBodyFittedRowStencilUncached(row);
    for col in bodyStencil {
        stencilSet.add(col);
        cellsPending.add(col);
    }

    while cellsPending.size > 0 {
        var nextIBMCells = new set(int);
        var cellsToProcess: [0..<cellsPending.size] int;
        var pendingIdx = 0;
        for elem in cellsPending {
            cellsToProcess[pendingIdx] = elem;
            pendingIdx += 1;
        }
        cellsPending.clear();

        for elem in cellsToProcess {
            if cellsScanned.contains(elem) then
                continue;
            cellsScanned.add(elem);

            const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
            const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];
            for faceIdx in faceStart..faceEnd {
                const face = sd.mesh_.elem2edge_[faceIdx];
                if !sd.wallFaceSet_.contains(face) then
                    continue;

                const ibmIndex = sd.getIBMEntryForWallFace(face);
                if ibmIndex <= 0 then
                    continue;

                for stencilIndex in 1..sd.ibmInterpStencilSize_[ibmIndex] {
                    const cellIndex = sd.ibmInterpStencilCellIndex_[ibmIndex, stencilIndex];
                    if cellIndex <= sd.nelemDomain_ && !stencilSet.contains(cellIndex) {
                        stencilSet.add(cellIndex);
                        nextIBMCells.add(cellIndex);
                    }
                }
            }
        }

        for cellIndex in nextIBMCells {
            const cellStencil = this.buildBodyFittedRowStencilUncached(cellIndex);
            for col in cellStencil {
                stencilSet.add(col);
                if !cellsScanned.contains(col) then
                    cellsPending.add(col);
            }
        }
    }

    const dom = {0..<stencilSet.size};
    var stencil: [dom] int;
    var idx = 0;
    for elem in stencilSet {
        stencil[idx] = elem;
        idx += 1;
    }
    sort(stencil);
    return stencil;
}

proc temporalDiscretization.enforceIBMConsistentGammaForCurrentPhi(maxIts: int = 12,
                                                                   tol: real(64) = 1.0e-14) {
    const sd = this.borrowIBMSpatialDisc();
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    const n = this.spatialDisc_.nelemDomain_: c_int;
    var phiGlobal: [phiDom] real(64);

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    var gamma = this.spatialDisc_.circulation_;
    for gammaIt in 1..maxIts {
        const residual = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gamma), n);
        if abs(residual) < tol then break;

        var dgamma = 0.0;
        __enzyme_autodiff(c_ptrTo(kuttaResidualForADIBM): c_ptr(void),
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

proc temporalDiscretization.computeIBMADReducedExactJacobian() {
    this.A_petsc.zeroEntries();

    const sd = this.borrowIBMSpatialDisc();
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    const n = this.spatialDisc_.nelemDomain_: c_int;
    var phiGlobal: [phiDom] real(64);
    var gammaGlobal = this.spatialDisc_.circulation_;
    var kuttaDphi: [phiDom] real(64) = 0.0;
    var kuttaDgamma = 0.0;
    const kuttaStencil = this.buildKuttaStencil();

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    __enzyme_autodiff(c_ptrTo(kuttaResidualForADIBM): c_ptr(void),
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
        halt("IBM reduced exact Jacobian failed: dR_gamma/dGamma is too small");

    for row in 1..this.spatialDisc_.nelemDomain_ {
        if AD_JACOBIAN_PROGRESS_FREQ > 0 &&
           (row == 1 || row % AD_JACOBIAN_PROGRESS_FREQ == 0 || row == this.spatialDisc_.nelemDomain_) then
            writeln("IBM AD reduced exact Jacobian assembly row ", row, " / ", this.spatialDisc_.nelemDomain_);

        var rowDphi: [phiDom] real(64) = 0.0;
        var rowDgamma = 0.0;
        __enzyme_autodiff(c_ptrTo(residualRowForADIBM): c_ptr(void),
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

proc temporalDiscretization.runIBMADRowCheck() {
    if AD_ROW < 0 {
        this.runIBMADKuttaCheck();
        return;
    }

    const row = this.pickIBMADRow();
    writeln("Running IBM AD row Jacobian check for row ", row);

    const stencil = this.buildRowStencil(row);
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    var phiGlobal: [phiDom] real(64);
    var dphiGlobal: [phiDom] real(64) = 0.0;
    var dgammaGlobal = 0.0;
    var gammaGlobal = this.spatialDisc_.circulation_;

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    const sd = this.borrowIBMSpatialDisc();
    const baseResidual = residualRowForADIBM(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                             this.spatialDisc_.nelemDomain_: c_int);
    __enzyme_autodiff(c_ptrTo(residualRowForADIBM): c_ptr(void),
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

    var kuttaDphi: [phiDom] real(64) = 0.0;
    var kuttaDgamma = 0.0;
    __enzyme_autodiff(c_ptrTo(kuttaResidualForADIBM): c_ptr(void),
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

    writeln("IBM AD row base residual = ", baseResidual, " stencil size = ", stencil.size);

    var maxADFD = 0.0;
    var maxMatrixAD = 0.0;
    var printed = 0;

    for col in stencil {
        const ad = dphiGlobal[col - 1] - (dgammaGlobal / kuttaDgamma) * kuttaDphi[col - 1];
        const matrix = this.A_petsc.get(row - 1, col - 1);

        const saved = phiGlobal[col - 1];
        phiGlobal[col - 1] = saved + AD_FD_EPS;
        var gammaPlus = gammaGlobal;
        for gammaIt in 1..12 {
            const rGamma = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaPlus),
                                                 this.spatialDisc_.nelemDomain_: c_int);
            if abs(rGamma) < 1.0e-14 then break;
            var dGamma = 0.0;
            __enzyme_autodiff(c_ptrTo(kuttaResidualForADIBM): c_ptr(void),
                              enzyme_const,
                              sd,
                              enzyme_const,
                              c_ptrTo(phiGlobal[0]),
                              enzyme_dup,
                              c_ptrTo(gammaPlus),
                              c_ptrTo(dGamma),
                              enzyme_const,
                              this.spatialDisc_.nelemDomain_: c_int);
            if abs(dGamma) < 1.0e-14 then break;
            gammaPlus -= rGamma / dGamma;
        }
        const rPlus = residualRowForADIBM(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaPlus),
                                          this.spatialDisc_.nelemDomain_: c_int);

        phiGlobal[col - 1] = saved - AD_FD_EPS;
        var gammaMinus = gammaGlobal;
        for gammaIt in 1..12 {
            const rGamma = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaMinus),
                                                 this.spatialDisc_.nelemDomain_: c_int);
            if abs(rGamma) < 1.0e-14 then break;
            var dGamma = 0.0;
            __enzyme_autodiff(c_ptrTo(kuttaResidualForADIBM): c_ptr(void),
                              enzyme_const,
                              sd,
                              enzyme_const,
                              c_ptrTo(phiGlobal[0]),
                              enzyme_dup,
                              c_ptrTo(gammaMinus),
                              c_ptrTo(dGamma),
                              enzyme_const,
                              this.spatialDisc_.nelemDomain_: c_int);
            if abs(dGamma) < 1.0e-14 then break;
            gammaMinus -= rGamma / dGamma;
        }
        const rMinus = residualRowForADIBM(sd, row, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaMinus),
                                           this.spatialDisc_.nelemDomain_: c_int);
        phiGlobal[col - 1] = saved;
        const fd = (rPlus - rMinus) / (2.0 * AD_FD_EPS);

        maxADFD = max(maxADFD, abs(ad - fd));
        maxMatrixAD = max(maxMatrixAD, abs(matrix - ad));

        if (abs(ad) > AD_ROW_PRINT_TOL || abs(fd) > AD_ROW_PRINT_TOL ||
            abs(matrix) > AD_ROW_PRINT_TOL) &&
           printed < AD_ROW_MAX_PRINT {
            writeln("  col ", col,
                    " matrix=", matrix,
                    " ad=", ad,
                    " fd=", fd,
                    " |ad-fd|=", abs(ad - fd),
                    " |matrix-ad|=", abs(matrix - ad));
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

    writeln("IBM AD row check summary: max |AD-FD| = ", maxADFD,
            " max |matrix-AD| = ", maxMatrixAD,
            " d/dGamma = ", dgammaGlobal,
            " outside-stencil nonzeros = ", outsideCount,
            " outside-stencil max = ", outsideMax);
}

proc temporalDiscretization.runIBMADKuttaCheck() {
    writeln("Running IBM AD Kutta Jacobian check");

    const stencil = this.buildKuttaStencil();
    const phiDom = {0..<this.spatialDisc_.nelemDomain_};
    var phiGlobal: [phiDom] real(64);
    var dphiGlobal: [phiDom] real(64) = 0.0;
    var gammaGlobal = this.spatialDisc_.circulation_;
    var dgammaGlobal = 0.0;

    forall elem in 1..this.spatialDisc_.nelemDomain_ do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    const sd = this.borrowIBMSpatialDisc();
    const baseResidual = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                               this.spatialDisc_.nelemDomain_: c_int);
    __enzyme_autodiff(c_ptrTo(kuttaResidualForADIBM): c_ptr(void),
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

    writeln("IBM AD kutta base residual = ", baseResidual, " stencil size = ", stencil.size);

    var maxADFD = 0.0;
    var printed = 0;
    for col in stencil {
        const saved = phiGlobal[col - 1];
        phiGlobal[col - 1] = saved + AD_FD_EPS;
        const rPlus = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                            this.spatialDisc_.nelemDomain_: c_int);
        phiGlobal[col - 1] = saved - AD_FD_EPS;
        const rMinus = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                             this.spatialDisc_.nelemDomain_: c_int);
        phiGlobal[col - 1] = saved;

        const fd = (rPlus - rMinus) / (2.0 * AD_FD_EPS);
        const ad = dphiGlobal[col - 1];
        maxADFD = max(maxADFD, abs(ad - fd));

        if (abs(ad) > AD_ROW_PRINT_TOL || abs(fd) > AD_ROW_PRINT_TOL) &&
           printed < AD_ROW_MAX_PRINT {
            writeln("  kutta col ", col,
                    " ad=", ad,
                    " fd=", fd,
                    " |ad-fd|=", abs(ad - fd));
            printed += 1;
        }
    }

    const savedGamma = gammaGlobal;
    gammaGlobal = savedGamma + AD_FD_EPS;
    const rPlus = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                        this.spatialDisc_.nelemDomain_: c_int);
    gammaGlobal = savedGamma - AD_FD_EPS;
    const rMinus = kuttaResidualForADIBM(sd, c_ptrTo(phiGlobal[0]), c_ptrTo(gammaGlobal),
                                         this.spatialDisc_.nelemDomain_: c_int);
    gammaGlobal = savedGamma;
    const fdGamma = (rPlus - rMinus) / (2.0 * AD_FD_EPS);

    writeln("IBM AD kutta check summary: max |AD-FD| = ", maxADFD,
            " d/dGamma ad = ", dgammaGlobal,
            " d/dGamma fd = ", fdGamma,
            " |ad-fd|_gamma = ", abs(dgammaGlobal - fdGamma));
}
}
