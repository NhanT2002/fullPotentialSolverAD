module temporalDiscretizationAnalyticalIBM {
use CTypes;
use Math;
use temporalDiscretization;
use temporalDiscretizationADIBM;
use spatialDiscretization;

inline proc ibmGhostPhiForFaceWithDerivativeIBM(sd: borrowed spatialDiscretizationIBM,
                                                phiPtr: c_ptr(real(64)),
                                                face: int,
                                                interiorElem: int,
                                                ghostElem: int,
                                                seedCol: int,
                                                wrtGamma: bool): (real(64), real(64)) {
    if isWallFace(sd, face) {
        const ibmIndex = sd.getIBMEntryForWallFace(face);
        if ibmIndex <= 0 then
            halt("Missing IBM entry for wall face ", face,
                 " while building the analytical IBM exact Jacobian");

        var ghostPhi = 0.0;
        var dGhostPhi = 0.0;
        for stencilIndex in 1..sd.ibmInterpStencilSize_[ibmIndex] {
            const cellIndex = sd.ibmInterpStencilCellIndex_[ibmIndex, stencilIndex];
            const weight = sd.ibmInterpStencilWeight_[ibmIndex, stencilIndex];
            ghostPhi += weight * phiFromPtr(phiPtr, cellIndex);
            dGhostPhi += weight * phiSeedDerivative(cellIndex, seedCol, wrtGamma);
        }
        return (ghostPhi, dGhostPhi);
    }

    return ghostPhiForFaceWithDerivative(sd, phiPtr, face, interiorElem, ghostElem, seedCol, wrtGamma);
}

proc computeVelocityForCellFromPtrWithDerivativeIBM(sd: borrowed spatialDiscretizationIBM,
                                                    phiPtr: c_ptr(real(64)),
                                                    gammaPtr: c_ptr(real(64)),
                                                    elem: int,
                                                    seedCol: int,
                                                    wrtGamma: bool): (real(64), real(64), real(64), real(64)) {
    const phiI = phiFromPtr(phiPtr, elem);
    const dphiI = phiSeedDerivative(elem, seedCol, wrtGamma);
    const gamma = gammaFromPtr(gammaPtr);
    const dgamma = gammaSeedDerivative(wrtGamma);

    var gx = 0.0;
    var gy = 0.0;
    var dgx = 0.0;
    var dgy = 0.0;
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];

    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        const elem1 = sd.mesh_.edge2elem_[1, face];
        const elem2 = sd.mesh_.edge2elem_[2, face];
        const neighbor = if elem1 == elem then elem2 else elem1;

        var phiJ: real(64);
        var dphiJ: real(64);
        if neighbor <= sd.nelemDomain_ {
            phiJ = phiFromPtr(phiPtr, neighbor);
            dphiJ = phiSeedDerivative(neighbor, seedCol, wrtGamma);

            const elemKuttaType = sd.kuttaCell_[elem];
            const neighborKuttaType = sd.kuttaCell_[neighbor];
            if elemKuttaType == 1 && neighborKuttaType == -1 {
                phiJ += gamma;
                dphiJ += dgamma;
            } else if elemKuttaType == -1 && neighborKuttaType == 1 {
                phiJ -= gamma;
                dphiJ -= dgamma;
            }
        } else {
            (phiJ, dphiJ) = ibmGhostPhiForFaceWithDerivativeIBM(sd, phiPtr, face, elem, neighbor,
                                                                seedCol, wrtGamma);
        }

        const dphi = phiJ - phiI;
        const ddphi = dphiJ - dphiI;
        if elem == elem1 {
            gx += sd.lsGradQR_!.wxFinal1_[face] * dphi;
            gy += sd.lsGradQR_!.wyFinal1_[face] * dphi;
            dgx += sd.lsGradQR_!.wxFinal1_[face] * ddphi;
            dgy += sd.lsGradQR_!.wyFinal1_[face] * ddphi;
        } else {
            gx += sd.lsGradQR_!.wxFinal2_[face] * dphi;
            gy += sd.lsGradQR_!.wyFinal2_[face] * dphi;
            dgx += sd.lsGradQR_!.wxFinal2_[face] * ddphi;
            dgy += sd.lsGradQR_!.wyFinal2_[face] * ddphi;
        }
    }

    return (gx, gy, dgx, dgy);
}

proc computeGhostVelocityForFaceWithDerivativeIBM(sd: borrowed spatialDiscretizationIBM,
                                                  phiPtr: c_ptr(real(64)),
                                                  gammaPtr: c_ptr(real(64)),
                                                  face: int,
                                                  interiorElem: int,
                                                  seedCol: int,
                                                  wrtGamma: bool): (real(64), real(64), real(64), real(64)) {
    if isWallFace(sd, face) {
        const ibmIndex = sd.getIBMEntryForWallFace(face);
        if ibmIndex <= 0 then
            halt("Missing IBM entry for wall face ", face,
                 " while building the analytical IBM exact Jacobian");

        const beta = sd.ibmWallToImageSignedDistance_[ibmIndex];
        if abs(beta) < 1.0e-14 then
            halt("IBM analytical reduced exact Jacobian encountered near-zero image distance ",
                 "for wall face ", face, " and IBM entry ", ibmIndex);

        const normalX = sd.ibmNormalX_[ibmIndex];
        const normalY = sd.ibmNormalY_[ibmIndex];
        const normalScale = sd.ibmWallToGhostSignedDistance_[ibmIndex] / beta - 1.0;
        var imageVelocityX = 0.0;
        var imageVelocityY = 0.0;
        var dImageVelocityX = 0.0;
        var dImageVelocityY = 0.0;

        for stencilIndex in 1..sd.ibmInterpStencilSize_[ibmIndex] {
            const cellIndex = sd.ibmInterpStencilCellIndex_[ibmIndex, stencilIndex];
            const weight = sd.ibmInterpStencilWeight_[ibmIndex, stencilIndex];
            const (uCell, vCell, duCell, dvCell) =
                computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, cellIndex,
                                                               seedCol, wrtGamma);
            imageVelocityX += weight * uCell;
            imageVelocityY += weight * vCell;
            dImageVelocityX += weight * duCell;
            dImageVelocityY += weight * dvCell;
        }

        const normalVelImage = imageVelocityX * normalX + imageVelocityY * normalY;
        const dNormalVelImage = dImageVelocityX * normalX + dImageVelocityY * normalY;

        return (imageVelocityX + normalScale * normalVelImage * normalX,
                imageVelocityY + normalScale * normalVelImage * normalY,
                dImageVelocityX + normalScale * dNormalVelImage * normalX,
                dImageVelocityY + normalScale * dNormalVelImage * normalY);
    }

    const (uInt, vInt, duInt, dvInt) =
        computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, interiorElem,
                                                       seedCol, wrtGamma);
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

    return (2.0 * uFaceBC - uInt,
            2.0 * vFaceBC - vInt,
            -duInt,
            -dvInt);
}

proc computeRhoMuForElemOnFaceWithDerivativeIBM(sd: borrowed spatialDiscretizationIBM,
                                                phiPtr: c_ptr(real(64)),
                                                gammaPtr: c_ptr(real(64)),
                                                face: int,
                                                elem: int,
                                                seedCol: int,
                                                wrtGamma: bool): (real(64), real(64), real(64), real(64)) {
    if elem <= sd.nelemDomain_ {
        const (u, v, du, dv) =
            computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, elem, seedCol, wrtGamma);
        const (rho, mu, drho, dmu) = computeDensityAndMuFromVelocityWithDerivative(sd, u, v, du, dv);
        return (rho, mu, drho, dmu);
    }

    const elem1 = sd.mesh_.edge2elem_[1, face];
    const elem2 = sd.mesh_.edge2elem_[2, face];
    const interiorElem = if elem1 <= sd.nelemDomain_ then elem1 else elem2;

    if isWallFace(sd, face) {
        const (u, v, du, dv) =
            computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, interiorElem,
                                                           seedCol, wrtGamma);
        const (rho, mu, drho, dmu) = computeDensityAndMuFromVelocityWithDerivative(sd, u, v, du, dv);
        return (rho, mu, drho, dmu);
    }

    const (uGhost, vGhost, duGhost, dvGhost) =
        computeGhostVelocityForFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, interiorElem,
                                                     seedCol, wrtGamma);
    const (rhoGhost, muGhost, drhoGhost, dmuGhost) =
        computeDensityAndMuFromVelocityWithDerivative(sd, uGhost, vGhost, duGhost, dvGhost);
    return (rhoGhost, muGhost, drhoGhost, dmuGhost);
}

proc exactResidualFaceDerivativeAnalyticalIBM(sd: borrowed spatialDiscretizationIBM,
                                              elem: int,
                                              face: int,
                                              phiPtr: c_ptr(real(64)),
                                              gammaPtr: c_ptr(real(64)),
                                              seedCol: int,
                                              wrtGamma: bool): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    const dgamma = gammaSeedDerivative(wrtGamma);
    const elem1 = sd.mesh_.edge2elem_[1, face];
    const elem2 = sd.mesh_.edge2elem_[2, face];
    const sign = if elem1 == elem then 1.0 else -1.0;

    var u1, v1, du1, dv1, u2, v2, du2, dv2: real(64);
    if elem1 <= sd.nelemDomain_ {
        (u1, v1, du1, dv1) =
            computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, elem1, seedCol, wrtGamma);
    } else {
        (u1, v1, du1, dv1) =
            computeGhostVelocityForFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, elem2,
                                                         seedCol, wrtGamma);
    }

    if elem2 <= sd.nelemDomain_ {
        (u2, v2, du2, dv2) =
            computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, elem2, seedCol, wrtGamma);
    } else {
        (u2, v2, du2, dv2) =
            computeGhostVelocityForFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, elem1,
                                                         seedCol, wrtGamma);
    }

    const uAvg = sd.weights1_[face] * u1 + sd.weights2_[face] * u2;
    const vAvg = sd.weights1_[face] * v1 + sd.weights2_[face] * v2;
    const duAvg = sd.weights1_[face] * du1 + sd.weights2_[face] * du2;
    const dvAvg = sd.weights1_[face] * dv1 + sd.weights2_[face] * dv2;

    var phi1, dphi1, phi2, dphi2: real(64);
    if elem1 <= sd.nelemDomain_ {
        phi1 = phiFromPtr(phiPtr, elem1);
        dphi1 = phiSeedDerivative(elem1, seedCol, wrtGamma);
    } else {
        (phi1, dphi1) = ibmGhostPhiForFaceWithDerivativeIBM(sd, phiPtr, face, elem2, elem1,
                                                            seedCol, wrtGamma);
    }

    if elem2 <= sd.nelemDomain_ {
        phi2 = phiFromPtr(phiPtr, elem2);
        dphi2 = phiSeedDerivative(elem2, seedCol, wrtGamma);
    } else {
        (phi2, dphi2) = ibmGhostPhiForFaceWithDerivativeIBM(sd, phiPtr, face, elem1, elem2,
                                                            seedCol, wrtGamma);
    }

    if elem1 <= sd.nelemDomain_ && elem2 <= sd.nelemDomain_ {
        const kuttaType1 = sd.kuttaCell_[elem1];
        const kuttaType2 = sd.kuttaCell_[elem2];
        if kuttaType1 == 1 && kuttaType2 == -1 {
            phi2 += gamma;
            dphi2 += dgamma;
        } else if kuttaType1 == -1 && kuttaType2 == 1 {
            phi2 -= gamma;
            dphi2 -= dgamma;
        }
    }

    const dPhidl = (phi2 - phi1) * sd.invL_IJ_[face];
    const ddPhidl = (dphi2 - dphi1) * sd.invL_IJ_[face];
    const vDotT = uAvg * sd.t_IJ_x_[face] + vAvg * sd.t_IJ_y_[face];
    const dvDotT = duAvg * sd.t_IJ_x_[face] + dvAvg * sd.t_IJ_y_[face];
    const delta = vDotT - dPhidl;
    const ddelta = dvDotT - ddPhidl;
    const uFace = uAvg - delta * sd.corrCoeffX_[face];
    const vFace = vAvg - delta * sd.corrCoeffY_[face];
    const duFace = duAvg - ddelta * sd.corrCoeffX_[face];
    const dvFace = dvAvg - ddelta * sd.corrCoeffY_[face];
    const (rhoIsen, muFaceUnused, drhoIsen, dmuFaceUnused) =
        computeDensityAndMuFromVelocityWithDerivative(sd, uFace, vFace, duFace, dvFace);

    const nx = sd.faceNormalX_[face];
    const ny = sd.faceNormalY_[face];
    const vDotN = uFace * nx + vFace * ny;
    const dvDotN = duFace * nx + dvFace * ny;

    const (rho1Blend, mu1Blend, drho1Blend, dmu1Blend) =
        computeRhoMuForElemOnFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, elem1, seedCol, wrtGamma);
    const (rho2Blend, mu2Blend, drho2Blend, dmu2Blend) =
        computeRhoMuForElemOnFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, elem2, seedCol, wrtGamma);
    const upwindWeight = if vDotN >= 0.0 then 1.0 else 0.0;
    const rhoUpwind = upwindWeight * rho1Blend + (1.0 - upwindWeight) * rho2Blend;
    const muUpwind = upwindWeight * mu1Blend + (1.0 - upwindWeight) * mu2Blend;
    const drhoUpwind = upwindWeight * drho1Blend + (1.0 - upwindWeight) * drho2Blend;
    const dmuUpwind = upwindWeight * dmu1Blend + (1.0 - upwindWeight) * dmu2Blend;

    var rhoFace = rhoIsen;
    var drhoFace = drhoIsen;
    if muUpwind > 0.0 {
        rhoFace = rhoIsen - muUpwind * (rhoIsen - rhoUpwind);
        drhoFace = (1.0 - muUpwind) * drhoIsen + muUpwind * drhoUpwind +
                   (rhoUpwind - rhoIsen) * dmuUpwind;
    }

    return sign * sd.faceArea_[face] * (drhoFace * vDotN + rhoFace * dvDotN);
}

proc approximateResidualFaceDerivativeAnalyticalIBM(sd: borrowed spatialDiscretizationIBM,
                                                    elem: int,
                                                    face: int,
                                                    phiPtr: c_ptr(real(64)),
                                                    gammaPtr: c_ptr(real(64)),
                                                    seedCol: int,
                                                    wrtGamma: bool): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    const dgamma = gammaSeedDerivative(wrtGamma);
    const elem1 = sd.mesh_.edge2elem_[1, face];
    const elem2 = sd.mesh_.edge2elem_[2, face];
    const sign = if elem1 == elem then 1.0 else -1.0;

    var u1, v1, du1, dv1, u2, v2, du2, dv2: real(64);
    if elem1 <= sd.nelemDomain_ {
        (u1, v1, du1, dv1) =
            computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, elem1, seedCol, wrtGamma);
    } else {
        (u1, v1, du1, dv1) =
            computeGhostVelocityForFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, elem2,
                                                         seedCol, wrtGamma);
    }

    if elem2 <= sd.nelemDomain_ {
        (u2, v2, du2, dv2) =
            computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, elem2, seedCol, wrtGamma);
    } else {
        (u2, v2, du2, dv2) =
            computeGhostVelocityForFaceWithDerivativeIBM(sd, phiPtr, gammaPtr, face, elem1,
                                                         seedCol, wrtGamma);
    }

    const uAvg = sd.weights1_[face] * u1 + sd.weights2_[face] * u2;
    const vAvg = sd.weights1_[face] * v1 + sd.weights2_[face] * v2;
    const duAvg = sd.weights1_[face] * du1 + sd.weights2_[face] * du2;
    const dvAvg = sd.weights1_[face] * dv1 + sd.weights2_[face] * dv2;

    var phi1, dphi1, phi2, dphi2: real(64);
    if elem1 <= sd.nelemDomain_ {
        phi1 = phiFromPtr(phiPtr, elem1);
        dphi1 = phiSeedDerivative(elem1, seedCol, wrtGamma);
    } else {
        (phi1, dphi1) = ibmGhostPhiForFaceWithDerivativeIBM(sd, phiPtr, face, elem2, elem1,
                                                            seedCol, wrtGamma);
    }

    if elem2 <= sd.nelemDomain_ {
        phi2 = phiFromPtr(phiPtr, elem2);
        dphi2 = phiSeedDerivative(elem2, seedCol, wrtGamma);
    } else {
        (phi2, dphi2) = ibmGhostPhiForFaceWithDerivativeIBM(sd, phiPtr, face, elem1, elem2,
                                                            seedCol, wrtGamma);
    }

    if elem1 <= sd.nelemDomain_ && elem2 <= sd.nelemDomain_ {
        const kuttaType1 = sd.kuttaCell_[elem1];
        const kuttaType2 = sd.kuttaCell_[elem2];
        if kuttaType1 == 1 && kuttaType2 == -1 {
            phi2 += gamma;
            dphi2 += dgamma;
        } else if kuttaType1 == -1 && kuttaType2 == 1 {
            phi2 -= gamma;
            dphi2 -= dgamma;
        }
    }

    const dPhidl = (phi2 - phi1) * sd.invL_IJ_[face];
    const ddPhidl = (dphi2 - dphi1) * sd.invL_IJ_[face];
    const vDotT = uAvg * sd.t_IJ_x_[face] + vAvg * sd.t_IJ_y_[face];
    const dvDotT = duAvg * sd.t_IJ_x_[face] + dvAvg * sd.t_IJ_y_[face];
    const delta = vDotT - dPhidl;
    const ddelta = dvDotT - ddPhidl;
    const duFace = duAvg - ddelta * sd.corrCoeffX_[face];
    const dvFace = dvAvg - ddelta * sd.corrCoeffY_[face];
    const dvDotN = duFace * sd.faceNormalX_[face] + dvFace * sd.faceNormalY_[face];
    const rhoFace = sd.rhoFace_[face];

    return sign * sd.faceArea_[face] * rhoFace * dvDotN;
}

proc residualRowDerivativeAnalyticalIBM(sd: borrowed spatialDiscretizationIBM,
                                        elem: int,
                                        phiPtr: c_ptr(real(64)),
                                        gammaPtr: c_ptr(real(64)),
                                        seedCol: int,
                                        wrtGamma: bool,
                                        n: c_int): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    const dgamma = gammaSeedDerivative(wrtGamma);
    var dres = 0.0;
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];

    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        dres += exactResidualFaceDerivativeAnalyticalIBM(sd, elem, face, phiPtr, gammaPtr,
                                                         seedCol, wrtGamma);
    }

    return dres * sd.res_scale_;
}

proc kuttaResidualDerivativeAnalyticalIBM(sd: borrowed spatialDiscretizationIBM,
                                          phiPtr: c_ptr(real(64)),
                                          gammaPtr: c_ptr(real(64)),
                                          seedCol: int,
                                          wrtGamma: bool,
                                          n: c_int): real(64) {
    const dgamma = gammaSeedDerivative(wrtGamma);
    const upperElem = sd.upperTEelem_;
    const lowerElem = sd.lowerTEelem_;
    const (uUpper, vUpper, duUpper, dvUpper) =
        computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, upperElem, seedCol, wrtGamma);
    const (uLower, vLower, duLower, dvLower) =
        computeVelocityForCellFromPtrWithDerivativeIBM(sd, phiPtr, gammaPtr, lowerElem, seedCol, wrtGamma);

    const dphiUpper = phiSeedDerivative(upperElem, seedCol, wrtGamma) +
                      (duUpper * sd.deltaSupperTEx_ + dvUpper * sd.deltaSupperTEy_);
    const dphiLower = phiSeedDerivative(lowerElem, seedCol, wrtGamma) +
                      (duLower * sd.deltaSlowerTEx_ + dvLower * sd.deltaSlowerTEy_);

    return (dgamma - (dphiUpper - dphiLower)) * sd.res_scale_;
}

proc approximateResidualRowDerivativeAnalyticalIBM(sd: borrowed spatialDiscretizationIBM,
                                                   elem: int,
                                                   phiPtr: c_ptr(real(64)),
                                                   gammaPtr: c_ptr(real(64)),
                                                   seedCol: int,
                                                   wrtGamma: bool,
                                                   n: c_int): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    const dgamma = gammaSeedDerivative(wrtGamma);
    var dres = 0.0;
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];
    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        dres += approximateResidualFaceDerivativeAnalyticalIBM(sd, elem, face, phiPtr, gammaPtr,
                                                               seedCol, wrtGamma);
    }

    return dres * sd.res_scale_;
}

inline proc isIBMWallAdjacentRow(sd: borrowed spatialDiscretizationIBM, elem: int): bool {
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];
    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        if isWallFace(sd, face) then
            return true;
    }
    return false;
}

proc temporalDiscretization.computeIBMApproximateAnalyticalJacobian() {
    this.A_petsc.zeroEntries();

    const sd = this.borrowIBMSpatialDisc();
    const nelem = this.spatialDisc_.nelemDomain_;
    const phiDom = {0..<nelem};
    const n = nelem: c_int;
    var phiGlobal: [phiDom] real(64);
    var gammaGlobal = this.spatialDisc_.circulation_;

    forall elem in 1..nelem do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    for row in 1..nelem {
        const rowStencil = this.buildRowStencil(row);
        const useExactRow = isIBMWallAdjacentRow(sd, row);
        var diag = 0.0;

        for col in rowStencil {
            const value =
                if useExactRow then
                    residualRowDerivativeAnalyticalIBM(sd, row, c_ptrTo(phiGlobal[0]),
                                                       c_ptrTo(gammaGlobal), col, false, n)
                else
                    approximateResidualRowDerivativeAnalyticalIBM(sd, row, c_ptrTo(phiGlobal[0]),
                                                                  c_ptrTo(gammaGlobal), col, false, n);
            if col == row then
                diag = value / this.spatialDisc_.res_scale_;

            if abs(value) > AD_ROW_PRINT_TOL || col == row then
                this.addJacobianEntry(row - 1, col - 1, value,
                                      "computeIBMApproximateAnalyticalJacobian row entry");
        }

        const dRes_dGamma =
            if useExactRow then
                residualRowDerivativeAnalyticalIBM(sd, row, c_ptrTo(phiGlobal[0]),
                                                   c_ptrTo(gammaGlobal), 0, true, n)
            else
                approximateResidualRowDerivativeAnalyticalIBM(sd, row, c_ptrTo(phiGlobal[0]),
                                                              c_ptrTo(gammaGlobal), 0, true, n);
        if dRes_dGamma != 0.0 {
            var upperTEinfluence = 0;
            var lowerTEinfluence = 0;
            if this.getReducedGammaSupportColumns(row, upperTEinfluence, lowerTEinfluence) {
                this.addJacobianEntry(row - 1, upperTEinfluence - 1, dRes_dGamma,
                                      "computeIBMApproximateAnalyticalJacobian reduced-gamma upper");
                this.addJacobianEntry(row - 1, lowerTEinfluence - 1, -dRes_dGamma,
                                      "computeIBMApproximateAnalyticalJacobian reduced-gamma lower");
            }
        }

        this.Jij_[row] = diag;
    }

    forall elem in 1..this.spatialDisc_.nelemDomain_ {
        const faces = this.spatialDisc_.mesh_.elem2edge_[
            this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1 ..
            this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1]
        ];

        for face in faces {
            const machFace = this.spatialDisc_.machFace_[face];
            if machFace >= this.spatialDisc_.inputs_.MACH_C_ {
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
                if elem1 <= this.spatialDisc_.nelemDomain_ && elem2 <= this.spatialDisc_.nelemDomain_ {
                    const upwindElem = this.spatialDisc_.upwindElem_[face];
                    const downwindElem = this.spatialDisc_.downwindElem_[face];
                    if downwindElem == elem {
                        const increase = this.inputs_.BETA_ * this.spatialDisc_.velMagFace_[face] *
                                         this.spatialDisc_.invL_IJ_[face] * this.spatialDisc_.res_scale_;
                        const diagTerm = this.Jij_[elem];
                        if diagTerm >= 0.0 {
                            this.addJacobianEntry(elem - 1, elem - 1, increase,
                                                  "computeIBMApproximateAnalyticalJacobian beta diag+");
                            this.addJacobianEntry(elem - 1, upwindElem - 1, -increase,
                                                  "computeIBMApproximateAnalyticalJacobian beta offdiag-");
                        } else {
                            this.addJacobianEntry(elem - 1, elem - 1, -increase,
                                                  "computeIBMApproximateAnalyticalJacobian beta diag-");
                            this.addJacobianEntry(elem - 1, upwindElem - 1, increase,
                                                  "computeIBMApproximateAnalyticalJacobian beta offdiag+");
                        }
                    }
                }
            }
        }
    }

    this.A_petsc.assemblyComplete();
}

proc temporalDiscretization.computeIBMAnalyticalReducedExactJacobian() {
    this.A_petsc.zeroEntries();

    const sd = this.borrowIBMSpatialDisc();
    const nelem = this.spatialDisc_.nelemDomain_;
    const phiDom = {0..<nelem};
    const n = nelem: c_int;
    var phiGlobal: [phiDom] real(64);
    var gammaGlobal = this.spatialDisc_.circulation_;
    var kuttaDphi: [phiDom] real(64) = 0.0;
    const kuttaStencil = this.buildKuttaStencil();

    forall elem in 1..nelem do
        phiGlobal[elem - 1] = this.spatialDisc_.phi_[elem];

    const kuttaDgamma = kuttaResidualDerivativeAnalyticalIBM(sd, c_ptrTo(phiGlobal[0]),
                                                             c_ptrTo(gammaGlobal), 0, true, n);
    if abs(kuttaDgamma) < 1.0e-14 then
        halt("IBM analytical reduced exact Jacobian failed: dR_gamma/dGamma is too small");

    for col in kuttaStencil do
        kuttaDphi[col - 1] = kuttaResidualDerivativeAnalyticalIBM(sd, c_ptrTo(phiGlobal[0]),
                                                                  c_ptrTo(gammaGlobal), col, false, n);

    for row in 1..nelem {
        if AD_JACOBIAN_PROGRESS_FREQ > 0 &&
           (row == 1 || row % AD_JACOBIAN_PROGRESS_FREQ == 0 || row == nelem) then
            writeln("IBM analytical reduced exact Jacobian assembly row ", row, " / ", nelem);

        const rowDgamma = residualRowDerivativeAnalyticalIBM(sd, row, c_ptrTo(phiGlobal[0]),
                                                             c_ptrTo(gammaGlobal), 0, true, n);
        const gammaFactor = rowDgamma / kuttaDgamma;
        const mergedStart = this.rowMergedOffsets[row];
        const mergedStop = this.rowMergedOffsets[row + 1];

        for mergedIdx in mergedStart..<mergedStop {
            const col = this.rowMergedCols[mergedIdx];
            const rowValue =
                if this.rowMergedRowIdx[mergedIdx] >= 0 then
                    residualRowDerivativeAnalyticalIBM(sd, row, c_ptrTo(phiGlobal[0]),
                                                       c_ptrTo(gammaGlobal), col, false, n)
                else
                    0.0;
            const kuttaValue =
                if this.rowMergedKuttaIdx[mergedIdx] >= 0 then
                    kuttaDphi[col - 1]
                else
                    0.0;
            const reducedValue = rowValue - gammaFactor * kuttaValue;
            if abs(reducedValue) > AD_ROW_PRINT_TOL || col == row then
                this.addJacobianEntry(row - 1, col - 1, reducedValue,
                                      "computeIBMAnalyticalReducedExactJacobian merged row entry");
        }
    }

    this.A_petsc.assemblyComplete();
}
}
