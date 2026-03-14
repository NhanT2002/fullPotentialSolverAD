module temporalDiscretization
{
use mesh;
use writeCGNS;
use Math;
use linearAlgebra;
use Time;
import input.potentialInputs;
use spatialDiscretization;
use PETSCapi;
use C_PETSC;
use petsc;
use CTypes;
use List;
use Set;
use gmres;
use Sort;
use IO;
// Keep nonlinear iteration control here and split the Jacobian assembly
// implementations into dedicated extension-method modules.
use temporalDiscretizationAD;
use temporalDiscretizationAnalyticalApprox;
use temporalDiscretizationAnalyticalExact;

config const CHECK_AD_ROW = false;
config const AD_ROW = 0;
config const AD_ROW_CHECK_ONLY = false;
config const AD_FD_EPS = 1.0e-7;
config const AD_ROW_PRINT_TOL = 1.0e-9;
config const AD_ROW_MAX_PRINT = 40;
config const AD_JACOBIAN_PROGRESS_FREQ = 0;

extern {
    int enzyme_dup;
    int enzyme_const;
    void __enzyme_autodiff(void*, ...);
}

inline proc phiFromPtr(phiPtr: c_ptr(real(64)), elem: int): real(64) {
    return phiPtr[elem - 1];
}

inline proc gammaFromPtr(gammaPtr: c_ptr(real(64))): real(64) {
    return gammaPtr[0];
}

inline proc phiSeedDerivative(elem: int, seedCol: int, wrtGamma: bool): real(64) {
    if wrtGamma then
        return 0.0;
    return if elem == seedCol then 1.0 else 0.0;
}

inline proc gammaSeedDerivative(wrtGamma: bool): real(64) {
    return if wrtGamma then 1.0 else 0.0;
}

inline proc isReducedExactJacobianType(jacobianType: string): bool {
    return jacobianType == "ad_reduced_exact" || jacobianType == "analytical_reduced_exact";
}

proc isWallFace(sd: borrowed spatialDiscretization, face: int): bool {
    return sd.wallFaceSet_.contains(face);
}

proc ghostPhiForFace(sd: borrowed spatialDiscretization,
                     phiPtr: c_ptr(real(64)),
                     face: int,
                     interiorElem: int,
                     ghostElem: int): real(64) {
    if isWallFace(sd, face) then
        return phiFromPtr(phiPtr, interiorElem);

    const x = sd.elemCentroidX_[ghostElem];
    const y = sd.elemCentroidY_[ghostElem];
    if sd.farfieldIsCylinder_ {
        const R = sd.inputs_.CYLINDER_RADIUS_;
        const r2 = x * x + y * y;
        const r = sqrt(r2);
        const theta = atan2(y, x);
        return sd.inputs_.VEL_INF_ * (r + R * R / r) * cos(theta);
    }

    return sd.inputs_.U_INF_ * x + sd.inputs_.V_INF_ * y;
}

proc ghostPhiForFaceWithDerivative(sd: borrowed spatialDiscretization,
                                   phiPtr: c_ptr(real(64)),
                                   face: int,
                                   interiorElem: int,
                                   ghostElem: int,
                                   seedCol: int,
                                   wrtGamma: bool): (real(64), real(64)) {
    if isWallFace(sd, face) then
        return (phiFromPtr(phiPtr, interiorElem),
                phiSeedDerivative(interiorElem, seedCol, wrtGamma));

    const x = sd.elemCentroidX_[ghostElem];
    const y = sd.elemCentroidY_[ghostElem];
    if sd.farfieldIsCylinder_ {
        const R = sd.inputs_.CYLINDER_RADIUS_;
        const r2 = x * x + y * y;
        const r = sqrt(r2);
        const theta = atan2(y, x);
        return (sd.inputs_.VEL_INF_ * (r + R * R / r) * cos(theta), 0.0);
    }

    return (sd.inputs_.U_INF_ * x + sd.inputs_.V_INF_ * y, 0.0);
}

proc computeVelocityForCellFromPtr(sd: borrowed spatialDiscretization,
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
            phiJ = ghostPhiForFace(sd, phiPtr, face, elem, neighbor);
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

proc computeVelocityForCellFromPtrWithDerivative(sd: borrowed spatialDiscretization,
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
            (phiJ, dphiJ) = ghostPhiForFaceWithDerivative(sd, phiPtr, face, elem, neighbor,
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

proc computeDensityAndMuFromVelocity(sd: borrowed spatialDiscretization,
                                     u: real(64),
                                     v: real(64)): (real(64), real(64), real(64)) {
    const rho = (1.0 + sd.gamma_minus_one_over_two_ * sd.inputs_.MACH_ * sd.inputs_.MACH_ *
                (1.0 - u * u - v * v)) ** sd.one_over_gamma_minus_one_;
    const mach = sd.mach(u, v, rho);
    const M2 = mach * mach;
    const Mc2 = sd.inputs_.MACH_C_ * sd.inputs_.MACH_C_;
    const excessMach2 = max(0.0, M2 - Mc2);
    const mu = sd.inputs_.MU_C_ * excessMach2;
    return (rho, mach, mu);
}

proc computeDensityAndMuFromVelocityWithDerivative(sd: borrowed spatialDiscretization,
                                                   u: real(64),
                                                   v: real(64),
                                                   du: real(64),
                                                   dv: real(64)): (real(64), real(64), real(64), real(64)) {
    const machInf2 = sd.inputs_.MACH_ * sd.inputs_.MACH_;
    const a = sd.gamma_minus_one_over_two_ * machInf2;
    const B = 1.0 + a * (1.0 - u * u - v * v);
    const rho = B ** sd.one_over_gamma_minus_one_;
    const dB = -2.0 * a * (u * du + v * dv);
    const drho = sd.one_over_gamma_minus_one_ * (B ** (sd.one_over_gamma_minus_one_ - 1.0)) * dB;

    const vel2 = u * u + v * v;
    const rhoPow = rho ** (1.0 - sd.inputs_.GAMMA_);
    const M2 = machInf2 * vel2 * rhoPow;
    var dM2 = machInf2 * rhoPow * (2.0 * u * du + 2.0 * v * dv);
    dM2 += machInf2 * vel2 * (1.0 - sd.inputs_.GAMMA_) * (rho ** (-sd.inputs_.GAMMA_)) * drho;

    const Mc2 = sd.inputs_.MACH_C_ * sd.inputs_.MACH_C_;
    const mu = if M2 > Mc2 then sd.inputs_.MU_C_ * (M2 - Mc2) else 0.0;
    const dmu = if M2 > Mc2 then sd.inputs_.MU_C_ * dM2 else 0.0;

    return (rho, mu, drho, dmu);
}

private proc computeGhostVelocityForFace(sd: borrowed spatialDiscretization,
                                         phiPtr: c_ptr(real(64)),
                                         gammaPtr: c_ptr(real(64)),
                                         face: int,
                                         interiorElem: int): (real(64), real(64)) {
    const (uInt, vInt) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, interiorElem);

    if isWallFace(sd, face) {
        const nx = sd.faceNormalX_[face];
        const ny = sd.faceNormalY_[face];
        const vDotN = uInt * nx + vInt * ny;
        return (uInt - 2.0 * vDotN * nx, vInt - 2.0 * vDotN * ny);
    }

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

private proc computeGhostVelocityForFaceWithDerivative(sd: borrowed spatialDiscretization,
                                                       phiPtr: c_ptr(real(64)),
                                                       gammaPtr: c_ptr(real(64)),
                                                       face: int,
                                                       interiorElem: int,
                                                       seedCol: int,
                                                       wrtGamma: bool): (real(64), real(64), real(64), real(64)) {
    const (uInt, vInt, duInt, dvInt) =
        computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, interiorElem,
                                                    seedCol, wrtGamma);

    if isWallFace(sd, face) {
        const nx = sd.faceNormalX_[face];
        const ny = sd.faceNormalY_[face];
        const vDotN = uInt * nx + vInt * ny;
        const dvDotN = duInt * nx + dvInt * ny;
        return (uInt - 2.0 * vDotN * nx,
                vInt - 2.0 * vDotN * ny,
                duInt - 2.0 * dvDotN * nx,
                dvInt - 2.0 * dvDotN * ny);
    }

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

private proc computeRhoMuForElemOnFace(sd: borrowed spatialDiscretization,
                                       phiPtr: c_ptr(real(64)),
                                       gammaPtr: c_ptr(real(64)),
                                       face: int,
                                       elem: int): (real(64), real(64)) {
    if elem <= sd.nelemDomain_ {
        const (u, v) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, elem);
        const (rho, mach, mu) = computeDensityAndMuFromVelocity(sd, u, v);
        return (rho, mu);
    }

    const elem1 = sd.mesh_.edge2elem_[1, face];
    const elem2 = sd.mesh_.edge2elem_[2, face];
    const interiorElem = if elem1 <= sd.nelemDomain_ then elem1 else elem2;

    if isWallFace(sd, face) {
        const (u, v) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, interiorElem);
        const (rho, mach, mu) = computeDensityAndMuFromVelocity(sd, u, v);
        return (rho, mu);
    }

    const (uGhost, vGhost) = computeGhostVelocityForFace(sd, phiPtr, gammaPtr, face, interiorElem);
    const (rhoGhost, machGhost, muGhost) = computeDensityAndMuFromVelocity(sd, uGhost, vGhost);
    return (rhoGhost, muGhost);
}

private proc computeRhoMuForElemOnFaceWithDerivative(sd: borrowed spatialDiscretization,
                                                     phiPtr: c_ptr(real(64)),
                                                     gammaPtr: c_ptr(real(64)),
                                                     face: int,
                                                     elem: int,
                                                     seedCol: int,
                                                     wrtGamma: bool): (real(64), real(64), real(64), real(64)) {
    if elem <= sd.nelemDomain_ {
        const (u, v, du, dv) =
            computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, elem, seedCol, wrtGamma);
        const (rho, mu, drho, dmu) = computeDensityAndMuFromVelocityWithDerivative(sd, u, v, du, dv);
        return (rho, mu, drho, dmu);
    }

    const elem1 = sd.mesh_.edge2elem_[1, face];
    const elem2 = sd.mesh_.edge2elem_[2, face];
    const interiorElem = if elem1 <= sd.nelemDomain_ then elem1 else elem2;

    if isWallFace(sd, face) {
        const (u, v, du, dv) =
            computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, interiorElem,
                                                        seedCol, wrtGamma);
        const (rho, mu, drho, dmu) = computeDensityAndMuFromVelocityWithDerivative(sd, u, v, du, dv);
        return (rho, mu, drho, dmu);
    }

    const (uGhost, vGhost, duGhost, dvGhost) =
        computeGhostVelocityForFaceWithDerivative(sd, phiPtr, gammaPtr, face, interiorElem,
                                                  seedCol, wrtGamma);
    const (rhoGhost, muGhost, drhoGhost, dmuGhost) =
        computeDensityAndMuFromVelocityWithDerivative(sd, uGhost, vGhost, duGhost, dvGhost);
    return (rhoGhost, muGhost, drhoGhost, dmuGhost);
}

proc residualRowForAD(sd: borrowed spatialDiscretization,
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
        if elem1 <= sd.nelemDomain_ {
            (u1, v1) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, elem1);
        } else {
            (u1, v1) = computeGhostVelocityForFace(sd, phiPtr, gammaPtr, face, elem2);
        }

        if elem2 <= sd.nelemDomain_ {
            (u2, v2) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, elem2);
        } else {
            (u2, v2) = computeGhostVelocityForFace(sd, phiPtr, gammaPtr, face, elem1);
        }

        const uAvg = sd.weights1_[face] * u1 + sd.weights2_[face] * u2;
        const vAvg = sd.weights1_[face] * v1 + sd.weights2_[face] * v2;

        var phi1: real(64);
        var phi2: real(64);
        if elem1 <= sd.nelemDomain_ then
            phi1 = phiFromPtr(phiPtr, elem1);
        else
            phi1 = ghostPhiForFace(sd, phiPtr, face, elem2, elem1);

        if elem2 <= sd.nelemDomain_ then
            phi2 = phiFromPtr(phiPtr, elem2);
        else
            phi2 = ghostPhiForFace(sd, phiPtr, face, elem1, elem2);

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

        const (rho1Blend, mu1Blend) = computeRhoMuForElemOnFace(sd, phiPtr, gammaPtr, face, elem1);
        const (rho2Blend, mu2Blend) = computeRhoMuForElemOnFace(sd, phiPtr, gammaPtr, face, elem2);
        const upwindWeight = if vDotN >= 0.0 then 1.0 else 0.0;
        const rhoUpwind = upwindWeight * rho1Blend + (1.0 - upwindWeight) * rho2Blend;
        const muUpwind = upwindWeight * mu1Blend + (1.0 - upwindWeight) * mu2Blend;
        if muUpwind > 0.0 then
            rhoFace = rhoFace - muUpwind * (rhoFace - rhoUpwind);

        res += sign * rhoFace * vDotN * sd.faceArea_[face];
    }

    return res * sd.res_scale_;
}

proc residualRowDerivativeAnalytical(sd: borrowed spatialDiscretization,
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
        const elem1 = sd.mesh_.edge2elem_[1, face];
        const elem2 = sd.mesh_.edge2elem_[2, face];
        const sign = if elem1 == elem then 1.0 else -1.0;

        var u1, v1, du1, dv1, u2, v2, du2, dv2: real(64);
        if elem1 <= sd.nelemDomain_ {
            (u1, v1, du1, dv1) =
                computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, elem1, seedCol, wrtGamma);
        } else {
            (u1, v1, du1, dv1) =
                computeGhostVelocityForFaceWithDerivative(sd, phiPtr, gammaPtr, face, elem2, seedCol, wrtGamma);
        }

        if elem2 <= sd.nelemDomain_ {
            (u2, v2, du2, dv2) =
                computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, elem2, seedCol, wrtGamma);
        } else {
            (u2, v2, du2, dv2) =
                computeGhostVelocityForFaceWithDerivative(sd, phiPtr, gammaPtr, face, elem1, seedCol, wrtGamma);
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
            (phi1, dphi1) = ghostPhiForFaceWithDerivative(sd, phiPtr, face, elem2, elem1,
                                                          seedCol, wrtGamma);
        }

        if elem2 <= sd.nelemDomain_ {
            phi2 = phiFromPtr(phiPtr, elem2);
            dphi2 = phiSeedDerivative(elem2, seedCol, wrtGamma);
        } else {
            (phi2, dphi2) = ghostPhiForFaceWithDerivative(sd, phiPtr, face, elem1, elem2,
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
            computeRhoMuForElemOnFaceWithDerivative(sd, phiPtr, gammaPtr, face, elem1, seedCol, wrtGamma);
        const (rho2Blend, mu2Blend, drho2Blend, dmu2Blend) =
            computeRhoMuForElemOnFaceWithDerivative(sd, phiPtr, gammaPtr, face, elem2, seedCol, wrtGamma);
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

        dres += sign * sd.faceArea_[face] * (drhoFace * vDotN + rhoFace * dvDotN);
    }

    return dres * sd.res_scale_;
}

proc kuttaResidualForAD(sd: borrowed spatialDiscretization,
                        phiPtr: c_ptr(real(64)),
                        gammaPtr: c_ptr(real(64)),
                        n: c_int): real(64) {
    const gamma = gammaFromPtr(gammaPtr);
    const upperElem = sd.upperTEelem_;
    const lowerElem = sd.lowerTEelem_;
    const (uUpper, vUpper) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, upperElem);
    const (uLower, vLower) = computeVelocityForCellFromPtr(sd, phiPtr, gammaPtr, lowerElem);

    const phiUpper = phiFromPtr(phiPtr, upperElem) +
                     (uUpper * sd.deltaSupperTEx_ + vUpper * sd.deltaSupperTEy_);
    const phiLower = phiFromPtr(phiPtr, lowerElem) +
                     (uLower * sd.deltaSlowerTEx_ + vLower * sd.deltaSlowerTEy_);
    const gammaComputed = phiUpper - phiLower;

    return (gamma - gammaComputed) * sd.res_scale_;
}

proc kuttaResidualDerivativeAnalytical(sd: borrowed spatialDiscretization,
                                               phiPtr: c_ptr(real(64)),
                                               gammaPtr: c_ptr(real(64)),
                                               seedCol: int,
                                               wrtGamma: bool,
                                               n: c_int): real(64) {
    const dgamma = gammaSeedDerivative(wrtGamma);
    const upperElem = sd.upperTEelem_;
    const lowerElem = sd.lowerTEelem_;
    const (uUpper, vUpper, duUpper, dvUpper) =
        computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, upperElem, seedCol, wrtGamma);
    const (uLower, vLower, duLower, dvLower) =
        computeVelocityForCellFromPtrWithDerivative(sd, phiPtr, gammaPtr, lowerElem, seedCol, wrtGamma);

    const dphiUpper = phiSeedDerivative(upperElem, seedCol, wrtGamma) +
                      (duUpper * sd.deltaSupperTEx_ + dvUpper * sd.deltaSupperTEy_);
    const dphiLower = phiSeedDerivative(lowerElem, seedCol, wrtGamma) +
                      (duLower * sd.deltaSlowerTEx_ + dvLower * sd.deltaSlowerTEy_);

    return (dgamma - (dphiUpper - dphiLower)) * sd.res_scale_;
}

proc solveConsistentGammaForPhi(sd: borrowed spatialDiscretization,
                                        phiPtr: c_ptr(real(64)),
                                        gammaInitial: real(64),
                                        n: c_int,
                                        maxIts: int = 12,
                                        tol: real(64) = 1.0e-14): real(64) {
    var gamma = gammaInitial;
    for gammaIt in 1..maxIts {
        const residual = kuttaResidualForAD(sd, phiPtr, c_ptrTo(gamma), n);
        if abs(residual) < tol then
            break;

        const dgamma = kuttaResidualDerivativeAnalytical(sd, phiPtr, c_ptrTo(gamma), 0, true, n);
        if abs(dgamma) < 1.0e-14 then
            break;
        gamma -= residual / dgamma;
    }
    return gamma;
}

proc findStencilIndex(stencil: [?D] int, col: int): int {
    var lo = D.low;
    var hi = D.high;

    while lo <= hi {
        const mid = (lo + hi) / 2;
        const value = stencil[mid];
        if value == col then
            return mid;
        if value < col then
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

proc accumulateCellVelocitySensitivityOnStencil(sd: borrowed spatialDiscretization,
                                                        elem: int,
                                                        stencil: [?D] int,
                                                        ref du: [D] real(64),
                                                        ref dv: [D] real(64),
                                                        ref duGamma: real(64),
                                                        ref dvGamma: real(64)) {
    const elemIdx = findStencilIndex(stencil, elem);
    const faceStart = sd.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = sd.mesh_.elem2edgeIndex_[elem + 1];

    for faceIdx in faceStart..faceEnd {
        const face = sd.mesh_.elem2edge_[faceIdx];
        const elem1 = sd.mesh_.edge2elem_[1, face];
        const elem2 = sd.mesh_.edge2elem_[2, face];
        const neighbor = if elem1 == elem then elem2 else elem1;

        var wx, wy: real(64);
        if elem == elem1 {
            wx = sd.lsGradQR_!.wxFinal1_[face];
            wy = sd.lsGradQR_!.wyFinal1_[face];
        } else {
            wx = sd.lsGradQR_!.wxFinal2_[face];
            wy = sd.lsGradQR_!.wyFinal2_[face];
        }

        if neighbor <= sd.nelemDomain_ {
            const neighborIdx = findStencilIndex(stencil, neighbor);
            if neighborIdx >= 0 {
                du[neighborIdx] += wx;
                dv[neighborIdx] += wy;
            }
            if elemIdx >= 0 {
                du[elemIdx] -= wx;
                dv[elemIdx] -= wy;
            }

            const elemKuttaType = sd.kuttaCell_[elem];
            const neighborKuttaType = sd.kuttaCell_[neighbor];
            if elemKuttaType == 1 && neighborKuttaType == -1 {
                duGamma += wx;
                dvGamma += wy;
            } else if elemKuttaType == -1 && neighborKuttaType == 1 {
                duGamma -= wx;
                dvGamma -= wy;
            }
        } else if !isWallFace(sd, face) {
            if elemIdx >= 0 {
                du[elemIdx] -= wx;
                dv[elemIdx] -= wy;
            }
        }
    }
}

proc accumulateDensityMuSensitivityFromVelocityOnStencil(sd: borrowed spatialDiscretization,
                                                                 u: real(64),
                                                                 v: real(64),
                                                                 const ref du: [?D] real(64),
                                                                 const ref dv: [D] real(64),
                                                                 duGamma: real(64),
                                                                 dvGamma: real(64),
                                                                 ref drho: [D] real(64),
                                                                 ref dmu: [D] real(64),
                                                                 ref drhoGamma: real(64),
                                                                 ref dmuGamma: real(64),
                                                                 out rho: real(64),
                                                                 out mu: real(64)) {
    const machInf2 = sd.inputs_.MACH_ * sd.inputs_.MACH_;
    const a = sd.gamma_minus_one_over_two_ * machInf2;
    const B = 1.0 + a * (1.0 - u * u - v * v);
    rho = B ** sd.one_over_gamma_minus_one_;
    const rhoFactor = -2.0 * a * sd.one_over_gamma_minus_one_ *
                      (B ** (sd.one_over_gamma_minus_one_ - 1.0));

    for idx in D do
        drho[idx] = rhoFactor * (u * du[idx] + v * dv[idx]);
    drhoGamma = rhoFactor * (u * duGamma + v * dvGamma);

    const vel2 = u * u + v * v;
    const rhoPow = rho ** (1.0 - sd.inputs_.GAMMA_);
    const M2 = machInf2 * vel2 * rhoPow;
    const Mc2 = sd.inputs_.MACH_C_ * sd.inputs_.MACH_C_;
    if M2 > Mc2 {
        mu = sd.inputs_.MU_C_ * (M2 - Mc2);
        const common = machInf2;
        for idx in D {
            const dM2 = common * rhoPow * (2.0 * u * du[idx] + 2.0 * v * dv[idx]) +
                        common * vel2 * (1.0 - sd.inputs_.GAMMA_) *
                        (rho ** (-sd.inputs_.GAMMA_)) * drho[idx];
            dmu[idx] = sd.inputs_.MU_C_ * dM2;
        }
        const dM2Gamma = common * rhoPow * (2.0 * u * duGamma + 2.0 * v * dvGamma) +
                         common * vel2 * (1.0 - sd.inputs_.GAMMA_) *
                         (rho ** (-sd.inputs_.GAMMA_)) * drhoGamma;
        dmuGamma = sd.inputs_.MU_C_ * dM2Gamma;
    } else {
        mu = 0.0;
        dmu = 0.0;
        dmuGamma = 0.0;
    }
}

proc computeAnalyticalKuttaDerivativesOnStencil(sd: borrowed spatialDiscretization,
                                                         stencil: [?D] int,
                                                         ref kuttaDphi: [D] real(64),
                                                         ref kuttaDgamma: real(64)) {
    const upperElem = sd.upperTEelem_;
    const lowerElem = sd.lowerTEelem_;
    var duUpper: [D] real(64) = 0.0;
    var dvUpper: [D] real(64) = 0.0;
    var duLower: [D] real(64) = 0.0;
    var dvLower: [D] real(64) = 0.0;
    var duUpperGamma = 0.0;
    var dvUpperGamma = 0.0;
    var duLowerGamma = 0.0;
    var dvLowerGamma = 0.0;

    accumulateCellVelocitySensitivityOnStencil(sd, upperElem, stencil,
                                               duUpper, dvUpper, duUpperGamma, dvUpperGamma);
    accumulateCellVelocitySensitivityOnStencil(sd, lowerElem, stencil,
                                               duLower, dvLower, duLowerGamma, dvLowerGamma);

    const upperIdx = findStencilIndex(stencil, upperElem);
    const lowerIdx = findStencilIndex(stencil, lowerElem);
    if upperIdx >= 0 then
        kuttaDphi[upperIdx] -= sd.res_scale_;
    if lowerIdx >= 0 then
        kuttaDphi[lowerIdx] += sd.res_scale_;

    for idx in D {
        kuttaDphi[idx] += sd.res_scale_ *
                          (-sd.deltaSupperTEx_ * duUpper[idx] - sd.deltaSupperTEy_ * dvUpper[idx] +
                            sd.deltaSlowerTEx_ * duLower[idx] + sd.deltaSlowerTEy_ * dvLower[idx]);
    }

    kuttaDgamma = sd.res_scale_ *
                  (1.0 - sd.deltaSupperTEx_ * duUpperGamma - sd.deltaSupperTEy_ * dvUpperGamma +
                   sd.deltaSlowerTEx_ * duLowerGamma + sd.deltaSlowerTEy_ * dvLowerGamma);
}

class temporalDiscretization {
    var spatialDisc_: shared spatialDiscretization;
    var inputs_: potentialInputs;
    var it_: int = 0;
    var t0_: real(64) = 0.0;
    var first_res_: real(64) = 1e12;
    
    var A_petsc : owned PETSCmatrix_c;
    var x_petsc : owned PETSCvector_c;
    var b_petsc : owned PETSCvector_c;
    var ksp : owned PETSCksp_c;

    var timeList_ = new list(real(64));
    var itList_ = new list(int);
    var resList_ = new list(real(64));
    var clList_ = new list(real(64));
    var cdList_ = new list(real(64));
    var cmList_ = new list(real(64));
    var circulationList_ = new list(real(64));

    // Gradient sensitivity to circulation: ∂(∇φ)/∂Γ for each cell
    // These capture how each cell's gradient depends on Γ through wake-crossing faces
    var gradSensitivity_dom: domain(1) = {1..0};
    var dgradX_dGamma_: [gradSensitivity_dom] real(64);
    var dgradY_dGamma_: [gradSensitivity_dom] real(64);
    var Jij_: [gradSensitivity_dom] real(64);
    var rowStencilOffsetDom: domain(1) = {1..0};
    var rowStencilOffsets: [rowStencilOffsetDom] int;
    var rowStencilDataDom: domain(1) = {0..<0};
    var rowStencilData: [rowStencilDataDom] int;
    var rowFaceOffsetDom: domain(1) = {1..0};
    var rowFaceOffsets: [rowFaceOffsetDom] int;
    var rowFaceDataDom: domain(1) = {0..<0};
    var rowFaceData: [rowFaceDataDom] int;
    var rowFaceMinusIdx: [rowFaceDataDom] int;
    var rowFacePlusIdx: [rowFaceDataDom] int;
    var rowMergedOffsetDom: domain(1) = {1..0};
    var rowMergedOffsets: [rowMergedOffsetDom] int;
    var rowMergedDataDom: domain(1) = {0..<0};
    var rowMergedCols: [rowMergedDataDom] int;
    var rowMergedRowIdx: [rowMergedDataDom] int;
    var rowMergedKuttaIdx: [rowMergedDataDom] int;
    var cellVelTemplateOffsetDom: domain(1) = {1..0};
    var cellVelTemplateOffsets: [cellVelTemplateOffsetDom] int;
    var cellVelTemplateDataDom: domain(1) = {0..<0};
    var cellVelTemplateCols: [cellVelTemplateDataDom] int;
    var cellVelTemplateUx: [cellVelTemplateDataDom] real(64);
    var cellVelTemplateUy: [cellVelTemplateDataDom] real(64);
    var cellVelTemplateGammaDom: domain(1) = {1..0};
    var cellVelTemplateGammaUx: [cellVelTemplateGammaDom] real(64);
    var cellVelTemplateGammaUy: [cellVelTemplateGammaDom] real(64);
    var kuttaStencilCacheDom: domain(1) = {0..<0};
    var kuttaStencilCache: [kuttaStencilCacheDom] int;

    proc init(spatialDisc: shared spatialDiscretization, ref inputs: potentialInputs) {
        writeln("Initializing temporal discretization...");
        this.spatialDisc_ = spatialDisc;
        this.inputs_ = inputs;

        const reducedExactSystem = isReducedExactJacobianType(this.inputs_.JACOBIAN_TYPE_);
        const systemSize = spatialDisc.nelemDomain_;
        const M = systemSize;
        const N = systemSize;
        
        this.A_petsc = new owned PETSCmatrix_c(PETSC_COMM_SELF, "seqaij", M, M, N, N);
        this.x_petsc = new owned PETSCvector_c(PETSC_COMM_SELF, N, N, 0.0, "seq");
        this.b_petsc = new owned PETSCvector_c(PETSC_COMM_SELF, N, N, 0.0, "seq");

        var nnz : [0..M-1] PetscInt;
        const baseNNZ = 4 * (this.spatialDisc_.mesh_.elem2edge_[
            this.spatialDisc_.mesh_.elem2edgeIndex_[1] + 1 ..
            this.spatialDisc_.mesh_.elem2edgeIndex_[1 + 1]
        ].size + 1);
        if reducedExactSystem {
            // The reduced phi-only exact Jacobian is a Schur complement:
            // J_red = dR/dphi - dR/dGamma * (dK/dphi)/(dK/dGamma).
            // Each row keeps its local stencil plus the Kutta row support, so
            // it needs a noticeably wider preallocation than the standard local Jacobian.
            nnz = max(64, baseNNZ + 16);
        } else {
            nnz = baseNNZ;
        }
        A_petsc.preAllocate(nnz);

        this.ksp = new owned PETSCksp_c(PETSC_COMM_SELF, "gmres");
        this.ksp.setTolerances(inputs.GMRES_RTOL_, inputs.GMRES_ATOL_, inputs.GMRES_DTOL_, inputs.GMRES_MAXIT_);
        this.ksp.GMRESSetRestart(inputs.GMRES_RESTART_);
        this.ksp.GMRESSetPreAllocateVectors();
        if this.inputs_.GMRES_PRECON_ == "jacobi" {
            writeln("Using Jacobi preconditioner for GMRES");
            this.ksp.setPreconditioner("jacobi");
        } else if this.inputs_.GMRES_PRECON_ == "ilu" {
            writeln("Using ILU preconditioner for GMRES");
            this.ksp.setPreconditioner("ilu");
        } else if this.inputs_.GMRES_PRECON_ == "lu" {
            writeln("Using lu preconditioner for GMRES");
            this.ksp.setPreconditioner("lu");
        } else if this.inputs_.GMRES_PRECON_ == "asm" {
            writeln("Using asm preconditioner for GMRES");
            this.ksp.setPreconditioner("asm");
        } else if this.inputs_.GMRES_PRECON_ == "gasm" {
            writeln("Using gasm preconditioner for GMRES");
            this.ksp.setPreconditioner("gasm");
        } else if this.inputs_.GMRES_PRECON_ == "bjacobi" {
            writeln("Using bjacobi preconditioner for GMRES");
            this.ksp.setPreconditioner("bjacobi");
        } else if this.inputs_.GMRES_PRECON_ == "none" {
            writeln("Using no preconditioner for GMRES");
            this.ksp.setPreconditioner("none");
        } else {
            writeln("No preconditioner for GMRES");
            this.ksp.setPreconditioner("none");
        }
        
        // Initialize gradient sensitivity arrays
        this.gradSensitivity_dom = {1..spatialDisc.nelemDomain_};
    }

    proc initializeJacobian() {
        if isReducedExactJacobianType(this.inputs_.JACOBIAN_TYPE_) {
            const kuttaStencil = this.buildKuttaStencil();
            for elem in 1..this.spatialDisc_.nelemDomain_ {
                var stencilSet = new set(int);
                const rowStencil = this.buildRowStencil(elem);
                for col in rowStencil do stencilSet.add(col);
                for col in kuttaStencil do stencilSet.add(col);
                for col in stencilSet do
                    this.A_petsc.set(elem - 1, col - 1, 0.0);
            }

            this.A_petsc.assemblyComplete();
            return;
        }

        forall elem in 1..this.spatialDisc_.nelemDomain_ {
            this.A_petsc.set(elem-1, elem-1, 0.0);
            const faces = this.spatialDisc_.mesh_.elem2edge_[this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1 .. this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1]];
            var isWakeFace = false;
            for face in faces {
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
                const neighbor = if elem1 == elem then elem2 else elem1;
                if neighbor <= this.spatialDisc_.nelemDomain_ {
                    this.A_petsc.set(elem-1, neighbor-1, 0.0);
                    if this.dgradX_dGamma_[neighbor] != 0.0 || this.dgradY_dGamma_[neighbor] != 0.0 {
                        isWakeFace = true;
                    }
                }
                if this.dgradX_dGamma_[elem] != 0.0 || this.dgradY_dGamma_[elem] != 0.0 {
                    isWakeFace = true;
                }
            }
            if isWakeFace {
                const wakeFaceIndex = this.spatialDisc_.wakeFaceIndexInfluenceOnElem_[elem];
                const upperTE_influences = this.spatialDisc_.wakeFaceUpper_[wakeFaceIndex];
                const lowerTE_influences = this.spatialDisc_.wakeFaceLower_[wakeFaceIndex];

                this.A_petsc.set(elem-1, upperTE_influences - 1, 0.0);
                this.A_petsc.set(elem-1, lowerTE_influences - 1, 0.0);
            }
        }

        this.A_petsc.assemblyComplete();
        // this.A_petsc.matView();
    }

    proc computeJacobian() {
        // Dispatch to the selected Jacobian implementation so the solve loop
        // stays separate from the assembly details.
        if this.inputs_.JACOBIAN_TYPE_ == "ad_reduced_exact" {
            this.computeADReducedExactJacobian();
            return;
        }
        if this.inputs_.JACOBIAN_TYPE_ == "analytical_reduced_exact" {
            this.computeAnalyticalReducedExactJacobian();
            return;
        }
        this.computeApproximateAnalyticalJacobian();
    }

    proc initialize() {
        this.spatialDisc_.initializeMetrics();
        this.spatialDisc_.initializeKuttaCells();
        this.spatialDisc_.initializeSolution();
        this.spatialDisc_.run();
        if isReducedExactJacobianType(this.inputs_.JACOBIAN_TYPE_) {
            this.enforceConsistentGammaForCurrentPhi();
            this.spatialDisc_.run();
        }
        this.computeGradientSensitivity();
        this.initializeStencilCache();
        this.initializeCellVelocityTemplateCache();
        this.initializeJacobian();
        
        this.computeJacobian();

        if this.inputs_.START_FILENAME_ != "" {
            writeln("Initializing solution from file: ", this.inputs_.START_FILENAME_);
            const (xElem, yElem, rho, phi, it, time, res, cl, cd, cm, circulation, wakeGamma) = readSolution(this.inputs_.START_FILENAME_);
            for i in it.domain {
                this.timeList_.pushBack(time[i]);
                this.itList_.pushBack(it[i]);
                this.resList_.pushBack(res[i]);
                this.clList_.pushBack(cl[i]);
                this.cdList_.pushBack(cd[i]);
                this.cmList_.pushBack(cm[i]);
                this.circulationList_.pushBack(circulation[i]);
            }
            this.it_ = it.last;
            this.t0_ = time.last;
            this.first_res_ = res.first;
        }

        if CHECK_AD_ROW {
            this.runADRowCheck();
        }
    }

    proc buildRowStencilUncached(row: int) {
        var stencilSet = new set(int);
        var firstRing = new set(int);

        stencilSet.add(row);
        firstRing.add(row);

        const rowNeighbors = this.spatialDisc_.mesh_.esuel_[
            this.spatialDisc_.mesh_.esuelIndex_[row] + 1 ..
            this.spatialDisc_.mesh_.esuelIndex_[row + 1]];
        for neighbor in rowNeighbors {
            if neighbor <= this.spatialDisc_.nelemDomain_ {
                stencilSet.add(neighbor);
                firstRing.add(neighbor);
            }
        }

        for elem in firstRing {
            const neighbors = this.spatialDisc_.mesh_.esuel_[
                this.spatialDisc_.mesh_.esuelIndex_[elem] + 1 ..
                this.spatialDisc_.mesh_.esuelIndex_[elem + 1]];
            for neighbor in neighbors do
                if neighbor <= this.spatialDisc_.nelemDomain_ then
                    stencilSet.add(neighbor);
        }

        const wakeInfluence = this.spatialDisc_.wakeFaceIndexInfluenceOnElem_[row];
        if wakeInfluence > 0 {
            stencilSet.add(this.spatialDisc_.wakeFaceUpper_[wakeInfluence]);
            stencilSet.add(this.spatialDisc_.wakeFaceLower_[wakeInfluence]);
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

    proc buildRowStencil(row: int) {
        if this.rowStencilOffsets.size > 0 {
            const start = this.rowStencilOffsets[row];
            const stop = this.rowStencilOffsets[row + 1];
            const dom = {0..<(stop - start)};
            var stencil: [dom] int;
            if stop > start then
                stencil = this.rowStencilData[start..<stop];
            return stencil;
        }

        return this.buildRowStencilUncached(row);
    }

    proc buildKuttaStencilUncached() {
        var stencilSet = new set(int);
        const upperStencil = this.buildRowStencilUncached(this.spatialDisc_.upperTEelem_);
        const lowerStencil = this.buildRowStencilUncached(this.spatialDisc_.lowerTEelem_);

        for col in upperStencil do stencilSet.add(col);
        for col in lowerStencil do stencilSet.add(col);

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

    proc buildKuttaStencil() {
        if this.kuttaStencilCache.size > 0 {
            const dom = {0..<this.kuttaStencilCache.size};
            var stencil: [dom] int = this.kuttaStencilCache;
            return stencil;
        }

        return this.buildKuttaStencilUncached();
    }

    proc initializeStencilCache() {
        const nelem = this.spatialDisc_.nelemDomain_;
        this.rowStencilOffsetDom = {1..nelem + 1};
        this.rowStencilOffsets = 0;
        this.rowFaceOffsetDom = {1..nelem + 1};
        this.rowFaceOffsets = 0;

        var totalSize = 0;
        var totalFaceEntries = 0;
        for row in 1..nelem {
            const faceStart = this.spatialDisc_.mesh_.elem2edgeIndex_[row] + 1;
            const faceEnd = this.spatialDisc_.mesh_.elem2edgeIndex_[row + 1];
            totalFaceEntries += faceEnd - faceStart + 1;

            this.rowStencilOffsets[row] = totalSize;
            totalSize += this.buildRowStencilUncached(row).size;
        }
        this.rowStencilOffsets[nelem + 1] = totalSize;

        this.rowStencilDataDom = {0..<totalSize};
        this.rowStencilData = 0;
        for row in 1..nelem {
            const stencil = this.buildRowStencilUncached(row);
            const start = this.rowStencilOffsets[row];
            for idx in stencil.domain do
                this.rowStencilData[start + idx] = stencil[idx];
        }

        const kuttaStencil = this.buildKuttaStencilUncached();
        this.kuttaStencilCacheDom = {0..<kuttaStencil.size};
        this.kuttaStencilCache = kuttaStencil;

        this.rowFaceDataDom = {0..<totalFaceEntries};
        this.rowFaceData = 0;
        this.rowFaceMinusIdx = -1;
        this.rowFacePlusIdx = -1;

        var totalMergedSize = 0;
        this.rowMergedOffsetDom = {1..nelem + 1};
        this.rowMergedOffsets = 0;
        for row in 1..nelem {
            const rowStencil = this.buildRowStencil(row);
            var rowIdx = rowStencil.domain.low;
            var kuttaIdx = kuttaStencil.domain.low;
            var mergedCount = 0;

            while rowIdx <= rowStencil.domain.high || kuttaIdx <= kuttaStencil.domain.high {
                if kuttaIdx > kuttaStencil.domain.high ||
                   (rowIdx <= rowStencil.domain.high && rowStencil[rowIdx] < kuttaStencil[kuttaIdx]) {
                    rowIdx += 1;
                } else if rowIdx > rowStencil.domain.high ||
                          kuttaStencil[kuttaIdx] < rowStencil[rowIdx] {
                    kuttaIdx += 1;
                } else {
                    rowIdx += 1;
                    kuttaIdx += 1;
                }
                mergedCount += 1;
            }

            this.rowMergedOffsets[row] = totalMergedSize;
            totalMergedSize += mergedCount;
        }
        this.rowMergedOffsets[nelem + 1] = totalMergedSize;

        this.rowMergedDataDom = {0..<totalMergedSize};
        this.rowMergedCols = 0;
        this.rowMergedRowIdx = -1;
        this.rowMergedKuttaIdx = -1;

        var faceOffset = 0;
        for row in 1..nelem {
            const rowStencil = this.buildRowStencil(row);
            const faceStart = this.spatialDisc_.mesh_.elem2edgeIndex_[row] + 1;
            const faceEnd = this.spatialDisc_.mesh_.elem2edgeIndex_[row + 1];

            this.rowFaceOffsets[row] = faceOffset;
            for faceIdx in faceStart..faceEnd {
                const face = this.spatialDisc_.mesh_.elem2edge_[faceIdx];
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];

                this.rowFaceData[faceOffset] = face;
                this.rowFaceMinusIdx[faceOffset] =
                    if elem1 <= this.spatialDisc_.nelemDomain_ then
                        findStencilIndex(rowStencil, elem1)
                    else if isWallFace(this.spatialDisc_.borrow(), face) then
                        findStencilIndex(rowStencil, elem2)
                    else
                        -1;
                this.rowFacePlusIdx[faceOffset] =
                    if elem2 <= this.spatialDisc_.nelemDomain_ then
                        findStencilIndex(rowStencil, elem2)
                    else if isWallFace(this.spatialDisc_.borrow(), face) then
                        findStencilIndex(rowStencil, elem1)
                    else
                        -1;
                faceOffset += 1;
            }

            const mergedStart = this.rowMergedOffsets[row];
            var mergedOffset = mergedStart;
            var rowIdx = rowStencil.domain.low;
            var kuttaIdx = kuttaStencil.domain.low;

            while rowIdx <= rowStencil.domain.high || kuttaIdx <= kuttaStencil.domain.high {
                if kuttaIdx > kuttaStencil.domain.high ||
                   (rowIdx <= rowStencil.domain.high && rowStencil[rowIdx] < kuttaStencil[kuttaIdx]) {
                    this.rowMergedCols[mergedOffset] = rowStencil[rowIdx];
                    this.rowMergedRowIdx[mergedOffset] = rowIdx;
                    rowIdx += 1;
                } else if rowIdx > rowStencil.domain.high ||
                          kuttaStencil[kuttaIdx] < rowStencil[rowIdx] {
                    this.rowMergedCols[mergedOffset] = kuttaStencil[kuttaIdx];
                    this.rowMergedKuttaIdx[mergedOffset] = kuttaIdx;
                    kuttaIdx += 1;
                } else {
                    this.rowMergedCols[mergedOffset] = rowStencil[rowIdx];
                    this.rowMergedRowIdx[mergedOffset] = rowIdx;
                    this.rowMergedKuttaIdx[mergedOffset] = kuttaIdx;
                    rowIdx += 1;
                    kuttaIdx += 1;
                }
                mergedOffset += 1;
            }
        }
        this.rowFaceOffsets[nelem + 1] = faceOffset;
    }

    proc initializeCellVelocityTemplateCache() {
        const nelem = this.spatialDisc_.nelemDomain_;
        this.cellVelTemplateOffsetDom = {1..nelem + 1};
        this.cellVelTemplateOffsets = 0;
        this.cellVelTemplateGammaDom = {1..nelem};
        this.cellVelTemplateGammaUx = 0.0;
        this.cellVelTemplateGammaUy = 0.0;

        var totalSize = 0;
        for elem in 1..nelem {
            var supportDom: domain(int);
            var uxCoeff: [supportDom] real(64);
            var uyCoeff: [supportDom] real(64);
            var gammaUx = 0.0;
            var gammaUy = 0.0;
            const faceStart = this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1;
            const faceEnd = this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1];

            for faceIdx in faceStart..faceEnd {
                const face = this.spatialDisc_.mesh_.elem2edge_[faceIdx];
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
                const neighbor = if elem1 == elem then elem2 else elem1;

                var wx, wy: real(64);
                if elem == elem1 {
                    wx = this.spatialDisc_.lsGradQR_!.wxFinal1_[face];
                    wy = this.spatialDisc_.lsGradQR_!.wyFinal1_[face];
                } else {
                    wx = this.spatialDisc_.lsGradQR_!.wxFinal2_[face];
                    wy = this.spatialDisc_.lsGradQR_!.wyFinal2_[face];
                }

                if neighbor <= nelem {
                    if !supportDom.contains(neighbor) then supportDom += neighbor;
                    uxCoeff[neighbor] += wx;
                    uyCoeff[neighbor] += wy;

                    if !supportDom.contains(elem) then supportDom += elem;
                    uxCoeff[elem] -= wx;
                    uyCoeff[elem] -= wy;

                    const elemKuttaType = this.spatialDisc_.kuttaCell_[elem];
                    const neighborKuttaType = this.spatialDisc_.kuttaCell_[neighbor];
                    if elemKuttaType == 1 && neighborKuttaType == -1 {
                        gammaUx += wx;
                        gammaUy += wy;
                    } else if elemKuttaType == -1 && neighborKuttaType == 1 {
                        gammaUx -= wx;
                        gammaUy -= wy;
                    }
                } else if !isWallFace(this.spatialDisc_.borrow(), face) {
                    if !supportDom.contains(elem) then supportDom += elem;
                    uxCoeff[elem] -= wx;
                    uyCoeff[elem] -= wy;
                }
            }

            this.cellVelTemplateOffsets[elem] = totalSize;
            totalSize += supportDom.size;
            this.cellVelTemplateGammaUx[elem] = gammaUx;
            this.cellVelTemplateGammaUy[elem] = gammaUy;
        }
        this.cellVelTemplateOffsets[nelem + 1] = totalSize;

        this.cellVelTemplateDataDom = {0..<totalSize};
        this.cellVelTemplateCols = 0;
        this.cellVelTemplateUx = 0.0;
        this.cellVelTemplateUy = 0.0;

        for elem in 1..nelem {
            var supportDom: domain(int);
            var uxCoeff: [supportDom] real(64);
            var uyCoeff: [supportDom] real(64);
            const faceStart = this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1;
            const faceEnd = this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1];

            for faceIdx in faceStart..faceEnd {
                const face = this.spatialDisc_.mesh_.elem2edge_[faceIdx];
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
                const neighbor = if elem1 == elem then elem2 else elem1;

                var wx, wy: real(64);
                if elem == elem1 {
                    wx = this.spatialDisc_.lsGradQR_!.wxFinal1_[face];
                    wy = this.spatialDisc_.lsGradQR_!.wyFinal1_[face];
                } else {
                    wx = this.spatialDisc_.lsGradQR_!.wxFinal2_[face];
                    wy = this.spatialDisc_.lsGradQR_!.wyFinal2_[face];
                }

                if neighbor <= nelem {
                    if !supportDom.contains(neighbor) then supportDom += neighbor;
                    uxCoeff[neighbor] += wx;
                    uyCoeff[neighbor] += wy;

                    if !supportDom.contains(elem) then supportDom += elem;
                    uxCoeff[elem] -= wx;
                    uyCoeff[elem] -= wy;
                } else if !isWallFace(this.spatialDisc_.borrow(), face) {
                    if !supportDom.contains(elem) then supportDom += elem;
                    uxCoeff[elem] -= wx;
                    uyCoeff[elem] -= wy;
                }
            }

            const start = this.cellVelTemplateOffsets[elem];
            const dom = {0..<supportDom.size};
            var cols: [dom] int;
            var idx = 0;
            for col in supportDom {
                cols[idx] = col;
                idx += 1;
            }
            sort(cols);
            for localIdx in dom {
                const col = cols[localIdx];
                const dataIdx = start + localIdx;
                this.cellVelTemplateCols[dataIdx] = col;
                this.cellVelTemplateUx[dataIdx] = uxCoeff[col];
                this.cellVelTemplateUy[dataIdx] = uyCoeff[col];
            }
        }
    }

    proc solve() {
        if CHECK_AD_ROW && AD_ROW_CHECK_ONLY {
            writeln("AD row check requested; skipping nonlinear solve.");
            return;
        }

        var normalized_res: real(64) = 1e12;
        var res : real(64) = 1e12;        // Current residual (absolute)
        var res_prev : real(64) = 1e12;  // Previous iteration residual for line search
        var omega : real(64) = this.inputs_.OMEGA_;  // Current relaxation factor
        var time: stopwatch;
        const sufficientDecrease = this.inputs_.SUFFICIENT_DECREASE_;

        // Initial residual
        res = RMSE(this.spatialDisc_.res_, this.spatialDisc_.elemVolume_);
        if this.inputs_.START_FILENAME_ == "" {
            this.first_res_ = res;
        }
        normalized_res = res / this.first_res_;
        const res_wall = RMSE(this.spatialDisc_.res_[this.spatialDisc_.wall_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.wall_dom]);
        const res_farfield = RMSE(this.spatialDisc_.res_[this.spatialDisc_.farfield_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.farfield_dom]);
        const res_fluid = RMSE(this.spatialDisc_.res_[this.spatialDisc_.fluid_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.fluid_dom]);
        const res_wake = RMSE(this.spatialDisc_.res_[this.spatialDisc_.wake_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.wake_dom]);
        const res_shock = RMSE(this.spatialDisc_.res_[this.spatialDisc_.shock_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.shock_dom]);
        const (Cl, Cd, Cm) = this.spatialDisc_.computeAerodynamicCoefficients();
        const elapsed = this.t0_;
        writeln(" Time: ", elapsed, " It: ", this.it_,
                " res: ", res, " norm res: ", normalized_res, " kutta res: ", this.spatialDisc_.kutta_res_,
                " res wall: ", res_wall, " res farfield: ", res_farfield, " res fluid: ", res_fluid, " res wake: ", res_wake, " res shock: ", res_shock,
                " Cl: ", Cl, " Cd: ", Cd, " Cm: ", Cm, " Circulation: ", this.spatialDisc_.circulation_);
        this.timeList_.pushBack(elapsed);
        this.itList_.pushBack(this.it_);
        this.resList_.pushBack(res);
        this.clList_.pushBack(Cl);
        this.cdList_.pushBack(Cd);
        this.cmList_.pushBack(Cm);
        this.circulationList_.pushBack(this.spatialDisc_.circulation_);

        while ((normalized_res > this.inputs_.CONV_TOL_ && res > this.inputs_.CONV_ATOL_) && this.it_ < this.inputs_.IT_MAX_ && isNan(normalized_res) == false) {
            this.it_ += 1;
            time.start();
            const reducedExactMode = isReducedExactJacobianType(this.inputs_.JACOBIAN_TYPE_);
            res_prev = res;
            var jacobianTimer: stopwatch;
            var rhsTimer: stopwatch;
            var linearSolveTimer: stopwatch;
            var lineSearchTimer: stopwatch;
            var postProcessTimer: stopwatch;

            jacobianTimer.start();
            this.computeJacobian();
            jacobianTimer.stop();
            
            // Always use full Newton step in RHS (omega will be applied to the solution update)
            rhsTimer.start();
            forall elem in 1..this.spatialDisc_.nelemDomain_ {
                this.b_petsc.set(elem-1, -this.spatialDisc_.res_[elem]);
            }
            this.b_petsc.assemblyComplete();
            rhsTimer.stop();

            var its: int;
            var reason: int;

            // === PETSC GMRES ===
            linearSolveTimer.start();
            const (petscIts, petscReason) = GMRES(this.ksp, this.A_petsc, this.b_petsc, this.x_petsc);
            linearSolveTimer.stop();
            its = petscIts;
            reason = petscReason;
            
            var lineSearchIts = 0;
            if this.inputs_.OMEGA_RAMP_ {
                const rampEvery = max(1, this.inputs_.OMEGA_RAMP_IT_);
                const rampStage = (this.it_ - 1) / rampEvery;
                var omegaRamp = this.inputs_.OMEGA_START_;
                for 1..rampStage do
                    omegaRamp *= this.inputs_.OMEGA_RAMP_FACTOR_;
                omega = min(omegaRamp, this.inputs_.OMEGA_RAMP_MAX_);
                if abs(omega - this.inputs_.lastOmega_) > 1.0e-14 {
                    writeln("OMEGA ramp increased to ", omega, " at iteration ", this.it_);
                    this.inputs_.lastOmega_ = omega;
                }
            } else if this.inputs_.ADAPTIVE_OMEGA_ {
                if normalized_res <= this.inputs_.OMEGA_THRESHOLD_ {
                    omega = this.inputs_.OMEGA_FINAL_;
                    if !this.inputs_.omegaAdapted_ {
                        writeln("Adaptive OMEGA switched to ", omega,
                                " at normalized residual ", normalized_res);
                        this.inputs_.omegaAdapted_ = true;
                        this.inputs_.lastOmega_ = omega;
                    }
                } else {
                    omega = this.inputs_.OMEGA_START_;
                    this.inputs_.lastOmega_ = omega;
                }
            } else {
                omega = this.inputs_.OMEGA_;
                this.inputs_.lastOmega_ = omega;
            }

            const phiDom = {1..this.spatialDisc_.nelemDomain_};
            var phiSaved: [phiDom] real(64);
            var phiBest: [phiDom] real(64);
            forall elem in phiDom {
                phiSaved[elem] = this.spatialDisc_.phi_[elem];
                phiBest[elem] = this.spatialDisc_.phi_[elem];
            }
            const gammaSaved = this.spatialDisc_.circulation_;
            var gammaBest = gammaSaved;
            var bestRes = max(real(64));
            var accepted = false;
            const maxLineSearch = if this.inputs_.LINE_SEARCH_ then max(0, this.inputs_.MAX_LINE_SEARCH_) else 0;

            lineSearchTimer.start();
            for ls in 0..maxLineSearch {
                lineSearchIts = ls;
                if ls > 0 then
                    omega *= 0.5;

                forall elem in phiDom do
                    this.spatialDisc_.phi_[elem] = phiSaved[elem] + omega * this.x_petsc.get(elem-1);
                this.spatialDisc_.circulation_ = gammaSaved;

                if reducedExactMode {
                    this.enforceConsistentGammaForCurrentPhi();
                } else {
                    const phi_upper = this.spatialDisc_.phi_[this.spatialDisc_.upperTEelem_] +
                                      (this.spatialDisc_.uu_[this.spatialDisc_.upperTEelem_] * this.spatialDisc_.deltaSupperTEx_ +
                                       this.spatialDisc_.vv_[this.spatialDisc_.upperTEelem_] * this.spatialDisc_.deltaSupperTEy_);
                    const phi_lower = this.spatialDisc_.phi_[this.spatialDisc_.lowerTEelem_] +
                                      (this.spatialDisc_.uu_[this.spatialDisc_.lowerTEelem_] * this.spatialDisc_.deltaSlowerTEx_ +
                                       this.spatialDisc_.vv_[this.spatialDisc_.lowerTEelem_] * this.spatialDisc_.deltaSlowerTEy_);
                    const gamma_computed = phi_upper - phi_lower;
                    this.spatialDisc_.circulation_ = gamma_computed;
                }

                this.spatialDisc_.run();
                const trialRes = RMSE(this.spatialDisc_.res_, this.spatialDisc_.elemVolume_);

                if !isNan(trialRes) && trialRes < bestRes {
                    bestRes = trialRes;
                    gammaBest = this.spatialDisc_.circulation_;
                    forall elem in phiDom do
                        phiBest[elem] = this.spatialDisc_.phi_[elem];
                }

                const acceptTrial = !this.inputs_.LINE_SEARCH_ ||
                                    (!isNan(trialRes) &&
                                     trialRes <= sufficientDecrease * res_prev);
                if acceptTrial {
                    res = trialRes;
                    accepted = true;
                    break;
                }
            }

            if !accepted {
                if bestRes < max(real(64)) {
                    forall elem in phiDom do
                        this.spatialDisc_.phi_[elem] = phiBest[elem];
                    this.spatialDisc_.circulation_ = gammaBest;
                    this.spatialDisc_.run();
                    res = bestRes;
                } else {
                    forall elem in phiDom do
                        this.spatialDisc_.phi_[elem] = phiSaved[elem];
                    this.spatialDisc_.circulation_ = gammaSaved;
                    this.spatialDisc_.run();
                    res = res_prev;
                    omega = 0.0;
                }
            }
            lineSearchTimer.stop();

            postProcessTimer.start();
            normalized_res = res / this.first_res_;

            const res_wall = RMSE(this.spatialDisc_.res_[this.spatialDisc_.wall_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.wall_dom]);
            const res_farfield = RMSE(this.spatialDisc_.res_[this.spatialDisc_.farfield_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.farfield_dom]);
            const res_fluid = RMSE(this.spatialDisc_.res_[this.spatialDisc_.fluid_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.fluid_dom]);
            const res_wake = RMSE(this.spatialDisc_.res_[this.spatialDisc_.wake_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.wake_dom]);
            const res_shock = RMSE(this.spatialDisc_.res_[this.spatialDisc_.shock_dom], this.spatialDisc_.elemVolume_[this.spatialDisc_.shock_dom]);

            const (Cl, Cd, Cm) = this.spatialDisc_.computeAerodynamicCoefficients();
            time.stop();
            postProcessTimer.stop();
            const iterTime = time.elapsed();
            const elapsed = iterTime + this.t0_;

            writeln(" Time: ", elapsed, " It: ", this.it_,
                    " res: ", res, " norm res: ", normalized_res, " kutta res: ", this.spatialDisc_.kutta_res_,
                    " res wall: ", res_wall, " res farfield: ", res_farfield, " res fluid: ", res_fluid, " res wake: ", res_wake, " res shock: ", res_shock,
                    " Cl: ", Cl, " Cd: ", Cd, " Cm: ", Cm, " Circulation: ", this.spatialDisc_.circulation_,
                    " GMRES its: ", its, " reason: ", reason, " omega: ", omega, " ls its: ", lineSearchIts);
            if this.inputs_.PROFILE_ITERATION_TIMINGS_ {
                writeln("    Timing dt: ", iterTime,
                        " jac: ", jacobianTimer.elapsed(),
                        " rhs: ", rhsTimer.elapsed(),
                        " linear: ", linearSolveTimer.elapsed(),
                        " ls/update: ", lineSearchTimer.elapsed(),
                        " post: ", postProcessTimer.elapsed());
            }

            this.timeList_.pushBack(elapsed);
            this.itList_.pushBack(this.it_);
            this.resList_.pushBack(res);
            this.clList_.pushBack(Cl);
            this.cdList_.pushBack(Cd);
            this.cmList_.pushBack(Cm);
            this.circulationList_.pushBack(this.spatialDisc_.circulation_);
            
            if this.it_ % this.inputs_.CGNS_OUTPUT_FREQ_ == 0 {
                this.spatialDisc_.writeSolution(this.timeList_, this.itList_, this.resList_, this.clList_, this.cdList_, this.cmList_, this.circulationList_);
            }
        }

        this.spatialDisc_.writeSolution(this.timeList_, this.itList_, this.resList_, this.clList_, this.cdList_, this.cmList_, this.circulationList_);

    }
}

}
