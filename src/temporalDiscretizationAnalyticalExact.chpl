module temporalDiscretizationAnalyticalExact {
use temporalDiscretization;

// Hand-coded reduced exact Jacobian assembly.

inline proc addActiveLocalIndex(idx: int,
                                ref activeMask: [?D] bool,
                                ref activeList: [D] int,
                                ref activeCount: int) {
    if idx >= 0 && !activeMask[idx] {
        activeMask[idx] = true;
        activeList[activeCount] = idx;
        activeCount += 1;
    }
}

proc temporalDiscretization.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem: int,
                                                                                   stencil: [?D] int,
                                                                                   ref du: [D] real(64),
                                                                                   ref dv: [D] real(64),
                                                                                   ref activeMask: [D] bool,
                                                                                   ref activeList: [D] int,
                                                                                   ref activeCount: int,
                                                                                   ref duGamma: real(64),
                                                                                   ref dvGamma: real(64)) {
    const start = this.cellVelTemplateOffsets[elem];
    const stop = this.cellVelTemplateOffsets[elem + 1];

    for dataIdx in start..<stop {
        const col = this.cellVelTemplateCols[dataIdx];
        const localIdx = findStencilIndex(stencil, col);
        if localIdx >= 0 {
            du[localIdx] += this.cellVelTemplateUx[dataIdx];
            dv[localIdx] += this.cellVelTemplateUy[dataIdx];
            addActiveLocalIndex(localIdx, activeMask, activeList, activeCount);
        }
    }

    duGamma += this.cellVelTemplateGammaUx[elem];
    dvGamma += this.cellVelTemplateGammaUy[elem];
}

proc accumulateDensityMuSensitivityOnActiveIndices(sd: borrowed spatialDiscretization,
                                                   u: real(64),
                                                   v: real(64),
                                                   const ref du: [?D] real(64),
                                                   const ref dv: [D] real(64),
                                                   const ref activeList: [D] int,
                                                   activeCount: int,
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

    for activePos in 0..<activeCount {
        const idx = activeList[activePos];
        drho[idx] = rhoFactor * (u * du[idx] + v * dv[idx]);
    }
    drhoGamma = rhoFactor * (u * duGamma + v * dvGamma);

    const vel2 = u * u + v * v;
    const rhoPow = rho ** (1.0 - sd.inputs_.GAMMA_);
    const M2 = machInf2 * vel2 * rhoPow;
    const Mc2 = sd.inputs_.MACH_C_ * sd.inputs_.MACH_C_;
    if M2 > Mc2 {
        mu = sd.inputs_.MU_C_ * (M2 - Mc2);
        const common = machInf2;
        for activePos in 0..<activeCount {
            const idx = activeList[activePos];
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
        dmuGamma = 0.0;
    }
}

proc temporalDiscretization.computeAnalyticalReducedExactJacobian() {
    this.A_petsc.zeroEntries();

    const sd = this.spatialDisc_.borrow();
    const kuttaStencil = this.buildKuttaStencil();
    const kuttaDom = {0..<kuttaStencil.size};
    var kuttaDphi: [kuttaDom] real(64) = 0.0;
    var kuttaDgamma = 0.0;
    computeAnalyticalKuttaDerivativesOnStencil(sd, kuttaStencil, kuttaDphi, kuttaDgamma);

    if abs(kuttaDgamma) < 1.0e-14 then
        halt("Analytical reduced exact Jacobian failed: dR_gamma/dGamma is too small");

    const rows = 1..this.spatialDisc_.nelemDomain_;
    var mergedVals: [this.rowMergedDataDom] real(64) = 0.0;

    forall row in rows with (ref mergedVals) {
        const rowStencil = this.buildRowStencil(row);
        const rowDom = {0..<rowStencil.size};
        const rowFaceStart = this.rowFaceOffsets[row];
        const rowFaceStop = this.rowFaceOffsets[row + 1];
        var rowDphi: [rowDom] real(64) = 0.0;
        var rowDgamma = 0.0;
        var du1: [rowDom] real(64) = 0.0;
        var dv1: [rowDom] real(64) = 0.0;
        var du2: [rowDom] real(64) = 0.0;
        var dv2: [rowDom] real(64) = 0.0;
        var duAvg: [rowDom] real(64) = 0.0;
        var dvAvg: [rowDom] real(64) = 0.0;
        var dDeltaPhi: [rowDom] real(64) = 0.0;
        var duFaceCoeff: [rowDom] real(64) = 0.0;
        var dvFaceCoeff: [rowDom] real(64) = 0.0;
        var drhoIsen: [rowDom] real(64) = 0.0;
        var drho1Blend: [rowDom] real(64) = 0.0;
        var dmu1Blend: [rowDom] real(64) = 0.0;
        var drho2Blend: [rowDom] real(64) = 0.0;
        var dmu2Blend: [rowDom] real(64) = 0.0;
        var drhoFace: [rowDom] real(64) = 0.0;
        var duTmp: [rowDom] real(64) = 0.0;
        var dvTmp: [rowDom] real(64) = 0.0;
        var active1Mask: [rowDom] bool = false;
        var active2Mask: [rowDom] bool = false;
        var activeTmpMask: [rowDom] bool = false;
        var faceActiveMask: [rowDom] bool = false;
        var active1List: [rowDom] int = 0;
        var active2List: [rowDom] int = 0;
        var activeTmpList: [rowDom] int = 0;
        var faceActiveList: [rowDom] int = 0;
        var active1Count = 0;
        var active2Count = 0;
        var activeTmpCount = 0;
        var faceActiveCount = 0;

        for rowFaceIdx in rowFaceStart..<rowFaceStop {
            for activePos in 0..<active1Count {
                const idx = active1List[activePos];
                du1[idx] = 0.0;
                dv1[idx] = 0.0;
                active1Mask[idx] = false;
            }
            active1Count = 0;

            for activePos in 0..<active2Count {
                const idx = active2List[activePos];
                du2[idx] = 0.0;
                dv2[idx] = 0.0;
                active2Mask[idx] = false;
            }
            active2Count = 0;

            for activePos in 0..<activeTmpCount {
                const idx = activeTmpList[activePos];
                duTmp[idx] = 0.0;
                dvTmp[idx] = 0.0;
                activeTmpMask[idx] = false;
            }
            activeTmpCount = 0;

            for activePos in 0..<faceActiveCount {
                const idx = faceActiveList[activePos];
                duAvg[idx] = 0.0;
                dvAvg[idx] = 0.0;
                dDeltaPhi[idx] = 0.0;
                duFaceCoeff[idx] = 0.0;
                dvFaceCoeff[idx] = 0.0;
                drhoIsen[idx] = 0.0;
                drho1Blend[idx] = 0.0;
                dmu1Blend[idx] = 0.0;
                drho2Blend[idx] = 0.0;
                dmu2Blend[idx] = 0.0;
                drhoFace[idx] = 0.0;
                faceActiveMask[idx] = false;
            }
            faceActiveCount = 0;

            const face = this.rowFaceData[rowFaceIdx];
            const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
            const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
            const sign = if elem1 == row then 1.0 else -1.0;
            const nx = this.spatialDisc_.faceNormalX_[face];
            const ny = this.spatialDisc_.faceNormalY_[face];
            const area = this.spatialDisc_.faceArea_[face];
            const qFace = this.spatialDisc_.uFace_[face] * nx + this.spatialDisc_.vFace_[face] * ny;
            const rhoFace = this.spatialDisc_.rhoFace_[face];
            const rhoIsen = this.spatialDisc_.rhoIsenFace_[face];
            const uFace = this.spatialDisc_.uFace_[face];
            const vFace = this.spatialDisc_.vFace_[face];
            var du1Gamma = 0.0;
            var dv1Gamma = 0.0;
            var du2Gamma = 0.0;
            var dv2Gamma = 0.0;

            if elem1 <= sd.nelemDomain_ {
                this.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem1, rowStencil,
                                                                            du1, dv1,
                                                                            active1Mask, active1List, active1Count,
                                                                            du1Gamma, dv1Gamma);
            } else if isWallFace(sd, face) {
                this.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem2, rowStencil,
                                                                            du1, dv1,
                                                                            active1Mask, active1List, active1Count,
                                                                            du1Gamma, dv1Gamma);
                for activePos in 0..<active1Count {
                    const idx = active1List[activePos];
                    const dvDotN = du1[idx] * nx + dv1[idx] * ny;
                    du1[idx] = du1[idx] - 2.0 * dvDotN * nx;
                    dv1[idx] = dv1[idx] - 2.0 * dvDotN * ny;
                }
                const dvDotNGamma = du1Gamma * nx + dv1Gamma * ny;
                du1Gamma = du1Gamma - 2.0 * dvDotNGamma * nx;
                dv1Gamma = dv1Gamma - 2.0 * dvDotNGamma * ny;
            }

            if elem2 <= sd.nelemDomain_ {
                this.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem2, rowStencil,
                                                                            du2, dv2,
                                                                            active2Mask, active2List, active2Count,
                                                                            du2Gamma, dv2Gamma);
            } else if isWallFace(sd, face) {
                this.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem1, rowStencil,
                                                                            du2, dv2,
                                                                            active2Mask, active2List, active2Count,
                                                                            du2Gamma, dv2Gamma);
                for activePos in 0..<active2Count {
                    const idx = active2List[activePos];
                    const dvDotN = du2[idx] * nx + dv2[idx] * ny;
                    du2[idx] = du2[idx] - 2.0 * dvDotN * nx;
                    dv2[idx] = dv2[idx] - 2.0 * dvDotN * ny;
                }
                const dvDotNGamma = du2Gamma * nx + dv2Gamma * ny;
                du2Gamma = du2Gamma - 2.0 * dvDotNGamma * nx;
                dv2Gamma = dv2Gamma - 2.0 * dvDotNGamma * ny;
            }

            const weight1 = sd.weights1_[face];
            const weight2 = sd.weights2_[face];
            for activePos in 0..<active1Count do
                addActiveLocalIndex(active1List[activePos], faceActiveMask, faceActiveList, faceActiveCount);
            for activePos in 0..<active2Count do
                addActiveLocalIndex(active2List[activePos], faceActiveMask, faceActiveList, faceActiveCount);

            const idxMinus = this.rowFaceMinusIdx[rowFaceIdx];
            const idxPlus = this.rowFacePlusIdx[rowFaceIdx];
            addActiveLocalIndex(idxMinus, faceActiveMask, faceActiveList, faceActiveCount);
            addActiveLocalIndex(idxPlus, faceActiveMask, faceActiveList, faceActiveCount);

            for activePos in 0..<faceActiveCount {
                const idx = faceActiveList[activePos];
                duAvg[idx] = weight1 * du1[idx] + weight2 * du2[idx];
                dvAvg[idx] = weight1 * dv1[idx] + weight2 * dv2[idx];
            }
            const duAvgGamma = weight1 * du1Gamma + weight2 * du2Gamma;
            const dvAvgGamma = weight1 * dv1Gamma + weight2 * dv2Gamma;

            var dDeltaPhiGamma = 0.0;
            if elem1 <= sd.nelemDomain_ && elem2 <= sd.nelemDomain_ {
                const kuttaType1 = sd.kuttaCell_[elem1];
                const kuttaType2 = sd.kuttaCell_[elem2];
                if kuttaType1 == 1 && kuttaType2 == -1 then
                    dDeltaPhiGamma = 1.0;
                else if kuttaType1 == -1 && kuttaType2 == 1 then
                    dDeltaPhiGamma = -1.0;
            }
            if idxMinus >= 0 then dDeltaPhi[idxMinus] -= 1.0;
            if idxPlus >= 0 then dDeltaPhi[idxPlus] += 1.0;

            for activePos in 0..<faceActiveCount {
                const idx = faceActiveList[activePos];
                const ddelta = duAvg[idx] * sd.t_IJ_x_[face] + dvAvg[idx] * sd.t_IJ_y_[face] -
                               dDeltaPhi[idx] * sd.invL_IJ_[face];
                duFaceCoeff[idx] = duAvg[idx] - ddelta * sd.corrCoeffX_[face];
                dvFaceCoeff[idx] = dvAvg[idx] - ddelta * sd.corrCoeffY_[face];
            }
            const ddeltaGamma = duAvgGamma * sd.t_IJ_x_[face] + dvAvgGamma * sd.t_IJ_y_[face] -
                                dDeltaPhiGamma * sd.invL_IJ_[face];
            const duFaceGamma = duAvgGamma - ddeltaGamma * sd.corrCoeffX_[face];
            const dvFaceGamma = dvAvgGamma - ddeltaGamma * sd.corrCoeffY_[face];

            const machInf2 = sd.inputs_.MACH_ * sd.inputs_.MACH_;
            const a = sd.gamma_minus_one_over_two_ * machInf2;
            const BFace = 1.0 + a * (1.0 - uFace * uFace - vFace * vFace);
            const rhoFaceFactor = -2.0 * a * sd.one_over_gamma_minus_one_ *
                                  (BFace ** (sd.one_over_gamma_minus_one_ - 1.0));
            for activePos in 0..<faceActiveCount {
                const idx = faceActiveList[activePos];
                drhoIsen[idx] = rhoFaceFactor * (uFace * duFaceCoeff[idx] + vFace * dvFaceCoeff[idx]);
            }
            const drhoIsenGamma = rhoFaceFactor * (uFace * duFaceGamma + vFace * dvFaceGamma);

            var rho1Blend, mu1Blend, rho2Blend, mu2Blend: real(64);
            var drho1BlendGamma, dmu1BlendGamma, drho2BlendGamma, dmu2BlendGamma: real(64);

            if elem1 <= sd.nelemDomain_ {
                rho1Blend = sd.rhorho_[elem1];
                mu1Blend = sd.mumu_[elem1];
                accumulateDensityMuSensitivityOnActiveIndices(sd, sd.uu_[elem1], sd.vv_[elem1],
                                                              du1, dv1, active1List, active1Count,
                                                              du1Gamma, dv1Gamma,
                                                              drho1Blend, dmu1Blend,
                                                              drho1BlendGamma, dmu1BlendGamma,
                                                              rho1Blend, mu1Blend);
            } else if isWallFace(sd, face) {
                rho1Blend = sd.rhorho_[elem2];
                mu1Blend = sd.mumu_[elem2];
                var duTmpGamma = 0.0;
                var dvTmpGamma = 0.0;
                this.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem2, rowStencil,
                                                                            duTmp, dvTmp,
                                                                            activeTmpMask, activeTmpList, activeTmpCount,
                                                                            duTmpGamma, dvTmpGamma);
                accumulateDensityMuSensitivityOnActiveIndices(sd, sd.uu_[elem2], sd.vv_[elem2],
                                                              duTmp, dvTmp, activeTmpList, activeTmpCount,
                                                              duTmpGamma, dvTmpGamma,
                                                              drho1Blend, dmu1Blend,
                                                              drho1BlendGamma, dmu1BlendGamma,
                                                              rho1Blend, mu1Blend);
                for activePos in 0..<activeTmpCount do
                    addActiveLocalIndex(activeTmpList[activePos], faceActiveMask, faceActiveList, faceActiveCount);
            } else {
                rho1Blend = 0.0;
                mu1Blend = 0.0;
                drho1BlendGamma = 0.0;
                dmu1BlendGamma = 0.0;
            }

            if elem2 <= sd.nelemDomain_ {
                rho2Blend = sd.rhorho_[elem2];
                mu2Blend = sd.mumu_[elem2];
                accumulateDensityMuSensitivityOnActiveIndices(sd, sd.uu_[elem2], sd.vv_[elem2],
                                                              du2, dv2, active2List, active2Count,
                                                              du2Gamma, dv2Gamma,
                                                              drho2Blend, dmu2Blend,
                                                              drho2BlendGamma, dmu2BlendGamma,
                                                              rho2Blend, mu2Blend);
            } else if isWallFace(sd, face) {
                rho2Blend = sd.rhorho_[elem1];
                mu2Blend = sd.mumu_[elem1];
                var duTmpGamma = 0.0;
                var dvTmpGamma = 0.0;
                this.accumulateCellVelocitySensitivityFromTemplateOnStencil(elem1, rowStencil,
                                                                            duTmp, dvTmp,
                                                                            activeTmpMask, activeTmpList, activeTmpCount,
                                                                            duTmpGamma, dvTmpGamma);
                accumulateDensityMuSensitivityOnActiveIndices(sd, sd.uu_[elem1], sd.vv_[elem1],
                                                              duTmp, dvTmp, activeTmpList, activeTmpCount,
                                                              duTmpGamma, dvTmpGamma,
                                                              drho2Blend, dmu2Blend,
                                                              drho2BlendGamma, dmu2BlendGamma,
                                                              rho2Blend, mu2Blend);
                for activePos in 0..<activeTmpCount do
                    addActiveLocalIndex(activeTmpList[activePos], faceActiveMask, faceActiveList, faceActiveCount);
            } else {
                rho2Blend = 0.0;
                mu2Blend = 0.0;
                drho2BlendGamma = 0.0;
                dmu2BlendGamma = 0.0;
            }

            const useElem1Upwind = qFace >= 0.0;
            const rhoUpwind = if useElem1Upwind then rho1Blend else rho2Blend;
            const muUpwind = if useElem1Upwind then mu1Blend else mu2Blend;
            const drhoUpwindGamma = if useElem1Upwind then drho1BlendGamma else drho2BlendGamma;
            const dmuUpwindGamma = if useElem1Upwind then dmu1BlendGamma else dmu2BlendGamma;

            var drhoFaceGamma = drhoIsenGamma;
            if muUpwind > 0.0 {
                for activePos in 0..<faceActiveCount {
                    const idx = faceActiveList[activePos];
                    const drhoUpwind = if useElem1Upwind then drho1Blend[idx] else drho2Blend[idx];
                    const dmuUpwind = if useElem1Upwind then dmu1Blend[idx] else dmu2Blend[idx];
                    drhoFace[idx] = (1.0 - muUpwind) * drhoIsen[idx] + muUpwind * drhoUpwind +
                                    (rhoUpwind - rhoIsen) * dmuUpwind;
                }
                drhoFaceGamma = (1.0 - muUpwind) * drhoIsenGamma + muUpwind * drhoUpwindGamma +
                                (rhoUpwind - rhoIsen) * dmuUpwindGamma;
            } else {
                for activePos in 0..<faceActiveCount {
                    const idx = faceActiveList[activePos];
                    drhoFace[idx] = drhoIsen[idx];
                }
            }

            for activePos in 0..<faceActiveCount {
                const idx = faceActiveList[activePos];
                const dq = duFaceCoeff[idx] * nx + dvFaceCoeff[idx] * ny;
                rowDphi[idx] += sign * area * (drhoFace[idx] * qFace + rhoFace * dq) * sd.res_scale_;
            }
            const dqGamma = duFaceGamma * nx + dvFaceGamma * ny;
            rowDgamma += sign * area * (drhoFaceGamma * qFace + rhoFace * dqGamma) * sd.res_scale_;
        }

        const mergedStart = this.rowMergedOffsets[row];
        const mergedStop = this.rowMergedOffsets[row + 1];
        const gammaFactor = rowDgamma / kuttaDgamma;
        for mergedIdx in mergedStart..<mergedStop {
            const rowIdx = this.rowMergedRowIdx[mergedIdx];
            const kuttaIdx = this.rowMergedKuttaIdx[mergedIdx];
            const rowValue = if rowIdx >= 0 then rowDphi[rowIdx] else 0.0;
            const kuttaValue = if kuttaIdx >= 0 then kuttaDphi[kuttaIdx] else 0.0;
            mergedVals[mergedIdx] = rowValue - gammaFactor * kuttaValue;
        }
    }

    if AD_JACOBIAN_PROGRESS_FREQ > 0 {
        writeln("Analytical reduced exact Jacobian assembly completed in parallel for ",
                this.spatialDisc_.nelemDomain_, " rows");
    }

    for row in rows {
        const mergedStart = this.rowMergedOffsets[row];
        const mergedStop = this.rowMergedOffsets[row + 1];
        for mergedIdx in mergedStart..<mergedStop {
            const col = this.rowMergedCols[mergedIdx];
            const val = mergedVals[mergedIdx];
            if abs(val) > AD_ROW_PRINT_TOL || col == row then
                this.A_petsc.add(row - 1, col - 1, val);
        }
    }

    this.A_petsc.assemblyComplete();
}
}
