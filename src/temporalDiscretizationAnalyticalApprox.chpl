module temporalDiscretizationAnalyticalApprox {
use temporalDiscretization;

// Original approximate analytical Jacobian path.

proc temporalDiscretization.computeGradientSensitivity() {
    this.dgradX_dGamma_ = 0.0;
    this.dgradY_dGamma_ = 0.0;

    forall elem in 1..this.spatialDisc_.nelemDomain_ {
        const kuttaType_elem = this.spatialDisc_.kuttaCell_[elem];

        if kuttaType_elem == 1 || kuttaType_elem == -1 {
            const faces = this.spatialDisc_.mesh_.elem2edge_[
                this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1 ..
                this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1]];

            var dgx = 0.0, dgy = 0.0;

            for face in faces {
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
                const neighbor = if elem1 == elem then elem2 else elem1;

                const kuttaType_neighbor = this.spatialDisc_.kuttaCell_[neighbor];
                const isWakeCrossing = (kuttaType_elem == 1 && kuttaType_neighbor == -1) ||
                                       (kuttaType_elem == -1 && kuttaType_neighbor == 1);

                if isWakeCrossing {
                    var wx, wy: real(64);
                    if elem == elem1 {
                        wx = this.spatialDisc_.lsGradQR_!.wxFinal1_[face];
                        wy = this.spatialDisc_.lsGradQR_!.wyFinal1_[face];
                    } else {
                        wx = this.spatialDisc_.lsGradQR_!.wxFinal2_[face];
                        wy = this.spatialDisc_.lsGradQR_!.wyFinal2_[face];
                    }

                    const gammaSgn = if kuttaType_elem == 1 then 1.0 else -1.0;
                    dgx += wx * gammaSgn;
                    dgy += wy * gammaSgn;
                }
            }

            this.dgradX_dGamma_[elem] = dgx;
            this.dgradY_dGamma_[elem] = dgy;
        }
    }
}

proc temporalDiscretization.computeCirculationSensitivityForRow(elem: int): real(64) {
    var dRes_dGamma = 0.0;
    const faceStart = this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1;
    const faceEnd = this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1];

    for faceIdx in faceStart..faceEnd {
        const face = this.spatialDisc_.mesh_.elem2edge_[faceIdx];
        const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
        const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
        const neighbor = if elem1 == elem then elem2 else elem1;
        const sign = if elem1 == elem then 1.0 else -1.0;

        const nx = this.spatialDisc_.faceNormalX_[face];
        const ny = this.spatialDisc_.faceNormalY_[face];
        const area = this.spatialDisc_.faceArea_[face];
        const rhoFace = this.spatialDisc_.rhoFace_[face];

        const t_x = this.spatialDisc_.t_IJ_x_[face];
        const t_y = this.spatialDisc_.t_IJ_y_[face];
        const invL = this.spatialDisc_.invL_IJ_[face];
        const nDotT = nx * t_x + ny * t_y;
        const k = 1.0 / nDotT;
        const mx = nx - k * t_x;
        const my = ny - k * t_y;
        const directCoeff = k * invL;

        var faceWeightElem: real(64);
        var faceWeightNeighbor: real(64);
        if elem1 == elem {
            faceWeightElem = this.spatialDisc_.weights1_[face];
            faceWeightNeighbor = this.spatialDisc_.weights2_[face];
        } else {
            faceWeightElem = this.spatialDisc_.weights2_[face];
            faceWeightNeighbor = this.spatialDisc_.weights1_[face];
        }

        const isInteriorFace = neighbor <= this.spatialDisc_.nelemDomain_;
        const wallFace = this.spatialDisc_.wallFaceSet_.contains(face);
        const kuttaType_elem = this.spatialDisc_.kuttaCell_[elem];
        const kuttaType_neighbor = this.spatialDisc_.kuttaCell_[neighbor];
        const isWakeCrossingFace = (kuttaType_elem == 1 && kuttaType_neighbor == -1) ||
                                   (kuttaType_elem == -1 && kuttaType_neighbor == 1);

        if isInteriorFace {
            const dgradX_elem = this.dgradX_dGamma_[elem];
            const dgradY_elem = this.dgradY_dGamma_[elem];
            var dFlux_dGamma = faceWeightElem * (dgradX_elem * mx + dgradY_elem * my);

            const dgradX_neighbor = this.dgradX_dGamma_[neighbor];
            const dgradY_neighbor = this.dgradY_dGamma_[neighbor];
            dFlux_dGamma += faceWeightNeighbor * (dgradX_neighbor * mx + dgradY_neighbor * my);
            dFlux_dGamma *= sign * area * rhoFace;

            if isWakeCrossingFace {
                const gammaSgn = if kuttaType_elem == 1 then 1.0 else -1.0;
                dFlux_dGamma += directCoeff * area * rhoFace * gammaSgn;
            }

            dRes_dGamma += dFlux_dGamma;
        } else if wallFace {
            const nDotM = nx * mx + ny * my;
            const mWallX = mx - nDotM * nx;
            const mWallY = my - nDotM * ny;
            const dgradX_elem = this.dgradX_dGamma_[elem];
            const dgradY_elem = this.dgradY_dGamma_[elem];
            dRes_dGamma += (dgradX_elem * mWallX + dgradY_elem * mWallY) * sign * area * rhoFace;
        }
    }

    return dRes_dGamma;
}

proc temporalDiscretization.computeApproximateAnalyticalJacobian() {
    this.A_petsc.zeroEntries();

    forall elem in 1..this.spatialDisc_.nelemDomain_ {
        const faces = this.spatialDisc_.mesh_.elem2edge_[
            this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1 ..
            this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1]];

        var diag = 0.0;
        var dRes_dGamma = 0.0;

        for face in faces {
            const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
            const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
            const neighbor = if elem1 == elem then elem2 else elem1;
            const sign = if elem1 == elem then 1.0 else -1.0;

            const nx = this.spatialDisc_.faceNormalX_[face];
            const ny = this.spatialDisc_.faceNormalY_[face];
            const area = this.spatialDisc_.faceArea_[face];
            const rhoFace = this.spatialDisc_.rhoFace_[face];

            const t_x = this.spatialDisc_.t_IJ_x_[face];
            const t_y = this.spatialDisc_.t_IJ_y_[face];
            const invL = this.spatialDisc_.invL_IJ_[face];
            const nDotT = nx * t_x + ny * t_y;
            const k = 1.0 / nDotT;
            const mx = nx - k * t_x;
            const my = ny - k * t_y;
            const directCoeff = k * invL;

            const sumWx_elem = this.spatialDisc_.lsGradQR_!.sumWx_[elem];
            const sumWy_elem = this.spatialDisc_.lsGradQR_!.sumWy_[elem];
            var wx_elemToNeighbor: real(64);
            var wy_elemToNeighbor: real(64);
            var wx_neighborToElem: real(64);
            var wy_neighborToElem: real(64);
            var faceWeightElem: real(64);
            var faceWeightNeighbor: real(64);

            if elem1 == elem {
                wx_elemToNeighbor = this.spatialDisc_.lsGradQR_!.wxFinal1_[face];
                wy_elemToNeighbor = this.spatialDisc_.lsGradQR_!.wyFinal1_[face];
                wx_neighborToElem = this.spatialDisc_.lsGradQR_!.wxFinal2_[face];
                wy_neighborToElem = this.spatialDisc_.lsGradQR_!.wyFinal2_[face];
                faceWeightElem = this.spatialDisc_.weights1_[face];
                faceWeightNeighbor = this.spatialDisc_.weights2_[face];
            } else {
                wx_elemToNeighbor = this.spatialDisc_.lsGradQR_!.wxFinal2_[face];
                wy_elemToNeighbor = this.spatialDisc_.lsGradQR_!.wyFinal2_[face];
                wx_neighborToElem = this.spatialDisc_.lsGradQR_!.wxFinal1_[face];
                wy_neighborToElem = this.spatialDisc_.lsGradQR_!.wyFinal1_[face];
                faceWeightElem = this.spatialDisc_.weights2_[face];
                faceWeightNeighbor = this.spatialDisc_.weights1_[face];
            }

            const isInteriorFace = neighbor <= this.spatialDisc_.nelemDomain_;
            const wallFace = this.spatialDisc_.wallFaceSet_.contains(face);
            const kuttaType_elem = this.spatialDisc_.kuttaCell_[elem];
            const kuttaType_neighbor = this.spatialDisc_.kuttaCell_[neighbor];
            const isWakeCrossingFace = (kuttaType_elem == 1 && kuttaType_neighbor == -1) ||
                                       (kuttaType_elem == -1 && kuttaType_neighbor == 1);

            if isInteriorFace {
                var face_diag = faceWeightElem * (-sumWx_elem * mx - sumWy_elem * my);
                face_diag += faceWeightNeighbor * (wx_neighborToElem * mx + wy_neighborToElem * my);
                diag += sign * face_diag * area * rhoFace;
                diag += -directCoeff * area * rhoFace;

                const sumWx_neighbor = this.spatialDisc_.lsGradQR_!.sumWx_[neighbor];
                const sumWy_neighbor = this.spatialDisc_.lsGradQR_!.sumWy_[neighbor];

                var offdiag = faceWeightElem * (wx_elemToNeighbor * mx + wy_elemToNeighbor * my);
                offdiag += faceWeightNeighbor * (-sumWx_neighbor * mx - sumWy_neighbor * my);
                offdiag *= sign * area * rhoFace;
                offdiag += directCoeff * area * rhoFace;

                this.A_petsc.add(elem - 1, neighbor - 1, offdiag * this.spatialDisc_.res_scale_);

                const dgradX_elem = this.dgradX_dGamma_[elem];
                const dgradY_elem = this.dgradY_dGamma_[elem];
                var dFlux_dGamma = faceWeightElem * (dgradX_elem * mx + dgradY_elem * my);

                const dgradX_neighbor = this.dgradX_dGamma_[neighbor];
                const dgradY_neighbor = this.dgradY_dGamma_[neighbor];
                dFlux_dGamma += faceWeightNeighbor * (dgradX_neighbor * mx + dgradY_neighbor * my);
                dFlux_dGamma *= sign * area * rhoFace;

                if isWakeCrossingFace {
                    var gammaSgn = 0.0;
                    if kuttaType_elem == 1 && kuttaType_neighbor == -1 {
                        gammaSgn = 1.0;
                    } else {
                        gammaSgn = -1.0;
                    }
                    dFlux_dGamma += directCoeff * area * rhoFace * gammaSgn;
                }

                dRes_dGamma += dFlux_dGamma;
            } else if wallFace {
                const nDotM = nx * mx + ny * my;
                const mWallX = mx - nDotM * nx;
                const mWallY = my - nDotM * ny;

                var face_diag = (-sumWx_elem + wx_elemToNeighbor) * mWallX +
                                (-sumWy_elem + wy_elemToNeighbor) * mWallY;
                diag += sign * face_diag * area * rhoFace;

                const dgradX_elem = this.dgradX_dGamma_[elem];
                const dgradY_elem = this.dgradY_dGamma_[elem];
                const dFlux_dGamma_wall =
                    (dgradX_elem * mWallX + dgradY_elem * mWallY) * sign * area * rhoFace;
                dRes_dGamma += dFlux_dGamma_wall;
            } else {
                diag -= directCoeff * area * rhoFace;
            }
        }

        this.A_petsc.add(elem - 1, elem - 1, diag * this.spatialDisc_.res_scale_);
        if dRes_dGamma != 0.0 {
            const dRes_dPhi_upper = dRes_dGamma;
            const dRes_dPhi_lower = -dRes_dGamma;

            const wakeFaceIndex = this.spatialDisc_.wakeFaceIndexInfluenceOnElem_[elem];
            const upperTE_influences = this.spatialDisc_.wakeFaceUpper_[wakeFaceIndex];
            const lowerTE_influences = this.spatialDisc_.wakeFaceLower_[wakeFaceIndex];

            this.A_petsc.add(elem - 1, upperTE_influences - 1, dRes_dPhi_upper);
            this.A_petsc.add(elem - 1, lowerTE_influences - 1, dRes_dPhi_lower);
        }

        // Store diag for possible use in upwinding
        this.Jij_[elem] = diag;
    }

    
    // === BETA-BASED UPWIND AUGMENTATION (element-centric for parallelization) ===
    // Loop over elements instead of faces to avoid race conditions.
    // Each element checks its faces to see if it's the downwind cell of a supersonic face.
    // Reuses upwindElem_ computed during artificial density calculation.
    forall elem in 1..this.spatialDisc_.nelemDomain_ {
        const faces = this.spatialDisc_.mesh_.elem2edge_[
            this.spatialDisc_.mesh_.elem2edgeIndex_[elem] + 1 ..
            this.spatialDisc_.mesh_.elem2edgeIndex_[elem + 1]];
            
        for face in faces {
            const machFace = this.spatialDisc_.machFace_[face];
            
            if machFace >= this.spatialDisc_.inputs_.MACH_C_ {
                const elem1 = this.spatialDisc_.mesh_.edge2elem_[1, face];
                const elem2 = this.spatialDisc_.mesh_.edge2elem_[2, face];
                
                // Both cells must be interior (not ghost cells)
                if elem1 <= this.spatialDisc_.nelemDomain_ && elem2 <= this.spatialDisc_.nelemDomain_ {
                    // Reuse upwind/downwind elements computed during artificial density
                    const upwindElem = this.spatialDisc_.upwindElem_[face];
                    const downwindElem = this.spatialDisc_.downwindElem_[face];
                    
                    // Only process if this element is the downwind cell
                    // This ensures each matrix entry is only written by one task
                    if downwindElem == elem {
                        // Use precomputed invL_IJ_ (inverse of cell centroid distance)
                        const increase = this.inputs_.BETA_ * this.spatialDisc_.velMagFace_[face] 
                        * this.spatialDisc_.invL_IJ_[face] * this.spatialDisc_.res_scale_;
                        
                        // Increase absolute value of diagonal term for downwind element
                        const diagTerm = this.Jij_[elem];
                        if diagTerm >= 0.0 {
                            // Increase diagonal and decrease off-diagonal
                            this.A_petsc.add(elem-1, elem-1, increase);
                            this.A_petsc.add(elem-1, upwindElem-1, -increase);
                        } else {
                            // Decrease diagonal and increase off-diagonal
                            this.A_petsc.add(elem-1, elem-1, -increase);
                            this.A_petsc.add(elem-1, upwindElem-1, increase);
                        }
                    }
                }
            }
        }
    }

    this.A_petsc.assemblyComplete();
}
}
