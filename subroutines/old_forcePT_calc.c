/* COMPUTATION OF FORCES BETWEEN NON-CONNECTED BEADS AND BEADS-SUBSTRATE_ATOMS(with neighbor-list) */
void ComputeForcesPT () {
        VecR dr, dv;            // distance and velocity difference between i- and j- beads (VECTORS)
        real scal, fDispRand, fDisp, fRand, theta, uVal, omega_r, omega_d;      // dissipative and random forces variables
        real fcVal, rr, rri, rri3;
        int j1, j2, n;
        int N0, j, N_j1, N_j2, j_min, j_max, top;       // slabs profile variable
        real dTop, dBot, xx, yy, zz, zi;                // distances for pressure tensor variables (dTop - upper particle, dBot - bottom particle) 
        real secTermX, secTermY, secTermZ, secTermXY;   // pressure tensor calculating

	VecR invWid, rs;
	VecI cc;
	int c, hSize;
	real uValPP, uValPS;	// polymer-polymer and polymer-surface energies

        DO_MOL V_ZERO (mol[n].ra);
        uSum = 0.;
        U_LJ_sub = 0.;
        virSum = 0.;
	
	uValPP = uValPS = 0.;

	V_ZERO (force1);
	V_ZERO (force2);
	V_ZERO (forceT);

	/* set all substrate Unit - polymer liquid interactions to zero */
	for (n = -1; n < totSUnitNum; n++) {
		enLJSUnit[n] = 0.;
	}
	
	if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepConVal == 0) {
		hSize = V_PROD (sizeVelGrid);		// number of grid cells in the system

		/* setting initial snapshots values to zero */
		for (j = 0; j < NHIST_E; j++) {
			for (n = 0; n < hSize; n++) snapEnGrid[j][n] = 0.;
		}
	}
	
        for (n = 0; n < nebrTabLen; n++) {
                j1 = nebrTab[2 * n];
                j2 = nebrTab[2 * n + 1];
                if ((mol[j1].inChain == -1) && (mol[j2].inChain == -1)) {
                }
                else {
                        V_SUB (dr, mol[j1].r, mol[j2].r);
                        V_WRAP_ALL (dr);
                        rr = V_LEN_SQ (dr);
                        if (rr < rrCut) {
                        /* defining slab contains j1 */
                                if (mol[j1].r.z < 0.)
                                        N_j1 = halfSl - (int) (-mol[j1].r.z / hSlab) - 1;
                                else
                                        N_j1 = halfSl + (int) (mol[j1].r.z / hSlab);
                        /* defining slab contains j2 */
                                if (mol[j2].r.z < 0.)
                                        N_j2 = halfSl - (int) (-mol[j2].r.z / hSlab) - 1;
                                else
                                        N_j2 = halfSl + (int) (mol[j2].r.z / hSlab);
                        /* x, y, z-distances between two beads (squared) */
                                xx = SQR(dr.x);
                                yy = SQR(dr.y);
                                zz = SQR(dr.z);
                                zi = 1. / sqrt (zz);
                        /* min and max slabs, containing interacting particles */
                                if (N_j1 > N_j2) {
                                        j_max = N_j1;
                                        j_min = N_j2;
                                        top = 1;
                                        /* defining distances to desired slab boundaries */
                                        if (j_max - j_min < halfSl) {
                                                dTop = mol[j1].r.z - hSlab * (j_max - halfSl);
                                                dBot = hSlab * (j_min - halfSl + 1) - mol[j2].r.z;
                                        } else {
                                                dTop = hSlab * (j_max - halfSl + 1) - mol[j1].r.z;
                                                dBot = mol[j2].r.z + hSlab * (halfSl - j_min);
                                        }
                                } else if (N_j1 < N_j2) {
                                        j_max = N_j2;
                                        j_min = N_j1;
                                        top = 2;
                                        /* defining distances to desired slab boundaries */
                                        if (j_max - j_min < halfSl) {
                                                dTop = mol[j2].r.z - hSlab * (j_max - halfSl);
                                                dBot = hSlab * (j_min - halfSl + 1) - mol[j1].r.z;
                                        } else {
                                                dTop = hSlab * (j_max - halfSl + 1) - mol[j2].r.z;
                                                dBot = mol[j1].r.z + hSlab * (halfSl - j_min);
                                        }
                                } else if (N_j1 == N_j2) {
                                        j_max = N_j1;
                                        j_min = N_j1;
                                        top = 0;
                                }

                                /* Lennard-Jones forces */
                                rri = 1. / rr;
                                rri3 = CUBE (rri);
                                if ((mol[j1].inChain != -1) && (mol[j2].inChain != -1)) {       // interaction of nonbonded chain particles
                                        fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                                        uVal = 4. * rri3 * (rri3 - 1.) - U0;
					uValPP = uVal;
					uValPS = 0.;
                                        V_V_S_ADD (mol[j1].ra, fcVal, dr);
                                        V_V_S_ADD (mol[j2].ra, -fcVal, dr);
                                } else if ((mol[j1].inChain == -1) && (mol[j1].inSUnit != -2) && (mol[j2].inChain != -1)) {
                                        fcVal = 48. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 0.5) * rri; // interaction with the substrate
                                        uVal = 4. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 1.) - U0_surf;
                                        U_LJ_sub += uVal;
					uValPP = 0.;
					uValPS = uVal;
					if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepSUEn == 0) {
						enLJSUnit[mol[j1].inSUnit] += uVal;
					}
                                        V_V_S_ADD (mol[j2].ra, -fcVal, dr);
					// controlling forces acting on substrate atoms
					if (mol[j1].r.z > centUCinZ[0]) {
						V_V_S_ADD (force1, fcVal, dr);
					} else {
						V_V_S_ADD (force2, fcVal, dr);
					}
					V_V_S_ADD (forceT, fcVal, dr);
				} else if ((mol[j1].inChain == -1) && (mol[j1].inSUnit == -2) && (mol[j2].inChain != -1)) {
					if (problem == 0) {
						fcVal = 4. * 48. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur * rri; // rep interaction with flatTopW
						uVal = 4. * 4. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur;
					} else {
						fcVal = 48. * epsTopW * rri3 * ss3_sur * (rri3 * ss3_sur - 0.5) * rri; // rep-attr interaction with flatTopW
						uVal = 4. * epsTopW * rri3 * ss3_sur * (rri3 * ss3_sur - 1.);
					}
					U_LJ_sub += uVal;
					uValPP = 0.;
					uValPS = uVal;
					if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepSUEn == 0) {
						enLJSUnit[mol[j1].inSUnit] += uVal;
					}
					V_V_S_ADD (mol[j2].ra, -fcVal, dr);
				} else if ((mol[j1].inChain != -1) && (mol[j2].inChain == -1) && (mol[j2].inSUnit != -2)) {
                                        fcVal = 48. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 0.5) * rri; // interaction with the substrate
                                        uVal = 4. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 1.) - U0_surf;
                                        U_LJ_sub += uVal;
					uValPP = 0.;
					uValPS = uVal;
					if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepSUEn == 0) {
						enLJSUnit[mol[j2].inSUnit] += uVal;
					}
                                        V_V_S_ADD (mol[j1].ra, fcVal, dr);
					// controlling forces acting on substrate atoms
					if (mol[j2].r.z > centUCinZ[0]) {
						V_V_S_ADD (force1, -fcVal, dr);
					} else {
						V_V_S_ADD (force2, -fcVal, dr);
					}
					V_V_S_ADD (forceT, -fcVal, dr);
                                } else if ((mol[j1].inChain != -1) && (mol[j2].inChain == -1) && (mol[j2].inSUnit == -2)) {
					if (problem == 0) {
						fcVal = 4. * 48. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur * rri; // rep interaction with flatTopW
						uVal = 4. * 4. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur;
					} else {
						fcVal = 48. * epsTopW * rri3 * ss3_sur * (rri3 * ss3_sur - 0.5) * rri; // rep-attr interaction with flatTopW
						uVal = 4. * epsTopW * rri3 * ss3_sur * (rri3 * ss3_sur - 1.);
					}
					U_LJ_sub += uVal;
					uValPP = 0.;
					uValPS = uVal;
					if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepSUEn == 0) {
						enLJSUnit[mol[j2].inSUnit] += uVal;
					}
					V_V_S_ADD (mol[j1].ra, fcVal, dr);
				} else {
				}

				if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepConVal == 0) {
					V_DIV (invWid, sizeVelGrid, region);
		
					/* sorting interaction energy for j1 */
					V_S_ADD (rs, mol[j1].r, 0.5, region);	//shifts coordinates to "positive" region (from 0 to L)
					V_MUL (cc, rs, invWid);
					c = V_LINEAR (cc, sizeVelGrid);
					snapEnGrid[1][c] += 0.5 * uValPP;
					snapEnGrid[2][c] += 0.5 * uValPS;
//					snapEnGrid[3][c] += 0.5 * produced by thermostat
					
					/* sorting interaction energy for j2 */
					V_S_ADD (rs, mol[j2].r, 0.5, region);	//shifts coordinates to "positive" region (from 0 to L)
					V_MUL (cc, rs, invWid);
					c = V_LINEAR (cc, sizeVelGrid);
					snapEnGrid[1][c] += 0.5 * uValPP;
					snapEnGrid[2][c] += 0.5 * uValPS;
//					snapEnGrid[3][c] += 0.5 * produced by thermostat
				}

                                uSum += uVal;
                                virSum += fcVal * rr;

                                /* calculating pressure tensor components */
                                if ((stepCount > stepEquil) && ((stepCount - stepEquil) % stepConVal == 0)) {
                                        secTermX = xx * fcVal;
                                        secTermY = yy * fcVal;
                                        secTermZ = zz * fcVal;
                                        secTermXY = secTermX + secTermY;

					/* interacting particles are in the same slab */
                                        if (top == 0) {
                                                histPx[j_max] += secTermX;
                                                histPy[j_max] += secTermY;

						// j1 particle is the substrate atom
                                                if (mol[j1].inChain == -1) { 		// if mol[j1] is an atom of the substrate
							histPnW[j_max] += secTermZ;
//                                                      histPnW[j_max] -= secTermZ * mol[j2].r.z / dr.z;
							histPtW[j_max] += secTermX * mol[j2].r.x / dr.x + secTermY * mol[j2].r.y / dr.y;
                                                } 
						// j2 particle is the substrate atom
						else if (mol[j2].inChain == -1) {	// if mol[j2] is an atom of the substrate
							histPnW[j_max] += secTermZ;
//							histPnW[j_max] += secTermZ * mol[j1].r.z / dr.z;
							histPtW[j_max] += secTermX * mol[j1].r.x / dr.x + secTermY * mol[j1].r.y / dr.y;
                                                }
						// neither j1 nor j2 are substrate atoms
						 else {				// the interaction happens in the polymer solution (not with the substrate atom)
                                                        histPn[j_max] += secTermZ;
                                                        histPt[j_max] += secTermXY;
                                                }
                                        }
					/* interacting particles are in different slabs */
					else {
						// treating slab with maximal number
						histPx[j_max] += secTermX * dTop * zi;
						histPy[j_max] += secTermY * dTop * zi;
						if (mol[j1].inChain == -1) {
							histPnW[j_max] += secTermZ * dTop * zi;
//							histPnW[j_max] += (topSubAtZ + .5 * region.z) * secTermZ * dTop * zi / dr.z;		// corr number 1
							histPtW[j_max] += (secTermX * mol[j2].r.x / dr.x + secTermY * mol[j2].r.y / dr.y) * dTop * zi;
							if (mol[j1].r.z < centUCinZ[0]) {
//								histPnW[j_max] -= secTermZ * zi * zi * dTop * 0.5 * gapSurf.z;	// corr number 2
							}
						} else if (mol[j2].inChain == -1) {
							histPnW[j_max] += secTermZ * dTop * zi;
//							histPnW[j_max] -= (topSubAtZ + .5 * region.z) * secTermZ * dTop * zi / dr.z;		// corr number 1
							histPtW[j_max] -= (secTermX * mol[j1].r.x / dr.x + secTermY * mol[j1].r.y / dr.y) * dTop * zi;
							if (mol[j2].r.z < centUCinZ[0]) {
//								histPnW[j_max] += secTermZ * zi * zi * dTop * 0.5 * gapSurf.z;	// corr number 2
							}
						} else {
							histPn[j_max] += secTermZ * dTop * zi;
							histPt[j_max] += secTermXY * dTop * zi;
						}

						// treating slab with minimal number
						histPx[j_min] += secTermX * dBot * zi;
						histPy[j_min] += secTermY * dBot * zi;
						if (mol[j1].inChain == -1) {
							histPnW[j_min] += secTermZ * dBot * zi;
//							histPnW[j_min] += (topSubAtZ + .5 * region.z) * secTermZ * dBot * zi / dr.z;	// corr number 1
							histPtW[j_min] += (secTermX * mol[j2].r.x / dr.x + secTermY * mol[j2].r.y / dr.y) * dBot * zi;
							if (mol[j1].r.z < centUCinZ[0]) {
								histPnW[j_min] -= secTermZ * zi * dBot;	// corr number 2
							}
						} else if (mol[j2].inChain == -1) {
							histPnW[j_min] += secTermZ * dBot * zi;
//							histPnW[j_min] -= (topSubAtZ + .5 * region.z) * secTermZ * dTop * zi / dr.z;	// corr number 1
							histPtW[j_min] -= (secTermX * mol[j1].r.x / dr.x + secTermY * mol[j1].r.y / dr.y) * dBot * zi;
							if (mol[j2].r.z < centUCinZ[0]) {
								histPnW[j_min] -= secTermZ * zi * dBot;	// corr number 2
							}
						} else {
							histPn[j_min] += secTermZ * dBot * zi;
							histPt[j_min] += secTermXY * dBot * zi;
						}
						
						// treating slabs between min and max numbers
						if (((j_max - j_min) < halfSl) && ((j_max - j_min) > 1)) {
							for (j = j_min + 1; j < j_max; j++) {
								histPx[j] += secTermX * hSlab * zi;
								histPy[j] += secTermY * hSlab * zi;
								if (mol[j1].inChain == -1) {
									histPnW[j] += secTermZ * hSlab * zi;
//									histPnW[j] += (topSubAtZ + .5 * region.z) * secTermZ * hSlab * zi / dr.z;	// corr number 1
									histPtW[j] += (secTermX * mol[j2].r.x / dr.x + secTermY * mol[j2].r.y / dr.y) * hSlab * zi;
									if ((mol[j1].r.z < centUCinZ[0]) && (j < N_high_surf)) {
										histPnW[j] -= secTermZ * zi * hSlab;		// corr number 2
									} else if ((mol[j1].r.z < centUCinZ[0]) && (j == N_high_surf)) {
										histPnW[j] -= secTermZ * zi * dTop_high_surf;	// corr number 2
									} else {
									}
								} else if (mol[j2].inChain == -1) {
									histPnW[j] += secTermZ * hSlab * zi;
//									histPnW[j] -= (topSubAtZ + .5 * region.z) * secTermZ * hSlab * zi / dr.z;	// corr number 1
									histPtW[j] -= (secTermX * mol[j1].r.x / dr.x + secTermY * mol[j1].r.y / dr.y) * hSlab * zi;
									if ((mol[j2].r.z < centUCinZ[0]) && (j < N_high_surf)) {
										histPnW[j] -= secTermZ * zi * hSlab;		// corr number 2
									} else if ((mol[j2].r.z < centUCinZ[0]) && (j == N_high_surf)) {
										histPnW[j] -= secTermZ * zi * dTop_high_surf;	// corr number 2
									} else {
									}
								} else {
									histPn[j] += secTermZ * hSlab * zi;
									histPt[j] += secTermXY * hSlab * zi;
								}
							}
						}
                                                // when the box has a periodicity in z-direction
/*                                              if (((j_max - j_min) > halfSl) && ((j_max - j_min) < (numSlabs - 1))) {
                                                        if (j_min > 0) {        
                                                                for (j = 0; j < j_min; j++) {
                                                                        histPt[j] += secTermXY * hSlab * zi;
                                                                        histPx[j] += secTermX * hSlab * zi;
                                                                        histPy[j] += secTermY * hSlab * zi;
                                                                        histPn[j] += secTermZ * hSlab * zi;
                                                                }
                                                        }
                                                        if (j_max < numSlabs - 1) {
                                                                for (j = j_max + 1; j < numSlabs; j++) {
                                                                        histPt[j] += secTermXY * hSlab * zi;
                                                                        histPx[j] += secTermX * hSlab * zi;
                                                                        histPy[j] += secTermY * hSlab * zi;
                                                                        histPn[j] += secTermZ * hSlab * zi;
                                                                }
                                                        }
                                                }*/
                                        }
                                }
                                /* dissipative forces */
                                omega_r = sqrt(rri) - 1. / rCut;                // omega_r / abs (Rij)
                                V_SUB (dv, mol[j1].rv, mol[j2].rv);
                                scal = V_DOT (dr, dv);
                                /* random forces */
                                zeta = sqrt (12. * 2. * temperature * gamma_d / deltaT);
                                theta = RandR () - 0.5;
                                fDispRand = (zeta * theta - gamma_d * scal * omega_r) * omega_r;
                                if ((mol[j1].inChain != -1) && (mol[j2].inChain != -1)) {
                                        V_V_S_ADD (mol[j1].ra, fDispRand, dr);
                                        V_V_S_ADD (mol[j2].ra, -fDispRand, dr);
                                } else {
//                                      if (mol[j1].inChain == -1) {
//                                              V_V_S_ADD (mol[j2].ra, -fDispRand, dr);
//                                      } else {
//                                              V_V_S_ADD (mol[j1].ra, fDispRand, dr);
//                                      }
                                }
                        }
                }
        }
}

