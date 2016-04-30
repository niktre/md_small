/* COMPUTATION OF FORCES BETWEEN CONNECTED BEADS IN CHAIN */
void ComputeChainBondForcesPT () {
   VecR dr, dv;            // distance and velocity difference between i- and j- beads (VECTORS)
   real fcVal, fTot, rr, rri, rri3, uVal, rrno;
   int i, j1, j2, n;
   
   real scal, fDispRand, theta, omega_r;    // DPD
   real zeta = sqrt (24. * temperature * gamma_d / deltaT);
   int j, N_j1, N_j2, j_min, j_max, top;       // slabs profile variable
   real dTop, dBot, xx, yy, zz, zi;                // distances for pressure tensor variables (dTop - upper particle, dBot - bottom particle)
   real secTermX, secTermY, secTermZ, secTermXY;   // pressure tensor calculating
   
   fTot = 0.;
   
   for (n = 0; n < nChain; n++) {
      for (i = 0; i < chainLen - 1; i++) {
         j1 = n * chainLen + i;
         j2 = j1 + 1;
         V_SUB (dr, mol[j1].r, mol[j2].r);
         V_WRAP_ALL (dr);
         rr = V_LEN_SQ (dr);
         if (rr < rrCut) {
            /* defining slab contains j1 */
            if (mol[j1].r.z < 0.) {
               N_j1 = halfSl - (int) (-mol[j1].r.z / hSlab) - 1;
            }
            else {
               N_j1 = halfSl + (int) (mol[j1].r.z / hSlab);
            }
            /* defining slab contains j2 */
            if (mol[j2].r.z < 0.) {
               N_j2 = halfSl - (int) (-mol[j2].r.z / hSlab) - 1;
            }
            else {
               N_j2 = halfSl + (int) (mol[j2].r.z / hSlab);
            }
            /* z-distance between two beads (square) */
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
            fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
            uVal = 4. * rri3 * (rri3 - 1.) - U0;
            V_V_S_ADD (mol[j1].ra, fcVal, dr);
            V_V_S_ADD (mol[j2].ra, -fcVal, dr);
            
            uSum += uVal;
            //			if ((stepCount > stepEquil) && ((stepCount - stepEquil) % stepConVal == 0)) {
            fTot = fcVal;
            //			}
            /* dissipative forces */
            omega_r = sqrt(rri) - 1. / rCut;        // %omega_r% divided by %|Rij|%
            V_SUB (dv, mol[j1].rv, mol[j2].rv);
            scal = V_DOT (dr, dv);
            /* random forces */
            theta = RandR () - 0.5;                 // %theta% is in range [-0.5 ; 0.5]
            fDispRand = (zeta * theta - gamma_d * scal * omega_r) * omega_r;
            V_V_S_ADD (mol[j1].ra, fDispRand, dr);
            V_V_S_ADD (mol[j2].ra, -fDispRand, dr);
         } else {
            ErrExit (ERR_BOND_SNAPPED);
         }
         /* FENE forces */
         if (rr < R2) {
            rrno = 1. - rr / R2;
            fcVal = -k / rrno;
            uVal = (-0.5) * k * R2 * log(rrno);
            V_V_S_ADD (mol[j1].ra, fcVal, dr);
            V_V_S_ADD (mol[j2].ra, -fcVal, dr);
            
            fTot += fcVal;
            uSum += uVal;
            virSum += fTot * rr;
            
            /* normal to (x,y) pressure */
            if ((stepCount > stepEquil) && ((stepCount - stepEquil) % stepConVal == 0)) {
               
               secTermX = xx * fTot;
               secTermY = yy * fTot;
               secTermZ = zz * fTot;
               secTermXY = secTermX + secTermY;
               
               if (top == 0) {
                  histPt[j_max] += secTermXY;
                  histPx[j_max] += secTermX;
                  histPy[j_max] += secTermY;
                  histPn[j_max] += secTermZ;
               } else {
                  histPt[j_max] += secTermXY * dTop * zi;
                  histPx[j_max] += secTermX * dTop * zi;
                  histPy[j_max] += secTermY * dTop * zi;
                  histPn[j_max] += secTermZ * dTop * zi;
                  histPt[j_min] += secTermXY * dBot * zi;
                  histPx[j_min] += secTermX * dBot * zi;
                  histPy[j_min] += secTermY * dBot * zi;
                  histPn[j_min] += secTermZ * dBot * zi;
                  
                  if (((j_max - j_min) < halfSl) && ((j_max - j_min) > 1)) {
                     for (j = j_min + 1; j < j_max; j++) {
                        histPt[j] += secTermXY * hSlab * zi;
                        histPx[j] += secTermX * hSlab * zi;
                        histPy[j] += secTermY * hSlab * zi;
                        histPn[j] += secTermZ * hSlab * zi;
                     }
                  }
                  
                  /*                                      if (((j_max - j_min) > halfSl) && ((j_max - j_min) < (numSlabs - 1))) {
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
         }
      }
   }
}

