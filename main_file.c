//////////////////////////////////////////////////////////////////
//**************************************************************//
//                                                              //
//      Author: Nikita Tretyakov                                //
//      PhD-student, Institut fuer Theoretische Physik          //
//      Georg-August-Universitaet Goettingen                    //
//                                                              //
//      PostDoc, Max-Planck-Institute for Polymer Research      //
//      Mainz                                                   //
//                                                              //
//      tretyakov@mpip-mainz.mpg.de       niktre@gmail.com      //
//                                                              //
//      Former Supervisor:        Prof. Dr. Marcus Mueller      //
//                                                              //
//                                Start Date 29 April 2009      //
//      Version 4.00              Modified   30 Juli  2011      //
//                                                              //
//**************************************************************//
//////////////////////////////////////////////////////////////////

#define NDIM	3
#define NHIST_C	1
#define NHIST_V	7
#include "in_functions/in_mddefs.h"
#include "in_functions/external_files.h"
#include "in_functions/files.h"
#include "main_file.h"

int main (int argc, char **argv) {
   MakeParamFile (stdout, stdout, stdout, argc, argv);
   GetNameList (argc, argv);
   PrintNameList (stdout);
   SetParams ();
   SetupJob ();
   moreCycles = 1;
   while (moreCycles) {
      SingleStep ();
      if (stepCount >= stepLimit) moreCycles = 0;
   }
}

/* MD SINGLE STEP IMPLEMENTATION */
void SingleStep () {
   ++stepCount;
   
   timeNow = stepCount * deltaT;
   LeapfrogStep (1);
   ApplyBoundaryCond ();
   if (nebrNow == 1) {
      nebrNow = 0;
      dispHi = 0.;
      BuildNebrList ();
   }
   
   /* do not calculate pressure tensor for droplets (useless) */
   if (problem == 2) {
      ComputeForces ();
      ComputeChainBondForces ();
   } else {
      ComputeForcesPT ();
      ComputeChainBondForcesPT ();
   }
   
   /* if we do not have flows, put a repulsive wall on the top */
   if (problem != 1 && flatTopW == 0) {
      ComputeWallForces ();
   } else {
   }
   
   /* if we have wall oscillations */
   if (problem == 2 && (stepCount > stepEquil) && ampSchwing != .0) {
      ComputeWandSchwingungen ();
   } else {
   }
   
   /* if we want to change epsWall during the run */
   if ( (stepCount > stepEquil) && (epsChemWall != .0) && ((stepCount - stepEquil) % (int)(periodChem / deltaT) == 0) ) {
      epsChemStore = epsWall;
      epsWall = epsChemWall;
      epsChemWall = epsChemStore;
      U0_surf = 4. * epsWall * rriCut3 * ss3_sur * (rriCut3 * ss3_sur - 1.); // shifting of LJ-potential between liquid and substrate
   }
   
   /* if there is a body force add it after equilibration */
   if ((stepCount > stepEquil) && (gravField != 0.)) {
      ComputeExternalForces ();
   } else {
   }
   
   /* if we want to change the width of the channel "on the fly" */
   if (problem == 1 && shiftWidth != 0. && (stepCount <= (int) (0.5*stepEquil))) {
      ModifyWidth ();
   }
   
   LeapfrogStep (2);
   EvalProps ();
   AccumProps (1);

   if ((VTF_stepPrint != 0) && (stepCount % VTF_stepPrint == 0)) {
      WrapChainsIntoCell ();
      PrintVTF (stdout);
   }
   if (stepCount % stepAvg == 0) {
      AccumProps (2);
      PrintSummary (stdout, problem);
      AccumProps (0);
      MakeFilename (stepCount / stepAvg, storeVelFor);
      CVFoutput (stdout, storeVelFor, stdout);		// function to store coord, vel, forces
   }
   if ((stepCount > stepEquil) && (stepCount - stepEquil) % stepConVal == 0) {
      countDenPn += stepConVal;
      EvalDenDist (1);
      if (problem != 2) {
         EvalVelDist (1);
      }
      EvalControlValues (controlFile, 1);
      CoordGridAverage (1);
      VelGridAverage (1);
      PrintPosVelCM(stdout);
      addCMStep = 1;
      if (countDenPn == limitConVal) {
         EvalDenDist (2);
         PrintDenDist (stdout);
         
         EvalPrDist (2);
         PrintPrDist (stdout);
         EvalPrDist (0);
         EvalDenDist (0);
         
         if (problem != 2) {
            EvalVelDist (2);
            PrintVelDist (stdout);
            EvalVelDist (0);
         } else {
         }
         
         EvalControlValues (controlFile, 2);
         PrintControlFile (stdout);
         addControlStep = 1;
         EvalControlValues (controlFile, 0);
         
         CoordGridAverage (2);
         VelGridAverage (2);
         EvalGridProfile ();
         PrintGridProfile (stdout, stdout, stdout, stdout);
         CoordGridAverage (0);
         VelGridAverage (0);
         
         countDenPn = 0;
      }
   }
}

void SetupJob () {
   AllocArrays ();
   InitRand (randSeed);
   
   if (restart == 0 && getIConf == 0) {		// start COMPLETELY NEW trajectory
      stepCount = 0;
      InitCoords (problem);
      InitSubstrate (problem);
      CopyAbsWrapChains ();
      AssignToChain ();
      PrintInitConf (stdout, stdout);
      InitVels ();
      InitAccels ();
      AccumProps (0);
      addControlStep = 0;	// stage of adding info to control-file
      addCMStep = 0;
   } else if (restart == 0 && getIConf == 1) {	// start with a file providing coords, vels, accelerations
      stepCount = 0;
      ReadCoordsVelsAccels (stdout, stdout, stdout, stdout);
      PrintInitConf (stdout, stdout);
      AccumProps (0);
      addControlStep = 0;	// stage of adding info to control-file
      addCMStep = 0;
   } else if (restart == 1 && oldEndStep != 0) {	// just CONTINUE EXISTING trajectory
      stepCount = oldEndStep;
      if (problem == 1 && shiftWidth != 0. && (stepCount >= (int) (0.5*stepEquil))) {
         region.z += shiftWidth;
      } else {
      }
      ReadCoordsVelsAccels (stdout, stdout, stdout, stdout);
      AccumProps (0);
      addControlStep = 0;
      addCMStep = 0;
   }
   EvalControlValues (controlFile, 0);
   EvalDenDist (0);
   EvalPrDist (0);
   if (problem != 2) {
      EvalVelDist (0);
   }
   CoordGridAverage (0);
   VelGridAverage (0);

   countDenPn = 0;
   nebrNow = 1;
   countConVal = 0;
}

/* SET SYSTEM PARAMETERS, PROVIDING SIMULATION BOX AND SYSTEM */
void SetParams () {
   int temp_bin;
   
   rCut = 2. * pow (2., 1. / 6.);		// cut-off radius
   zCutWall = 3.0;				// cut-off distance for wall-forces
   nChain = 2 * V_PROD (initUchain);	// 2 chains in bcc lattice
   if (nChain == 2) nChain = 1;		// correcting 1-chain bug in bcc-lattice
   nPolyMol = nChain * chainLen;		// total number of beads in chains
   
   /* defining region size */
   if (oldSub == 1) {
      density = 0.788013429;
      V_SET (initUcell, 33775, 214500, 51952);
   } else {
      V_SET (initUcell, 2111.584723, 351.9307871, 573.0043);
   }
   V_S_COPY (region, pow (nPolyMol / (density * V_PROD(initUcell)), 1./3.), initUcell);	// dimensions of one cell with 2 chains
   V_SCALE (region, pow (1. / V_PROD(initUchain), 1./3.));
   V_MUL (region, region, initUchain);
   
   /* specifying surface and vapor space of the box */
   surfCell.z = nSurfLayer;					// replication of surface-cell in z-direction
   if (oldSub == 1) surfCell.y = (int)(region.y / surf_dist);	// replication of surface-cell in y-direction
   else surfCell.y = (int)(region.y / (surf_dist * sqrt(3.)));	// replication of surface-cell in y-direction
   surfCell.x = (int)(region.x / (surf_dist * sqrt(3.)));		// replication of surface-cell in x-direction
   
   if (flatTopW == 1) {
      nSurf = 4 * V_PROD (surfCell) + 6 * surfCell.x * surfCell.y;	// preliminary number of surface atoms for bottom surface
   } else if (flatTopW == 0 && problem == 1) {
      nSurf = 2 * 4 * V_PROD (surfCell);				// preliminary number of surface atoms for bottom surface
   } else if (flatTopW == 0 && problem != 1) {
      nSurf = 4 * V_PROD (surfCell);					// preliminary number of surface atoms for bottom surface
   }
   
   if (surfCell.z == 1 && oldSub == 0) {
      nSurf += 2 * V_PROD (surfCell);
   }
   
   /* substrate unit variables
    * unter Bedingung dass alle Oberflächigen Units sind gleich(!!!) (Oberfläche, Symmetrie, usw.)
    */
   // betrachten x-Richtung Struktur der Oberfläche
   if (surfDel.x == 0 && slope == 0) {
      totSUnitNum = surfCell.y;
      effSUArea = region.x;
   } else {
      totSUnitNum = (surfCell.x * surfCell.y);
      effSUArea = (surfSave.x + 0.5) * region.x / surfCell.x;
   }
   // betrachten y-Richtung Struktur der Oberfläche
   if (surfDel.y == 0) {
      totSUnitNum /= ( (surfSave.x + surfDel.x + 1) * (surfCell.y) );
      effSUArea *= region.y;
   } else {
      totSUnitNum /= ( (surfSave.x + surfDel.x + 1) * (surfSave.y + surfDel.y + 1) );
      effSUArea *= (surfSave.y + 0.5) * region.y / surfCell.y;
   }
   if (problem == 1) {
      totSUnitNum *= 2;
      effSUArea *= 2.;
   }
   
   nMol = nChain * chainLen + nSurf;	// total number of atoms in the system
   if (problem == 1 && flatTopW == 0 && oldSub == 0) { // WORKS
      addSpace = nSurfLayer * 2. * surf_dist * 1.01 + addWidth;
   } else if (problem == 1 && flatTopW == 1 && oldSub == 0 && slope == 0) { // WORKS
      addSpace = (nSurfLayer + 1.5) * surf_dist * 1.01 + addWidth;
   } else if (problem == 1 && flatTopW == 1 && oldSub == 0 && slope == 1) { //
      addSpace = (nSurfLayer + 1 + 1.5) * surf_dist * 1.01 + addWidth;
   } else if (problem == 1 && flatTopW == 0 && oldSub == 1 ) { //
      addSpace = nSurfLayer * 2. * surf_dist * 1.01 * sqrt(3.) + addWidth;
   } else if (problem == 1 && flatTopW == 1 && oldSub == 1 ) { //
      addSpace = (nSurfLayer + 1.5) * surf_dist * 1.01 * sqrt(3.) + addWidth;
   } else {
      //		addSpace = 125.; // height of vapor space for sh_dynamic 19712 and 9856 chains
      addSpace = 75.; // height of vapor space for sh_dynamic 4928 and 2112 chains
      //		addSpace = 50.; // height of vapor space for sh_dynamic 352 and 1056 chains and other drops
   }
   //	region.z += addSpace;			// add vapor space to the box
   
   /* compensate the third plane of atoms in case of surfCell.z == 1 for newSub */
   if (oldSub == 0 && surfCell.z == 1) {
      addSpace = 0.5 * 1.59;
      //		region.z += addSpace;			// add vapor space to the box
   }
   
   /* add space to fit integer number of hSlab's*/
   addSpace += hSlab * (int) (region.z / hSlab) - region.z;
   region.z += addSpace;			// add vapor space to the box
   
   /* general simulation constants */
   numSlabs = (int) (region.z / hSlab);	// number of slabs in the system
   numSlabsX = (int) (region.x / hSlab) + 1;	// number of slabs in the system
   halfSl = (int) (numSlabs / 2.);		// number of slabs in positve or negative region
   halfSlX = (int) (numSlabsX / 2.);	// number of slabs in positve or negative region
   rrCut = SQR (rCut);			// squared cut-off distance
   rriCut = 1. / SQR (rCut);		// inverted squared cut-off distance
   rriCut3 = CUBE (rriCut);		// inverted sixth power of cut-off distance
   s3_sur = CUBE (sigma_sur);		// cube of sigma
   ss3_sur = SQR (s3_sur);			// six-power of sigma
   zzCutWall = SQR (zCutWall);		// squared cut-off distance for wall-forces
   U0 = 4. * rriCut3 * (rriCut3 - 1.);	// shifting of LJ-potential in liquid
   U0_surf = 4. * epsWall * rriCut3 * ss3_sur * (rriCut3 * ss3_sur - 1.); // shifting of LJ-potential between liquid and substrate
   R2 = SQR (1.5);				// squared cut-off distance for FENE-forces
   k = 30.;				// FENE-constant of interaction
   V_S_COPY (cells, 1. / (rCut + rNebrShell), region);
   nebrTabMax = nebrTabFac * nMol;
   omegaSchwing = 2. * M_PI / periodSchwing;
   
   /* set coordinate and velocity grids */
   temp_bin = (int) (region.x / sizeCoordBin);
   if (temp_bin % 2 == 0)
      sizeCoordGrid.x = temp_bin;
   else
      sizeCoordGrid.x = temp_bin + 1;
   temp_bin = (int) (region.y / sizeCoordBin);
   if (temp_bin % 2 == 0)
      sizeCoordGrid.y = temp_bin;
   else
      sizeCoordGrid.y = temp_bin + 1;
   temp_bin = (int) (region.z / sizeCoordBin);
   if (temp_bin % 2 == 0)
      sizeCoordGrid.z = temp_bin;
   else
      sizeCoordGrid.z = temp_bin + 1;
   
   temp_bin = (int) (region.x / sizeVelBin);
   if (temp_bin % 2 == 0)
      sizeVelGrid.x = temp_bin;
   else
      sizeVelGrid.x = temp_bin + 1;
   sizeVelGrid.y = 4;
   temp_bin = (int) (region.z / sizeVelBin);
   if (temp_bin % 2 == 0)
      sizeVelGrid.z = temp_bin;
   else
      sizeVelGrid.z = temp_bin + 1;
}

/* DYNAMICAL VARIABLES */
void AllocArrays () {
   long k;
   
   ALLOC_MEM (mol, nMol, Mol);
   ALLOC_MEM (ghost, nSurf, Ghost);
   ALLOC_MEM (centUCinZ, 2 * surfCell.z - 1, real);
   ALLOC_MEM (histDen, numSlabs, real);
   ALLOC_MEM (snapDen, numSlabs, real);
   ALLOC_MEM (histAvDen, numSlabs, real);
   ALLOC_MEM (histDenX, numSlabsX, real);
   ALLOC_MEM (histPn, numSlabs, real);
   ALLOC_MEM (histPnW, numSlabs, real);
   ALLOC_MEM (histPt, numSlabs, real);
   ALLOC_MEM (histPtW, numSlabs, real);
   ALLOC_MEM (histPx, numSlabs, real);
   ALLOC_MEM (histPy, numSlabs, real);
   ALLOC_MEM (snapVx, numSlabs, real);
   ALLOC_MEM (snapTx, numSlabs, real);
   ALLOC_MEM (snapT, numSlabs, real);
   ALLOC_MEM (histVx, numSlabs, real);
   ALLOC_MEM (histTx, numSlabs, real);
   ALLOC_MEM (histT, numSlabs, real);
   ALLOC_MEM (histAvVx, numSlabs, real);
   ALLOC_MEM (histAvTx, numSlabs, real);
   ALLOC_MEM (histAvT, numSlabs, real);
   ALLOC_MEM (profileT, sizeVelGrid.z, real);
   ALLOC_MEM (profileVx, sizeVelGrid.z, real);
   ALLOC_MEM (profileZX, sizeCoordGrid.x, real);
   ALLOC_MEM (profile_botZX, sizeCoordGrid.x, real);
   
   ALLOC_MEM (part_y_p_coord, sizeCoordGrid.y, int);
   ALLOC_MEM (part_y_m_coord, sizeCoordGrid.y, int);
   ALLOC_MEM (rx_p_coord, sizeCoordGrid.y, real);
   ALLOC_MEM (rx_m_coord, sizeCoordGrid.y, real);
   
   ALLOC_MEM (ocCellXY, sizeVelGrid.z, int);
   ALLOC_MEM (cellList, V_PROD (cells) + nMol, int);
   ALLOC_MEM (nebrTab, 2 * nebrTabMax, int);
   ALLOC_MEM2 (den2D, sizeCoordGrid.x, sizeCoordGrid.z, real);
   ALLOC_MEM2 (denMap, sizeCoordGrid.x, sizeCoordGrid.z, real);
   ALLOC_MEM2 (histCoordGrid, NHIST_C, V_PROD (sizeCoordGrid), real);
   ALLOC_MEM2 (snapCoordGrid, NHIST_C, V_PROD (sizeCoordGrid), real);
   ALLOC_MEM2 (histVelGrid, NHIST_V, V_PROD (sizeVelGrid), real);
   ALLOC_MEM2 (snapVelGrid, NHIST_V, V_PROD (sizeVelGrid), real);
   ALLOC_MEM2 (vel2D_x, sizeVelGrid.x, sizeVelGrid.z, real);
   ALLOC_MEM2 (vel2D_z, sizeVelGrid.x, sizeVelGrid.z, real);
   ALLOC_MEM2 (vel2D_xx, sizeVelGrid.x, sizeVelGrid.z, real);
   ALLOC_MEM2 (vel2D_zz, sizeVelGrid.x, sizeVelGrid.z, real);
   ALLOC_MEM2 (dropVel2D, sizeVelGrid.x, sizeVelGrid.z, real);
}

/* NEIGHBOR LIST */
void BuildNebrList () {
   VecR dr, invWid, rs, shift;
   VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
   real rrNebr;
   int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
   
   rrNebr = SQR (rCut + rNebrShell);
   V_DIV (invWid, cells, region);
   for (n = nMol; n < nMol + V_PROD (cells); n++) cellList[n] = -1;
   DO_MOL {
      V_S_ADD (rs, mol[n].r, 0.5, region);
      V_MUL (cc, rs, invWid);
      c = V_LINEAR (cc, cells) + nMol;
      cellList[n] = cellList[c];
      cellList[c] = n;
   }
   nebrTabLen = 0;
   for (m1z = 0; m1z < cells.z; m1z++) {
      for (m1y = 0; m1y < cells.y; m1y++) {
         for (m1x = 0; m1x < cells.x; m1x++) {
            V_SET (m1v, m1x, m1y, m1z);
            m1 = V_LINEAR (m1v, cells) + nMol;
            for (offset = 0; offset < N_OFFSET; offset++) {
               V_ADD (m2v, m1v, vOff[offset]);
               V_ZERO (shift);
               V_CELL_WRAP (x);
               V_CELL_WRAP (y);
               if (m2v.z < 0 || m2v.z >= cells.z) continue;		// cancellation of periodicity in z-direction
               //			V_CELL_WRAP_ALL ();
               m2 = V_LINEAR (m2v, cells) + nMol;
               DO_CELL (j1, m1) {
                  DO_CELL (j2, m2) {
                     if ( (m1 != m2 || j2 < j1) &&
                        ((mol[j1].inChain == -1 && mol[j2].inChain != -1) ||	// if j1 is surface atom and j2 is not
                        (mol[j2].inChain == -1 && mol[j1].inChain != -1) ||		// or if j2 is surface atom and j1 is not
                        ((mol[j1].inChain != -1 && mol[j2].inChain != -1) &&	// or j1 and j2 not from the surface both combined with:
                        (mol[j1].inChain != mol[j2].inChain ||			// (j1 and j2 from different chains
                        abs (j1 - j2) > 1) )) ) {					// or j1 and j2 separated by at least 1 bead)
                           V_SUB (dr, mol[j1].r, mol[j2].r);
                           V_V_SUB (dr, shift);
                           if (V_LEN_SQ (dr) < rrNebr) {
                              if (nebrTabLen >= nebrTabMax) {
                                 ErrExit (ERR_TOO_MANY_NEBRS);
                              }
                              nebrTab[2 * nebrTabLen] = j1;
                              nebrTab[2 * nebrTabLen + 1] = j2;
                              ++nebrTabLen;
                           }
                        }
                  }
               }
            }
         }
      }
   }
}

/* COMPUTATION OF FORCES BETWEEN NON-CONNECTED BEADS AND BEADS-SUBSTRATE_ATOMS(with neighbor-list) */
void ComputeForces () {
   VecR dr, dv;		// distance and velocity difference between i- and j- beads (VECTORS)
   real scal, fDispRand, theta, uVal, omega_r;	// dissipative and random forces variables
   real zeta = sqrt (24. * temperature * gamma_d / deltaT);
   real fcVal, rr, rri, rri3;
   int j1, j2, n;
   
   DO_MOL V_ZERO (mol[n].ra);
   uSum = 0.;
   U_LJ_sub = 0.;
   virSum = 0.;
   
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
            /* Lennard-Jones forces */
            rri = 1. / rr;
            rri3 = CUBE (rri);
            if ((mol[j1].inChain != -1) && (mol[j2].inChain != -1)) {	// interaction of nonbonded chain particles
               fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
               uVal = 4. * rri3 * (rri3 - 1.) - U0;
               V_V_S_ADD (mol[j1].ra, fcVal, dr);
               V_V_S_ADD (mol[j2].ra, -fcVal, dr);

            } else if ((mol[j1].inChain == -1) && (mol[j1].inSUnit != -2) && (mol[j2].inChain != -1)) {
               fcVal = 48. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 0.5) * rri; // interaction with the substrate
               uVal = 4. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 1.) - U0_surf;
               U_LJ_sub += uVal;
               V_V_S_ADD (mol[j2].ra, -fcVal, dr);

            } else if ((mol[j1].inChain == -1) && (mol[j1].inSUnit == -2) && (mol[j2].inChain != -1)) {
               fcVal = 48. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur * rri; // repulsive interaction with the flatTopW
               uVal = 4. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur;
               U_LJ_sub += uVal;
               V_V_S_ADD (mol[j2].ra, -fcVal, dr);

            } else if ((mol[j1].inChain != -1) && (mol[j2].inChain == -1) && (mol[j2].inSUnit != -2)) {
               fcVal = 48. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 0.5) * rri; // interaction with the substrate
               uVal = 4. * epsWall * rri3 * ss3_sur * (rri3 * ss3_sur - 1.) - U0_surf;
               U_LJ_sub += uVal;
               V_V_S_ADD (mol[j1].ra, fcVal, dr);

            } else if ((mol[j1].inChain != -1) && (mol[j2].inChain == -1) && (mol[j2].inSUnit == -2)) {
               fcVal = 48. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur * rri; // repulsive interaction with the flatTopW
               uVal = 4. * epsTopW * rri3 * ss3_sur * rri3 * ss3_sur;
               U_LJ_sub += uVal;
               V_V_S_ADD (mol[j1].ra, fcVal, dr);

            } else {
            }
            
            uSum += uVal;
            virSum += fcVal * rr;
            
            /* dissipative forces */
            omega_r = sqrt(rri) - 1. / rCut;		// omega_r / abs (Rij)
            V_SUB (dv, mol[j1].rv, mol[j2].rv);
            scal = V_DOT (dr, dv);
            /* random forces */
            theta = RandR () - 0.5;
            fDispRand = (zeta * theta - gamma_d * scal * omega_r) * omega_r;

				if (mol[j1].inChain != -1) {
					V_V_S_ADD (mol[j1].ra, fDispRand, dr);
				} else if (mol[j2].inChain != -1) {
					V_V_S_ADD (mol[j2].ra, -fDispRand, dr);
				}
         }
      }
   }
}

/* COMPUTATION OF FORCES BETWEEN CONNECTED BEADS IN CHAIN */
void ComputeChainBondForces () {
   VecR dr, dv;		// distance and velocity difference between i- and j- beads (VECTORS)
   real fcVal, fTot, rr, rri, rri3, uVal, rrno;
   int i, j1, j2, n;
   
   real scal, fDispRand, theta, omega_r;	// DPD
   real zeta = sqrt (24. * temperature * gamma_d / deltaT);
   
   fTot = 0.;
   
   for (n = 0; n < nChain; n++) {
      for (i = 0; i < chainLen - 1; i++) {
         j1 = n * chainLen + i;
         j2 = j1 + 1;
         V_SUB (dr, mol[j1].r, mol[j2].r);
         V_WRAP_ALL (dr);
         rr = V_LEN_SQ (dr);
         if (rr < rrCut) {
            /* Lennard-Jones forces */
            rri = 1. / rr;
            rri3 = CUBE (rri);
            fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
            uVal = 4. * rri3 * (rri3 - 1.) - U0;
            V_V_S_ADD (mol[j1].ra, fcVal, dr);
            V_V_S_ADD (mol[j2].ra, -fcVal, dr);
            
            uSum += uVal;
            //		        if ((stepCount > stepEquil) && ((stepCount - stepEquil) % stepConVal == 0)) {
            fTot = fcVal;
            //			}
            /* dissipative forces */
            omega_r = sqrt(rri) - 1. / rCut;	// %omega_r% divided by %|Rij|%
            V_SUB (dv, mol[j1].rv, mol[j2].rv);
            scal = V_DOT (dr, dv);
            /* random forces */
            theta = RandR () - 0.5;			// %theta% is in range [-0.5 ; 0.5]
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
         }
      }
   }
}

/* WALL FORCES ON THE TOP FACE OF THE BOX */
void ComputeWallForces () {
   int n, j, N_j;
   real zz, zzi, zzi3, fwVal, wVal;
   real z_dist, zi, dTop;
   real secTermZ;	// pressure tensor calculating
   
   DO_POLY_MOL {
      if (mol[n].r.z >= 0) {
         z_dist = region.z / 2. - mol[n].r.z;
         N_j = halfSl + (int) (mol[n].r.z / hSlab);
         zz = SQR (z_dist);
         zi = 1. / z_dist;
         if (zz <= zzCutWall) {
            zzi = 1. / zz;
            zzi3 = CUBE(zzi);
            fwVal = 4. * 9. * zzi3 * zzi3 * z_dist;
            wVal = zzi3 * zzi / z_dist;
            mol[n].ra.z -= fwVal * z_dist;
            uSum += wVal;
            virSum += fwVal * zz;
            if (stepCount > stepEquil && (stepCount - stepEquil) % stepConVal == 0) {
               secTermZ = zz * fwVal;
               dTop = z_dist - hSlab * (numSlabs - 1 - N_j);
               for (j = N_j; j < numSlabs; j++) {
                  if (j > N_j) histPnW[j] += secTermZ * hSlab * zi;
                  else histPnW[j] += secTermZ * dTop * zi;
               }
            }
         }
      }
   }
}

void ComputeWandSchwingungen () {
   int at;
   real phaseTopW, phaseBotW;
   
   DO_SURF {
      if (mol[nPolyMol + at].r.z > 0. && flatTopW == 0) {
         phaseTopW = omegaSchwing * deltaT * stepEquil;
         mol[nPolyMol + at].rv.z = - ampSchwing * omegaSchwing * sin (omegaSchwing * timeNow - phaseTopW);
      } else if (mol[nPolyMol + at].r.z > 0. && flatTopW == 1) {
      } else if (mol[nPolyMol + at].r.z < 0.) {
         phaseBotW = omegaSchwing * deltaT * stepEquil;
         mol[nPolyMol + at].rv.z = ampSchwing * omegaSchwing * sin (omegaSchwing * timeNow - phaseBotW);
      } else {
      }
   }
}

void ComputeExternalForces () {
   int n;
   
   DO_POLY_MOL mol[n].ra.x += gravField;
}

/* LEAP-FROG METHOD FOR INTEGRATING NEWTON'S EQUATIONS*/
void LeapfrogStep (int part) {
   int n;
   
   if (part == 1) {
      DO_MOL {
         V_V_S_ADD (mol[n].rv, 0.5 * deltaT, mol[n].ra);
         V_V_S_ADD (mol[n].r, deltaT, mol[n].rv);
         V_V_S_ADD (mol[n].r_abs, deltaT, mol[n].rv);
      }
   }
   else {
      DO_MOL V_V_S_ADD (mol[n].rv, 0.5 * deltaT, mol[n].ra);
   }
}

/* BOUNDARY CONDITIONS - WRAP CHAINS */
void ApplyBoundaryCond () {
   int n;
   
   if (shearRate == 0.) {
      DO_POLY_MOL V_WRAP_ALL (mol[n].r);
   } else {
      DO_MOL V_WRAP_ALL (mol[n].r);
   }
}

/* MAKE INITIAL CONFORMATION */
void InitCoords (int SysType) {
   VecR c, gap;
   real bx, bz, z_defl, film_exist_z;
   int j, m, n, tempCh;
   int nx, ny, nz;
   int nxLoc, nyLoc, nzLoc; 	// local number of chain-construction cell (depends on the problem)
   
   /* constructing chains */
   bx = 1.01 * cos (M_PI / 4.);
   bz = 0.90 * sin (M_PI / 4.);
   film_exist_z = region.z - addSpace;
   z_defl = (film_exist_z / initUchain.z - bz) / 2.;
   n = 0;
   tempCh = 0;
   gap.x = region.x / initUchain.x;	// sets x-size of chain unit cell
   gap.y = region.y / initUchain.y;	// sets y-size of chain unit cell
   gap.z = film_exist_z / initUchain.z;	// sets z-size of chain unit cell
   
   if (SysType == 2) {
      nxLoc = (int) (initUchain.x / 2.);
      nyLoc = initUchain.y;
      nzLoc = 2 * initUchain.z;
   } else {
      nxLoc = initUchain.x;
      nyLoc = initUchain.y;
      nzLoc = initUchain.z;
   }
   for (nz = 0; nz < nzLoc; nz++) {
      for (ny = 0; ny < nyLoc; ny++) {
         for (nx = 0; nx < nxLoc; nx++) {
            //		if ( ((chainRed != 0) && ((nx + ny + nz + 1) % chainRed != 0)) || (chainRed == 0) ) {	// do we want to substract some chains?
            //		if ( ((chainRed != 0) && (nx >= (int)(.5 * nxLoc / chainRed) && nx < nxLoc - (int)(.5 * nxLoc / chainRed) )) || (chainRed == 0) ) {	// do we want to substract some chains?
												// value chainRed == 0 means we don't want to substract!
            if ( tempCh < chainRed || chainRed == 0) {
               V_SET (c, nx, ny, nz);
               V_MUL (c, c, gap);
               if (SysType == 2) c.x -= 0.25 * region.x;	//if a droplet, reduce the x size 2 times and put it on the top
               else c.x -= 0.5 * region.x;
               c.y -= 0.5 * region.y;
               if (oldSub == 1) {
                  c.z -= 0.5 * region.z - surfCell.z * sqrt(3.) * surf_dist + 0.07937 - bz;
               } else if (oldSub == 0 && surfCell.z == 1) {
                  c.z -= 0.5 * region.z - surfCell.z * surf_dist - surf_dist * 0.5 - bz;
               } else if (oldSub == 0 && surfCell.z != 1) {
                  c.z -= 0.5 * region.z - surfCell.z * surf_dist - bz;
               } else {
               }
               
               if (SysType == 2) c.z += surf_dist * 3.5;	// lift initial droplet by half unit cell
               if (SysType == 1 && addWidth != 0.0) c.z += .5 * addWidth; // shift the liquid when modifying the channel width
               
               for (j = 0; j < 2; j++) {
                  for (m = 0; m < chainLen; m++) {
                     if (z_defl >= 0 && oldSub == 0) {
                        V_SET (mol[n].r, m * bx, 0., (m % 2 - 0.5) * bz);
                     } else if (z_defl >= 0 && oldSub == 1){
                        V_SET (mol[n].r, 0.,  m * bx, (m % 2 - 0.5) * bz);
                     } else {
                        V_SET (mol[n].r, 0., (m % 2 - .5) * bx, (m + 4) * bz);
                     }
                     mol[n].r.x = mol[n].r.x + 0.5 * j * gap.x;
                     mol[n].r.y = mol[n].r.y + 0.5 * j * gap.y;
                     mol[n].r.z = mol[n].r.z + 0.5 * j * gap.z;
                     V_V_ADD (mol[n].r, c);
                     mol[n].inSUnit = -1;		// doesn't belong to the unit of a substrate
                     ++n;
                  }
               }
               tempCh += 2;
            } else {
            }
         }
      }
   }
   nPolyMol = n;
   nChain = (int) nPolyMol / chainLen;
}

void InitSubstrate (int SysType) {
   VecR cSurf;
   int j, n, at;
   int nx, ny, nz;
   int SUnitNum, MinSUnitNum, MaxSUnitNum, botSUnitNum;
   
   n = nPolyMol;
   
   SUnitNum = 0;
   MinSUnitNum = 0;
   MaxSUnitNum = 0;
   botSUnitNum = 0;
   N_high_surf = 0;
   dTop_high_surf = -region.z / 2.;
   topSubAtZ = -region.z / 2.;
   
   /* constructing surface */
   switch (SysType) {
      case 0:
         if (oldSub == 1) V_SET (gapSurf, sqrt(3.) * surf_dist, surf_dist, sqrt(3.) * surf_dist);
         else V_SET  (gapSurf, sqrt(3.) * surf_dist, sqrt(3.) * surf_dist, surf_dist);
         //			V_SET  (gapSurf, surf_dist, sqrt(3.) * surf_dist, sqrt(3.) * surf_dist);
         break;
      case 1:
         if (oldSub == 1) V_SET (gapSurf, sqrt(3.) * surf_dist, surf_dist, sqrt(3.) * surf_dist);
         else V_SET  (gapSurf, sqrt(3.) * surf_dist, sqrt(3.) * surf_dist, surf_dist);
         //			V_SET  (gapSurf, surf_dist, sqrt(3.) * surf_dist, sqrt(3.) * surf_dist);
         break;
      case 2:
         if (oldSub == 1) V_SET (gapSurf, sqrt(3.) * surf_dist, surf_dist, sqrt(3.) * surf_dist);
         else V_SET  (gapSurf, sqrt(3.) * surf_dist, sqrt(3.) * surf_dist, surf_dist);
         break;
   }
   
   /* bottom substrate */
   for (nz = 0; nz < surfCell.z; nz++) {
      /* set min and max numbers of surface units to zero for every new layer of the substrate */
      MinSUnitNum = 0;
      MaxSUnitNum = 0;
      for (ny = 0; ny < surfCell.y; ny++) {
         /* woodoo with assigning numbers to surface units */
         if ( (ny > 0) && (ny % (surfSave.y + surfDel.y + 1) >= (surfSave.y + 1)) && ((ny - 1) % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) { 		// if we are on the latest row of substrate Units
            MinSUnitNum = MaxSUnitNum;
         } else if ( (ny > 0) && (ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) && ((ny - 1) % (surfSave.y + surfDel.y + 1) >= (surfSave.y + 1)) ) {	// if we are on the first row of substrate Units
            SUnitNum = MinSUnitNum;
         } else if ( (ny > 0) && (ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) && ((ny - 1) % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) { 	// if we are already inside the row of substrate Units
            SUnitNum = MinSUnitNum;
         } else if (ny == 0) {
            SUnitNum = 0;
         } else {
         }
         for (nx = 0; nx < surfCell.x; nx++) {
            /* small part of woodoo, relating only current y-layer of the cells */
            if ( (nx > 0) && (nx % (surfSave.x + surfDel.x + 1) >= (surfSave.x + 1)) && ((nx - 1) % (surfSave.x + surfDel.x + 1) < (surfSave.x + 1)) && (ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) ++SUnitNum;
            // account for the case with sloped ridges, but surfDel.x = 0
            if (nx > 0 && slope == 1 && surfDel.x == 0 && (nx % (surfSave.x + 1)) == 0) ++SUnitNum;
            // construction of the bottom layer of the substrate
            if (nz == 0) {
               V_SET (cSurf, nx + 0.1, ny + 0.1, nz + 0.01);
               V_MUL (cSurf, cSurf, gapSurf);
               V_V_S_ADD (cSurf, -0.5, region);
               if (surfCell.z == 1 && oldSub == 0) {
                  for (j = 0; j < 6; j++) {
                     mol[n].r = cSurf;
                     if (j < 3) {
                        if (j != 0) mol[n].r.z += 0.5 * gapSurf.z;
                        if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                        if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                     } else if (j == 4) {
                        mol[n].r.z += gapSurf.z;
                        mol[n].r.x += 0.5 * gapSurf.x;
                        mol[n].r.y += 0.5 * gapSurf.y;
                     } else if (j == 5) {
                        mol[n].r.z += gapSurf.z;
                     } else {
                     }
                     N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                     dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                     topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                     mol[n].inSUnit = -1;		// doesn't belong to the unit of a substrate (this layer is totally flat)
                     ++n;
                  }
               } else {
                  for (j = 0; j < 4; j++) {
                     mol[n].r = cSurf;
                     if (j != 3) {
                        if (j != 0) mol[n].r.z += 0.5 * gapSurf.z;
                        if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                        if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                     }
                     N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                     dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                     topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                     mol[n].inSUnit = -1;		// doesn't belong to the unit of a substrate (this layer is totally flat)
                     ++n;
                  }
               }
            } else {
               // construction of the upper layers of the substrate
               if ( (surfDel.x == 0) || ( nx % (surfSave.x + surfDel.x + 1) < (surfSave.x + 1)) ) {
                  if ( (surfDel.y == 0) || ( ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) {
                     V_SET (cSurf, nx + 0.1, ny + 0.1, nz + 0.01);
                     V_MUL (cSurf, cSurf, gapSurf);
                     V_V_S_ADD (cSurf, -0.5, region);
                     
                     // the ridge of the pillar in x-dir
                     if ( (nx % (surfSave.x + surfDel.x + 1) == surfSave.x) && (surfDel.x != 0) ) {
                        if ( (ny % (surfSave.y + surfDel.y + 1) == surfSave.y) && (surfDel.y != 0) ) {
                           // the corner of the pillar
                           mol[n].r = cSurf;
                           N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                           dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                           topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                           mol[n].inSUnit = SUnitNum;
                           ++n;
                        } else {
                           // x-ridge of the pillar
                           if ( (slope == 1) && (nx % (surfSave.x + surfDel.x + 1) > (surfCell.z - 1 - nz)) ){
                           } else {
                              for (j = 0; j < 2; j++) {
                                 mol[n].r = cSurf;
                                 if (j != 1) {
                                    mol[n].r.y += 0.5 * gapSurf.y;
                                    mol[n].r.z += 0.5 * gapSurf.z;
                                 }
                                 N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                                 dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                                 topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                                 mol[n].inSUnit = SUnitNum;
                                 ++n;
                              }
                           }
                        }
                     } else if ( (ny % (surfSave.y + surfDel.y + 1) == surfSave.y) && (nx % (surfSave.x + surfDel.x + 1) != surfSave.x) && (surfDel.y != 0) ) {
                        // the ridge of the pillar in y-dir
                        for (j = 0; j < 2; j++) {
                           mol[n].r = cSurf;
                           if (j != 1) {
                              mol[n].r.x += 0.5 * gapSurf.x;
                              mol[n].r.z += 0.5 * gapSurf.z;
                           }
                           N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                           dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                           topSubAtZ = MAX (topSubAtZ, mol[n].r.z);
                           mol[n].inSUnit = SUnitNum;
                           ++n;
                        }
                     } else {
                        // the body of the pillar
                        if ( (slope == 1) && (nx % (surfSave.x + surfDel.x + 1) == (surfCell.z - 1 - nz)) ){
                           for (j = 0; j < 3; j++) {
                              mol[n].r = cSurf;
                              if (j != 2) {
                                 mol[n].r.y += 0.5 * gapSurf.y;
                                 if (j != 0) mol[n].r.z += 0.5 * gapSurf.z;
                                 if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                              }
                              N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                              dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                              topSubAtZ = MAX (topSubAtZ, mol[n].r.z);
                              mol[n].inSUnit = SUnitNum;
                              ++n;
                           }
                        } else if ( (slope == 1) && (nx % (surfSave.x + surfDel.x + 1) > (surfCell.z - 1 - nz)) ){
                        } else {
                           for (j = 0; j < 4; j++) {
                              mol[n].r = cSurf;
                              if (j != 3) {
                                 if (j != 0) mol[n].r.z += 0.5 * gapSurf.z;
                                 if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                                 if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                              }
                              N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                              dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                              topSubAtZ = MAX (topSubAtZ, mol[n].r.z);
                              mol[n].inSUnit = SUnitNum;
                              ++n;
                           }
                        }
                     }
                     
                  }
               }
            }
         }
         MaxSUnitNum = MAX (MaxSUnitNum, SUnitNum);
      }
      centUCinZ[nz] = cSurf.z + 0.25 * gapSurf.z;
   }
   botSUnitNum = MaxSUnitNum;

   /* top substrate */
   atTopWall = 0;		// is needed for setting up a flow in channel with walls of different topology
   if (problem == 1 && flatTopW == 0) {
      for (nz = 0; nz < surfCell.z; nz++) {
         /* set min and max numbers of surface units to zero for every new layer of the substrate */
         MinSUnitNum = botSUnitNum;
         MaxSUnitNum = botSUnitNum;
         for (ny = 0; ny < surfCell.y; ny++) {
            /* woodoo with assigning numbers to surface units */
            if ( (ny > 0) && (ny % (surfSave.y + surfDel.y + 1) >= (surfSave.y + 1)) && ((ny - 1) % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) { 		// if we are on the latest row of substrate Units
               MinSUnitNum = MaxSUnitNum;
            } else if ( (ny > 0) && (ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) && ((ny - 1) % (surfSave.y + surfDel.y + 1) >= (surfSave.y + 1)) ) {	// if we are on the first row of substrate Units
               SUnitNum = MinSUnitNum;
            } else if ( (ny > 0) && (ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) && ((ny - 1) % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) { 	// if we are already inside the row of substrate Units
               SUnitNum = MinSUnitNum;
            } else if (ny == 0) {
               SUnitNum = botSUnitNum;
            } else {
            }
            for (nx = 0; nx < surfCell.x; nx++) {
               /* small part of woodoo, relating only current y-layer of the cells */
               if ( (nx > 0) && (nx % (surfSave.x + surfDel.x + 1) >= (surfSave.x + 1)) && ((nx - 1) % (surfSave.x + surfDel.x + 1) < (surfSave.x + 1)) && (ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) ++SUnitNum;
               // account for the case with sloped ridges, but surfDel.x = 0
               if (nx > 0 && slope == 1 && surfDel.x == 0 && (nx % (surfSave.x + 1)) == 0) ++SUnitNum;
               // construction of the bottom layer of the substrate
               if (nz == 0) {
                  V_SET (cSurf, nx + 0.1, ny + 0.1, -nz - 0.01);
                  V_MUL (cSurf, cSurf, gapSurf);
                  cSurf.x -= 0.5 * region.x;
                  cSurf.y -= 0.5 * region.y;
                  cSurf.z += 0.5 * region.z;
                  for (j = 0; j < 4; j++) {
                     mol[n].r = cSurf;
                     if (j != 3) {
                        if (j != 0) mol[n].r.z -= 0.5 * gapSurf.z;
                        if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                        if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                     }
                     //				N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                     //				dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                     //				topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                     mol[n].inSUnit = -1;		// doesn't belong to the unit of a substrate (this layer is totally flat)
                     ++n;
                  }
               } else if (nz > 0) {
                  // construction of the upper layers of the substrate
                  if ( (surfDel.x == 0) || ( nx % (surfSave.x + surfDel.x + 1) < (surfSave.x + 1)) ) {
                     if ( (surfDel.y == 0) || ( ny % (surfSave.y + surfDel.y + 1) < (surfSave.y + 1)) ) {
                        V_SET (cSurf, nx + 0.1, ny + 0.1, -nz - 0.01);
                        V_MUL (cSurf, cSurf, gapSurf);
                        cSurf.x -= 0.5 * region.x;
                        cSurf.y -= 0.5 * region.y;
                        cSurf.z += 0.5 * region.z;
                        
                        // the ridge of the pillar in x-dir
                        if ( (nx % (surfSave.x + surfDel.x + 1) == surfSave.x) && (surfDel.x != 0) ) {
                           if ( (ny % (surfSave.y + surfDel.y + 1) == surfSave.y) && (surfDel.y != 0) ) {
                              // the corner of the pillar
                              mol[n].r = cSurf;
                              //						N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                              //						dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                              //						topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                              mol[n].inSUnit = SUnitNum;
                              ++n;
                           } else {
                              // x-ridge of the pillar
                              if ( (slope == 1) && (nx % (surfSave.x + surfDel.x + 1) > (surfCell.z - 1 - nz)) ){
                              } else {
                                 for (j = 0; j < 2; j++) {
                                    mol[n].r = cSurf;
                                    if (j != 1) {
                                       mol[n].r.y += 0.5 * gapSurf.y;
                                       mol[n].r.z -= 0.5 * gapSurf.z;
                                    }
                                    //								N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                                    //								dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                                    //								topSubAtZ = MAX(topSubAtZ, mol[n].r.z);
                                    mol[n].inSUnit = SUnitNum;
                                    ++n;
                                 }
                              }
                           }
                        } else if ( (ny % (surfSave.y + surfDel.y + 1) == surfSave.y) && (nx % (surfSave.x + surfDel.x + 1) != surfSave.x) && (surfDel.y != 0) ) {
                           // the ridge of the pillar in y-dir
                           for (j = 0; j < 2; j++) {
                              mol[n].r = cSurf;
                              if (j != 1) {
                                 mol[n].r.x += 0.5 * gapSurf.x;
                                 mol[n].r.z -= 0.5 * gapSurf.z;
                              }
                              //						N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                              //						dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                              //						topSubAtZ = MAX (topSubAtZ, mol[n].r.z);
                              mol[n].inSUnit = SUnitNum;
                              ++n;
                           }
                        } else {
                           // the body of the pillar
                           if ( (slope == 1) && (nx % (surfSave.x + surfDel.x + 1) == (surfCell.z - 1 - nz)) ){
                              for (j = 0; j < 3; j++) {
                                 mol[n].r = cSurf;
                                 if (j != 2) {
                                    mol[n].r.y += 0.5 * gapSurf.y;
                                    if (j != 0) mol[n].r.z -= 0.5 * gapSurf.z;
                                    if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                                 }
                                 //								N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                                 //								dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                                 //								topSubAtZ = MAX (topSubAtZ, mol[n].r.z);
                                 mol[n].inSUnit = SUnitNum;
                                 ++n;
                              }
                           } else if ( (slope == 1) && (nx % (surfSave.x + surfDel.x + 1) > (surfCell.z - 1 - nz)) ){
                           } else {
                              for (j = 0; j < 4; j++) {
                                 mol[n].r = cSurf;
                                 if (j != 3) {
                                    if (j != 0) mol[n].r.z -= 0.5 * gapSurf.z;
                                    if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                                    if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                                 }
                                 //								N_high_surf = MAX (N_high_surf, halfSl - (int) (-mol[n].r.z / hSlab) - 1);      // highest slab with surface atom
                                 //								dTop_high_surf = MAX (dTop_high_surf, mol[n].r.z + hSlab * (halfSl - N_high_surf));     // part of line connecting surface atoms in slab contained highest surf atom
                                 //								topSubAtZ = MAX (topSubAtZ, mol[n].r.z);
                                 mol[n].inSUnit = SUnitNum;
                                 ++n;
                              }
                           }
                        }
                        
                     }
                  }
               } else {
               }
            }
         }
      }
   }
   
   /* check the need to create a flat wall at the top of the box */
   if (flatTopW == 1 && slope == 0) {
      // if we have only one layer of the surface surface lattice (i.e. this is just repulsive wall)
      for (nz = 0; nz < 1; nz++) {
         for (ny = 0; ny < surfCell.y; ny++) {
            for (nx = 0; nx < surfCell.x; nx++) {
               V_SET (cSurf, nx + 0.1, ny + 0.1, -nz - 0.01);
               V_MUL (cSurf, cSurf, gapSurf);
               cSurf.x -= 0.5 * region.x;
               cSurf.y -= 0.5 * region.y;
               cSurf.z += 0.5 * region.z;
               for (j = 0; j < 6; j++) {
                  mol[n].r = cSurf;
                  if (j != 3 && j != 4 && j != 5) {
                     //				if (j != 3) {
                     if (j != 0) mol[n].r.z -= 0.5 * gapSurf.z;
                     if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                     if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                  } else if (j == 4) {
                     mol[n].r.z -= gapSurf.z;
                     mol[n].r.x += 0.5 * gapSurf.x;
                     mol[n].r.y += 0.5 * gapSurf.y;
                  } else if (j == 5) {
                     mol[n].r.z -= gapSurf.z;
                  }
                  mol[n].inSUnit = -2;		// doesn't belong to the unit of a substrate (this layer is totally flat)
                  ++n;
               }
            }
         }
      }
   } else if (flatTopW == 1 && slope == 1) {
      for (nz = 0; nz < 2; nz++) {
         for (ny = 0; ny < surfCell.y; ny++) {
            for (nx = 0; nx < surfCell.x; nx++) {
               V_SET (cSurf, nx + 0.1, ny + 0.1, -nz - 0.01);
               V_MUL (cSurf, cSurf, gapSurf);
               cSurf.x -= 0.5 * region.x;
               cSurf.y -= 0.5 * region.y;
               cSurf.z += 0.5 * region.z;
               //			for (j = 0; j < 6; j++) {
               for (j = 0; j < 4; j++) {
                  mol[n].r = cSurf;
                  //				if (j != 3 && j != 4 && j != 5) {
                  if (j != 3) {
                     if (j != 0) mol[n].r.z -= 0.5 * gapSurf.z;
                     if (j != 1) mol[n].r.x += 0.5 * gapSurf.x;
                     if (j != 2) mol[n].r.y += 0.5 * gapSurf.y;
                     //				} else if (j == 4) {
                     //					mol[n].r.z -= gapSurf.z;
                     //					mol[n].r.x += 0.5 * gapSurf.x;
                     //					mol[n].r.y += 0.5 * gapSurf.y;
                     //				} else if (j == 5) {
                     //					mol[n].r.z -= gapSurf.z;
                  }
                  mol[n].inSUnit = -1;		// it does belong to the unit of a substrate!!
                  ++n;
                  ++atTopWall;
               }
            }
         }
      }
   }
   
   nMol = n;
   nSurf = nMol - nPolyMol;
   atBotWall = nSurf - atTopWall;
   
   // copy initial positions of substrate atoms to ghost centers
   DO_SURF V_COPY (ghost[at].r, mol[nPolyMol + at].r);
}

void CopyAbsWrapChains () {
   int n;
   
   /* absolute coordinates for visualization and transport coefficients */
   DO_MOL V_COPY (mol[n].r_abs, mol[n].r);
   DO_MOL V_COPY (mol[n].r_vtf, mol[n].r);
   WrapChainsIntoCell ();
   ApplyBoundaryCond ();
}

/* SET INITIAL VELOCITIES FOR BEADS TO RANDOM NUMBERS (based on initial temperature) */
void InitVels () {
   int n, atWall;
   
   V_ZERO (vSum);
   velMag = sqrt (NDIM * (1. - 1. / nPolyMol) * temperature);
   DO_MOL V_ZERO (mol[n].rv);
   DO_POLY_MOL {
      VRand (&mol[n].rv);
      V_SCALE (mol[n].rv, velMag);
      V_V_ADD (vSum, mol[n].rv);
   }
   
   DO_POLY_MOL V_V_S_ADD (mol[n].rv, - 1. / nPolyMol, vSum);
   
   /* set initial velocities to substrate atoms */
   atWall = (int) (nSurf / 2);
   
   if (shearRate == 0.) {
   } else if (atTopWall == 0) {
      for (n = nPolyMol; n < nPolyMol + atWall; n++) mol[n].rv.x = -shearRate / 2.;
      for (n = nPolyMol + atWall; n < nMol; n++) mol[n].rv.x = shearRate / 2.;
   } else if (atTopWall != 0) {
      for (n = nPolyMol; n < nPolyMol + atBotWall; n++) mol[n].rv.x = shearRate / 2.;
      for (n = nPolyMol + atBotWall; n < nMol; n++) mol[n].rv.x = -shearRate / 2.;
   }
}

void ModifyWidth () {
   int at;
   
   DO_SURF {
      if (mol[nPolyMol + at].r.z > 0. && flatTopW == 0) {
         mol[nPolyMol + at].rv.z = shiftWidth / (stepEquil * deltaT);
      } else if (mol[nPolyMol + at].r.z > 0. && flatTopW == 1 && slope == 0) {
      } else if (mol[nPolyMol + at].r.z > 0. && flatTopW == 1 && slope == 1) {
         mol[nPolyMol + at].rv.z = - shiftWidth / (stepEquil * deltaT);
      } else if (mol[nPolyMol + at].r.z < 0.) {
         mol[nPolyMol + at].rv.z = - shiftWidth / (stepEquil * deltaT);
      } else {
      }
   }
   
   region.z += 2. * shiftWidth / stepEquil;
   
   if (stepCount == (int) (0.5*stepEquil)) DO_SURF mol[nPolyMol + at].rv.z = 0.;
}

/* SET INITIAL ACCELERATIONS FOR BEADS TO ZERO */
void InitAccels () {
   int n;
   
   DO_MOL V_ZERO (mol[n].ra);
}

/* ASSIGNING ATOMS TO CHAINS */
void AssignToChain () {
   int i, j, n;
   
   n = 0;
   for (i = 0; i < nChain; i++) {
      for (j = 0; j < chainLen; j++) {
         mol[n].inChain = i;
         ++n;
      }
   }
   for (n = nChain * chainLen; n < nMol; n++) mol[n].inChain = -1;
}

/* EVALUATION OF PROPERTIES */
void EvalProps () {
   real vv, vvMax;
   int n;
   
   V_ZERO (vSum);
   vvSum = 0.;
   vvMax = 0.;
   DO_POLY_MOL {
      V_V_ADD (vSum, mol[n].rv);
      vv = V_LEN_SQ (mol[n].rv);
      vvSum += vv;
      vvMax = MAX (vvMax, vv);
   }
   dispHi += sqrt (vvMax) * deltaT;
   if (dispHi > 0.5 * rNebrShell) nebrNow = 1;
   kinEnergy.val = 0.5 * vvSum / nPolyMol;
   totEnergy.val = kinEnergy.val + uSum / nPolyMol;
   pressure.val = density * (vvSum + virSum) / (nPolyMol * NDIM);	// not real pressure (doesn't work for inhomogeneous systems)
   
}

void AccumProps (int icode) {
   if (icode == 0) {
      PROP_ZERO (totEnergy);
      PROP_ZERO (kinEnergy);
      PROP_ZERO (pressure);
   } 
   else if (icode == 1) {
      PROP_ACCUM (totEnergy);
      PROP_ACCUM (kinEnergy);
      PROP_ACCUM (pressure);
   } 
   else if (icode == 2) {
      PROP_AVG (totEnergy, stepAvg);
      PROP_AVG (kinEnergy, stepAvg);
      PROP_AVG (pressure, stepAvg);
   }
}

void PrintSummary (FILE *fp, int SysType) {
   if (stepCount == stepAvg) {
      fprintf (fp,
               "%7d %8.3f %8.5f  %7.4f %6.4f  %7.4f %6.4f  %7.4f %6.4f  %7.5f ",
               stepCount, timeNow, V_C_SUM (vSum) / nPolyMol, PROP_EST (totEnergy),
               PROP_EST (kinEnergy), PROP_EST (pressure), U0_surf);
      switch (SysType) {
         case 0:
            fprintf (fp, "   film\n");
            break;
         case 1:
            fprintf (fp, "   flow\n");
            break;
         case 2:
            fprintf (fp, "droplet\n");
            break;
      }
   } else {
      fprintf (fp,
               "%7d %8.3f %8.5f  %7.4f %6.4f  %7.4f %6.4f  %7.4f %6.4f  %7.5f\n",
               stepCount, timeNow, V_C_SUM (vSum) / nPolyMol, PROP_EST (totEnergy),
               PROP_EST (kinEnergy), PROP_EST (pressure), U0_surf);
   }
   fflush (fp);
}

/* GRID CALCULATIONS: DENSITY, VELOCITIES etc. */
void EvalGridProfile () {
   int k, n, i, j, enter_droplet;
   real halfDen;
   
   halfDen = density / 2.;				// profile's line criteria
   
   /* set histogram values calcilated for coordinate grid to 0 */
   for (i = 0; i < sizeCoordGrid.x; i++) {
      for (k = 0; k < sizeCoordGrid.z; k++) {
         den2D[i][k] = 0.;
         denMap[i][k] = 0.;
      }
      profileZX[i] = 0.;
      profile_botZX[i] = 0.;
   }

   /* collect grid histograms for coordinate grid */
   for (i = 0; i < sizeCoordGrid.x; i++) {
      enter_droplet = 0;
      for (k = 0; k < sizeCoordGrid.z; k++) {
         /* collecting densities along y-axis [making 2D density map out of 3D] */
         for (j = 0; j < sizeCoordGrid.y; j++) {
            n = k * sizeCoordGrid.x * sizeCoordGrid.y + j * sizeCoordGrid.x + i;
            den2D[i][k] += histCoordGrid[0][n];
         }
         den2D[i][k] /= (V_PROD(region) / (sizeCoordGrid.x * sizeCoordGrid.z) );
         denMap[i][k] = den2D[i][k];
         
         /* applying profile's line criteria (\rho_liquid + \rho_vapor)/2 */
         if ( (enter_droplet == 0) && (den2D[i][k] >= halfDen) ) {
            profile_botZX[i] = (k + 0.5 - sizeCoordGrid.z / 2.) * region.z / sizeCoordGrid.z;
            enter_droplet = 1;
         } else if ( (enter_droplet == 1) && (profile_botZX[i] < 0.) && 	// avoiding all "flying" chains in the region z > 0
                    (den2D[i][k-1] >= halfDen) && (den2D[i][k] < halfDen) ) {
            profileZX[i] = (k + 0.5 - sizeCoordGrid.z / 2.) * region.z / sizeCoordGrid.z;
         }
      }
   }
   
   
   /* set histogram values calcilated for velocity grid to 0 */
   for (i = 0; i < sizeVelGrid.x; i++) {
      for (k = 0; k < sizeVelGrid.z; k++) {
         vel2D_x[i][k] = 0.;
         vel2D_z[i][k] = 0.;
         vel2D_xx[i][k] = 0.;
         vel2D_zz[i][k] = 0.;
         dropVel2D[i][k] = 0.;
      }
   }
   
   /* collect grid histograms for velocity and energy grids */
   for (i = 0; i < sizeVelGrid.x; i++) {
      for (k = 0; k < sizeVelGrid.z; k++) {
         /* collecting velocities along y-axis [making 2D map out of 3D] */
         for (j = 0; j < sizeVelGrid.y; j++) {
            n = k * sizeVelGrid.x * sizeVelGrid.y + j * sizeVelGrid.x + i;
            //				nMolVel2D[i][k] += histVelGrid[0][n];
            vel2D_x[i][k] += histVelGrid[2][n];
            vel2D_z[i][k] += histVelGrid[4][n];
            vel2D_xx[i][k] += histVelGrid[5][n];
            vel2D_zz[i][k] += histVelGrid[6][n];
            dropVel2D[i][k] += histVelGrid[1][n];
			}
         
         vel2D_x[i][k] /= sizeVelGrid.y;
         vel2D_z[i][k] /= sizeVelGrid.y;
         vel2D_xx[i][k] /= sizeVelGrid.y;
         vel2D_zz[i][k] /= sizeVelGrid.y;
         dropVel2D[i][k] /= sizeVelGrid.y;
         
      }
   }
}

void EvalControlValues (int make_controlFile, int opCode) {
   if (make_controlFile == 1) {
      VecR dr_ree;
      VecR ijDist;
      int i, j;
      
      if (opCode == 0) {
         PROP_ZERO (ReeAv);
         PROP_ZERO (ReeAv_x);
         PROP_ZERO (ReeAv_y);
         PROP_ZERO (ReeAv_z);
         bondLength = 0.;
      } else if (opCode == 1) {
         for (i = 0; i < nPolyMol; i += chainLen) {
            V_ZERO (dr_ree);
            for (j = chainLen - 1; j > 0; j--) {
               V_SUB (ijDist, mol[i + j].r, mol[i + j - 1].r);
               V_WRAP_ALL (ijDist);
               bondLength += V_LEN_SQ (ijDist);
               V_V_ADD (dr_ree, ijDist);
            }
            ReeAv.val = V_LEN (dr_ree);
            ReeAv_x.val = sqrt(SQR(dr_ree.x));
            ReeAv_y.val = sqrt(SQR(dr_ree.y));
            ReeAv_z.val = sqrt(SQR(dr_ree.z));
            PROP_ACCUM (ReeAv);
            PROP_ACCUM (ReeAv_x);
            PROP_ACCUM (ReeAv_y);
            PROP_ACCUM (ReeAv_z);
         }
      } else if (opCode == 2) {
         bondLength /= (nChain * (chainLen - 1));	// averaging over number of bonds
         bondLength /= (limitConVal / stepConVal);	// averaging over number of measurements
         PROP_AVG (ReeAv, (nChain * limitConVal / stepConVal));
         PROP_AVG (ReeAv_x, (nChain * limitConVal / stepConVal));
         PROP_AVG (ReeAv_y, (nChain * limitConVal / stepConVal));
         PROP_AVG (ReeAv_z, (nChain * limitConVal / stepConVal));
         
      }
   }
}

void WrapChainsIntoCell () {
   int i, j, n;
   VecR rMean, ijDist;
   
   DO_MOL V_COPY (mol[n].r_vtf, mol[n].r);
   
   /* COLLECTING CHAIN BEADS TOGETHER */	
   for (i = 0; i < nPolyMol; i += chainLen) {
      for (j = chainLen - 1; j > 0; j--) {
         V_SUB (ijDist, mol[i + j].r_vtf, mol[i + j - 1].r_vtf);
         
         if (ijDist.x >= 0.5 * region.x) {
            mol[i + j - 1].r_vtf.x += region.x;
         } else if (ijDist.x < -0.5 * region.x) {
            mol[i + j - 1].r_vtf.x -= region.x;
         }
         
         if (ijDist.y >= 0.5 * region.y) {
            mol[i + j - 1].r_vtf.y += region.y;
         } else if (ijDist.y < -0.5 * region.y) {
            mol[i + j - 1].r_vtf.y -= region.y;
         }
      }
   }
   
   /* TRANSLATION OF CHAINS DEPENDING ON THEIR CENTER OF MASS */
   for (i = 0; i < nPolyMol; i += chainLen) {
      V_ZERO (rMean);
      for (j = chainLen - 1; j >= 0; j--) {
         V_V_ADD (rMean, mol[i + j].r_vtf);
      }
      V_SCALE (rMean, 1. / chainLen);
      if (rMean.x >= 0.5 * region.x) {
         for (j = chainLen - 1; j >= 0; j--) mol[i + j].r_vtf.x -= region.x;
      }
      else if (rMean.x < -0.5 * region.x) {
         for (j = chainLen - 1; j >= 0; j--) mol[i + j].r_vtf.x += region.x;
      }
      if (rMean.y >= 0.5 * region.y) {
         for (j = chainLen - 1; j >= 0; j--) mol[i + j].r_vtf.y -= region.y;
      }
      else if (rMean.y < -0.5 * region.y) {
         for (j = chainLen - 1; j >= 0; j--) mol[i + j].r_vtf.y += region.y;
      }
      // box doesn't have a periodicity in z-direction, so
      /*		if (rMean.z >= 0.5 * region.z) {
       for (j = chainLen - 1; j >= 0; j--) mol[i + j].r_vtf.z -= region.z;
       }
       else if (rMean.z < -0.5 * region.z) {
       for (j = chainLen - 1; j >= 0; j--) mol[i + j].r_vtf.z += region.z;
       }*/
   }
}

#include "subroutines/forcePT_calc.c"
#include "subroutines/chainForcePT_calc.c"
#include "subroutines/den_profile.c"
#include "subroutines/pr_profile.c"
#include "subroutines/grid_profiles.c"
#include "subroutines/vel_profile_slab_geometry.c"
#include "subroutines/itoa.c"
#include "subroutines/reverse.c"

#include "in_functions/out_file.c"
#include "in_functions/in_rand.c"
#include "in_functions/in_errexit.c"
#include "in_functions/in_namelist.c"
