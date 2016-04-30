void CoordGridAverage (int opCode) {
        VecR invWid, rs;
        VecI cc;
        int c, hSize, i, j, k, l, n, m_cm, half_x, half_y;
        real vv;
	int part_p_coord, part_m_coord;
	real rx_av_p_coord, rx_av_m_coord;

        half_x = (int) (sizeCoordGrid.x / 2);    // the number of the cell in the middle of region.x
        half_y = (int) (sizeCoordGrid.y / 2);    // the number of the cell in the middle of region.y
        hSize = V_PROD (sizeCoordGrid);          // number of grid cells in the system

        /* SETTING INITIAL HISTOGRAMS VALUES TO ZERO */
        if (opCode == 0) {
                for (j = 0; j < NHIST_C; j++) {
                        for (n = 0; n < hSize; n++) {
                                histCoordGrid[j][n] = 0.;
                        }
                }
        } else if (opCode == 1) {
                /* setting initial snapshots values to zero */
                for (j = 0; j < NHIST_C; j++) {
                        for (n = 0; n < hSize; n++) snapCoordGrid[j][n] = 0.;
                }

/* if we measure center of mass position in every slab in y-direction, uncomment */
                /* setting two parameters to zero */
/*		for (j = 0; j < sizeCoordGrid.y; j++) {
                        part_y_p_coord[j] = 0;
                        rx_p_coord[j] = 0.;
                        part_y_m_coord[j] = 0;
                        rx_m_coord[j] = 0.;
		}
*/

		/* global CM parameters */
		V_SET_ALL (globCM, 0.);
		part_p_coord = part_m_coord = 0;
		rx_av_p_coord = rx_av_m_coord = 0.;

                V_DIV (invWid, sizeCoordGrid, region);
                DO_POLY_MOL {
                        V_S_ADD (rs, mol[n].r, 0.5, region);    //shifts coordinates to "positive" region (from 0 to L)
                        V_MUL (cc, rs, invWid);
                        c = V_LINEAR (cc, sizeCoordGrid);
                        ++ snapCoordGrid[0][c];
/*			if (sizeCoordGrid.y == 1) {
				j = 0;
			} else {
                        	if (mol[n].r.y < 0.)
                        	        j = half_y - (int) (-mol[n].r.y * sizeCoordGrid.y / region.y) - 1;
                        	else
                        	        j = half_y + (int) (mol[n].r.y * sizeCoordGrid.y / region.y);
			}
*/			/* y and z cannot intersect the borders of the simulation region, x will be looked individually(!!!) */
			globCM.y += mol[n].r.y;
			globCM.z += mol[n].r.z;
			/* droplet can intersect x-border. calculate the CoM in every half of region.x */
			if (mol[n].r.x < 0.) {
				++ part_m_coord;
	                        rx_av_m_coord += mol[n].r.x;
//	                        ++ part_y_m_coord[j];
//	                        rx_m_coord[j] += mol[n].r.x;
			} else {
				++ part_p_coord;
	                        rx_av_p_coord += mol[n].r.x;
//	                        ++ part_y_p_coord[j];
//	                        rx_p_coord[j] += mol[n].r.x;
			}
                }

                /* making snapshot averages in coord-based-histograms */
                for (n = 0; n < hSize; n++) {
                        if (snapCoordGrid[0][n] > 0. && NHIST_C > 1) {
                                for (j = 1; j < NHIST_C; j++) snapCoordGrid[j][n] /= snapCoordGrid[0][n];
                        }
                }

		/* add rx from "plus" and "minus" regions to get real rx of the center of mass, scaling come next */
		if (part_m_coord > 0) rx_av_m_coord /= part_m_coord;
		if (part_p_coord > 0) rx_av_p_coord /= part_p_coord;
			
		if ((rx_av_p_coord - rx_av_m_coord) < region.x / 2.) {
			globCM.x = (part_p_coord*rx_av_p_coord + part_m_coord*rx_av_m_coord);
		} else {
			if ((rx_av_p_coord + rx_av_m_coord) < 0.) {
				globCM.x = (part_p_coord*rx_av_p_coord + part_m_coord*(rx_av_m_coord + region.x));
			} else {
				globCM.x = (part_p_coord*(rx_av_p_coord - region.x) + part_m_coord*rx_av_m_coord);
			}
		}

		/* scale CM parameters */
		V_SCALE (globCM, 1. / nPolyMol);

		CentrifyHistogram(sizeCoordGrid, half_x, NHIST_C, snapCoordGrid, histCoordGrid);

        } else if (opCode == 2) {
                for (j = 0; j < NHIST_C; j++) {
                        for (n = 0; n < hSize; n++) histCoordGrid[j][n] /= (limitConVal / stepConVal);
                }
        } /* hier haben wir histCoordGrid mit den Mittelwerten der Okkupazion der Zelle */
}

/* CENTRIFYING THE DIAGRAMM WITH RESPECT TO THE CENTER OF MASS */
void CentrifyHistogram (VecLI sizeHistGrid, int half_x, int NHIST, real **snapGrid, real **histGrid) {
	int i, j, k, l, n, m_cm;

		/* find out the layer where the center of mass in x-direction is */
                if (globCM.x < 0.)
                        m_cm = half_x - (int) (-globCM.x * sizeHistGrid.x / region.x) - 1;
                else
                        m_cm = half_x + (int) (globCM.x * sizeHistGrid.x / region.x);

                /* centrifying of the histograms */
                for (j = 0; j < sizeHistGrid.y; j++) {
/* unkommentieren, wenn Schwerpunkt in jeder y-Schichte ausgerechnet ist */
/*			rx_av = 0.;
			if (part_y_m[j] > 0) rx_m[j] /= part_y_m[j];
                        if (part_y_p[j] > 0) rx_p[j] /= part_y_p[j];
		
*/

			/* die Funktion wurde zur CoordGridAverage übergegeben */
			/* add rx from "plus" and "minus" regions to get real rx of the center of mass */	
/*			if ((rx_p[j] - rx_m[j]) <= region.x / 2.) {
				rx_av = (part_y_p[j] * rx_p[j] + part_y_m[j] * rx_m[j]) / (part_y_m[j] + part_y_p[j]);
			} else {
				if ((rx_p[j] + rx_m[j]) < 0.) {
					rx_av = (part_y_p[j] * rx_p[j] + part_y_m[j] * (rx_m[j] + region.x)) / (part_y_m[j] + part_y_p[j]);
				} else {
					rx_av = (part_y_p[j] * (rx_p[j] - region.x) + part_y_m[j] * rx_m[j]) / (part_y_m[j] + part_y_p[j]);
				}
			}
*/
                        if (m_cm > half_x) {
                                for (k = 0; k < sizeHistGrid.z; k++) {
                                        for (i = m_cm - half_x; i < sizeHistGrid.x; i++) {
                                                n = k * sizeHistGrid.x * sizeHistGrid.y + j * sizeHistGrid.x + i;
                                                for (l = 0; l < NHIST; l++) histGrid[l][n - m_cm + half_x] += snapGrid[l][n];
                                        }
                                        for (i = 0; i < m_cm - half_x; i++) {
                                                n = k * sizeHistGrid.x * sizeHistGrid.y + j * sizeHistGrid.x + i;
                                                for (l = 0; l < NHIST; l++) histGrid[l][sizeHistGrid.x + n - m_cm + half_x] += snapGrid[l][n];
                                        }
                                }
                        } else if (m_cm < half_x) {
                                for (k = 0; k < sizeHistGrid.z; k++) {
                                        for (i = 0; i < sizeHistGrid.x - half_x + m_cm; i++) {
                                                n = k * sizeHistGrid.x * sizeHistGrid.y + j * sizeHistGrid.x + i;
                                                for (l = 0; l < NHIST; l++) histGrid[l][n + half_x - m_cm] += snapGrid[l][n];
                                        }
                                        for (i = sizeHistGrid.x - half_x + m_cm; i < sizeHistGrid.x; i++) {
                                                n = k * sizeHistGrid.x * sizeHistGrid.y + j * sizeHistGrid.x + i;
                                                for (l = 0; l < NHIST; l++) histGrid[l][n - sizeHistGrid.x + half_x - m_cm] += snapGrid[l][n];
                                        }
                                }
                        } else if (m_cm == half_x) {
                                for (k = 0; k < sizeHistGrid.z; k++) {
                                        for (i = 0; i < sizeHistGrid.x; i++) {
                                                n = k * sizeHistGrid.x * sizeHistGrid.y + j * sizeHistGrid.x + i;
                                                for (l = 0; l < NHIST; l++) histGrid[l][n] += snapGrid[l][n];
                                        }
                                }
                        }
                }
}

void VelGridAverage (int opCode) {
        VecR invWid, rs, va;
        VecI cc, mitSlab;
        int c, hSize, i, j, k, l, n;
        real vv;

	V_SET (mitSlab, (int) (sizeVelGrid.x / 2), (int) (sizeVelGrid.y / 2), (int) (sizeVelGrid.z / 2)); // die mittlere Schichte
        hSize = V_PROD (sizeVelGrid);          // number of grid cells in the system

        /* SETTING INITIAL HISTOGRAMS VALUES TO ZERO */
        if (opCode == 0) {
		V_SET_ALL (fastCMvel, 0.);
                for (j = 0; j < NHIST_V; j++) {
                        for (n = 0; n < hSize; n++) {
                                histVelGrid[j][n] = 0.;
                        }
                }
        } else if (opCode == 1) {
                /* setting initial snapshots values to zero */
                for (j = 0; j < NHIST_V; j++) {
                        for (n = 0; n < hSize; n++) snapVelGrid[j][n] = 0.;
                }

                /* setting two parameters to zero */
/*                for (j = 0; j < sizeVelGrid.y; j++) {
                        part_y_p_vel[j] = 0;
                        rx_p_vel[j] = 0.;
                        part_y_m_vel[j] = 0;
                        rx_m_vel[j] = 0.;
                }
*/		
		/* set global CoM and it's velocity to zero and it's cell to -1 */
//		V_SET_ALL (VglobCM, 0.);
		V_SET_ALL (snapCMvel, 0.);
		V_SET_ALL (globCMslab, 0);
		n_globCM = -1;

                V_DIV (invWid, sizeVelGrid, region);
                DO_POLY_MOL {
                        V_S_ADD (rs, mol[n].r, 0.5, region);    //shifts coordinates to "positive" region (from 0 to L)
                        V_MUL (cc, rs, invWid);
                        c = V_LINEAR (cc, sizeVelGrid);
                        ++ snapVelGrid[0][c];
                        vv = V_LEN_SQ (mol[n].rv);
                        snapVelGrid[1][c] += vv;
                        snapVelGrid[2][c] += mol[n].rv.x;
                        snapVelGrid[3][c] += mol[n].rv.y;
                        snapVelGrid[4][c] += mol[n].rv.z;
                        snapVelGrid[5][c] += SQR(mol[n].rv.x);
                        snapVelGrid[6][c] += SQR(mol[n].rv.z);
/*			if (sizeVelGrid.y == 1) {
				j = 0;
			} else {
	                        if (mol[n].r.y < 0.)
	                                j = mitSlab.y - (int) (-mol[n].r.y * sizeVelGrid.y / region.y) - 1;
	                        else
					j = mitSlab.y + (int) (mol[n].r.y * sizeVelGrid.y / region.y);
			}
*/			/* droplet can intersect x-border. calculate the CoM in every half of region.x */
/*			if (mol[n].r.x < 0.) {
	                        ++ part_y_m_vel[j];
	                        rx_m_vel[j] += mol[n].r.x;
			} else {
	                        ++ part_y_p_vel[j];
	                        rx_p_vel[j] += mol[n].r.x;
			}
*/			V_V_ADD (snapCMvel, mol[n].rv);
                }
		
                /* making snapshot averages in velocity-histograms */
/*		for (n = 0; n < hSize; n++) {
                        if (snapVelGrid[0][n] > 0. && NHIST_V > 1) {
                                for (j = 1; j < NHIST_V; j++) snapVelGrid[j][n] /= snapVelGrid[0][n];
                        }
		 }
*/	
		/* add rx from "plus" and "minus" regions to get real rx of the center of mass, scaling come next */
/*		if (part_m_vel > 0) rx_av_m_vel /= part_m_vel;
		if (part_p_vel > 0) rx_av_p_vel /= part_p_vel;
			
		if ((rx_av_p_vel - rx_av_m_vel) < region.x / 2.) {
			globCM.x = (part_p_vel*rx_av_p_vel + part_m_vel*rx_av_m_vel);
		} else {
			if ((rx_av_p_vel + rx_av_m_vel) < 0.) {
				globCM.x = (part_p_vel*rx_av_p_vel + part_m_vel*(rx_av_m_vel + region.x));
			} else {
				globCM.x = (part_p_vel*(rx_av_p_vel - region.x) + part_m_vel*rx_av_m_vel);
			}
		}
*/
		/* scale CM parameters */
//		V_SCALE (globCM, 1. / nPolyMol);
		V_SCALE (snapCMvel, 1. / nPolyMol);
		V_V_ADD (fastCMvel, snapCMvel);
				
		/* finding global Centre of Mass velocity and substracting it from other grid cells*/
		FIND_SLAB (globCMslab.x, sizeVelGrid.x, globCM.x, region.x);
		if (sizeVelGrid.y == 1) globCMslab.y = 0;
		else FIND_SLAB (globCMslab.y, sizeVelGrid.y, globCM.y, region.y);
		FIND_SLAB (globCMslab.z, sizeVelGrid.z, globCM.z, region.z);

		/* straighten CoM from 3D representation to 1D */
		n_globCM = globCMslab.z * sizeVelGrid.x * sizeVelGrid.y + globCMslab.y * sizeVelGrid.x + globCMslab.x;
		if (n_globCM == -1) {
			ErrExit(ERR_CENT_O_MASS);
		} else {
//			V_SET (VglobCM, snapVelGrid[2][n_globCM], snapVelGrid[3][n_globCM], snapVelGrid[4][n_globCM]);
		}

		/* substracting Global Centre of Mass velocity from all other grid cells */
/*		for (n = 0; n < hSize; n++) {
*/			/* appplying substraction only to meaningful cells */
/*			if (snapVelGrid[0][n] > 0. && NHIST_V > 1) {
				snapVelGrid[2][n] -= fastCMvel.x;
				snapVelGrid[3][n] -= fastCMvel.y;
				snapVelGrid[4][n] -= fastCMvel.z;
			}
		}
*/
		/* centrifying diagram with respect to the local centre of mass of every y-layer*/
		CentrifyHistogram(sizeVelGrid, mitSlab.x, NHIST_V, snapVelGrid, histVelGrid);
	
        } else if (opCode == 2) {
		fastCMvel.x /= (limitConVal / stepConVal);
		fastCMvel.y /= (limitConVal / stepConVal);
		fastCMvel.z /= (limitConVal / stepConVal);

		for (n = 0; n < hSize; n++) {
//			for (n = 0; n < hSize; n++) histVelGrid[j][n] /= (limitConVal / stepConVal);
			if (histVelGrid[0][n] > 0.) {
				for (j = 1; j < NHIST_V; j++) histVelGrid[j][n] /= histVelGrid[0][n];
				V_SET (va, histVelGrid[2][n], histVelGrid[3][n], histVelGrid[4][n]);
				histVelGrid[1][n] = (histVelGrid[1][n] - V_LEN_SQ (va)) / NDIM; 
				histVelGrid[5][n] = histVelGrid[5][n] - SQR (histVelGrid[2][n]);
				histVelGrid[6][n] = histVelGrid[6][n] - SQR (histVelGrid[4][n]);
			}
                }
	}
}

void EnGridAverage (int opCode) {
        VecR invWid, rs;
        VecI cc, mitSlab;
        int c, hSize, j, n;
	
	V_SET (mitSlab, (int) (sizeVelGrid.x / 2), (int) (sizeVelGrid.y / 2), (int) (sizeVelGrid.z / 2)); // die mittlere Schichte
	hSize = V_PROD (sizeVelGrid);
	
        if (opCode == 0) {
                for (j = 0; j < NHIST_E; j++) {
                        for (n = 0; n < hSize; n++) {
                                histEnGrid[j][n] = 0.;
                        }
                }
	} else if (opCode == 1) {

                V_DIV (invWid, sizeVelGrid, region);
                DO_POLY_MOL {
                        V_S_ADD (rs, mol[n].r, 0.5, region);    //shifts coordinates to "positive" region (from 0 to L)
                        V_MUL (cc, rs, invWid);
                        c = V_LINEAR (cc, sizeVelGrid);
                        ++ snapEnGrid[0][c];
/*			if (sizeVelGrid.y == 1) {
				j = 0;
			} else {
	                        if (mol[n].r.y < 0.)
	                                j = mitSlab.y - (int) (-mol[n].r.y * sizeVelGrid.y / region.y) - 1;
	                        else
	                                j = mitSlab.y + (int) (mol[n].r.y * sizeVelGrid.y / region.y);
			}
*/			/* droplet can intersect x-border. calculate the CoM in every half of region.x */
/*			if (mol[n].r.x < 0.) {
	                        ++ part_y_m_en[j];
	                        rx_m_en[j] += mol[n].r.x;
			} else {
	                        ++ part_y_p_en[j];
	                        rx_p_en[j] += mol[n].r.x;
			}
*/		}

		/* making snapshot averages in enegry-based-histograms */
		for (n = 0; n < hSize; n++) {
                        if (snapEnGrid[0][n] > 0. && NHIST_E > 1) {
                                for (j = 1; j < NHIST_E; j++) snapEnGrid[j][n] /= snapEnGrid[0][n];
                        }
		}

		CentrifyHistogram(sizeVelGrid, mitSlab.x, NHIST_E, snapEnGrid, histEnGrid);

	} else if (opCode == 2) {
		for (j = 0; j < NHIST_E; j++) {
                        for (n = 0; n < hSize; n++) histEnGrid[j][n] /= (limitConVal / stepConVal);
                }
	}


}

void PrintPosVelCM (FILE *file_cm) {
	extern int addCMStep;
	extern int n_globCM;
	extern VecI globCMslab;
	
	extern char parameter_filepath[];
	extern char foldername[];
	extern char pos_vel_cm_filename[];
	extern char pos_vel_cm_fullname[];

	strcpy (pos_vel_cm_fullname,parameter_filepath);
	strcat (pos_vel_cm_fullname,foldername);
	strcat (pos_vel_cm_fullname,pos_vel_cm_filename);

	/* PRINTING CM-PARAMETERS */
	file_cm = fopen (pos_vel_cm_fullname, "a");
	if (file_cm != NULL) {
		if (stepCount == (stepEquil + stepConVal) ) {
			fprintf (file_cm, "#parameters of CM from grid cells method\n");
			fprintf (file_cm, "#   step    time   epsW    x_pos    y_pos    z_pos     vx       vy       vz       v^2  ");
			fprintf (file_cm, "globCMslab        trajType\n");
		}
		fprintf (file_cm, "%8d %8.4f %5.2f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %4d %4d %4d %7d\n", 
				stepCount, timeNow, epsWall, globCM.x, globCM.y, globCM.z, snapCMvel.x, snapCMvel.y, snapCMvel.z, V_LEN_SQ(snapCMvel),
				globCMslab.x, globCMslab.y, globCMslab.z, trajType);
	}

	fclose (file_cm);

}

/* PRINTING DENSITIES, DENSITY MAPS, BOUNDARIES AND VELOCITIES FOR GRID GEOMETRY (SUITABLE FOR DROPLETS) */
void PrintGridProfile (FILE *file_grid, FILE *file_denX_hist, FILE *file_denMap, FILE *file_dropVel2D) {
        real zVal, xVal;
        int n, k, j;
        extern char file_grid_fullname[];
        extern char file_denX_fullname[];
        extern char denMap_fullname[];
        extern char dropVel2D_fullname[];
	extern VecR fastCMvel; 

	int hSize;
	
	hSize = V_PROD (sizeVelGrid);

        /* PRINTING GRID VELOCITY PROFILE AT CONSTANT Z-LEVEL 
	*	(similar to slab geometry, works not so good for droplets
	*	[uncomment if needed]
        file_grid = fopen (file_grid_fullname, "a");
        fprintf (file_grid, "#Vx(z) profile after (%d) steps with cells method\n", stepCount);
        fprintf (file_grid, "# zVal       Vx        T\n");
        for (n = 0; n < sizeVelGrid.z; n++) {
                zVal = (n + 0.5 - sizeVelGrid.z / 2.) * region.z / sizeVelGrid.z;
                fprintf (file_grid, "%7.3f %9.5f %9.5f\n", zVal, profileVx[n], profileT[n]);
        }
        fflush (file_grid);
	*/

        /* PRINTING DROPLET'S BOUNDARIES IN X-DIRECTION */
        file_denX_hist = fopen (denX_hist_fullname, "a");
	if (file_denX_hist != NULL) {
	        fprintf (file_denX_hist, "#densityX profile (step %d) with cells method. Coord Grid is %ld %ld %ld \n",
					 stepCount, sizeCoordGrid.x, sizeCoordGrid.y, sizeCoordGrid.z);
	        fprintf (file_denX_hist, "#   x     z_bot    z_top\n");
       		for (n = 0; n < sizeCoordGrid.x; n++) {
                	xVal = (n + 0.5 - sizeCoordGrid.x / 2.) * region.x / sizeCoordGrid.x;
               		fprintf (file_denX_hist, "%8.4f %8.4f %8.4f\n", xVal, profile_botZX[n], profileZX[n]);
        	}
	}

        fclose (file_denX_hist);

        /* PRINTING 2D DENSITY MAP */
        file_denMap = fopen (denMap_fullname, "a");
	if (file_denMap != NULL) {
	        fprintf (file_denMap, "#density map 2D (step %d) with cells method\n", stepCount);
	        fprintf (file_denMap, "#   x        z       den  \n");
       		for (n = 0; n < sizeCoordGrid.x; n++) {
       	        	for (k = 0; k < sizeCoordGrid.z; k++) {
                        	xVal = (n + 0.5 - sizeCoordGrid.x / 2.) * region.x / sizeCoordGrid.x;
                        	zVal = (k + 0.5 - sizeCoordGrid.z / 2.) * region.z / sizeCoordGrid.z;
				if (denMap[n][k] > 0.) {
                        		fprintf (file_denMap, "%8.4f %8.4f %8.4f\n", xVal, zVal, denMap[n][k]);
				}
				else {
                        		fprintf (file_denMap, "%8.4f %8.4f %8.4f\n", xVal, zVal, 0.);
				}
                	}
        	}
	}

        fclose (file_denMap);

        /* PRINTING VELOCITIES VALUES INSIDE THE DROPLET and WORK OF THERMOSTAT AND POT ENERGY */
        file_dropVel2D = fopen (dropVel2D_fullname, "a");
	if (file_dropVel2D != NULL) {
	        fprintf (file_dropVel2D, "#velocity droplet map 2D (step %d) with cells method. Vel Grid is %ld %ld %ld; fastCMvel %8.4f %8.4f %8.4f \n", 
						stepCount, sizeVelGrid.x, sizeVelGrid.y, sizeVelGrid.z, fastCMvel.x, fastCMvel.y, fastCMvel.z);
		fprintf (file_dropVel2D, "#   x        z       vx       vz       vx^2     vz^2   kinEn2D  potEnPP  potEnPS  totEn2D  thermTaG\n");
	
		for (n = 0; n < sizeVelGrid.x; n++) {
			for (k = 0; k < sizeVelGrid.z; k++) {
				xVal = (n + 0.5 - sizeVelGrid.x / 2.) * region.x / sizeVelGrid.x;
				zVal = (k + 0.5 - sizeVelGrid.z / 2.) * region.z / sizeVelGrid.z;
				fprintf (file_dropVel2D, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
					xVal, zVal, vel2D_x[n][k], vel2D_z[n][k], vel2D_xx[n][k], vel2D_zz[n][k], 
					dropVel2D[n][k], potEnPP[n][k], potEnPS[n][k], dropVel2D[n][k] + potEnPP[n][k] + potEnPS[n][k], thermTaG2D[n][k]);
       	  	       }
		}
	}

        fclose (file_dropVel2D);
}

