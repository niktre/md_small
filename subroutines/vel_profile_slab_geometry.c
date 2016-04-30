/* VELOCITY PROFILES FOR SLAB GEOMETRY (SUITABLE FOR FLOW STUDIES) */
void EvalVelDist (int opCode) {
	extern int halfSl;
	int j, n;
	real vv, vvx;

	if (opCode == 0) {	
		for (j = 0; j < numSlabs; j++) {
			histAvDen[j] = 0.;
			histAvVx[j] = 0.;
			histAvTx[j] = 0.;
			histAvT[j] = 0.;

			histVx[j] = 0.;
			histTx[j] = 0.;
			histT[j] = 0.;
		}
	} else if (opCode == 1) {
		for (j = 0; j < numSlabs; j++) {
			snapDen[j] = 0.;
			snapVx[j] = 0.;
			snapTx[j] = 0.;
			snapT[j] = 0.;
		}

		DO_POLY_MOL {
			if (mol[n].r.z < 0.)
				j = halfSl - (int) (-mol[n].r.z / hSlab) - 1;
			else
				j = halfSl + (int) (mol[n].r.z / hSlab);
			
			++snapDen[j];
			++histAvDen[j];
			vvx = SQR (mol[n].rv.x);
			vv = V_LEN_SQ (mol[n].rv);

			snapVx[j] += mol[n].rv.x;
			snapTx[j] += vvx;
			snapT[j] += vv;

			histAvVx[j] += mol[n].rv.x;
			histAvTx[j] += vvx;
			histAvT[j] += vv;
		}

		/* add to statistics averaged Vx, Tx and T */
		for (j = 0; j < numSlabs; j++) {
			if (snapDen[j] > 0.) {
				histVx[j] += (snapVx[j] / snapDen[j]);
				histTx[j] += (snapTx[j] / snapDen[j]);
				histT[j] += (snapT[j] / snapDen[j]);
			}
		}

	} else if (opCode == 2) {
		for (j = 0; j < numSlabs; j++) {
			if (histAvDen[j] > 0.) {
				histAvVx[j] /= histAvDen[j];
				histAvTx[j] /= histAvDen[j];
				histAvT[j] /= histAvDen[j];
			}
			histVx[j] /= (limitConVal / stepConVal);
			histTx[j] /= (limitConVal / stepConVal);
			histT[j] /= (limitConVal / stepConVal);
		}
	}
}

/* PRINTING VELOCITY PROFILES FOR SLAB GEOMETRY (SUITABLE FOR FLOW STUDIES) */
void PrintVelDist (FILE *file_hist) {
        real vBin;
        int j;
        extern char vel_hist_fullname[80];

        file_hist = fopen (vel_hist_fullname, "a");

	if (file_hist != NULL) {
	        fprintf (file_hist, "#vdist (step %d) with slabs method\n", stepCount);
	        fprintf (file_hist, "#pnBin  MarcusVx MarcusVx^2 MarcusV^2     AvVx    AvVx^2     AvV^2\n");
	        for (j = 0; j < numSlabs; j++) {
	                vBin = (j + 0.5 - halfSl) * region.z / numSlabs;
	                fprintf (file_hist, "%7.3f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", 
					vBin, histAvVx[j], histAvTx[j], histAvT[j], histVx[j], histTx[j], histT[j]);
	        }
	}

        fclose (file_hist);
}
