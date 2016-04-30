void EvalDenDist (int opCode) {
	real vv;
	int j, m, n;

	if (opCode == 0) {
		for (j = 0; j < numSlabs; j++) histDen[j] = 0.;
		instTemp = 0.;
	} else if (opCode == 1) { 
	/* collecting data for density profile and instant temperature */
		DO_POLY_MOL {
			if (mol[n].r.z < 0)
				j = halfSl - (int) (-mol[n].r.z / hSlab) - 1;
			else
				j = halfSl + (int) (mol[n].r.z / hSlab);
			histDen[j]++;
			
			vv = V_LEN_SQ (mol[n].rv);
			instTemp += vv;
		}       
	} else if (opCode == 2) {
	/* averaging pressure profile values */
			for (j = 0; j < numSlabs; j++) histDen[j] /= (hSlab * region.x * region.y * limitConVal / stepConVal);
			instTemp /= (3. * nPolyMol * limitConVal / stepConVal);
	}
}

/* PRINTING DENSITY PROFILES */
void PrintDenDist (FILE *file_den_hist) {
        real denBin, denXBin;
        int n;
        extern char den_hist_fullname[];

        file_den_hist = fopen (den_hist_fullname, "a");
	
	if (file_den_hist != NULL) {
	        fprintf (file_den_hist, "#density profile (step %d) \n", stepCount);
	        for (n = 0; n < numSlabs; n++) {
	                denBin = (n + 0.5 - halfSl) * region.z / numSlabs;
	                fprintf (file_den_hist, "%8.3f %8.3f\n", denBin, histDen[n]);
	        }
	}
        fclose (file_den_hist);
}
