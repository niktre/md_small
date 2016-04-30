void EvalPrDist (int opCode) {
        int j;

        if (opCode == 0) {
                for (j = 0; j < numSlabs; j++) {
                        histPn[j] = 0.;
                        histPnW[j] = 0.;
                        histPt[j] = 0.;
                        histPtW[j] = 0.;
                        histPx[j] = 0.;
                        histPy[j] = 0.;
                }
        } else if (opCode == 2) {
                for (j = 0; j < numSlabs; j++) {
                        histPn[j] /= (hSlab * region.x * region.y * limitConVal / stepConVal);
                        histPnW[j] /= (hSlab * region.x * region.y * limitConVal / stepConVal);
                        histPx[j] /= (hSlab * region.x * region.y * limitConVal / stepConVal);
                        histPy[j] /= (hSlab * region.x * region.y * limitConVal / stepConVal);
                        histPt[j] /= (2. * hSlab * region.x * region.y * limitConVal / stepConVal);
                        histPtW[j] /= (2. * hSlab * region.x * region.y * limitConVal / stepConVal);
                }
        }
}

/* PRINTING PRESSURE TENSOR COMPONENTS AND CALCULATION OF SURFACE TENSION*/
void PrintPrDist (FILE *file_pn_hist) {
        real pnBin;
        int j;
        extern char pn_hist_fullname[];

        file_pn_hist = fopen (pn_hist_fullname, "a");

	if (file_pn_hist != NULL) {
	        fprintf (file_pn_hist, "#normal (x,y) pressure tensor profile (step %d) temp (%8.3f)\n", stepCount, instTemp);
	        fprintf (file_pn_hist, "#pnBin      first     secPn     walPn     secPx     secPy     secPt    secPtW     P_x       P_y       P_n       P_t        dP      surTen\n");
	        for (j = 0; j < numSlabs; j++) {
	                pnBin = (j + 0.5 - halfSl) * region.z / numSlabs;
	                fprintf (file_pn_hist, "%7.3f %9.3f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	                                        pnBin, histDen[j] * instTemp, histPn[j], histPnW[j], histPx[j], histPy[j], histPt[j], histPtW[j],
       	                                 histDen[j] * instTemp + histPx[j], histDen[j] * instTemp + histPy[j], histDen[j] * instTemp + histPn[j] + histPnW[j],
       	                                 histDen[j] * instTemp + histPt[j], histPn[j] + histPnW[j] - histPt[j], hSlab * (histPn[j] + histPnW[j] - histPt[j]));
        	}
	}

        fclose (file_pn_hist);
}
