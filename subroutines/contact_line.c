void EvalContactRegion (int opCode) {
	VecR invWid, rs;
	VecI cc;
	real crSize;
	int c, i, j, k, n;
	int enter_droplet;
	real halfDen;

	crSize = V_PROD (contRegGrid);
        
	if (opCode == 0) {
		for (i = 0; i < contRegGrid.x; i++) {
			for (j = 0; j < contRegGrid.y; j++) {
				contRegDen2D[i][j] = 0.;
				contRegDenMap[i][j] = 0.;
			}
		}
		for (n = 0; n < crSize; n++) {
			contRegDen[n] = 0.;
		}
	} else if (opCode == 1) {
		V_DIV (invWid, sizeCoordGrid, region);
		DO_POLY_MOL {
			V_S_ADD (rs, mol[n].r, 0.5, region);    //shifts coordinates to "positive" region (from 0 to L)
			V_MUL (cc, rs, invWid);
			c = V_LINEAR (cc, sizeCoordGrid);
			/* take a look at the contact region */
			if (cc.z >= k0) {
				c -= k0 * sizeCoordGrid.x * sizeCoordGrid.y;
				if (c < crSize) ++contRegDen[c];
			}
		}
	} else if (opCode == 2) {
		halfDen = density / 2.;                         // profile's line criteria

		for (n = 0; n < crSize; n++) contRegDen[n] /= (limitContReg / stepContReg);
	
	        for (j = 0; j < contRegGrid.y; j++) {
			enter_droplet = 0;
 			for (i = 0; i < contRegGrid.x; i++) {
				for (k = 0; k < contRegGrid.z; k++) {
					n = k * contRegGrid.x * contRegGrid.y + j * contRegGrid.x + i;
					contRegDen2D[i][j] += contRegDen[n];
				}
//				contRegDen2D[i][k] /= (region.x * region.y * contRegGrid.z * region.z / sizeCoordGrid.z) / (contRegGrid.x * contRegGrid.z);
				contRegDen2D[i][j] /= ( V_PROD(region) / (contRegGrid.x * contRegGrid.y) );
				contRegDenMap[i][j] = contRegDen2D[i][j];

				/* trying enter droplet in x-direction */
				if ( (enter_droplet == 0) && (contRegDen2D[i][j] >= halfDen) ) {
					cRegLeftProf[j] = (i + 0.5 - contRegGrid.x / 2.) * region.x / contRegGrid.x;
					enter_droplet = 1;
				} else if ( (enter_droplet == 1) && (contRegDen2D[i-1][j] >= halfDen) && 
				      (contRegDen2D[i][j] < halfDen) ) {
						cRegRightProf[j] = (i + 0.5 - contRegGrid.x / 2.) * region.x / contRegGrid.x;
				}
			}
		}	
	}
}
