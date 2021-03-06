void MakeFilename (int temp_file_num, int temp_store_vel_for) {
	char fileext[] = ".txt";                        // file extention (look in input_value.c)

	int size_s = 0;                                 // number of digits in the name of cvf_file
        int num_null;                                   // number of "0" simbols in the name of cvf_file
	extern char foldername[];			// name of the correspondent folder
	extern char parameter_filepath[];
        extern char coord_fullname[80];			// full name of the coordinate file
        extern char coordfile[];

	extern char den_hist_fullname[80];
	extern char denX_hist_fullname[80];
	extern char denMap_fullname[80];
	extern char dropVel2D_fullname[80];
	extern char pn_hist_fullname[80];
	extern char vel_hist_fullname[80];
	extern char file_grid_fullname[80];

	extern char den_hist_filename[];
	extern char denX_hist_filename[];
	extern char denMap_filename[];
	extern char dropVel2D_filename[];
	extern char pn_hist_filename[];
	extern char vel_hist_filename[];
	extern char file_grid_filename[];
	
        extern void itoa (int temp_i, char temp_s[]);          // function converts int to char
        extern void reverse (char temp_s[]);                   // function reverses char

        int tempLimit;                                  // temporary variable for calculating number of digits

        tempLimit = (int) stepLimit / stepAvg;		// variable that gives filenumber (depends on current step)

	if (tempLimit < 10) size_s++;

        while (tempLimit != 0) {                        // cycle for number of digit's calculatuion
	        size_s++;
        	tempLimit /= 10;
        }

        char s[size_s];                                 // char variable for digitals in output filename

        num_null = size_s;                              // number of digitals in output file

	/* creating a filename */
        itoa (temp_file_num, s);                        // convert int to string
	strcpy (coord_fullname, parameter_filepath);		// copy "parameter_filepath" to "coord_fullname"
	strcpy (den_hist_fullname, parameter_filepath);		// copy "parameter_filepath" to "den_hist_fullname"
	strcpy (denX_hist_fullname, parameter_filepath);	// copy "parameter_filepath" to "denX_hist_fullname"
	strcpy (denMap_fullname, parameter_filepath);		// copy "parameter_filepath" to "denMap_fullname"
	strcpy (dropVel2D_fullname, parameter_filepath);	// copy "parameter_filepath" to "dropVel2D_fullname"
	strcpy (pn_hist_fullname, parameter_filepath);		// copy "parameter_filepath" to "pn_hist_fullnsme"
	strcpy (vel_hist_fullname, parameter_filepath);		// copy "parameter_filepath" to "pn_hist_fullnsme"
	strcpy (file_grid_fullname, parameter_filepath);	// copy "parameter_filepath" to "pn_hist_fullnsme"

	strcat (coord_fullname,foldername);		// copy "foldername" to "coord_fullname"
	strcat (den_hist_fullname,foldername);		// copy "foldername" to "den_hist_fullname"
	strcat (denX_hist_fullname,foldername);		// copy "foldername" to "denX_hist_fullname"
	strcat (denMap_fullname,foldername);		// copy "foldername" to "denMap_fullname"
	strcat (dropVel2D_fullname,foldername);		// copy "foldername" to "dropVel2D_fullname"
	strcat (pn_hist_fullname,foldername);		// copy "foldername" to "pn_hist_fullnsme"
	strcat (vel_hist_fullname,foldername);		// copy "foldername" to "pn_hist_fullnsme"
	strcat (file_grid_fullname,foldername);		// copy "foldername" to "pn_hist_fullnsme"

        strcat (coord_fullname,coordfile);		// concatenate "coordfile" with "coord_fullname"
	strcat (den_hist_fullname,den_hist_filename);	// concatenate "den_hist_filename" with "den_hist_fullname"
	strcat (denX_hist_fullname,denX_hist_filename);	// concatenate "denX_hist_filename" with "denX_hist_fullname"
	strcat (denMap_fullname,denMap_filename);	// concatenate "denMap_filename" with "denMap_fullname"
	strcat (dropVel2D_fullname,dropVel2D_filename);	// concatenate "dropVel2D_filename" with "dropVel2D_fullname"
	strcat (pn_hist_fullname,pn_hist_filename);	// concatenate "pn_hist_filename" with "pn_hist_fullname"
	strcat (vel_hist_fullname,vel_hist_filename);	// concatenate "vel_hist_filename" with "vel_hist_fullname"
	strcat (file_grid_fullname,file_grid_filename);	// concatenate "file_grid_filename" with "file_grid_fullname"

        if (temp_store_vel_for == 1) {
                extern char vf_fullname[80];            // full name of the cvf_file
                extern char vffile[];
		strcpy (vf_fullname, parameter_filepath);	// copy "foldername" to "vf_fullname"
		strcat (vf_fullname,foldername);	// copy "foldername" to "vf_fullname"
                strcat (vf_fullname,vffile);           // concatenate "vffile" with "vf_fullname"
        }

/* define the number of zeros before the filenumber */
                if      (temp_file_num < 10)            // if filenumber from 1 to 9
                        num_null = num_null - 1;
                else if (temp_file_num < 100)           // if filenumber from 10 to 99
                        num_null = num_null - 2;
                else if (temp_file_num < 1000)          // if filenumber from 100 to 999
                        num_null = num_null - 3;
                else if (temp_file_num < 10000)         // if filenumber from 1000 to 9999
                        num_null = num_null - 4;
                else if (temp_file_num < 100000)        // if filenumber from 10000 to 99999
                        num_null = num_null - 5;
                else if (temp_file_num < 1000000)       // if filenumber from 100000 to 999999
                        num_null = num_null - 6;

/* write corresponding number of zeros to filename */
        if (temp_store_vel_for == 1) {
                while (num_null != 0) {
                        strcat(coord_fullname,"0");                   // concatenate with the number of zeros
                        strcat(vf_fullname,"0");
                        strcat(den_hist_fullname,"0");                   // concatenate with the number of zeros
                        strcat(denX_hist_fullname,"0");                   // concatenate with the number of zeros
			strcat(denMap_fullname,"0");			// concatenate with the number of output file
			strcat(dropVel2D_fullname,"0");			// concatenate with the number of output file
                        strcat(pn_hist_fullname,"0");                   // concatenate with the number of zeros
                        strcat(vel_hist_fullname,"0");                   // concatenate with the number of zeros
                        strcat(file_grid_fullname,"0");                   // concatenate with the number of zeros
                        --num_null;
                }
        }
        else {
                while (num_null != 0) {
                        strcat (coord_fullname,"0");			// concatenate with the number of zeros
                        strcat (den_hist_fullname,"0");                 // concatenate with the number of zeros
                        strcat (denX_hist_fullname,"0");                // concatenate with the number of zeros
                        strcat (denMap_fullname,"0");                 // concatenate with the number of output file
			strcat (dropVel2D_fullname,"0");			// concatenate with the number of output file
                        strcat (pn_hist_fullname,"0");                  // concatenate with the number of zeros
                        strcat (vel_hist_fullname,"0");                 // concatenate with the number of zeros
                        strcat (file_grid_fullname,"0");                // concatenate with the number of zeros
                        --num_null;
                }
        }
/***************************************************/
	strcat (coord_fullname,s);                            // concatenate with the number of output file
	strcat (coord_fullname,fileext);                      // concatenate with file extention. Result: the filename is created

        if (temp_store_vel_for == 1) {
                strcat (vf_fullname,s);
                strcat (vf_fullname,fileext);
        }

	strcat (den_hist_fullname,s);				// concatenate with the number of output file
	strcat (den_hist_fullname,".dat");			// concatenate with file extention. Result: the filename is created
	strcat (denX_hist_fullname,s);				// concatenate with the number of output file
	strcat (denX_hist_fullname,".dat");			// concatenate with file extention. Result: the filename is created
	strcat (denMap_fullname,s);				// concatenate with the number of output file
	strcat (denMap_fullname,".dat");			// concatenate with the number of output file
	strcat (dropVel2D_fullname,s);				// concatenate with the number of output file
	strcat (dropVel2D_fullname,".dat");			// concatenate with the number of output file
	strcat (pn_hist_fullname,s);				// concatenate with the number of output file
	strcat (pn_hist_fullname,".dat");			// concatenate with file extention. Result: the filename is created
	strcat (vel_hist_fullname,s);				// concatenate with the number of output file
	strcat (vel_hist_fullname,".dat");			// concatenate with file extention. Result: the filename is created
	strcat (file_grid_fullname,s);				// concatenate with the number of output file
	strcat (file_grid_fullname,".dat");			// concatenate with the number of output file
}

void CVFoutput (FILE *coord_file, int temp_store_vel_for, FILE *vf_file) {

        int n;

        extern char coord_fullname[80];
	extern char vf_fullname[80];

	// coordinate file
        coord_file = fopen (coord_fullname, "a");            // open file
	if (coord_file != NULL) {
		fprintf (coord_file, "# nMol %6d, box size (%9.4f x %9.4f x %9.4f) chains: %d\n", nMol, region.x, region.y, region.z, nChain);
		fprintf (coord_file, "# atom     x        y        z        ch   \n");
	        DO_MOL fprintf (coord_file,
	        	"%6d %9.4f %9.4f %9.4f %6d\n",
	                n + 1, mol[n].r.x, mol[n].r.y, mol[n].r.z, mol[n].inChain + 1);
	}
	fclose (coord_file);
	
	// velocity and force file
	vf_file = fopen (vf_fullname, "a");
	if (vf_file != NULL) {
//		fprintf (vf_file, "# atom    vx        vy        vz        ax        ay        az         ch   \n");
		fprintf (vf_file, "# atom    vx        vy        vz\n");
	        DO_MOL fprintf (vf_file,
//			"%6d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %6d\n",
			"%6d %9.4f %9.4f %9.4f\n",
//			n + 1, mol[n].rv.x, mol[n].rv.y, mol[n].rv.z, mol[n].ra.x, mol[n].ra.y, mol[n].ra.z, mol[n].inChain + 1);
       		        n + 1, mol[n].rv.x, mol[n].rv.y, mol[n].rv.z);
	}
        fclose (vf_file);
}

void PrintInitConf (FILE *init_file, FILE *vtf_file) {
	int at, n, i;

	extern char foldername[];
	extern char parameter_filepath[];
	extern char init_filename[];
	extern char init_fullname[];
	extern char vtf_filename[];
	extern char vtf_fullname[];

//	extern VecR globCM;

	strcpy (init_fullname,parameter_filepath);
	strcpy (vtf_fullname,parameter_filepath);
	strcat (init_fullname,foldername);		// copy "foldername" to "init_fullname"
	strcat (vtf_fullname,foldername);		// copy "foldername" to "vtf_fullname"
	strcat (init_fullname,init_filename);		// concatenate "init_filename" with "init_fullname"
	strcat (vtf_fullname,vtf_filename);		// concatenate "vtf_filename" with "vtf_fullname"

	// creating initial conformation file
	init_file = fopen (init_fullname, "a");
	fprintf (init_file, "# nMol %6d, box size (%9.4f x %9.4f x %9.4f)\n", nMol, region.x, region.y, region.z);
	fprintf (init_file, "# atom     x        y        z       vx       vy       vz        ch   \n");
        DO_MOL fprintf (init_file,
        	"%6d %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%6d\n",
                n + 1, mol[n].r.x, mol[n].r.y, mol[n].r.z, mol[n].rv.x, mol[n].rv.y, mol[n].rv.z, mol[n].inChain + 1);

	fclose (init_file);
	
	// creating vtf trajectory file
	vtf_file = fopen (vtf_fullname, "w");
	fprintf (vtf_file, "# STRUCTURE BLOCK\n");
	for (i = 0; i < nPolyMol; i += chainLen) {
		fprintf (vtf_file, "a %d:%d r 0.8 resid %d\nb %d::%d\n", i, i + chainLen - 1, mol[i].inChain, i, i + chainLen - 1);
	}
	// making a flag "SURFACE" and assigning to Units of the Substrate their numbers
	for (i = nPolyMol; i < nMol; i++) {
			if (mol[i].inSUnit == -1) {		// if it is a sole of the substrate
				fprintf (vtf_file, "a %d r 0.5 name SSole\n", i);
			} else if (mol[i].inSUnit == -2) {	// if it is a top wall 
				fprintf (vtf_file, "a %d r 0.5 name STopW\n", i);
			} else { 				// if it is normal unit of the substrate
				fprintf (vtf_file, "a %d r 0.5 name SUnit%d\n", i, mol[i].inSUnit);
			}
	}
	// making ghost atoms
	if (einstCryst == 1) {
		for (i = nMol; i < nMol + nSurf; i++) {
			fprintf (vtf_file, "a %d r 0.5 name Ghost\n", i);
		}
	}
//	fprintf (vtf_file, "a %d r 2.0 name CM\n", nMol);
//	fprintf (vtf_file, "a %d:%d r 0.5 name SURFACE\n", nPolyMol, nMol - 1);
	
	// creating timestep	
	fprintf (vtf_file, "# TIMESTEP BLOCK\n# Init Conf\n");
	fprintf (vtf_file, "timestep\npbc %9.4f %9.4f %9.4f\n", region.x, region.y, region.z);
        DO_MOL fprintf (vtf_file, " %9.4f %9.4f %9.4f\n", mol[n].r_vtf.x, mol[n].r_vtf.y, mol[n].r_vtf.z);
	// making ghost atoms
	if (einstCryst == 1) {
		DO_SURF fprintf (vtf_file, " %9.4f %9.4f %9.4f\n", ghost[at].r.x, ghost[at].r.y, ghost[at].r.z);
	}	
	// print CM
//	fprintf (vtf_file, " %9.4f %9.4f %9.4f\n", globCM.x, -region.y / 2. + .1, globCM.z);

	fclose (vtf_file);
   fflush (vtf_file);

}

void PrintVTF (FILE *vtf_file) {
	int at, n;
	extern int restart;
	extern VecR region;
//	extern VecR globCM;
	extern int stepCount;

	extern char foldername[];
        extern char parameter_filepath[];
	extern char vtf_filename[];
	extern char vtf_fullname[];

	if (restart == 0) {
	} else if (restart == 1) {
		strcpy (vtf_fullname,parameter_filepath);	// copy "parameter_filepat" to "vtf_fullname"
		strcat (vtf_fullname,foldername);		// copy "foldername" to "vtf_fullname"
		strcat (vtf_fullname,vtf_filename);		// concatenate "vtf_filename" with "vtf_fullname"
	}

	vtf_file = fopen (vtf_fullname, "a");

	if (vtf_file != NULL) {
		fprintf (vtf_file, "#step %d\n", stepCount);
		fprintf (vtf_file, "timestep\n");
		fprintf (vtf_file, "pbc %9.4f %9.4f %9.4f\n", region.x, region.y, region.z);

		DO_MOL fprintf (vtf_file, " %9.4f %9.4f %9.4f\n", mol[n].r_vtf.x, mol[n].r_vtf.y, mol[n].r_vtf.z);
		// making ghost atoms
		if (einstCryst == 1) {
			DO_SURF fprintf (vtf_file, " %9.4f %9.4f %9.4f\n", ghost[at].r.x, ghost[at].r.y, ghost[at].r.z);
		}	
		// print CM
//		fprintf (vtf_file, " %9.4f %9.4f %9.4f\n", globCM.x, -region.y / 2. + .1, globCM.z);
	}

	fclose (vtf_file);
}

void PrintControlFile (FILE *control_file) {
	extern int addControlStep, halfSl;
	extern char foldername[];
        extern char parameter_filepath[];
	extern char control_filename[];
	extern char control_fullname[];
	extern real bondLength;
	extern Prop ReeAv, ReeAv_x, ReeAv_y, ReeAv_z;

	strcpy (control_fullname,parameter_filepath);		// copy "foldername" to "control_fullname"
	strcat (control_fullname,foldername);		// copy "foldername" to "control_fullname"
	strcat (control_fullname,control_filename);	// concatenate "control_filename" with "control_fullname"

	control_file = fopen (control_fullname, "a");

	if (control_file != NULL) {
		if (addControlStep == 0) {
			fprintf (control_file, "     #step");
			fprintf (control_file, "   <bond>\n ");
		}

		fprintf (control_file, "  %8d", stepCount);
		fprintf (control_file, " %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n", 
				PROP_EST (ReeAv), PROP_EST (ReeAv_x), PROP_EST (ReeAv_y), PROP_EST (ReeAv_z), 
				sqrt(bondLength)); 
	}

        fclose (control_file);
}

void MakeParamFile (FILE *param_file, FILE *hist_param_file, FILE *init_file, int argc, char **argv) {
	
	char param_str[100];
	extern char foldername[];
        extern char parameter_filepath[];
        extern char param_filename[];
        extern char param_fullname[80];
        extern char hist_param_filename[];
        extern char hist_param_fullname[80];
	int argz = 1;
	int epsWall_flag = 0;
	int initUchain_flag = 0;
	int surfSave_flag = 0;
	int surfDel_flag = 0;
	int gravField_flag = 0;
	int shearRate_flag = 0;
	int chainRed_flag = 0;
	int addWidth_flag = 0;
	int shiftWidth_flag = 0;

	/* check the existence of the correct folder */
	if (strcmp(argv[argz], "-foldername") == 0) {
		strcpy(foldername,argv[++argz]);
		++argz;
	} else {
		ErrExit(ERR_NO_FOLDER);
	}

	strcpy (param_fullname, parameter_filepath);		// copy destination to "param_fullname"
	strcpy (hist_param_fullname, parameter_filepath);	// copy destination to "param_fullname"
	strcat (param_fullname,foldername);			// concatenate "foldername" with "param_fullname"
	strcat (hist_param_fullname,foldername);		// concatenate "foldername" with "param_fullname"
	strcat (param_fullname,param_filename);			// concatenate "param_filename" with "param_fullname"
	strcat (hist_param_fullname,hist_param_filename);	// concatenate "param_filename" with "param_fullname"

	init_file = fopen (initialParam_file, "r");
	param_file = fopen (param_fullname, "w");

	/* copy the general parameters from init_file to param_file */
	while (!feof(init_file)) {
		if (fgets (param_str, 100, init_file) == NULL) {
			continue;
		} else {
			fputs (param_str, param_file);
		}
	}
        fclose (init_file);

	/* adding the specific parameters to param_file */
	while ((argz < argc) && (argv[argz][0] == '-')) {	// while we have additional parameters, specified by '-'
		if (strcmp(argv[argz], "-epsWall") == 0) {	// if one of parameters is epsWall
			fputs ("epsWall \t", param_file);	// add "epsWall %f \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			epsWall_flag = 1;
		} else if (strcmp(argv[argz], "-initUchain") == 0) {	// if one of parameters is initUchain
			fputs ("initUchain \t", param_file);	// add "initUchain %d %d %d \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs (" ", param_file);	
			fputs (argv[++argz], param_file);	
			fputs (" ", param_file);	
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			initUchain_flag = 1;
		} else if (strcmp(argv[argz], "-surfSave") == 0) {	// if one of parameters is surfSave
			fputs ("surfSave \t", param_file);	// add "surfSave %d %d %d \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs (" ", param_file);	
			fputs (argv[++argz], param_file);	
			fputs (" ", param_file);	
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			surfSave_flag = 1;
		} else if (strcmp(argv[argz], "-surfDel") == 0) {	// if one of parameters is surfDel
			fputs ("surfDel \t", param_file);	// add "surfDel %d %d %d \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs (" ", param_file);	
			fputs (argv[++argz], param_file);	
			fputs (" ", param_file);	
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			surfDel_flag = 1;
		} else if (strcmp(argv[argz], "-gravField") == 0) {     // if one of parameters is gravField
			fputs ("gravField \t", param_file);     // add "gravField %f \n" to parameters file
			fputs (argv[++argz], param_file);
			fputs ("\n", param_file);
			gravField_flag = 1;
		} else if (strcmp(argv[argz], "-shearRate") == 0) {     // if one of parameters is shearRate
			fputs ("shearRate \t", param_file);     // add "shearRate %f \n" to parameters file
			fputs (argv[++argz], param_file);
			fputs ("\n", param_file);
			shearRate_flag = 1;
		} else if (strcmp(argv[argz], "-addWidth") == 0) {	// if one of parameters is addWidth
			fputs ("addWidth \t", param_file);	// add "addWidth %f \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			addWidth_flag = 1;
		} else if (strcmp(argv[argz], "-shiftWidth") == 0) {	// if one of parameters is shiftWidth
			fputs ("shiftWidth \t", param_file);	// add "shiftWidth %f \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			shiftWidth_flag = 1;
		} else if (strcmp(argv[argz], "-chainRed") == 0) {	// if one of parameters is chainRed
			fputs ("chainRed \t", param_file);	// add "chainRed %d \n" to parameters file
			fputs (argv[++argz], param_file);	
			fputs ("\n", param_file);
			chainRed_flag = 1;
		}
		argz++; 
	}
	/* if specific parameters weren't specified set them to the default values */ 
	if (epsWall_flag == 0) {
		fputs ("epsWall \t0.4\n", param_file);		// if there was no epsWall, make default epsWall = 1.0
	}
	if (initUchain_flag == 0) {
		fputs ("initUchain \t8 20 4\n", param_file);	// if there was no initUchain, make default initUchain = 38 6 1
	}
	if (surfSave_flag == 0) {
		fputs ("surfSave \t0 0 0\n", param_file);	// if there was no surfSave, make default surfSave = 0 0 0
	}
	if (surfDel_flag == 0) {
		fputs ("surfDel \t0 0 0\n", param_file);	// if there was no surfDel, make default surfDel = 0 0 0
	}
	if (gravField_flag == 0) {
		fputs ("gravField \t0.0\n", param_file);	// if there was no gravField, make default gravField = 0.0
	}
	if (shearRate_flag == 0) {
		fputs ("shearRate \t0.0\n", param_file);	// if there was no shearRate, make default shearRate = 0.0
	}
	if (chainRed_flag == 0) {
		fputs ("chainRed \t0\n", param_file);		// if there was no chainRed, make default chainRed = 0
	}
	if (addWidth_flag == 0) {
		fputs ("addWidth \t0.0\n", param_file);		// if there was no addWidth, make default addWidth = 0.0
	}
	if (shiftWidth_flag == 0) {
		fputs ("shiftWidth \t0.0\n", param_file);	// if there was no shiftWidth, make default shiftWidth = 0.0
	}
	fclose (param_file);

	/* creating history input values file */
	param_file = fopen (param_fullname, "r");
	hist_param_file = fopen (hist_param_fullname, "a");
	while (!feof(param_file)) {
		if (fgets (param_str, 100, param_file) == NULL)	{	// read the parameter-string from parameter file
			continue;
		} else {
			fputc ('#', hist_param_file);
			fputc (' ', hist_param_file);
			fputs (param_str, hist_param_file);
		}
	}
	fprintf (hist_param_file, "\n");
	fclose (param_file);
	fclose (hist_param_file);

}

void ReadCoordsVelsAccels (FILE *end_coord_file, FILE *end_coord_file_test, FILE *end_vf_file, FILE *end_vf_file_test) {
        int n, oldFileNum, chNumber;
	int numCoords;			// how many numbers were obtained in coord file?
	int maxChNum;			// maximal chain number (to know the number of molecules)
	int numVelsAccels;
	extern int oldEndStep, stepAvg, inChain, halfSl, problem;
	extern int stepLimit;
	extern real hSlab;
	real x,y,z, vx, vy, vz;
	int num;
	real minSubZ, maxSubZ;
	VecR vSum;

	char line[200];
	char fileext[] = ".txt";
	extern char parameter_filepath[];
        extern char end_coord_fullname[80];
	extern char end_coord_fullname_test[80];
	extern char end_vf_fullname[80];
	extern char end_vf_fullname_test[80];
	extern char init_coord_filename[];
	extern char init_vel_filename[];

	minSubZ = region.z / 2.;
	maxSubZ = -region.z / 2.;
	N_high_surf = 0;
	dTop_high_surf = -region.z / 2.;
	maxChNum = 0;
	topSubAtZ = -region.z / 2.;
	
	oldFileNum = (int) (oldEndStep / stepAvg);
	int numFilesExpect = (int) (stepLimit / stepAvg);	// number of files in the simulation we expect
	int numFileDigits = 0;					// number of digits in the filename

        strcpy (end_coord_fullname, parameter_filepath);	// copy "parameter_filepath" to "end_coord_fullname"
        strcpy (end_vf_fullname, parameter_filepath);		// copy "parameter_filepath" to "end_coord_fullname"
        strcat (end_coord_fullname, foldername);		// copy "foldername" to "end_coord_fullname"
	strcat (end_vf_fullname, foldername);			// copy "foldername" to "end_vf_fullname"

	if (oldFileNum == 0 && getIConf == 0) {			// if we start from another initial conformation
	        strcat (end_coord_fullname,start_conform);	// concatenate "coordfile" with "end_coord_fullname"
		strcat (end_vf_fullname,start_vf);		// concatenate "coordfile" with "end_coord_fullname"
	} else if (oldFileNum != 0 && getIConf == 0) {
		strcat (end_coord_fullname,coordfile);		// concatenate "coordfile" with "end_coord_fullname"
		strcat (end_vf_fullname,vffile);		// concatenate "coordfile" with "end_coord_fullname"
	} else {
	}

	/* CREATE NUMBERS USED IN FILENAME */
	/* calculate number of digits in total */
	while (numFilesExpect > 0) {
		++numFileDigits;
		numFilesExpect /= 10;
	}

	/* calculate number of meaningfull digits (i.e. in "0010" meaningfull are "10") */
	int tempOldFileNum = oldFileNum;
	int numMeanDigits = 0;					// number of meaningfull digits in filename
	
	while (tempOldFileNum > 0) {
		++numMeanDigits;
		tempOldFileNum /= 10;
	}

	char s[numMeanDigits];					// char with the size of meaningfull digits in filename

	int numPrefixNull = numFileDigits - numMeanDigits;	// number of "0" in front of the filenumber

	/* write the leading "0" into filenames */
	while (numPrefixNull != 0) {
		strcat (end_coord_fullname,"0");		// concatenate "coordfile" with "0"
		strcat (end_vf_fullname,"0");			// concatenate "coordfile" with "0"
		--numPrefixNull;
	}

	if (getIConf == 1) {
		strcat (end_coord_fullname,init_coord_filename);	// concatenate "coordfile" with "end_coord_fullname"
		strcat (end_vf_fullname,init_vel_filename);		// concatenate "coordfile" with "end_coord_fullname"
	}

	/* write filenumber into the filename */
	if (oldFileNum == 0 && getIConf == 0) {
	} else {
		itoa (oldFileNum,s);
		strcat (end_coord_fullname,s);
		strcat (end_vf_fullname,s);
	}

	if (getIConf != 1) {
	        strcat (end_coord_fullname,fileext);		// concatenate with file extention. Result: the filename is created
		strcat (end_vf_fullname,fileext);		// concatenate with file extention. Result: the filename is created
	}
	/* TESTING FILES*/
        strcpy (end_coord_fullname_test, end_coord_fullname);	// copy 
	strcpy (end_vf_fullname_test, end_vf_fullname);		// copy 
        strcat (end_coord_fullname_test,"t");			// concatenation 
	strcat (end_vf_fullname_test,"t");			// concatenation 

        // coordinate file
        end_coord_file = fopen (end_coord_fullname, "r");		// open file
	end_coord_file_test = fopen (end_coord_fullname_test, "a");	// open test file

	while (!feof(end_coord_file)) {
		fgets (line, 199, end_coord_file);
		if (line[0] == '#') {
			fprintf (end_coord_file_test, "%s", line);
			continue;
		} else {
			numCoords = sscanf (line, "%d%lf%lf%lf%d", &num, &x, &y, &z, &chNumber);
			mol[num-1].r.x = x;
			mol[num-1].r.y = y;
			mol[num-1].r.z = z;
			maxChNum = MAX (maxChNum, chNumber);
			mol[num-1].inChain = chNumber - 1;
			
			if (chNumber == 0 && problem == 1) {
				N_high_surf = MAX (N_high_surf, halfSl - (int) (-z / hSlab) - 1);       // highest slab with surface atom
				dTop_high_surf = MAX (dTop_high_surf, z + hSlab * (halfSl - N_high_surf));      // part of line connecting surface atoms in slab contained highest surf atom
				topSubAtZ = MAX (topSubAtZ, z);
//				minSubZ = MIN (minSubZ, mol[num-1].r.z);
//				maxSubZ = MAX (maxSubZ, mol[num-1].r.z); // maybe a bit unnecessary, but can define later number of substrate layers with help of it
			} else {
			}
		}

		if (numCoords !=5) {
			ErrExit(ERR_COORD_READ);
		}
	}

	/* resetting particle number constants */
	nChain = maxChNum;
	nPolyMol = maxChNum * chainLen;

	if (getIConf == 1 && problem == 1) {
		nMol = num;
		nSurf = nMol - nPolyMol;
	} else {
		InitSubstrate (problem);
	}


        /* absolute coordinates for visualization and transport coefficients */
        DO_MOL V_COPY (mol[n].r_abs, mol[n].r);
        DO_MOL V_COPY (mol[n].r_vtf, mol[n].r);

	DO_MOL fprintf (end_coord_file_test, "%6d%9.4f%9.4f%9.4f%6d\n",
			n + 1, mol[n].r.x, mol[n].r.y, mol[n].r.z, mol[n].inChain + 1);

	fclose (end_coord_file);
	fclose (end_coord_file_test);

        // vf file
	end_vf_file = fopen (end_vf_fullname, "r");			// open file
	end_vf_file_test = fopen (end_vf_fullname_test, "a");		// open test file

	while (!feof(end_vf_file)) {
		fgets (line, 199, end_vf_file);
		if (line[0] == '#') {
			fprintf (end_vf_file_test, "%s", line);
			continue;
		}
		numVelsAccels = sscanf (line, "%d %lf %lf %lf",
			 &num, &vx, &vy, &vz);
		mol[num-1].rv.x = vx;
		mol[num-1].rv.y = vy;
		mol[num-1].rv.z = vz;
		if (numVelsAccels !=4) {
			ErrExit(ERR_VELACC_READ);
		}
	}
	DO_MOL fprintf (end_vf_file_test, "%6d %9.4f %9.4f %9.4f\n",
        	        n + 1, mol[n].rv.x, mol[n].rv.y, mol[n].rv.z);

	fclose (end_vf_file);
	fclose (end_vf_file_test);

	/* make CM velocity equal to Zero */
	if (oldFileNum == 0) {
		V_SET_ALL (vSum, 0.);

		DO_POLY_MOL {
			V_V_ADD (vSum, mol[n].rv);
		}
		
		DO_POLY_MOL V_V_S_ADD (mol[n].rv, - 1. / nPolyMol, vSum);	
	}
}
	
