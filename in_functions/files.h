char subfoldername_1[] = "a/";
char coord_fullname[80];				// full coordinat file's name
char end_coord_fullname[80];				// full coordinat file's name
char end_coord_fullname_test[80];			// full coordinat file's name
char vf_fullname[80];					// full velocity-force file's name
char end_vf_fullname[80];				// full coordinat file's name
char end_vf_fullname_test[80];				// full coordinat file's name
char init_fullname[80];                                 // full initial conformation file's name
char control_fullname[80];                              // full control file's name
char en_mom_fullname[80];                               // full name of the file with energy and momentum components
char su_en_fullname[80];				// full name of the file with substrate unit - polymer liquid energy
char const_fullname[80];				// full name of the file with constraint forces
char free_en_fullname[80];				// full name of the file with constraint forces
char param_fullname[80];				// full name of the simulation parameter file
char hist_param_fullname[80];				// full name of the simulation parameter file history
char vtf_fullname[80];                                  // full vtf-file's name

char coordfile[] = "coord_out";                         // general coordinat file's name
char vffile[] = "vf_out";                               // general velocity-force file's name
char init_filename[] = "initial.dat";                   // initial configuration file's name
char start_conform[] = "coord_start";			// file with new starting point conformation
char start_vf[] = "vf_start";				// file with new starting point velocities and forces
char control_filename[] = "file_control.dat";           // file of control values
char en_mom_filename[] = "en_momentum.dat";             // file of energy and momentum in slab 0
char su_en_filename[] = "su_energy.dat";		// file of substrate unit - polymer liquid energy
char const_filename[] = "constraint.dat";		// file of constraint forces
char free_en_filename[] = "freeEnDiff.dat";		// file of differences in free energy
char param_filename[] = "input_values.in";		// file of simulation parameters
char hist_param_filename[] = "hist_input_values.in";	// history file of simulation parameters
char vtf_filename[] = "trajectory.vtf";                 // vtf-file for VMD
char init_coord_filename[] = "ICoord.dat";		// initial coordinates
char init_vel_filename[] = "IVel.dat";			// initial velocities 

// density, pressure tensor and velocity profile
char den_hist_filename[] = "den_hist";			// general density profile file's name
char denX_hist_filename[] = "xden_hist";		// general density profile file's name
char denXcr_hist_filename[] = "xden_contReg";		// general density profile file's name
char pn_hist_filename[] = "pn_hist";			// general pressure tensor profile file's name
char vel_hist_filename[] = "vel_hist";			// velocity profile file's name
char file_grid_filename[] = "grid_profile";
char pos_vel_cm_filename[] = "cmPosVel.dat";		// center of mass parameters file's name
char denMap_filename[] = "2D_den";
char contRegDenMap_filename[] = "2D_contReg";
char dropVel2D_filename[] = "dropVel";
char den_hist_fullname[80];				// full density profile file's name
char denX_hist_fullname[80];				// full density profile file's name
char denXcr_hist_fullname[80];				// full density profile file's name
char pn_hist_fullname[80];				// full pressure tensor profile file's name
char vel_hist_fullname[80];				// full velocity
char file_grid_fullname[80];
char pos_vel_cm_fullname[80];				// full center of mass parameters file's name
char denMap_fullname[80];
char contRegDenMap_fullname[80];
char dropVel2D_fullname[80];
