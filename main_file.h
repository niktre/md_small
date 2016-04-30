typedef struct {
   VecR r, rv, ra, r_abs, r_vtf;		// r-vect, v-vect and a-vect, r-absolute and r-vtf
   int inChain;				// chain identifier
   int inSUnit;				// identificator of belonging to a substrate unit (ridge, pillar, etc.)
} Mol;

typedef struct {
   VecR r;
} Ghost;

Ghost *ghost;
Mol *mol;

/* TYPE OF THE PROBLEM: "0" - film, "1" - flow, "2" - droplet */
int problem;

/* NUMBER OF TYPICAL ATOMS AND MONOMERS */
int nMol, nPolyMol, nSurf;    // number of molecules in the box, number of polymer monomers and of surface atoms

/* FILE VARIABLE */
char foldername[30];          // name of working folder

/* TYPICAL SIZES IN THE SYSTEM */
VecR region;                  // vect-size of the simulation box
VecR initUcell;					// vect-size of initial cell with chains (more interested in ratio)
VecR gapSurf;                 // size of surface unit cell
VecI surfCell;                // number of replications of unit surface cell
real addSpace;                // space adding above the film
real addWidth;                // control the width of the channel in FLOW problem
real shiftWidth;              // take an old channel and move the walls toward each other
int N_high_surf;              // number of slab where the highest substrate atome is
int oldSub;                   // do we use the substrate from old simulations (for surface tension calculation)
real dTop_high_surf;          // distance from the highest substrate atom to the top border of the it's slab
real *centUCinZ;              // z-coordinate of the center of cell layer

/* PARAMETERS OF POTENTIALS */
real rCut;                    // cut-off distance for potentials
real rrCut, rriCut, rriCut3;  // it's powers
real zCutWall;                // cut-off distance for repulsive wall interations
real zzCutWall;               // it's power
real U0, U0_surf;             // shifting of L-J-potential in poly-liquid and between poly-liquid and substrate
real epsWall, surf_dist;      // strength of polymer-substrate interaction and unitlength-parameter in surface cell
real sigma_sur, s3_sur, ss3_sur;		// sigma in polymer-substrate potential
real R2, k;                   // FENE-constants
real gamma_d;           // DPD-constants

/* GENERAL TRAJECTORY PARAMETERS */
int controlFile;				// whether we need the file with control-values (Ree, etc.)
int randSeed;					// initial random value (depends on clock when = 0)
int storeVelFor;				// whether we need to store velocities and forces
int stepConVal, limitConVal;			// analyse values every stepConVal steps till come to limitConVal step
int stepAvg, stepCount, stepEquil, stepLimit;	// parameters of the trajectory
int VTF_stepPrint;				// print configuration into VTF-file every VTF_stepPrint steps
int restart;					// continuation of the trajectory for longer run
int oldEndStep;					// last step in previous trajectory
int einstCryst;					// represent substrate as a einstein crystal
real ampSchwing, periodSchwing, omegaSchwing;	// oscillations of walls
real epsChemWall, periodChem;			// 2nd eps for a running droplet and period of chemical changing of the substrate
real epsChemStore;				// reservoir for keeping two epsilons
int getIConf;					// start from existing conformation

/* THERMOSTAT PARAMETERS */
real deltaT, density;				// timestep, initial density
real temperature;				// initial temperature in the box (sets initial velocities)

/* CALCULATION PROPERTIES OF THE SYSTEM*/
real timeNow;					// current "box"-time
VecR vSum;					// "impulse"-vect of the system
real U_LJ_sub, uSum, velMag, velMagSub, vvSum, virSum;	// LJ-potential for pol-liquid-substrate, OTHERS
Prop kinEnergy, totEnergy, pressure;	// kinEn, totEn, Polymer-Surface Energy and pressure in the whole cell
real instTemp;					// temperature in the box (averaged over time)

/* CHAIN PARAMETERS */
VecI initUchain;				// numbezer of chain-cell replications
int chainLen; 					// length of the chain
int nChain;					// number of chains
int chainRed;					// not construct every chainRed chain
Prop ReeAv, ReeAv_x, ReeAv_y, ReeAv_z;		// end-to-end distances
real bondLength;				// bond length in chain

/* SURFACE PARAMETERS */
int nSurfLayer;					// number of surface layers
VecI surfSave;					// number of surface cells to save in a layer
VecI surfDel;					// number of surface cells to delete from the layer
int slope;					// incline the substrate's consituents
int flatTopW;					// create a flat wall at the top of the simulation box
real epsTopW;					// value of eps strength for top wall
VecR force1, force2, forceT;			// forces acting on top, bottom and both layers of the substrate
int atTopWall, atBotWall;			// number of atoms in top and bottom walls
real topSubAtZ;					// coordinate of the highest substrate atom

/* NEIGHBOR-LIST*/
// cells division
VecI cells;					//
int *cellList;					//
// neighbor-list variables
real dispHi;					//
real rNebrShell;				//
int nebrTabLen;					//
int *nebrTab, nebrNow, nebrTabFac, nebrTabMax;	//

/* INCREMENTS FOR DIFFERENT LOOPS */
int moreCycles;					// continuation of Step-cycle
int addControlStep;				// trigger of opening control file for "append" info
int countConVal;				// control-file increment value

/* HISTOGRAMS */
// density, normal and transverse pressure histograms
real hSlab;					// height of the slab
real *histDen, *snapDen, *histDenX, *histAvDen;	// density-histograms in z- and x- directions and for snapshot
real *histPn, *histPnW, *histPt, *histPtW; 	// normal and normal wall, transverse and transverse wall
real *histPx, *histPy;				// x and y pressures
real *snapVx, *histVx;				// Vx - component
real *snapTx, *histTx;				// Tx - component
real *snapT, *histT;				// T - value
real *histAvVx;					// average for slab Vx - component
real *histAvTx;					// average for slab Tx - component
real *histAvT;					// average for slab T - value
int countDenPn;					// counter for Den P averaging
int halfSl, halfSlX, numSlabs, numSlabsX;	// half-number of slabs, number of slabs

real sizeCoordBin;				// number of coordinate grid cells in x-, y- and z-directions
real sizeVelBin;				// number of velocity grid cells in x-, y- and z-directions
VecLI sizeCoordGrid;				// number of coordinate grid cells in x-, y- and z-directions
VecLI sizeVelGrid;				// number of velocity grid cells in x-, y- and z-directions
real coordGrid;					// typical size of the grid for coordinates (~0.1 sigma)
int *part_y_p_coord, *part_y_m_coord;		// number of particles in every y-grid-plane for "+" and "-" region.x (for coord calc)
int *ocCellXY;					// occupied cells in XY-plane (with the same Z)
real **histCoordGrid, **histVelGrid;		// snapshot values centrified with respect to CM for every y-layer
real **snapCoordGrid, **snapVelGrid;		// snapshot values in coord and velocity grid
real *rx_p_coord, *rx_m_coord;			// x-coord of CM in coord grid for every y-layer
real *rx_p_vel, *rx_m_vel;
real *profileT, *profileVx;			// Temp and Vx profiles (z-layer) based on grid geometry (not so good for droplet!!!)
real *profileZX, *profile_botZX;		// lines defining top and bottom borders of a drop
real **den2D, **denMap;				// both are density maps (just one doesn't work)
real **vel2D_x, **vel2D_z;		// x- and z-velocities based on averaging over y-layers
real **vel2D_xx, **vel2D_zz, **dropVel2D;	// x-, z- and total energies based on averaging over y-layers

/* DROP MOTION ANALYSIS */
int addCMStep;					// indicator of CM-file existence
VecR globCM, VglobCM;				// coordinates and velocity of Global Centre of Mass

/* THERMOSTAT WORK ANALYSIS */
int *part_y_p_en, *part_y_m_en;
real *rx_p_en, *rx_m_en;
int totSUnitNum;				// total number of substrate units
real effSUArea;					// effective area of substrate Unit
VecI globCMslab;
int n_globCM;
real gravField, shearRate;
VecR fastCMvel, snapCMvel;

NameList nameList[] = {
   NAME_I (problem),		// type of the system under simulation
   NAME_I (flatTopW),		// creates a flat wall made of atoms at the top of the simulation box
   NAME_R (epsTopW),		// eps strength for top wall of the simulation box
   NAME_R (epsChemWall),		// 2nd eps for a running droplet
   NAME_R (periodChem),		// period of chemical changing of the substrate
   NAME_R (ampSchwing),		// amplitude of wall's oscillations
   NAME_R (periodSchwing),		// period of wall's oscillations
   NAME_I (chainLen),		// length of the chain
   NAME_I (chainRed),		// not constract every chainRed chain (redicung the thickness)
   NAME_I (controlFile),		// whether we need the file with control-values (Ree, etc.)
   NAME_R (deltaT),		// timestep
   NAME_R (density),		// initial density of polymer liquid
   NAME_R (gamma_d),		// dissipative constant in DPD
   NAME_I (getIConf),		// initial conformation from external file
   NAME_R (hSlab),			// height of the slab
   NAME_I (oldSub),		// do we use the substrate from old simulations (for surface tension calculation)
   NAME_I (limitConVal),
   NAME_I (nebrTabFac),
   NAME_I (oldEndStep),		// the end step from previous simulation run
   NAME_I (randSeed),		// random number initial value
   NAME_I (restart),		// shows whether current run is continuation of previous
   NAME_R (rNebrShell),		// radius of neighbor-shell
   NAME_R (sigma_sur),		// defines the sigma-constant in VdW interaction
   NAME_R (sizeCoordBin),		// size of coordinate histogram grid
   NAME_R (sizeVelBin),		// size of velocity histogram grid
   NAME_I (slope),			// make a slope in the substrate
   NAME_I (stepAvg),		// number of steps for averaging
   NAME_I (stepConVal),		// evaluate control values every stepConVal steps
   NAME_I (stepEquil),		// number of equilibration steps
   NAME_I (stepLimit),		// length of trajectory (in steps)
   NAME_I (storeVelFor),		// indicates whether we need to store velocities and forces for trajectory
   NAME_R (surf_dist),		// defines the density of surface atoms
   NAME_R (gravField),		// Poiseuille flow
   NAME_R (shearRate),		// Couette flow
   NAME_R (temperature),		// temperature in the system
   NAME_I (VTF_stepPrint),		// print configuration into VTF-file every VTF_stepPrint steps
   NAME_I (nSurfLayer),
   NAME_I (surfSave),
   NAME_I (surfDel),
   NAME_R (epsWall),		// strength of liquid-substrate interaction
   NAME_I (initUchain),		// replications of chains-unit
   NAME_R (addWidth),		// control the width of the channel in FLOW problem;
   NAME_R (shiftWidth),		// detailed control of the channel's width;
};