void AccumProps (int);

void AdjustTemp (void);

void AllocArrays (void);

void ApplyBoundaryCond (void);

void AssignToChain (void);

void BuildNebrList (void);

void *BuildNebrListT (void *);

void CentrifyHistogram (VecLI, int, int, real **, real **);

void ComputeChainBondForces (void);

void ComputeChainBondForcesPT (void);

void ComputeEinstCryst (void);

void ComputeExternalForces (void);

void ComputeForces (void);

void ComputeForcesPT (void);

void ComputeWallForces (void);

void ComputeWandSchwingungen (void);

void CoordGridAverage (int);

void CopyAbsWrapChains (void);

void CVFoutput (FILE *, int, FILE *);

void EnGridAverage(int);

void ErrExit (int);

void EvalProps (void);

void EvalControlValues (int, int);

void EvalDenDist (int);

void EvalGridProfile (void);

void EvalPrDist (int);

void EvalVelDist (int);

int  GetNameList (int, char **);

void InitAccels (void);

void InitCoords (int);

void InitHistograms (void);

void InitRand (int);

void InitState (void);

void InitSubstrate (int);

void InitVels (void);

void itoa (int, char *);          // function converts int to char

void LeapfrogStep (int);

void MakeFilename (int, int);

void MakeParamFile (FILE *, FILE *, FILE *, int, char **);

void ModifyWidth (void);

void PrintChainProps (FILE *);

void PrintConstraint (FILE *, FILE *, FILE *);

void PrintControlFile (FILE *);

void PrintDenDist (FILE *);

void PrintEnMom (FILE *);

void PrintGridProfile (FILE *, FILE *, FILE *, FILE *);

void PrintInitConf (FILE *, FILE *);

void PrintNameList (FILE *);

void PrintPosVelCM (FILE *);

void PrintPrDist (FILE *);

void PrintSUEn (FILE *);

void PrintSummary (FILE *, int);

void PrintVelDist (FILE *);

void PrintVTF (FILE *);

real RandR (void);

void ReadCoordsVelsAccels (FILE *, FILE *, FILE *, FILE *);

void reverse (char *);                   // function reverses char

void SetParams (void);

void SetupJob (void);

void SingleStep (void);

void SolveCubic (real *, real *);

void VRand (VecR *);

void VelGridAverage(int);

void WrapChainsIntoCell (void);
