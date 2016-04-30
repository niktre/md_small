#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <sys/time.h>

typedef double real;

#define ALLOC_MEM(a, n, t)  a = (t *) malloc ((n) * sizeof (t))

#define ALLOC_MEM2(a, n1, n2, t)  			\
	ALLOC_MEM (a, n1, t *);				\
	ALLOC_MEM (a[0], (n1) * (n2), t);		\
	for (k = 1; k < n1; k ++) a[k] = a[k - 1] + n2;

#include "in_vdefs.h"
#include "in_namelist.h"
#include "in_proto.h"

#define DO_MOL  for (n = 0; n < nMol; n ++)
#define DO_POLY_MOL  for (n = 0; n < nPolyMol; n++)
#define DO_SURF  for (at = 0; at < nSurf; at ++)
#define DO_CELL(j,m)  for (j = cellList[m]; j >= 0; j = cellList[j])

//#define V_WRAP(v, t)								\
   if (v.t >= .5 * region.t) do(v.t -= region.t); while(v.t >= .5 * region.t);	\
   else if (v.t < -.5 * region.t) do(v.t += region.t); while(v.t < -.5 * region.t);	
#define FIND_SLAB(slab, slabNum, pos, regSize)							\
	if (pos < 0.) slab = (int) (slabNum / 2.) - (int) (-pos * slabNum / regSize) - 1;	\
	else slab = (int) (slabNum / 2.) + (int) (pos * slabNum / regSize);

#define V_WRAP(v, t)						\
	if (v.t >= 0.5 * region.t) v.t -= region.t;		\
	else if (v.t < -0.5 * region.t) v.t += region.t		
#define V_SHIFT(v, t)						\
	if (v.t >= 0.5 * region.t) shift.t -= region.t		\
	else if (v.t < -0.5 * region.t) shift.t += region.t
#define V_SHIFT_WRAP(v, t)					\
	if (v.t >= 0.5 * region.t) {				\
		shift.t -= region.t;				\
		v.t -= region.t;				\
	} else if (v.t < -0.5 * region.t) {			\
	  shift.t += region.t;					\
	  v.t += region.t;					\
	}
#define V_CELL_WRAP(t)						\
	if (m2v.t >= cells.t) {					\
		m2v.t = 0;					\
		shift.t = region.t;				\
	} else if (m2v.t < 0) {					\
		m2v.t = cells.t - 1;				\
		shift.t = -region.t;				\
	}

#define V_WRAP_ALL(v)						\
	{V_WRAP (v, x);						\
	V_WRAP (v, y);}//					\
	V_WRAP (v, z);}
#define V_SHIFT_ALL(v)						\
	{V_SHIFT (v, x);					\
	V_SHIFT (v, y);}//					\
	V_SHIFT (v, z);}
#define V_CELL_WRAP_ALL()					\
	{V_CELL_WRAP (x);					\
	V_CELL_WRAP (y);}//					\
	V_CELL_WRAP (z);}

#define N_OFFSET	14

#define OFFSET_VALS						\
	{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1},	\
	{1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1},		\
	{-1,-1,1}, {0,-1,1}, {1,-1,1}}

typedef struct {
  real val, sum, sum2;
} Prop;

#define PROP_ZERO(v)  v.sum = v.sum2 = 0.
#define PROP_ACCUM(v)						\
	v.sum += v.val,						\
	v.sum2 += SQR (v.val)

#define PROP_AVG(v, n)						\
	v.sum /= n, 						\
	v.sum2 = sqrt (MAX (v.sum2 / n - SQR (v.sum), 0.))

#define PROP_EST(v)  v.sum, v.sum2

enum {ERR_NONE, ERR_BOND_SNAPPED, ERR_CHECKPT_READ, ERR_CHECKPT_WRITE,
   ERR_COPY_BUFF_FULL, ERR_EMPTY_EVPOOL, ERR_MSG_BUFF_FULL,
   ERR_OUTSIDE_REGION, ERR_SNAP_READ, ERR_SNAP_WRITE,
   ERR_SUBDIV_UNFIN, ERR_TOO_MANY_CELLS, ERR_TOO_MANY_COPIES,
   ERR_TOO_MANY_LAYERS, ERR_TOO_MANY_LEVELS, ERR_TOO_MANY_MOLS,
   ERR_TOO_MANY_MOVES, ERR_TOO_MANY_NEBRS, ERR_TOO_MANY_REPLICAS, ERR_NO_FOLDER,
   ERR_COORD_READ, ERR_VELACC_READ, ERR_CENT_O_MASS};

char *errorMsg[] = {"", "bond snapped", "read checkpoint data",
   "write checkpoint data", "copy buffer full", "empty event pool",
   "message buffer full", "outside region", "read snap data",
   "write snap data", "subdivision unfinished", "too many cells",
   "too many copied mols", "too many layers", "too many levels",
   "too many mols", "too many moved mols", "too many neighbors",
   "too many replicas","no folder name", "wrong input coordinates",
   "wrong input vels and accels", "wrong grid cell for CoM"};

