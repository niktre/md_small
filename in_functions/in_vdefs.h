#ifndef V_DEFS

#define V_DEFS

#ifndef NDIM
#define NDIM  3
#endif

#define SQR(x)     ((x) * (x))
#define CUBE(x)    ((x) * (x) * (x))

#define MIN(x1, x2)  (((x1) < (x2)) ? (x1) : (x2))
#define MAX(x1, x2)  (((x1) > (x2)) ? (x1) : (x2))
#define MIN3(x1, x2, x3)  					\
       (((x1) < (x2)) ? (((x1) < (x3)) ? (x1) : (x3)) :		\
			(((x2) < (x3)) ? (x2) : (x3))) 
#define MAX3(x1, x2, x3)					\
       (((x1) > (x2)) ? (((x1) > (x3)) ? (x1) : (x3)) :		\
			(((x2) > (x3)) ? (x2) : (x3))) 

typedef struct {real x, y;} VecR2;
typedef struct {real x, y, z;} VecR3;
typedef struct {int x, y;} VecI2;
typedef struct {int x, y, z;} VecI3;
typedef struct {long x, y;} VecLI2;
typedef struct {long x, y, z;} VecLI3;

#if NDIM == 3
typedef VecR3 VecR;
typedef VecI3 VecI;
typedef VecLI3 VecLI;
#define V_SET(v, sx, sy, sz)                                 \
   (v).x = sx,                                              \
   (v).y = sy,                                              \
   (v).z = sz
#define V_COPY(v1, v2)                                       \
   (v1).x = (v2).x,                                         \
   (v1).y = (v2).y,                                         \
   (v1).z = (v2).z
#define V_SCALE(v, s)                                        \
   (v).x *= s,                                              \
   (v).y *= s,                                              \
   (v).z *= s
#define V_S_COPY(v2, s1, v1)                                  \
   (v2).x = (s1) * (v1).x,                                  \
   (v2).y = (s1) * (v1).y,                                  \
   (v2).z = (s1) * (v1).z
#define V_ADD(v1, v2, v3)                                    \
   (v1).x = (v2).x + (v3).x,                                \
   (v1).y = (v2).y + (v3).y,                                \
   (v1).z = (v2).z + (v3).z
#define V_SUB(v1, v2, v3)                                    \
   (v1).x = (v2).x - (v3).x,                                \
   (v1).y = (v2).y - (v3).y,                                \
   (v1).z = (v2).z - (v3).z
#define V_MUL(v1, v2, v3)                                    \
   (v1).x = (v2).x * (v3).x,                                \
   (v1).y = (v2).y * (v3).y,                                \
   (v1).z = (v2).z * (v3).z
#define V_DIV(v1, v2, v3)                                    \
   (v1).x = (v2).x / (v3).x,                                \
   (v1).y = (v2).y / (v3).y,                                \
   (v1).z = (v2).z / (v3).z
#define V_S_ADD(v1, v2, s3, v3)                               \
   (v1).x = (v2).x + (s3) * (v3).x,                         \
   (v1).y = (v2).y + (s3) * (v3).y,                         \
   (v1).z = (v2).z + (s3) * (v3).z
#define V_S_S_ADD(v1, s2, v2, s3, v3)                          \
   (v1).x = (s2) * (v2).x + (s3) * (v3).x,                  \
   (v1).y = (s2) * (v2).y + (s3) * (v3).y,                  \
   (v1).z = (s2) * (v2).z + (s3) * (v3).z
#define V_DOT(v1, v2)                                        \
   ((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)

#define V_PROD(v)                                            \
   ((v).x * (v).y * (v).z)

#define V_LINEAR(p,s)						\
	(((p).z * (s).y + (p).y) * (s).x + (p).x)

#define V_SET_ALL(v, s)                                       \
   V_SET (v, s, s, s)

#define V_C_SUM(v)                                            \
   ((v).x + (v).y + (v).z)

#endif

#define V_ZERO(v)  V_SET_ALL (v, 0)
#define V_LEN_SQ(v)  V_DOT (v, v)
#define V_LEN(v)  sqrt (V_DOT (v, v))
#define V_V_ADD(v1, v2)  V_ADD (v1, v1, v2)
#define V_V_SUB(v1, v2)  V_SUB (v1, v1, v2)
#define V_V_S_ADD(v1, s2, v2) V_S_ADD (v1, v1, s2, v2)

#endif

