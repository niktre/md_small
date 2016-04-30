typedef enum {N_I, N_LI, N_R} VType;

#define NAME_I(x)	{#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NAME_LI(x)	{#x, &x, N_LI, sizeof (x) / sizeof (long)}
#define NAME_R(x)	{#x, &x, N_R, sizeof (x) / sizeof (real)}

typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	int vLen, vStatus;
} NameList;
