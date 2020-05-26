/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
#include "triangle.h"

void report(struct triangulateio*,int,int,int,int,int,int);
int tricall(int,REAL**,int*,int**,int**,int*,int**);
