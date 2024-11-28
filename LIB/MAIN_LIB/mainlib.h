#ifndef MAINLIB_H
#define MAINLIB_H

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>

#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef __APPLE__
    #include <sys/time.h> 
    #include <malloc/malloc.h>
#elif defined(__linux__)
    #include <time.h> 
    #include <malloc.h> 
#endif

#include <sys/stat.h>
#include <sys/types.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "GEN_LIB/gen_lib.h"
#include "MATH_LIB/math_lib.h"

#endif /* MAINLIB_H */