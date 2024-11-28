#ifndef math_lib_H
#define math_lib_H

#include "mainlib.h"

#define PI 3.14159265358979323846264338327

double drand48();
void srand48(long int seedvall);
double DRAND48();
double RANDOM(double min, double max, int precision);

long double GAUSSIAN_F(long double x, long double mean, long double var);
	
double* KDE1D(double *data, int comptMax, double xMin, double xMax, double dx, double hx) ;
void freeKDE1D(double* KDE);

double** KDE2D(double *datax, double *datay, int comptMax, double xMin, double xMax, double yMin, double yMax, double dx, double dy, double hx, double hy);
void freeKDE2D(double** KDE, double dx, double dy, double xMin, double xMax, double yMin, double yMax);

typedef struct {double **PCA_DATA;double PCA_EVA_DATA[2];double PCA_EVP_DATA[2][2];} PCA_STRUCT;
PCA_STRUCT* PCA(double *datax, double *datay, int comptMax);
void freePCA(PCA_STRUCT *structname, int comptMax);

#endif /* math_lib_H */