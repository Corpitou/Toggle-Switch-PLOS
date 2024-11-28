#include "math_lib.h"

#pragma GCC diagnostic ignored "-Wunused-variable"


//---------------  MATH FUNCTIONS -----------------// 

long double GAUSSIAN_F(long double x, long double mean, long double var)
{
	return exp(-pow((x-mean),2)/(2*var))/(sqrt(2*PI*var));
}

//---------------  RANDOM NUMBERS -----------------//

double DRAND48()
{
	srand48(time(NULL));
	return drand48();
}

double RANDOM(double min, double max, int precision)
{
    double factor = pow(10.0, (double)precision);
    double range = max - min;
    double random_number = (double)drand48() * range + min;
    random_number = round(random_number * factor) / factor;
    return random_number;
}


//---------------  DATA TRANSFORMATIONS -----------------//


double* KDE1D(double *data, int comptMax, double xMin, double xMax, double dx, double hx) 
{
    int compteur, j;
    double sumdata = 0;
    double avdata;
    double stdv1data = 0;
    double stdv2data;
	int nx = floor((xMax - xMin) / dx);
    double sumKDE[nx];
    double normal_termKDE = 0;
	double test1 = 0; 
    double xValue;
	
	double* KDE = (double*)calloc(nx,sizeof(double));
	
    for (compteur = 0; compteur <= comptMax; compteur++) 
	{
        sumdata += data[compteur];
    }

    avdata = sumdata / (double)comptMax;
	

    for (compteur = 0; compteur <= comptMax; compteur++) 
	{
        stdv1data += pow(data[compteur] - avdata, 2);
    }

    stdv2data = sqrt(stdv1data / (double)comptMax);
	
    if (hx < 0) 
	{
        hx = pow(pow((double)comptMax, (-1. / ((double)1 + 4))) * stdv2data, 2);
    }
	else{;}

    for (j = 0; j < nx; j++) 
	{
        xValue = (j * dx) - fabs(xMin);
        sumKDE[j] = 0;

        for (compteur = 0; compteur < comptMax; compteur++) 
		{
            sumKDE[j] += GAUSSIAN_F(xValue, data[compteur], hx);
        }

        KDE[j] = sumKDE[j] / (double)comptMax;
    }
	
	normal_termKDE = 0; 
	
	for (j = 0; j < nx ; j++)
	{
		normal_termKDE += KDE[j];		
	}
	
	for (j = 0; j < nx ; j++)
	{
		test1 += KDE[j]/normal_termKDE;
	}
	
	printf("KDE NORM = %f\n", test1);
	
	for (j = 0; j < nx ; j++)
	{
		KDE[j] /= 1;
	}

    return KDE;
}
void freeKDE1D(double* KDE) 
{
    free(KDE);
}

double** KDE2D(double *datax, double *datay, int comptMax, double xMin, double xMax, double yMin, double yMax, double dx, double dy, double hx, double hy)
{
    int compteur, j, k;
    double sumxdata = 0, sumydata = 0;
    double avxdata, avydata;
    double stdv1xdata = 0, stdv1ydata = 0;
    double stdv2xdata, stdv2ydata;
    int nx = floor((xMax - xMin) / dx);
    int ny = floor((yMax - yMin) / dy);
	double normal_termKDE = 0;
	double test1 = 0; 
    double sumKDE[nx][ny];
    double xValue, yValue;
	
	double** KDE = (double**)calloc(nx, sizeof(double*));
	for (j = 0; j < nx; j++)
	{
		KDE[j] = (double*)calloc(ny, sizeof(double));
	}
	
    for (compteur = 0; compteur < comptMax; compteur++)
    {
        sumxdata += datax[compteur];
        sumydata += datay[compteur];
    }

    avxdata = sumxdata / (double)comptMax;
    avydata = sumydata / (double)comptMax;
	

    for (compteur = 0; compteur < comptMax; compteur++)
    {
        stdv1xdata += pow(datax[compteur] - avxdata, 2);
        stdv1ydata += pow(datay[compteur] - avydata, 2);
    }

    stdv2xdata = sqrt(stdv1xdata / (double)comptMax);
    stdv2ydata = sqrt(stdv1ydata / (double)comptMax);
	
	
    if ((hx < 0) && (hy < 0))
    {
        hx = pow(pow((double)comptMax, (-1. / (2. + 4.))) * stdv2xdata, 2);
        hy = pow(pow((double)comptMax, (-1. / (2. + 4.))) * stdv2ydata, 2);
    }
	else{;}
	

    for (j = 0; j < nx; j++)
    {
        for (k = 0; k < ny; k++)
        {
            xValue = (j * dx) - fabs(xMin);
            yValue = (k * dy) - fabs(yMin);
            sumKDE[j][k] = 0;

            for (compteur = 0; compteur < comptMax; compteur++)
            {
                sumKDE[j][k] += GAUSSIAN_F(xValue, datax[compteur], hx) * GAUSSIAN_F(yValue, datay[compteur], hy);
            }

            KDE[j][k] = sumKDE[j][k] / (double)comptMax;
        }
    }
	
	
	normal_termKDE = 0; 
	
	for (j = 0; j < nx ; j++)
	{
		for (k = 0; k < ny; k++)
		{
			normal_termKDE += KDE[j][k];		
		}
	}
	
	for (j = 0; j < nx ; j++)
	{
		for (k = 0; k < ny; k++)
		{
			test1 += KDE[j][k]/normal_termKDE;
		}
	}
	
	printf("KDE NORM = %f\n", test1);
	
	for (j = 0; j < nx ; j++)
	{
		for (k = 0; k < ny; k++)
		{
			KDE[j][k] /= normal_termKDE;
		}
	}

    return KDE;
}

void freeKDE2D(double** KDE, double dx, double dy, double xMin, double xMax, double yMin, double yMax) 
{
    int nx = floor((xMax - xMin)/dx);
    int ny = floor((yMax - yMin)/dy);
	
	(void) ny;
	
    for (int j = 0; j < nx; j++) 
	{
        free(KDE[j]);
    }
    free(KDE);
}

PCA_STRUCT* PCA(double *datax, double *datay, int comptMax)
{
	
	PCA_STRUCT* res = calloc(1, sizeof(PCA_STRUCT));
	
	int compteur;
	
    double avxdata = 0.0, avydata = 0.0;
    for (compteur = 0; compteur < comptMax; compteur++) 
	{
        avxdata += datax[compteur] / (double)(comptMax);
        avydata += datay[compteur] / (double)(comptMax);
    }
    
    double stdv1xdata = 0.0, stdv1ydata = 0.0;
    for (compteur = 0; compteur < comptMax; compteur++) 
	{
        stdv1xdata += pow(datax[compteur] - avxdata, 2);
        stdv1ydata += pow(datay[compteur] - avydata, 2);
    }
    
    double stdv2xdata = sqrt(stdv1xdata / (double)(comptMax - 1));
    double stdv2ydata = sqrt(stdv1ydata / (double)(comptMax - 1));
	
	
    double datax_center[comptMax];
    double datay_center[comptMax];
    for (compteur = 0; compteur < comptMax; compteur++) 
	{
        datax_center[compteur] = (datax[compteur] - avxdata) / (stdv2xdata);
        datay_center[compteur] = (datay[compteur] - avydata) / (stdv2ydata);
    }
	
    double avxdata_center = 0.0, avydata_center = 0.0;
    for (compteur = 0; compteur < comptMax; compteur++) 
	{
        avxdata_center += datax_center[compteur] / (double)(comptMax);
        avydata_center += datay_center[compteur] / (double)(comptMax);
    }
	
	
    double covxx = 0.0, covyy = 0.0, covxy = 0.0, covyx = 0.0;
	
	(void) covyx;
	
    for (compteur = 0; compteur < comptMax; compteur++) 
	{
        covxx += ((datax_center[compteur] - avxdata_center) * (datax_center[compteur] - avxdata_center)) / ((double)(comptMax - 1));
        covyy += ((datay_center[compteur] - avydata_center) * (datay_center[compteur] - avydata_center)) / ((double)(comptMax - 1));
        covxy += ((datax_center[compteur] - avxdata_center) * (datay_center[compteur] - avydata_center)) / ((double)(comptMax - 1));
        covyx += ((datay_center[compteur] - avydata_center) * (datax_center[compteur] - avxdata_center)) / ((double)(comptMax - 1));
    }
	
	double covariance[2][2] = {{covxx, covxy}, {covxy, covyy}};
	
    gsl_matrix_view data_matrix = gsl_matrix_view_array(&covariance[0][0], 2, 2);
    
	gsl_vector *EVA_v = gsl_vector_alloc(2);
	gsl_matrix *EVP_m = gsl_matrix_alloc(2, 2);
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(2);
	
	gsl_eigen_symmv(&data_matrix.matrix, EVA_v, EVP_m, workspace);
	
	gsl_eigen_symmv_sort(EVA_v, EVP_m, GSL_EIGEN_SORT_VAL_DESC);
	// i = 0 => PC1
	// i = 1 => PC2
	//...
       
    res->PCA_EVA_DATA[0] = gsl_vector_get(EVA_v, 0);
    res->PCA_EVA_DATA[1] = gsl_vector_get(EVA_v, 1);
	
    for (int i = 0; i < 2; i++)
    {
        gsl_vector_view EVPi_v = gsl_matrix_column(EVP_m, i);

        if (gsl_vector_get(&EVPi_v.vector, 0) > 0)
        {
            gsl_vector_scale(&EVPi_v.vector, -1);
        }

        res->PCA_EVP_DATA[i][0] = gsl_vector_get(&EVPi_v.vector, 0);
        res->PCA_EVP_DATA[i][1] = gsl_vector_get(&EVPi_v.vector, 1);
    }
	
	gsl_matrix *centered_data = gsl_matrix_alloc(comptMax, 2);
	
	for (int i = 0; i < comptMax; i++) 
	{
	    gsl_matrix_set(centered_data, i, 0, datax_center[i]);
	    gsl_matrix_set(centered_data, i, 1, datay_center[i]);
	}
	
	gsl_vector_view EVP1_v = gsl_matrix_column(EVP_m, 0);
	gsl_vector_view EVP2_v = gsl_matrix_column(EVP_m, 1);
	
	gsl_matrix *projection_matrix = gsl_matrix_alloc(2, 2);
	gsl_matrix_set_col(projection_matrix, 0, &EVP1_v.vector);
	gsl_matrix_set_col(projection_matrix, 1, &EVP2_v.vector);
	
	gsl_matrix *transformed_data = gsl_matrix_alloc(comptMax, 2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, centered_data, projection_matrix, 0.0, transformed_data);
	
	res->PCA_DATA = (double**)calloc(comptMax, sizeof(double*));
	
	for (int i = 0; i < comptMax; i++) 
	{
	    res->PCA_DATA[i] = (double*)calloc(2, sizeof(double));
	}
	
    for (int i = 0; i < comptMax; i++) 
	{
        res->PCA_DATA[i][0] = gsl_matrix_get(transformed_data, i, 0);
        res->PCA_DATA[i][1] = gsl_matrix_get(transformed_data, i, 1);
    }
    
    gsl_matrix_free(projection_matrix);
    gsl_matrix_free(transformed_data);
    gsl_matrix_free(centered_data);
    
    gsl_vector_free(EVA_v);
    gsl_matrix_free(EVP_m);
    gsl_eigen_symmv_free(workspace);
    
    return res;
}

void freePCA(PCA_STRUCT *structname, int comptMax) 
{
    for (int i = 0; i < comptMax; i++) 
	{
        free(structname->PCA_DATA[i]);
    }
    free(structname->PCA_DATA);
    
    free(structname);
}