/*
 * --------------------------------------------------------------------------------
 * Project       : Gillespie Algorithm Simulation for Gene Regulatory Networks
 * Author        : ROBERT Corentin
 * Institution   : Université Libre de Bruxelles (ULB)
 * Date          : 08/10/2024
 * Language      : C
 * --------------------------------------------------------------------------------
 *
 * Description   :
 * This code implements the Gillespie Stochastic Simulation Algorithm to model
 * the dynamics of a gene regulatory network, focusing on the robustness of differentiation
 * processes and the influence of noise within a genetic regulatory system. The model
 * includes multiple reaction pathways, represented by rate equations that depend on 
 * several parameters and variables describing gene activation and degradation.
 *
 * --------------------------------------------------------------------------------
 */

#include "mainlib.h"

// ======================================================================================================== //
										/* VARIABLE DECLARATION */							
// ======================================================================================================== //

/*--------------------- Simulation Counters and State Variables ---------------------*/

int i = 0, j = 0, k = 0;                     // Iterators for loops.
int compteur = 0, comptMax = 0;              // General counters.
int comptMaxKDE = 0;                         // Maximum count for KDE analysis.
int choix = 0;                               // User choice variable.
int comptTraj = 0, comptTrajMax = 0;         // Trajectory count variables.
int choix_discr = 0;                         // Discretization choice variable.

/*------------------------ Time Tracking Variables ------------------------*/

double temps = 0, tMax = 0;                  // Current time and maximum time.
double dt = 0, dtImpr = 0;                   // Time steps for simulation and output.
double tempsImpr = 0;                        // Time for output printing.
double tempsDistr = 0;                       // Time for distribution sampling.
int tprintMax = 0;                           // Maximum time print intervals.
int t = 0;                                   // Discrete time step count.
double tSTAT = 0;                            // Stationary state threshold.

/*----------------------- Model State and Parameters ----------------------*/

int x = 0, y = 0;                                    // Molecular counts for two species.
double cix = 0, ciy = 0;                             // Initial concentrations.

double FI = 0, FA = 0;                               // Activation/inhibition rates.
double kappa = 0;                                    // Dissociation constant ratio.
double dI = 0, dA = 0, dKI = 0, dKA = 0, dD = 0;     // Perturbation terms.
int hl = 0;                                          // Hill coefficient.


/*----------------- Gillespie Algorithm Intermediate Variables ----------------*/

double r1 = 0, r2 = 0;                       // Random variables for reaction choices.
double w1 = 0., w2 = 0., w3 = 0.;            // Reaction propensity functions.
double w4 = 0., w5 = 0., w6 = 0.;            // Reaction propensity functions.
double a0 = 0.;                              // Total propensity function.
double Omega = 0.;                           // System size.


/*---------------- Kernel Density Estimation Variables ----------------*/

double test1 = 0;                            // Test variable for KDE analysis.
int nx = 0, xMin = 0, xMax = 0;              // KDE x-axis parameters.
double xValue = 0, dx = 0;                   // KDE x-axis calculations.
int ny = 0, yMin = 0, yMax = 0;              // KDE y-axis parameters.
double yValue = 0, dy = 0;                   // KDE y-axis calculations.

double sumydata = 0., avydata = 0.;          // KDE y-data aggregation and averages.
double avydata_center = 0.;                  // KDE y-data centered average.
double stdv1ydata = 0., stdv2ydata = 0.;     // KDE y-data variance measures.

double sumxdata = 0., avxdata = 0.;          // KDE x-data aggregation and averages.
double avxdata_center = 0.;                  // KDE x-data centered average.
double stdv1xdata = 0., stdv2xdata = 0.;     // KDE x-data variance measures.

double covxx = 0, covxy = 0;                 // Covariance terms for KDE.
double covyy = 0, covyx = 0;                 // Covariance terms for KDE.

double hx = 0., hy = 0.;                     // Bandwidths for KDE.
int d = 0;                                   // Dimensionality variable for KDE.

double normal_termKDE = 0;                   // KDE normalization factor.

bool A_STATE_COND = 0;                       // Condition for state A.
bool B_STATE_COND = 0;                       // Condition for state B.

double AS, BS, CS;                           // Proportion metrics for different states.
int comptST;                                 // Counter for steady-state.
int condA = 0, condB = 0;                    // State-specific condition flags.

double compteur_test = 0;                    // Counter for tests.
double step_IAcolor = 0;                     // Step size for proportion colors.
double param1 = 0, param2 = 0;               // Parameter sweeps for robustness.
double step_d = 0;                           // Step size for parameter sweeps.

void GILLESPIE_PT1()
{
	// Define conditions for A and B states based on concentration ratios and threshold factors.

	A_STATE_COND = (x / Omega) / (y / Omega) >= 10 && x / Omega >= FI * (1 + dI) / (2 * (1 + dD));
	B_STATE_COND = (y / Omega) / (x / Omega) >= 10 && y / Omega >= FI * (1 - dI) / (2 * (1 - dD));
	
	// Initialize the first random parameter for time step selection.
	
	r1 = drand48();
	
	// Initialize the second random parameter for reaction selection.
	
	r2 = drand48();
	
	// Calculate reaction propensities based on the current system state.
	
	//-------------- X -> X+1 ----------------//
	
	if (FI == 0 || dI == -1 || dKI == -1)
	{
		w1 = 0; 
	}
	else
	{
		w1 = (FI * (1 + dI) * Omega * pow(1 + dKI, hl)) / (pow(1 + dKI, hl) + pow((y / Omega), hl));
	}
	
	if (kappa == 0 || dKA == -1)
	{
		w3 = FA * (1 + dA) * Omega;
	}
	else
	{
		w3 = (FA * (1 + dA) * Omega * pow((x / Omega) ,hl)) / ((pow(kappa, hl) * pow(1 + dKA, hl)) + pow((x / Omega), hl));
	}
	
	//------------- X -> X-1 -----------------//
	
	w5 = (1 + dD) * x;
	
	//------------- Y -> Y+1 ----------------//
	
	if (FI == 0 || dI == 1 || dKI == 1)
	{
		w2 = 0; 
	}
	else 
	{
		w2 = (FI * (1 - dI) * Omega * pow(1 - dKI, hl)) / (pow(1 - dKI, hl) + pow((x / Omega), hl));
	}
	
	if (kappa == 0 || dKA == 1)
	{
		w4 = FA * Omega * (1 - dA);
	}
	else
	{
		w4 = (FA * (1 - dA) * Omega * pow((y / Omega), hl)) / ((pow(kappa, hl)*pow(1 - dKA, hl)) + pow((y / Omega),hl)) ;
	}
	
	//------------- Y -> Y-1 -----------------//
	
	w6 = (1 - dD) * y;
	
	// Compute the total propensity.
	
	a0 = w1 + w2 + w3 + w4 + w5 + w6;
}


void GILLESPIE_PT2()
{
	// Update time using the total propensity a0.
	
	temps -= (1/a0)*log(r1); 
	
	// Determine which reaction occurs based on the cumulative probabilities
	
	if (r2 >= 0 && r2 <= w1 / a0)  // Condition for reaction 1: X -> X+1
	{
		x = x + 1;
	}
	else if (r2 <= (w1 + w2) / a0)  // Condition for reaction 2: Y -> Y+1
	{
		y = y + 1;
	}
	else if (r2 <= (w1 + w2 + w3) / a0)  // Condition for reaction 3: X -> X+1 (alt)
	{
		x = x + 1;
	}
	else if (r2 <= (w1 + w2 + w3 + w4) / a0)  // Condition for reaction 4: Y -> Y+1 (alt)
	{
		y = y + 1;
	}
	else if (r2 <= (w1 + w2 + w3 + w4 + w5) / a0)  // Condition for reaction 5: X -> X-1
	{
		x = x - 1;
	}
	else if (r2 <= (w1 + w2 + w3 + w4 + w5 + w6) / a0)  // Condition for reaction 6: Y -> Y-1
	{
		y = y - 1;
	}
	else 
	{
		// No action needed; this branch should never occur given proper r2 values
	}
}

// ======================================================================================================== //
															/* TEMPORAL SERIES */							
// ======================================================================================================== //


void ALG_TEMPORALSERIES_CONC_PROP()
{	
	tprintMax = floor(tMax/dtImpr);
	
	//============================================================================================//
	
	double *A = NULL; double *B = NULL; double *C = NULL;
	double *avA = NULL; double *avB = NULL; double *avC = NULL;
	
	int n1 = 1; size_t size1[] = {tprintMax + 1};
	int n2 = 1; size_t size2[] = {tprintMax + 1};
	int n3 = 1; size_t size3[] = {tprintMax + 1};
	int n4 = 1; size_t size4[] = {tprintMax + 1};
	int n5 = 1; size_t size5[] = {tprintMax + 1};
	int n6 = 1; size_t size6[] = {tprintMax + 1};
	
	createNDimensionalArray((void***)&A, n1, size1, sizeof(long double));
	createNDimensionalArray((void***)&B, n2, size2, sizeof(long double));
	createNDimensionalArray((void***)&C, n3, size3, sizeof(long double));
	createNDimensionalArray((void***)&avA, n4, size4, sizeof(long double));
	createNDimensionalArray((void***)&avB, n5, size5, sizeof(long double));
	createNDimensionalArray((void***)&avC, n6, size6, sizeof(long double));
	
	if (A == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (B == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (C == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (avA == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (avB == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (avC == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	
	//============================================================================================//
	
	DRAND48();
	
	for (compteur = 0; compteur <= comptMax ; compteur ++)
	{
		x = cix * Omega;
		y = ciy * Omega;
		
		tempsImpr = 0;
		temps = 0;
		t = 0;
		
		printf("Count: %d/%d\n", compteur,comptMax);
		
		while(temps<=tMax)
		{
			GILLESPIE_PT1();
			
			if (temps >= tempsImpr)
			{
				if (temps > tempsImpr + 1)
				{
					printf("ERROR : Desynchronised times. => Try a bigger dtImpr\n");
					exit(0);
				}
				else{;}
				
				if (compteur <= comptST)
				{
					fprintf(fp[1],"%f\t %f\t %f\t %f\t %f\n",r1, r2, temps, x/Omega, y/Omega);
				}
				else{;}

				tempsImpr += dtImpr;
				t++;
			}
			else{;}
			
			if (A_STATE_COND)
			{
				A[t] = A[t] + 1;
				B[t] = B[t] + 0;
				C[t] = C[t] + 0;
			}
			else if (B_STATE_COND)
			{
				A[t] = A[t] + 0;
				B[t] = B[t] + 1;
				C[t] = C[t] + 0;
			}
			else
			{
				A[t] = A[t] + 0;
				B[t] = B[t] + 0;
				C[t] = C[t] + 1;
			}
		
			GILLESPIE_PT2();
		}
		
		if (compteur <= comptST)
		{
			fprintf(fp[1], "\n");
			fprintf(fp[2],"%f\t %f\n",x/Omega, y/Omega);
		}
	}
	
	for (t = 0 ; t < tprintMax ; t++)
	{
		avA[t] = (double)(A[t])/(comptMax);
		avB[t] = (double)(B[t])/(comptMax);
		avC[t] = (double)(C[t])/(comptMax);
		
		fprintf(fp[3],"%f %f %f %f\n", t*dtImpr, avA[t],avB[t],avC[t]);
	}
	
	freeNDimensionalArray((void*)A, n1, size1);
	freeNDimensionalArray((void*)B, n2, size2);
	freeNDimensionalArray((void*)C, n3, size3);
	freeNDimensionalArray((void*)avA, n4, size4);
	freeNDimensionalArray((void*)avB, n5, size5);
	freeNDimensionalArray((void*)avC, n6, size6);
}
void PROG_TEMPORALSERIES_CONC_PROP()
{
	DIR_FILE("DATA/PRINCIPAL","XY_TIMESERIES.csv",1,0);
	DIR_FILE("DATA/PRINCIPAL","STEADY_STATES.csv",2,0);
	DIR_FILE("DATA/PRINCIPAL","ABC_TIMESERIES.csv",3,0);
	
	ALG_TEMPORALSERIES_CONC_PROP();
	
	fclose(fp[1]);
	fclose(fp[2]);
	fclose(fp[3]);
}

void ALG_STEADY_STATES_POTLANDSCAPE() 
{
	double *dataxmy = NULL; double *dataxpy = NULL;
	double *PCA1_KDE = NULL; double *PCA2_KDE = NULL;
	int n1 = 1; size_t size1[] = {comptMax + 1};
	int n2 = 1; size_t size2[] = {comptMax + 1};
	int n3 = 1; size_t size3[] = {comptMax + 1};
	int n4 = 1; size_t size4[] = {comptMax + 1};
	createNDimensionalArray((void***)&dataxmy, n1, size1, sizeof(long double));
	createNDimensionalArray((void***)&dataxpy, n2, size2, sizeof(long double));
	createNDimensionalArray((void***)&PCA1_KDE, n3, size3, sizeof(long double));
	createNDimensionalArray((void***)&PCA2_KDE, n4, size4, sizeof(long double));
	if (dataxmy == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (dataxpy == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (PCA1_KDE == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (PCA2_KDE == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	
	DRAND48();
	
	for (compteur = 0; compteur <= comptMax ; compteur ++)
	{
		x = cix * Omega; 
		y = ciy * Omega;
		temps = 0; 
		
		if (compteur % 100 == 0)
		{
			printf (" Compteur: %d/%d\n", compteur,comptMax);
		}
		
		while(temps<=tMax)
		{
			GILLESPIE_PT1();
			GILLESPIE_PT2();
		}
		
		dataxmy[compteur] = (x/Omega) - (y/Omega); 
		dataxpy[compteur] = (x/Omega) + (y/Omega); 
		
		fprintf(fp[1],"%f\t %f\n",x/Omega, y/Omega);
		
	}
	
	xMin = -20, xMax = 20, dx = 0.05, hx = -1; 
	nx = floor((xMax - xMin)/dx);
	
	double* KDEPCA1 = KDE1D(dataxmy, comptMax, xMin, xMax, dx, hx);
	double* KDEPCA2 = KDE1D(dataxpy, comptMax, xMin, xMax, dx, hx);

	for (j = 0; j < nx; j++) 
	{
		xValue = (j*dx) - abs(xMin);
	
		fprintf(fp[2], "%e\t%e\t%e\n", xValue, KDEPCA1[j], KDEPCA2[j]);
	}
	
	freeKDE1D(KDEPCA1);
	freeKDE1D(KDEPCA2);
	freeNDimensionalArray((void*)dataxmy, n1, size1);
	freeNDimensionalArray((void*)dataxpy, n2, size2);
	freeNDimensionalArray((void*)PCA1_KDE, n3, size3); 
	freeNDimensionalArray((void*)PCA2_KDE, n4, size4);

}

void PROG_STEADY_STATE_POTLANDSCAPE()
{
	
	DIR_FILE("DATA/PRINCIPAL","STEADY_STATES_DISTR.csv",1,0);
	DIR_FILE("DATA/PRINCIPAL","STEADY_STATES_DISTR_KDE.csv",2,0);

	ALG_STEADY_STATES_POTLANDSCAPE();
	
	fclose(fp[1]);
	fclose(fp[2]);
}


// ======================================================================================================== //
						 									/* ROBUSTNESS (1-parameter) */								
// ======================================================================================================== //

void ALG_ROBUSTNESS_PROP_1P()
{
	DRAND48();
	
	for (param1 = -1; param1 <= 1.01; param1 += step_d)
	{
		dI = param1; 
			
		printf("P1 = %.2f\t FI = %.2f\t FA = %.2f\t K= %.2f\t O= %.0f\n", param1, FI, FA, kappa, Omega);
		
		for (compteur = 0; compteur < comptMax ; compteur ++)
		{
			x = 0;
			y = 0;
			
			temps = 0;
			condA = 0; 
			
			while(temps<=tMax)
			{
				GILLESPIE_PT1 ();
				GILLESPIE_PT2();
			}
			
			if ((x/Omega)/(y/Omega) >= 10 && x/Omega >= FI * (1+dI)/2)
			{
				AS = AS + 1;
				BS = BS + 0;
				CS = CS + 0;
			}
			else if ((y/Omega)/(x/Omega) >= 10 && y/Omega >= FI * (1-dI)/2)
			{
				AS = AS + 0;
				BS = BS + 1;
				CS = CS + 0;
			}
			else
			{
				AS = AS + 0;
				BS = BS + 0;
				CS = CS + 1;
			}
		}
		
		fprintf(fp[1], "%f\t %Lf\t %Lf\t %Lf\n", param1, (long double)(AS/comptMax), (long double)(BS/comptMax), (long double)(CS/comptMax));
		fflush(fp[1]);
		
		AS = BS = CS = 0;
	}
}

void PROG_ROBUSTNESS_PROP_1P()
{
	DIR_FILE("DATA/PROPORTIONS","PROPORTIONS_1P.csv",1,0);
	
	ALG_ROBUSTNESS_PROP_1P();
	
	fclose(fp[1]);
}

// ======================================================================================================== //
									    				/* ROBUSTNESS DIAGRAMS */
// ======================================================================================================== //

void ALG_ROBUSTNESS_DIAGRAM_FAK()
{
	
	DRAND48();
	
	compteur_test = floor(2./(step_d));
	
	//============================================================================================//
	
	double **Array = NULL;
	double *mean_param1 = NULL;
	double *var_param1 = NULL;
	
	int n1 = 2; size_t size1[] = {compteur_test + 1,6};
	int n2 = 1; size_t size2[] = {6};
	int n3 = 1; size_t size3[] = {6};
	
	createNDimensionalArray((void***)&Array, n1, size1, sizeof(long double));
	createNDimensionalArray((void***)&mean_param1, n2, size2, sizeof(long double));
	createNDimensionalArray((void***)&var_param1, n3, size3, sizeof(long double));
	
	if (Array == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (mean_param1 == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (var_param1 == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	
	//============================================================================================//
	
	for (FA = 0; FA <= 10; FA += step_IAcolor)
	{
		for (kappa = 0; kappa <= 10; kappa += step_IAcolor)
		{
			i = 0;
				
			for (param1 = -1; param1 <= 1.01; param1 += step_d)
			{
				dI = param1;
				
				printf("DIFI10_O100 %.2f\t %.2f\t (%d)\t k= %.2f\t %.0f\t %f\n", FI, FA, i, kappa, Omega, param1);
				
				for (compteur = 0; compteur < comptMax ; compteur ++)
				{
					x = 0; 
					y = 0;
					
					temps = 0; 
					
					while(temps<=tMax)
					{
						GILLESPIE_PT1 ();
						GILLESPIE_PT2 ();
					}
					
					if ((x/Omega)/(y/Omega) >= 10 && x/Omega >= FI*(1+dI)/2)
					{
						AS = AS + 1;
						BS = BS + 0;
						CS = CS + 0;
					}
					else if ((y/Omega)/(x/Omega) >= 10 && y/Omega >= FI*(1-dI)/2)
					{
						AS = AS + 0;
						BS = BS + 1;
						CS = CS + 0;
					}
					else  
					{
						AS = AS + 0;
						BS = BS + 0;
						CS = CS + 1;
					}
				}
				
				Array[i][0] = (long double)(AS)/(long double)(comptMax);
				Array[i][1] = (long double)(BS)/(long double)(comptMax);
				Array[i][2] = (long double)(CS)/(long double)(comptMax);
				Array[i][3] = fabs((double)(AS)/(double)(comptMax) - (double)(BS)/(double)(comptMax));
				Array[i][4] = param1;
				
				AS = BS = CS = 0;
				i++;
			}
			
			for (i = 0; i <= compteur_test; i++)
			{
				mean_param1[0] = 0;
				mean_param1[1] = 0;
				mean_param1[2] = 0;
				mean_param1[3] = 0;
				
				var_param1[0] = 0;
				var_param1[1] = 0;
				var_param1[2] = 0;
				var_param1[3] = 0;
			}
			
			for (i = 0; i <= compteur_test; i++)
			{
				mean_param1[0] += Array[i][0]/(compteur_test+1);
				mean_param1[1] += Array[i][1]/(compteur_test+1);
				mean_param1[2] += Array[i][2]/(compteur_test+1);
				mean_param1[3] += Array[i][3]/(compteur_test+1);
			}
			
			for (i = 0; i <= compteur_test; i++)
			{
				var_param1[0] += pow(Array[i][0] - mean_param1[0],2)/(compteur_test+1);
				var_param1[1] += pow(Array[i][1] - mean_param1[1],2)/(compteur_test+1);
				var_param1[2] += pow(Array[i][2] - mean_param1[2],2)/(compteur_test+1);
				var_param1[3] += pow(Array[i][3] - mean_param1[3],2)/(compteur_test+1);
			} 
			printf("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", FI, FA, kappa, mean_param1[0], mean_param1[1], mean_param1[2], mean_param1[3], var_param1[0], var_param1[1], var_param1[2], var_param1[3]);
			fprintf(fp[1],"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", FI, FA, kappa, mean_param1[0], mean_param1[1], mean_param1[2], mean_param1[3], var_param1[0], var_param1[1], var_param1[2], var_param1[3]);
			
			mean_param1[0] = 0; mean_param1[1] = 0; mean_param1[2] = 0; mean_param1[3] = 0;
			
			var_param1[0] = 0; var_param1[1] = 0; var_param1[2] = 0; var_param1[3] = 0;
		}
		
	}
	
	freeNDimensionalArray((void*)Array, n1, size1);
	freeNDimensionalArray((void*)mean_param1, n2, size2);
	freeNDimensionalArray((void*)var_param1, n3, size3);
}

void PROG_ROBUSTNESS_DIAGRAM_FAK()
{
	DIR_FILE("DATA/PROPORTIONS","ROBUSTNESS_DIAGRAM_FAK.csv",1,0);
	
	ALG_ROBUSTNESS_DIAGRAM_FAK();
	
	fclose(fp[1]);
	
}

// ======================================================================================================== //
														/* PROPORTION DIAGRAMS */
// ======================================================================================================== //

void ALG_PROPORTIONS_FAK()
{
	tprintMax = floor(tMax/dtImpr);
	
	DRAND48();
	
	for (FA = 0; FA <= 10; FA += step_IAcolor)
	{
		for (kappa = 0; kappa <= 10; kappa += step_IAcolor)
		{
			printf("%.2f\t %.2f\t k= %.2f\t %.0f\n", FI, FA, kappa, Omega);
			
			for (compteur = 0; compteur < comptMax ; compteur ++)
			{
				x = 0; 
				y = 0;
				
				temps = 0;  
				
				while(temps<=tMax)
				{
					GILLESPIE_PT1 ();
					GILLESPIE_PT2 ();
				}
				
				if ((x/Omega)/(y/Omega) >= 10 && x/Omega >= FI*(1+dI)/2)
				{
					AS = AS + 1;
					BS = BS + 0;
					CS = CS + 0;
				}
				else if ((y/Omega)/(x/Omega) >= 10 && y/Omega >= FI*(1-dI)/2)
				{
					AS = AS + 0;
					BS = BS + 1;
					CS = CS + 0;
				}
				else
				{
					AS = AS + 0;
					BS = BS + 0;
					CS = CS + 1;
				}
			}
			
			fprintf(fp[1], "%f\t %f\t %Lf\t %Lf\t %Lf\n", FA, kappa, (long double)(AS)/(long double)(comptMax), (long double)(BS)/(long double)(comptMax), (long double)(CS)/(long double)(comptMax));
			fflush(fp[1]);
			
			AS = BS = CS = 0; 
		}
	}
}

void PROG_PROPORTIONS_FAK()
{
	DIR_FILE("DATA/PROPORTIONS","PROPORTIONS_FAK.csv",1,0);
	
	ALG_PROPORTIONS_FAK();
	
	fclose(fp[1]);
}

// ======================================================================================================== //
										/* PROBABILITY DISTRIBUTION */
// ======================================================================================================== //

void ALG_PROBA_DISTR()
{
	//============================================================================================//
	
	double *datax = NULL; double *datay = NULL;
	double *PCA1_KDE = NULL; double *PCA2_KDE = NULL;
	int n1 = 1; size_t size1[] = {comptMax + 1};
	int n2 = 1; size_t size2[] = {comptMax + 1};
	int n6 = 1; size_t size6[] = {comptMax + 1};
	int n7 = 1; size_t size7[] = {comptMax + 1};
	createNDimensionalArray((void***)&datax, n1, size1, sizeof(long double));
	createNDimensionalArray((void***)&datay, n2, size2, sizeof(long double));
	createNDimensionalArray((void***)&PCA1_KDE, n6, size6, sizeof(long double));
	createNDimensionalArray((void***)&PCA2_KDE, n7, size7, sizeof(long double));
	if (datax == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (datay == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (PCA1_KDE == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	if (PCA2_KDE == NULL){fprintf(stderr, "Error allocating memory for the array\n");exit(1);}
	
	//============================================================================================//
	
	DRAND48();
	
	for (compteur = 0; compteur < comptMax ; compteur ++)
	{
		
		x = cix * Omega; 
		y = ciy * Omega;
		printf (" Compteur: %d/%d\n", compteur,comptMax);
			
		temps = 0.0;
			
		while(temps <= tMax)
		{
			GILLESPIE_PT1 ();
				
			if (temps >= tempsDistr)
			{
				datax[compteur] = x/Omega;
				datay[compteur] = y/Omega;
				fprintf(fp[1], "%e\t%e\t%e\n",temps,datax[compteur], datay[compteur]);
				break;
			}
			else{;}
			
			GILLESPIE_PT2 ();
		}
	}
	
	//---------------------------------------------------------------------------------------------//
												// KDE //
	//---------------------------------------------------------------------------------------------//
	
	xMin = 0, xMax = 20, yMin = 0, yMax = 20, dx = 0.1, dy = 0.1, hx = -1, hy = -1; 
	nx = floor((xMax - xMin)/dx);  ny = floor((yMax - yMin)/dy);
	double** KDE = KDE2D(datax, datay, comptMax, xMin, xMax, yMin, yMax, dx, dy, hx, hy);
	
	double maxvalKDE = 0;
	double testvalKDE = 0;  
	
	for (j = 0; j <= nx; j++) 
	{
		for (k = 0; k <= ny; k++) 
		{
			xValue = (j*dx) - abs(xMin);
			yValue = (k*dy) - abs(yMin);
			
			testvalKDE = KDE[j][k];
			
			if (maxvalKDE < testvalKDE)
			{
				maxvalKDE = testvalKDE;
			}
			 
			if (KDE[j][k] < 1E-10)
			{
				KDE[j][k] = 1E-10;
			}
			fprintf(fp[2], "%e\t%e\t%e\n",xValue, yValue, KDE[j][k]);
		}
	}
	
	freeKDE2D(KDE, dx, dy, xMin, xMax, yMin, yMax);
	
	//---------------------------------------------------------------------------------------------//
												// PCA //
	//---------------------------------------------------------------------------------------------//
	
	PCA_STRUCT* res = PCA(datax, datay, comptMax);
	
	for (i = 0; i < comptMax; i++)
	{
		PCA1_KDE[i] = res->PCA_DATA[i][0]; 
		PCA2_KDE[i] = res->PCA_DATA[i][1];
		
		fprintf(fp[3], "%e\t%e\n", res->PCA_DATA[i][0], res->PCA_DATA[i][1]);
	}
	
	printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", res->PCA_EVA_DATA[0], res->PCA_EVP_DATA[0][0], res->PCA_EVP_DATA[0][1], res->PCA_EVA_DATA[1], res->PCA_EVP_DATA[1][0], res->PCA_EVP_DATA[1][1]);
	fprintf(fp[5],"%e\t%e\t%e\t%e\t%e\t%e\n", res->PCA_EVA_DATA[0], 0., 0., res->PCA_EVA_DATA[1], 0., 0.);
	fprintf(fp[5],"%e\t%e\t%e\t%e\t%e\t%e\n", res->PCA_EVA_DATA[0], res->PCA_EVP_DATA[0][0], res->PCA_EVP_DATA[0][1], res->PCA_EVA_DATA[1], res->PCA_EVP_DATA[1][0], res->PCA_EVP_DATA[1][1]);
	
	// KDE DES HISTOGRAMMES DES DONNEES //
	
	xMin = -10, xMax = 10, dx = 0.01, hx = -1; 
	nx = floor((xMax - xMin)/dx);
	
	double* KDEPCA1 = KDE1D(PCA1_KDE, comptMax, xMin, xMax, dx, hx);
	double* KDEPCA2 = KDE1D(PCA2_KDE, comptMax, xMin, xMax, dx, hx);
	
	for (j = 0; j < nx; j++) 
	{
		xValue = (j*dx) - abs(xMin);
		
		fprintf(fp[4], "%e\t%e\t%e\n", xValue, KDEPCA1[j], KDEPCA2[j]);
	}
	
	freeKDE1D(KDEPCA1);
	freeKDE1D(KDEPCA2);
	freePCA(res,comptMax);
	
	freeNDimensionalArray((void*)datax, n1, size1);
	freeNDimensionalArray((void*)datay, n2, size2);
	freeNDimensionalArray((void*)PCA1_KDE, n6, size6); 
	freeNDimensionalArray((void*)PCA2_KDE, n7, size7);
}

void PROG_PROBA_DISTR()
{
	DIR_FILE("DATA/PRINCIPAL","PROBA_DISTR.csv",1,0);
	DIR_FILE("DATA/PRINCIPAL","PROBA_DISTR_KDE.csv",2,0);
	DIR_FILE("DATA/PRINCIPAL","PCA.csv",3,0);
	DIR_FILE("DATA/PRINCIPAL","KDEPCA.csv",4,0);
	DIR_FILE("DATA/PRINCIPAL","PCA_EIGENVALUE_VECTOR.csv",5,0);
	
	ALG_PROBA_DISTR();
	
	fclose(fp[1]);
	fclose(fp[2]);
	fclose(fp[3]);
	fclose(fp[4]);
	fclose(fp[5]);
}


// ======================================================================================================== //
											/* TRANSITIONS */
// ======================================================================================================== //


void ALG_TRANSITIONS_AB_FAK()
{
	
	DRAND48();
	int transcomptA = 0;
	int transcomptB = 0;
	int transA = 0, transB = 0;
	double PtransA = 0, PtransB = 0; 
	
	for (FA = 0; FA <= 10; FA += step_IAcolor)
	{
		for (kappa = 0; kappa <= 10; kappa += step_IAcolor)
		{
			printf("PAB SYM %.2f\t %.2f\t k= %.2f\t %.0f\n", FI, FA, kappa, Omega);
			
			for (compteur = 0; compteur < comptMax ; compteur ++)
			{
				x = 0;
				y = 0;
				temps = 0;
				transA = transB = 0;
				
				while(temps<=tMax)
				{
					GILLESPIE_PT1 ();
						
					if (temps >= tSTAT)
					{
						if ((x/Omega)/(y/Omega) >= 10 && x/Omega >= FI*(1+dI)/2)
						{
							if (transA == 0 && transB == 1)
							{
								transcomptA += 1;
							}
							transA = 1; transB = 0;
						}
						else if ((y/Omega)/(x/Omega) >= 10 && y/Omega >= FI*(1-dI)/2)
						{
							if (transA == 1 && transB == 0)
							{
								transcomptB += 1;
							}
							transA = 0; transB = 1;
						}
					}
					else{;}
						
					GILLESPIE_PT2 ();
					
				}
				
				PtransA += (double)(transcomptA) / ((double)(comptMax)*(double)(tMax-tSTAT));  
				PtransB += (double)(transcomptB) / ((double)(comptMax)*(double)(tMax-tSTAT));
			}
		
			fprintf(fp[1], "%e\t%e\t%e\t%e\n", FA, kappa, PtransA, PtransB);
		
			PtransA = PtransB = transcomptA = transcomptB = 0; 
		}
	}
}

void PROG_TRANSITIONS_AB_FAK()
{
	DIR_FILE("DATA/TRANSITIONS","AB_TRANSITIONS_FAK.csv",1,0);
	
	ALG_TRANSITIONS_AB_FAK();
	
	fclose(fp[1]);
	
}

void ALG_TRANSITIONS_NC_FAK()
{
	
	DRAND48();
	
	double cond1 = 0;
	double cond2 = 0;
	double AS1 = 0, AS2 = 0, BS1 = 0, BS2 = 0, CS1 = 0, CS2 = 0;
	

	for (FA = 0; FA <= 10; FA += step_IAcolor)
	{
		for (kappa = 0; kappa <= 10; kappa += step_IAcolor)
		{
			printf("DD %.2f\t %.2f\t k= %.2f\t %.0f\n", FI, FA, kappa, Omega);
			
			for (compteur = 0; compteur < comptMax ; compteur ++)
			{
				x = 0;
				y = 0;
				temps = 0;
				cond1 = cond2 = 0;
				
				while(temps<=tMax)
				{
					GILLESPIE_PT1 ();
						
					if (temps >= 20 && cond1 == 0)
					{
						if ((x/Omega)/(y/Omega) >= 10 && x/Omega >= FI*(1+dI)/2)
						{
							AS1 = AS1 + 1;
							BS1 = BS1 + 0;
							CS1 = CS1 + 0;
						}
						else if ((y/Omega)/(x/Omega) >= 10 && y/Omega >= FI*(1-dI)/2)
						{
							AS1 = AS1 + 0;
							BS1 = BS1 + 1;
							CS1 = CS1 + 0;
						}
						else  
						{
							AS1 = AS1 + 0;
							BS1 = BS1 + 0;
							CS1 = CS1 + 1;
						}
						
						cond1 = 1;
					}
					else if (temps >= 200 && cond2 == 0)
					{
						if ((x/Omega)/(y/Omega) >= 10 && x/Omega >= FI*(1+dI)/2)
						{
							AS2 = AS2 + 1;
							BS2 = BS2 + 0;
							CS2 = CS2 + 0;
						}
						else if ((y/Omega)/(x/Omega) >= 10 && y/Omega >= FI*(1-dI)/2)
						{
							AS2 = AS2 + 0;
							BS2 = BS2 + 1;
							CS2 = CS2 + 0;
						}
						else
						{
							AS2 = AS2 + 0;
							BS2 = BS2 + 0;
							CS2 = CS2 + 1;
						}
						cond2 = 1;
						break; 
					}
					else{;}
						
					GILLESPIE_PT2 ();
				}
			}
				
			fprintf(fp[1], "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", FI, FA, kappa, (double)(AS1)/comptMax, (double)(BS1)/comptMax, (double)(CS1)/comptMax, (double)(AS2)/comptMax, (double)(BS2)/comptMax, (double)(CS2)/comptMax); 
			
			AS1 = BS1 = CS1 = 0;
			AS2 = BS2 = CS2 = 0;
		}
	}
}

void PROG_TRANSITIONS_NC_FAK()
{
	DIR_FILE("DATA/TRANSITIONS","NC_TRANSITIONS_FAK.csv",1,0);
	
	ALG_TRANSITIONS_NC_FAK();
	
	fclose(fp[1]);
	
}



// ======================================================================================================== //
									/* PARAMETRES NUMERIQUES DES PROGRAMMES */													
// ======================================================================================================== //

void START_TEMPORALSERIES_CONC_PROP()
{
	comptMax = 100;
	comptST = 50; 
	
	tMax = 100;
	dtImpr = 0.1;

	
	PROG_TEMPORALSERIES_CONC_PROP();
}

void START_STEADY_STATE_POTLANDSCAPE()
{
	comptMax = 5000;
	
	tMax = 100;
	
	PROG_STEADY_STATE_POTLANDSCAPE();
}

//-----------------------------------------------//
			/* PROPORTIONS - TDIFF */
// ----------------------------------------------//

void START_ROBUSTNESS_PROP_1P()
{
	tMax = 100;
	tSTAT = 20; 
	comptMax = 200;
	step_d = 0.05;
	
	PROG_ROBUSTNESS_PROP_1P();
}

void START_ROBUSTNESS_DIAGRAM_FAK()
{
	tMax = 100;
	tSTAT = 20; 
	comptMax = 200;
	step_d = 0.05;
	step_IAcolor = 0.1;
	
	//GOOD compt = 200 (linéaire); step_d = 0.05 (linéaire); step_IAcolor = 0.10 (carré); tSTAT = 20 (linéaire) (ARTICLE)
	
	PROG_ROBUSTNESS_DIAGRAM_FAK();
}

void START_PROPORTIONS_FAK()
{
	tMax = 100;
	tSTAT = 20; 
	comptMax = 200;
	step_d = 0.05;
	step_IAcolor = 0.1;
	
	//GOOD compt = 200 (linéaire); step_d = 0.05 (linéaire); step_IAcolor = 0.10 (carré); tSTAT = 20 (linéaire) (ARTICLE)
	
	PROG_PROPORTIONS_FAK();
}


//-----------------------------------------------//
			/* PCA/POTENTIELS */
// ----------------------------------------------//

void START_PROBA_DISTR()
{
	tMax = 500;
	tempsDistr = 20;
	
	comptMax = 20000;
	
	PROG_PROBA_DISTR();
}

//-----------------------------------------------//
					/* TRANSITIONS */
// ----------------------------------------------//

void START_TRANSITIONS_AB_FAK()
{
	tMax = 500;
	tSTAT = 20; 
	comptMax = 200;
	step_IAcolor = 0.1;

	PROG_TRANSITIONS_AB_FAK();
	
	//GOODFAK : compt = 200 (linéaire); step_IAcolor = 0.1 (carré); "tSTAT = 20 ; tMax = 500
}

void START_TRANSITIONS_NC_FAK()
{
	tMax = 300;
	tempsImpr = 20; 
	comptMax = 1000;
	step_IAcolor = 0.1;
	
	//GOOD : compt = 1000 (linéaire); step_IAcolor = 0.1 (carré); "t_test" = 20&200
	
	PROG_TRANSITIONS_NC_FAK();
}

// ======================================================================================================== //
													/* MAIN */													
// ======================================================================================================== //

int main(void)
{
	//----------------------------------------------//
	
	cix = 0;
	ciy = 0;
	
	// CORE PARAMETERS //
	
	FI = 5;
	FA = 0;
	kappa = 1;
	hl = 4;
	
	// ASYMMETRY PARAMETERS //
	
	dI = 0.0;
	dA = 0.0;
	dKI = 0.0;
	dKA = 0.0;
	dD = 0.0;
	
	// SYSTEM SIZE //
	
	Omega = 10;
	
	//----------------------------------------------//
	
	// Choose the "START_..." program and call the function to start the simulation
	//START_... 
	//START_STEADY_STATE_POTLANDSCAPE();
	
	printf("Hello, it's working properly !\n");
	
	return 0;
}

//==========================================================================================================//
													/* END */													
//==========================================================================================================//









