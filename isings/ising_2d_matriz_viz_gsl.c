/*****************************************************************************
 *	 		         Ising Model 2D			             *
 *			         Pedro H Mendes	    		             *
 *								             *
 *   	     gcc -Wall ising_new.c -lgsl -lgslcblas -lm -static		     *
 *								             *
 ****************************************************************************/

/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>

/*****************************************************************************
 *                               DEFINITIONS                                 *
 ****************************************************************************/
#define				L	 		16 
#define				L2 	 		(L*L)
#define 			TRAN	 		100000 	//1e5
#define 			TMAX	 		1000000	//1e6
#define 			J	 		1.0

/*****************************************************************************
 *                             GLOBAL VARIABLES                              *
 ****************************************************************************/
int dE, M, ET;

/*****************************************************************************
 *                                 FUNCTIONS                                 *
 ****************************************************************************/
void initialize(gsl_rng *mt, double *boltz, int *spin, int **neigh,double TEMP);
void mc_routine(gsl_rng *mt, double *boltz, int *spin, int **neigh);

/*****************************************************************************
 *                                MAIN PROGRAM                               *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	clock_t t_i, t_f;
	t_i = clock();

	int i, mcs;
	unsigned long int seed;
	char Arq1[100];
	FILE *arq1;
	
	double TEMP, CPU_TIME;
	TEMP = atof(argv[1]);
	
	int *spin, **neigh;
	double *boltz;
	size_t size = L2*sizeof(int); 

	boltz = (double*)malloc(sizeof(double)*9);
	spin = (int*)malloc(size);
	neigh = (int**)malloc(size*sizeof(int*));
	
	for(i=0; i<L2; i++)
	{
		neigh[i] = (int *)malloc(5*sizeof(int));
	}

	srand(time(NULL));
	seed = rand();
	
	gsl_rng *mt;
	mt = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(mt, seed);

	initialize(mt, boltz, spin, neigh, TEMP);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		mc_routine(mt, boltz, spin, neigh);
	}

	sprintf(Arq1, "temp_T%.3lfL%d.dsf", TEMP, L);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#seed = %ld\n#MCS\tM\tET\n", seed);

	for(mcs=0; mcs<TMAX; mcs++)
	{
		mc_routine(mt, boltz, spin, neigh);
		fprintf(arq1, "%d\t%d\t%d\n", mcs, M, ET);
	}

	fclose(arq1);

	free(spin);
	free(neigh);
	free(boltz);

	gsl_rng_free(mt);

	t_f = clock();
	CPU_TIME = (double)(t_f - t_i)/CLOCKS_PER_SEC;

	printf("%lf\n", CPU_TIME);

	return 0;
}

/*****************************************************************************
 *                              INITIALIZATION                               *
 ****************************************************************************/
void initialize(gsl_rng *mt, double *boltz, int *spin, int **neigh, double TEMP)
{
	int i, j;

	boltz[4] = exp(-4.0*J/TEMP);
	boltz[8] = exp(-8.0*J/TEMP);

	for(i=0; i<L2; i++)
	{
		spin[i] = 2*(gsl_rng_uniform_int(mt,2)) - 1;
		
		for(j=0; j<5; j++)
		{
			neigh[i][j] = 0;
		}
	}

	ET = 0.0;
	M = 0.0;

	for(i=0; i<L2; i++)
	{
		neigh[i][0] = (i-L+L2)%L2;
		neigh[i][1] = (i+1)%L + (i/L)*L;
		neigh[i][2] = (i+L)%L2;
		neigh[i][3] = (i-1+L)%L + (i/L)*L;
		
		for(j=0; j<4; j++)
		{
			neigh[i][4] += spin[neigh[i][j]];
		}

		ET += spin[i]*neigh[i][4];
		M += spin[i];
	}

	ET = (ET*(-J))/2.0;

	return;
}

/*****************************************************************************
 *                   	      MONTE CARLO ROUTINE                            *
 ****************************************************************************/
void mc_routine(gsl_rng *mt, double *boltz, int *spin, int **neigh)
{
	int i, j, t;
	double r, prob;
	
	for(t=0; t<L2; t++)
	{
		i = gsl_rng_uniform_int(mt, L2);
	
		dE = 0;

		dE = 2*spin[i]*neigh[i][4];

		if(dE <= 0)
		{
			spin[i] *= -1;
			ET = ET + dE;
			M = M + 2*spin[i];
			
			for(j=0; j<4; j++)
			{
				neigh[neigh[i][j]][4] += 2*spin[i];
			}
		}
		else
		{
			prob = boltz[dE];
			r = gsl_rng_uniform(mt);
		
			if(r < prob)
			{
				spin[i] *= -1;
				ET = ET + dE;
				M = M + 2*spin[i];
	
				for(j=0; j<4; j++)
				{
					neigh[neigh[i][j]][4] += 2*spin[i];
				}
			}
		}
	}

	return;
}
