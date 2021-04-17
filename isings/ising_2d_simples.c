/*****************************************************************************
 *			        Ising Model 2D			             *
 *			        Pedro H Mendes 			             *
 ****************************************************************************/

/*****************************************************************************
 *                             	   INCLUDES                                  *
 ****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"mc.h"

/*****************************************************************************
 *                               DEFINITIONS                                 *
 ****************************************************************************/
#define 			L				64 
#define 			L2 	 			(L*L)
#define 			TRAN				100000 	//1e5
#define 			TMAX				1000000 //1e6
#define 			J				1.0

/*****************************************************************************
 *                           GLOBAL VARIABLES                                *
 ****************************************************************************/
int dE, M, ET;

/*****************************************************************************
 *                              FUNCTIONS                                    *
 ****************************************************************************/
void initialize(double *boltz, int *spin, int *neigh, double TEMP);
void mc_routine(double *boltz, int *spin, int *neigh);

/*****************************************************************************
 *                             MAIN PROGRAM                                  *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	clock_t t_i, t_f;
	t_i = clock();

	int mcs;
	char Arq1[100];
	FILE *arq1;
	
	double TEMP, CPU_TIME;
	TEMP = atof(argv[1]);
	
	int *spin, *neigh;
	double *boltz;
	size_t size = L2*sizeof(int); 

	spin = (int*)malloc(size);
	neigh = (int*)malloc(size);
	boltz = (double*)malloc(sizeof(double)*9);
	
	seed = start_randomic();

	initialize(boltz, spin, neigh, TEMP);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		mc_routine(boltz, spin, neigh);
	}

	sprintf(Arq1, "temp_T%.3lfL%d.dsf", TEMP, L);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#seed = %ld\n#MCS\tM\tET\n", seed);

	for(mcs=0; mcs<TMAX; mcs++)
	{
		mc_routine(boltz, spin, neigh);
		fprintf(arq1, "%d\t%d\t%d\n", mcs, M, ET);
	}

	fclose(arq1);

	free(spin);
	free(neigh);
	free(boltz);

	t_f = clock();
	CPU_TIME = (double)(t_f - t_i)/CLOCKS_PER_SEC;

	printf("%lf\n", CPU_TIME);

	return 0;
}

/*****************************************************************************
 *                             INITIALIZATION                                *
 ****************************************************************************/
void initialize(double *boltz, int *spin, int *neigh, double TEMP)
{
	int i;

	boltz[4] = exp(-4.0*J/TEMP);
	boltz[8] = exp(-8.0*J/TEMP);

	for(i=0; i<L2; i++)
	{
		if(FRANDOM < 0.5)
		{
			spin[i] = 1;
		}
		else
		{
			spin[i] = -1;
		}

		neigh[i] = 0;
	}

	ET = 0.0;
	M = 0.0;

	for(i=0; i<L2; i++)
	{
		neigh[i] += spin[(i-L+L2)%L2];
		neigh[i] += spin[(i+1)%L + (i/L)*L];
		neigh[i] += spin[(i+L)%L2];
		neigh[i] += spin[(i-1+L)%L + (i/L)*L];

		ET += spin[i]*neigh[i];
		M += spin[i];
	}

	ET = (ET*(-J))/2.0;

	return;
}

/*****************************************************************************
 *                     	      MONTE CARLO ROUTINE                            *
 ****************************************************************************/
void mc_routine(double *boltz, int *spin, int *neigh)
{
	int i, t;
	
	for(t=0; t<L2; t++)
	{
		i = FRANDOM*L2;

		dE = 0;
		dE = 2*spin[i]*neigh[i];

		if(dE <= 0 || FRANDOM < boltz[dE])
		{
			spin[i] *= -1;
			ET = ET + dE;
			M = M + 2*spin[i];

			neigh[(i-L+L2)%L2] += 2*spin[i];
			neigh[(i+1)%L + (i/L)*L] += 2*spin[i];
			neigh[(i+L)%L2] += 2*spin[i];
			neigh[(i-1+L)%L + (i/L)*L] += 2*spin[i];
		}
	}

	return;
}
