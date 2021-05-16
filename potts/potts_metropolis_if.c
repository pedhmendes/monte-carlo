/*****************************************************************************
 *                              Potts Model 2D                               *
 *                              Pedro H Mendes                               *
 *                                Metropolis				     *
 ****************************************************************************/

/*****************************************************************************
 *                            	   INCLUDES 			             *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mc.h"

/*****************************************************************************
 *                            	  DEFINITIONS                      	     *
 ****************************************************************************/
#define 			L           			64
#define 			L2          			(L*L)
#define 			J           			1.
#define 			KB			        1.
#define 			TRAN      			1//100000
#define 			TMAX     			1000000

/*****************************************************************************
 *                            GLOBAL VARIABLES                   	     *
 ****************************************************************************/
int M, ET;

/*****************************************************************************
 *                                FUNCTIONS 		                     *
 ****************************************************************************/
void initialize(int *spin, int **neigh, int _Q);
void sweep(int *spin, int **neigh, int _Q, double _TEMP);
void states(int *spin, int **neigh, int _Q);
void visualize(int *spin);

/*****************************************************************************
 *                         	 MAIN PROGRAM  				     *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	int *spin, **neigh;
	int i, mcs;
	int Q;	
	double TEMP;

	Q = atoi(argv[1]);
	TEMP = atof(argv[2]);

	spin = (int*)malloc(L2*sizeof(int));
	neigh = (int**)malloc(L2*sizeof(int*));
	
	for(i=0; i<L2; i++)
	{
		neigh[i] = (int *)malloc(4*sizeof(int));
	}

	seed = start_randomic();

	initialize(spin, neigh, Q);
	states(spin, neigh, Q);

	for(mcs=0; mcs<TRAN; mcs++)
	{
      		sweep(spin, neigh, Q, TEMP);	
    	}

	char Arq[100];
	FILE *arq;
	sprintf(Arq, "mp_Q%dL%dT%.3lf.dsf", Q, L, TEMP);
	arq = fopen(Arq, "w");
	fprintf(arq, "#SEED: %ld\n#MCS\tM\tET\n", seed);

    	
	for(mcs=0; mcs<TMAX; mcs++) 
	{
	  	sweep(spin, neigh, Q, TEMP);
	  	states(spin, neigh, Q);
#ifdef VIEW	 
	       	if(mcs % 10 == 0)
		{
			printf("set title 'T = %d MCS'\n", mcs);
			visualize(spin);
		}
#endif
		fprintf(arq, "%d\t%d\t%d\n", mcs, M, ET);

    	}

	fclose(arq);

	return 0;
}

/*****************************************************************************
 *                        	INITIALIZATION  			     *
 ****************************************************************************/
void initialize(int *spin, int **neigh, int _Q)
{
	int i;

     	for(i=0; i<L2; i++)
       	{
	   	spin[i]=FRANDOM*_Q;
     	}
	
     	for(i=0; i<L2; i++)
       	{
		neigh[i][0] = (i-L+L2)%L2;
                neigh[i][1] = (i+1)%L + (i/L)*L;
                neigh[i][2] = (i+L)%L2;
                neigh[i][3] = (i-1+L)%L + (i/L)*L;
	}
   
	return;
 }
/*****************************************************************************
 *                            MONTE CARLO ROUTINE                            *
 ****************************************************************************/
void sweep(int *spin, int **neigh, int _Q, double _TEMP)
{
	for (int i=0; i<L2; i++) 
	{
    		int site = FRANDOM*L2;
		int E1=0, E2 = 0;
  		int new_state = FRANDOM*_Q;
  		
		for(int j=0; j<4; j++)
		{
  		 	if(spin[site]==spin[neigh[site][j]])
			{
			  	E1-=J;
		    	}
      		}
	  	
		for (int j=0; j<4; j++)
		{
		    	if(new_state==spin[neigh[site][j]])
			{
			  	E2-=J;
		    	}
      		}
		
		double dE = E2 - E1;
		
		if(dE <= 0)
		{
		    	spin[site] = new_state;
      		}
		
		else
       		{
		    	double r = FRANDOM;
			double ALPHA = exp(-dE/(KB*_TEMP));  
		
			if(r<=ALPHA) 
			{     
			       	spin[site] = new_state;
		    	}
  		}
  	}

      	return;
}
/*****************************************************************************
 *                            	   STATES     		                     *
 ****************************************************************************/
void states(int *spin, int **neigh, int _Q)
{
	int i;

	ET=0;
  	M=0;

  	if(_Q!=2)
	{
		for (i=0; i<L2; i++)
		{
		  	if(spin[i]==spin[neigh[i][0]])
			{
				ET-=J;
			}
		  	if(spin[i]==spin[neigh[i][1]])
			{
				ET-=J;
			}
		  	if(spin[i]==spin[neigh[i][2]])
			{	
				ET-=J;
			}
		  	if(spin[i]==spin[neigh[i][3]])
			{
				ET-=J;
			}
		  
			M += spin[i];
	    	}
  	}

  	else
	{
	    	for(i=0; i<L2; i++)
		{
		  	if(spin[i]==spin[neigh[i][0]])
			{
				ET-=J;
			}
		  	if(spin[i]==spin[neigh[i][1]])
			{
				ET-=J;
			}
		  	if(spin[i]==spin[neigh[i][2]])
			{	
				ET-=J;
			}
		  	if(spin[i]==spin[neigh[i][3]])
			{
				ET-=J;
			}		
		  	
			if(spin[i]==0)
			{
				M -= 1;  
			}
		  	if(spin[i]==1)
			{
				M += 1;  
			}
    		}
  	}

	return;
}
/*****************************************************************************
 *                           GNUPLOT VISUALIZATION                           *
 ****************************************************************************/
void visualize(int *spin)
{
  	int i;
	printf("unset xtics\nunset ytics\n");	
	printf("set size square\n");
	printf("set xrange [0:%d]\nset yrange [0:%d]\n",L-1,L-1);
    	printf("pl '-' matrix w image\n");
	for(i=L2-1; i>=0; i--) 
	{
	    	printf("%d ", spin[i]);
		if(i%L == 0)
		{
		       	printf("\n");
		}
      	}
	
      	printf("e\n\n");
}

