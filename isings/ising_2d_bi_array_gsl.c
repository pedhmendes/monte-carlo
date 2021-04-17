////////////////////////
/////MODELO DE ISING////
//////PEDRO MENDES//////
////////////////////////

//REQUER BIBLIOTECA GSL ----> LER NA PAGINA

//gcc -Wall ising_2d_full.c -lgsl -lgslcblas -lm -static

//opções de compilacao
//-DSEED : cada execucao do programa usa uma seed diferente
//-DSHOT : faz screenshot inicial e final do sistema
//-DMEAN : calcula suscepitibilidade magnetica e calor especifico

/*******INCLUDES*******/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>

/*******DEFINICOES*******/
#define L 20
#define L2 (L*L)
#define transient 100000
#define tmax 1000000

/*******VARIAVEIS GLOBAIS*******/
int spin[L][L];
int neigh[L][L];
int dE, J=1, M, ET;
int h_M[2*L2], h_E[4*L2];
double T=1.0, boltz[9];

/*******FUNCOES*******/
void initialize(gsl_rng *mt);
void neighbour(int i, int j);
void update_neighbour(int i, int j);
void Boltz(void);
int E_flip(int i, int j);
void mc_routine(gsl_rng *mt);
void print_func(FILE *file);
void simulate(void);

/*******MAIN*******/
int main()
{
	simulate();

	return 0;
}
	
void simulate(void)
{
/*******VARIAVEIS LOCAIS*******/
	int mcs, k;
	unsigned long int seed;
	char Arq1[100], Arq2[100], Arq3[100];
	FILE *arq1, *arq2, *arq3;

/*******RNG MERSENNE TWISTER LIBGLS*******/
	srand(time(NULL));
#ifdef SEED
	seed = 0;
#else
	seed = rand();
#endif
	gsl_rng *mt;
	mt = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(mt, seed);

/*******INICIALIZA O SISTEMA*******/
	initialize(mt);

/*******SNAPSHOT INICIAL*******/
#ifdef SHOT
	char Arq4[100];
	FILE *arq4;

	sprintf(Arq4, "shot_i_T%.2lfL%d.dat", T, L);
	arq4 = fopen(Arq4, "w");
	print_func(arq4);
	fclose(arq4);
#endif

/*******TEMPO DE EQUILIBRIO*******/	
	for(mcs=0; mcs<transient; mcs++)
	{
		mc_routine(mt);
	}

/*******ZERA VARIAVEIS DE MEDIA*******/	
#ifdef MEAN
	double sum_M, sum_M2;
	double sum_E, sum_E2;
	double X, C;
	sum_M = 0; sum_M2 = 0;
	sum_E = 0; sum_E2 = 0;
#endif

/*******COMECA AS MEDIDAS*******/	
	sprintf(Arq1, "temp_T%.2lfL%d.dat", T, L);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#seed = %ld\n#MCS\tM\tET\n", seed);

	for(mcs=0; mcs<tmax; mcs++)
	{
		mc_routine(mt);
	
		fprintf(arq1, "%d\t%d\t%d\n", mcs, M, ET);

		h_M[M + L2] += 1;
		h_E[ET + 2*L2] += 1;

#ifdef MEAN
		sum_M += fabs(1.0*M);
		sum_E += 1.0*ET;
		sum_M2 += pow((1.0*M),2);
		sum_E2 += pow((1.0*ET),2);
#endif
	}

/*******SUSCETIBILIDADE MAGNETICA E CALOR ESPECIFICO*******/
#ifdef MEAN
	sum_M = sum_M/tmax;
	sum_E = sum_E/tmax;
	sum_M = pow(sum_M,2);
	sum_E = pow(sum_E,2);

	sum_M2 = sum_M2/tmax;
	sum_E2 = sum_E2/tmax;
	
	X = (sum_M2 - sum_M)/(L2*T);
	C = (sum_E2 - sum_E)/(L2*T*T);
	
	printf("#T\tX\tC\n%.2lf\t%.14lf\t%.14lf\n", T, X, C);
#endif

/*******IMPRIME HISTOGRAMAS*******/
	sprintf(Arq2, "hM_T%.2lfL%d.dat", T, L);
	sprintf(Arq3, "hE_T%.2lfL%d.dat", T, L);
	arq2 = fopen(Arq2, "w");
	arq3 = fopen(Arq3, "w");
	fprintf(arq2, "#seed = %ld\n#i\th_M[i]\n", seed);
	fprintf(arq3, "#seed = %ld\n#i\th_E[i]\n", seed);

	for(k=0; k<2*L2; k+=2)
	{
		fprintf(arq2, "%d\t%d\n", (k-L2), h_M[k]);
	}

	for(k=0; k<4*L2; k+=4)
	{
		fprintf(arq3, "%d\t%d\n", (k-2*L2), h_E[k]);
	}

/*******FECHA ARQUIVOS E LIMPA MEMORIA*******/	
	fclose(arq1);
	fclose(arq2);
	fclose(arq3);

	gsl_rng_free(mt);

/*******SNAPSHOT FINAL*******/
#ifdef SHOT
	char Arq5[100];
	FILE *arq5;

	sprintf(Arq5, "shot_f_T%.2lfL%d.dat", T, L);
	arq4 = fopen(Arq5, "w");
	print_func(arq5);
	fclose(arq5);
#endif
	
/*******TERMINA O PROGRAMA*******/
	return;
}

void initialize(gsl_rng *mt)
{
	int i, j;

/*******INICIALIZA PROBABILIDADES*******/
	Boltz();

/*******ESTADOS ALEATORIOS DOS SPINS*******/
	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			spin[i][j] = 2*(gsl_rng_uniform_int(mt,2)) -1;
		}
	}

/*******MATRIZ DE VIZINBOS*******/
	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			neighbour(i, j);
		}
	}

/*******ENERGIA INICIAL*******/
	ET = 0;

	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			ET = ET + spin[i][j]*(neigh[i][j]);
		}
	}

	ET = (ET*(-J))/2.0;

/*******MAGNETIZACAO INICIAL*******/
	M = 0;

	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			M = M + spin[i][j];
		}
	}

/*******ZERA OS HISTOGRAMAS*******/
	for(i=0; i<2*L2; i++)
	{
		h_M[i] = 0;
	}

	for(i=0; i<4*L2; i++)
	{
		h_E[i] = 0;
	}

	return;
}

void neighbour(int i, int j)
{
/*******ATUALIZA UMA POSICAO DA MATRIZ DOS VIZINHOS*******/
	neigh[i][j] = 0;
	neigh[i][j] += spin[(i == 0 ? (L-1) : i - 1)][j];
	neigh[i][j] += spin[(i == (L-1) ? 0 : i + 1)][j];
	neigh[i][j] += spin[i][(j == 0 ? (L-1) : j - 1)];
	neigh[i][j] += spin[i][(j == (L-1) ? 0 : j + 1)];

	return;
}

void update_neighbour(int i, int j)
{
/*******ATUALIZA OS VIZINHOS DO SITIO FLIPADO*******/
	neigh[(i == 0 ? (L-1) : i - 1)][j] += 2*spin[i][j];
	neigh[(i == (L-1) ? 0 : i + 1)][j] += 2*spin[i][j];
	neigh[i][(j == 0 ? (L-1) : j - 1)] += 2*spin[i][j];
	neigh[i][(j == (L-1) ? 0 : j + 1)] += 2*spin[i][j];

	return;
}

void Boltz(void)
{
/*******PROBABILIDADES*******/
	boltz[4] = exp(-4.0/T);
	boltz[8] = exp(-8.0/T);

	return;
}

int E_flip(int i, int j)
{
/*******VARIACAO DA ENERGIA*******/
	dE = 2.0*J*(spin[i][j]*neigh[i][j]);

	return dE;
}

void mc_routine(gsl_rng *mt)
{
	int i, j, t;
	double r, prob;

	for(t=0; t<L2; t++)
	{
/*******POSICAO ALEATORIA*******/
		i = gsl_rng_uniform_int(mt, L);
		j = gsl_rng_uniform_int(mt, L);
	
/*******VARIACAO DE ENERGIA*******/
		dE = 0;
		dE = E_flip(i, j);

/*******METROPOLIS*******/
		if(dE <= 0)
		{
			spin[i][j] = -spin[i][j];
			ET = ET + dE;
			M = M + 2*(spin[i][j]);

			update_neighbour(i, j);
		}
		else
		{
			prob = boltz[dE];
			r = gsl_rng_uniform(mt);
		
			if(r < prob)
			{
				spin[i][j] = -spin[i][j];
				ET = ET + dE;
				M = M + 2*(spin[i][j]);

				update_neighbour(i, j);
			}
		}
	
	}
}

void print_func(FILE *file)
{
/*******PRINTAR MATRIZES*******/
	int a, b;

	for(a=0; a<L; a++)
	{
		for(b=0; b<L; b++)
		{
			fprintf(file,"%d ", spin[a][b]);
		}
		fprintf(file,"\n");
	}

	return;
}
