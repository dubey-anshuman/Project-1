/*Finite time lyapunov exponent using "QR decompsition", but in this case not really QR because the state space is just one dimensional*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <mpi.h>

/*Global variables & constants*/
#define N_FRR 100		/*number of fourier coefficients considered for simulation*/
#define sigma 0.03
#define eta 0.01
#define max 500		/*number of total time steps*/
#define N 1000		/*total number of particles*/

double complex aa[2*N_FRR+1];		/*array to store fourier coefficients*/

FILE *output;				/* internal file name */

/*function prototypes*/
double complex gaussrand();
void randvec();
double f_derive(double xx);
double ff(double xx);

/*main*/

int main(int argc, char **argv) 
	{
		
		
		
		/*Variable declaration*/
		int i,m;
		time_t t;
		int my_id, root_process, ierr, num_procs;
		MPI_Status status;
		int N_each;
		/*initialization*/
		srand((unsigned) time(&t));
		ierr = MPI_Init(&argc, &argv);
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
		ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

		N_each = N/4;
		
		static double r[N],r_init[N], ftle[N];
		
		static double r_each[N/4],r_init_each[N/4], ftle_each[N/4];
		
		clock_t begin = clock();
		if (my_id == 0)
			{
				MPI_Scatter(r, N_each, MPI_DOUBLE, r_each, N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatter(r_init, N_each, MPI_DOUBLE, r_init_each, N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatter(ftle, N_each, MPI_DOUBLE, ftle_each, N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
		for (m = 0; m < N_each; m++)
			{
				r_each[m] = ((double)rand() / RAND_MAX);		/* clear array */
				r_init_each[m] = r_each[m];
				ftle_each[m] = 1.0;
			}
		MPI_Gather( r_each, N_each, MPI_DOUBLE, r, N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather( r_init_each, N_each, MPI_DOUBLE,r_init, N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather( ftle_each, N_each, MPI_DOUBLE, ftle, N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		/*program content*/
		if (my_id == 0)
		{
		for (i = 1; i < max; i++)
			{
				randvec();
				for (m = 0; m < N; m++)
					{
						r[m] += ff(r[m]);	/* numbers between */
						ftle[m] *= fabs(1.0+f_derive(r[m]));
					}
			}
			
		for (m = 0; m < N; m++)
			{
				if (ftle[m] < 1e-20) 
					{
						ftle[m] *= 1e22;
						ftle[m] = (log10(ftle[m])-22.0)/max;
					}
				else if (ftle[m] < 1e-10) 
					{
						ftle[m] *= 1e12;
						ftle[m] = (log10(ftle[m])-12.0)/max;
					}
				else if (ftle[m] > 1e10) 
					{
						ftle[m] *= 1e-12;
						ftle[m] = (log10(ftle[m])+12.0)/max;
					}
				else if (ftle[m] > 1e20) 
					{
						ftle[m] *= 1e-22;
						ftle[m] = (log10(ftle[m])+22.0)/max;
					}
				else
					{
						ftle[m] = log10(ftle[m])/max;
					}
			}
			clock_t end = clock();
			double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			printf("%e\n", time_spent);
		}
					
		
		ierr = MPI_Finalize(); /*check position of MPI_Finalize*/

		output= fopen ("ftle.dat", "w");	/* external file name */
		for (m = 0; m < N; m++)	fprintf(output,"%e\t%e\n", r_init[m],ftle[m]);
		fclose (output);
		
		return 0;
	}


/*function definitions*/
double complex gaussrand()
{
	double u1,u2,v1,v2,S;
	double complex X;
	do
	{
	u1 = (double) rand() / RAND_MAX; 
	u2 = (double) rand() / RAND_MAX;

	v1 = 2*u1-1;
	v2 = 2*u2-1;

	S = v1*v1 + v2*v2;
	} while (S >= 1 || S == 0);
	
	X = v1*sqrt(-2.0*log(S)/S) + (v2*sqrt(-2.0*log(S)/S))*I;
	X /= sqrt(2.0); //normalize X to have variance = 1
	return X;
}

void randvec()
{
	int i;
	for (i = N_FRR+1; i <= 2*N_FRR+1; i++)
	{
		aa[i] = gaussrand();
		aa[2*N_FRR-i] = conj(aa[i]);
	}
	aa[N_FRR] = gaussrand();
}

double f_derive(double xx)
{
	double complex f;
	double i_dbl, kk;
	int i;
	f=0+I*0;
	for (i = -N_FRR; i <=N_FRR; i++)
	{	
		i_dbl = (double) i;
		kk = 2*M_PI*i_dbl;
	f += aa[i+N_FRR]*kk*(I*cos(kk*xx) - sin(kk*xx))*exp(-eta*eta*M_PI*M_PI*i_dbl*i_dbl); /*derivative of the function f*/
	}
	f *= sigma*sqrt(sqrt(2*M_PI)*eta);
	return creal(f);
}

double ff(double xx)
{
	double complex f;
	double i_dbl, kk;
	int i;
	f=0+I*0;
	for (i = -N_FRR; i <=N_FRR; i++)
	{	
		i_dbl = (double) i;
		kk = 2*M_PI*i_dbl;
		f += aa[i+N_FRR]*(cos(kk*xx) + I*sin(kk*xx))*exp(-eta*eta*M_PI*M_PI*(i_dbl*i_dbl));
	}
	f *= sigma*sqrt(sqrt(2*M_PI)*eta);
	return creal(f);
}
