#include <sys/time.h>
#include <stdio.h>
#define FALSE 0
#define TRUE 1
typedef int          bool;
typedef short int          DataType;  
typedef double (change_stats_func_t)(DataType **, int, int, int, int );
  
double IsingChangeField(DataType **Data, int x, int y, int Lx, int Ly);
double IsingChangeInteraction(DataType **Data, int x, int y, int LX, int Ly);
double MetropolisSampler(DataType ** Data, int Lx, int Ly,  double * parameters,int Nparameters,
   change_stats_func_t *change_stats_funcs[], long  double * SumChangeStats, int sampler_m, double *Vari, bool CD);

void algorithm_S(DataType **Data, int Lx, int Ly, double *parameters, int Nparameters,
                 change_stats_func_t *change_stats_funcs[],
                 int M1, int sampler_m,  double K0,
                 double Dmean[], FILE * theta_outfile);
void EE_algorithm(DataType **Data, int Lx, int Ly, double *theta, int Nparameters,
                   change_stats_func_t *change_stats_funcs[],
                   int m_steps,  double c2, int Mouer, int m2,
                   double D0[], double p2, double c1,FILE * theta_outfilei,
		    FILE *dzA_outfile, FILE* mean_sd_file, long double * SumChangeStats,
		    double *Res, int t0);

 void Deriviative( int DS, double eps, double *Dervi, DataType **Data, int Lx, int Ly,
                   int  Nparameters, change_stats_func_t *change_stats_funcs[],
		   double *parameters, int sampler_m);
void  SeedR64(int time) ;
double  dR64(void) ;
int iR64 (void)  ;


static inline void init_prng(int tasknum)
{
SeedR64(tasknum);
}

static inline double urand(void)
{
return dR64();
}

static inline int int_urand(int n)
{
  return iR64() % n; 
}

 
void *safe_malloc(size_t size);

double mean_and_sd(double values[], int nvalues, double *sd);

 int readMatrix(char *infile, int Lmax, int *Lx, int *Ly, DataType ** Data );
