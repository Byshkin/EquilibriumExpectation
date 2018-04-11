#include <math.h>
#include "headers.h"
#include "stdlib.h"
#include "stdio.h"

 
void Deriviative( int DS, double eps, double *Dervi, DataType **Data, int Lx, int Ly,
		  int  Nparameters, change_stats_func_t *change_stats_funcs[], double *parameters,
		  int sampler_m)
{
int i,l,ind;
double dx;
double *SumChangeStats = (double *)safe_malloc(Nparameters*sizeof(double));
//printf("\nForward deriviatives:\n");
for(ind=0;ind<Nparameters;ind++)
{
	parameters[ind]+=eps;
        dx=0;
	for(i=0;i<DS;i++)
	{
	int acceptance_rate =MetropolisSampler(Data, Lx, Ly, parameters, Nparameters,
	change_stats_funcs, SumChangeStats, sampler_m, NULL, TRUE);
	dx+=SumChangeStats[ind];
	}
	parameters[ind]-=eps;
        double der=(dx)/eps/DS;
 //       printf("%g  ",der);
        Dervi[ind]=der;    
}
//printf("\nBackward deriviatives:\n");
for(ind=0;ind<Nparameters;ind++)
{
	parameters[ind]-=eps;
        dx=0;
	for(i=0;i<DS;i++)
	{
	int acceptance_rate =MetropolisSampler(Data, Lx, Ly, parameters, Nparameters,
	change_stats_funcs, SumChangeStats, sampler_m, NULL, FALSE);
	dx+=SumChangeStats[ind];
	}	
	parameters[ind]+=eps;
        double der=-(dx)/eps/DS;
//        printf("%g  ",der);
	Dervi[ind]+=der;
        Dervi[ind]*=0.5;
}
printf("\n");
free(SumChangeStats);
}


