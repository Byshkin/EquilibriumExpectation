//This software is free, open source and is distributed under the GPL-2 license
//Maksym Byshkin, April 2018
#include <math.h>
#include "headers.h"
#include "stdlib.h"
#include "stdio.h"
//Ising model on squere Lattice
//compute change of statistics when value of one variable is increased from -1 to 1 
double IsingChangeField(DataType **Data, int x, int y, int Lx, int Ly)
{
//Ising statistics t1=Sum(Data[i])
//if variable value is increased from -1 to 1 then t1 statistics is increased by 2
//printf("LX=%i",Lx);
return 2.0;
}

double IsingChangeInteraction(DataType **Data, int x, int y, int Lx, int Ly)
{
//t2=Sum(Data[i]*Data[neignors[i]])
int var=Data[x][y];
if(var!=-1) printf("Data[%i][%i] %i\n", x,y, var);
//on squere lattice coordinate (x,y) has 4 neigbors:
// (x+1,y), (x-1,y), (x,y-1), (x,y+1)
//unless (x,y) in on the boundary
double res=0;
if(x-1>=0) res+=2*Data[x-1][y];
if(y-1>=0) res+=2*Data[x][y-1];
if(x+1<Lx) res+=2*Data[x+1][y];
if(y+1<Ly) res+=2*Data[x][y+1];
return res;
}


double MetropolisSampler(DataType ** Data, int Lx, int  Ly,  double * parameters,int Nparameters,
    change_stats_func_t *change_stats_funcs[], long double * SumChangeStats, int sampler_m, double *Vari, boolean CD)
 {

  int Naccepted = 0;   
  double energy_change;
  int k,l;
  boolean   ONE,OLD;
  double *changestats = (double *)safe_malloc(Nparameters*sizeof(double));
  double *changestatsM;//changes of statistics to compute variance
  if (Vari) changestatsM = (double *) safe_malloc(Nparameters*sampler_m*sizeof(double));
  int x,y;
  if(CD) for (l = 0; l < Nparameters; l++)  SumChangeStats[l] = 0;

  for (k = 0; k < sampler_m; k++) {
    //select a random variable from the data  
    x = int_urand(Lx);
    y = int_urand(Ly);
    //Safe the varibale value
    OLD=Data[x][y];
    ONE = (Data[x][y]==1);
    //compute the  energy change when variable value is increased
    if (ONE)  Data[x][y]=-1;
    energy_change = 0;
    for (l = 0; l < Nparameters; l++) { 
      changestats[l] = (*change_stats_funcs[l])(Data, x, y, Lx, Ly);
      energy_change += parameters[l] * (ONE ? -1 : 1) * changestats[l];
    }
   //compute Metropolis-Hastings acceptance probability
   double AccProb= exp(energy_change);
   boolean ACC=(urand()<AccProb);
    if (ACC) { //move accepted
      Naccepted++;      
      if (!ONE)     Data[x][y]=1;
      for(l=0;l<Nparameters;l++)    SumChangeStats[l] +=(ONE ? -1 : 1)* changestats[l];
    } else { // move not accepted
      Data[x][y]=OLD;
    }
   if(CD) Data[x][y]=OLD;//restore the old state for CD
//   printf("val%i",Data[x][y] );
   if(Vari) {//safe changes of statitics to compute variance of statistics
	for(l=0;l<Nparameters;l++) changestatsM[k*Nparameters+l];//changestatsM[k*sampler_m+l]=0;
	if (ACC) for(l=0;l<Nparameters;l++) changestatsM[k*Nparameters+l]=(ONE ? -1 : 1)*changestats[l]; 
	}  
  }
if(Vari) //compute variance of statistics
{
	for(l=0;l<Nparameters;l++)
	{
	double mean=0;
	for (k = 0; k < sampler_m; k++) mean+=changestatsM[k*Nparameters+l];
	mean/=sampler_m;
	double Var=0;
	for (k = 0; k < sampler_m; k++) Var+=(changestatsM[k*Nparameters+l]-mean)*
					(changestatsM[k*Nparameters+l]-mean);
	Vari[l]=Var/sampler_m;
	}
free(changestatsM);
}
  free(changestats);
  return (double)Naccepted / sampler_m;
}

