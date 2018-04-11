//A demonstration of the Equilibrium Expectation algorithm
//for Monte Carlo Maximum Likelihood estimation
//Copyright (C) <2018>  <Maksym Byshkin>
//This program is free software: you can redistribute it and/or modify it
//under the terms of the GPL-2 license <https://www.gnu.org/licenses/>

//The Equilibrium Expectation algorithm is fully described in Byshkin et al,
//Fast Maximum Likelihood estimation via Equilibrium Expectation for Large Network Data. 
//Scientific Reports, under review, 2018, preprint https://arxiv.org/abs/1802.10311

#include "headers.h"
#include "stdlib.h"
#include <stdio.h>
#include <string.h>
//Maximum size of data matrix
const int Lmax=15000;
//Number of model parameters to estimate
const int Npar=1;
//Constants of Contrastive Divergance algorithm
const double K_CD= 0.01;
const int CDsteps=5000;
const int m_steps=1000;
//Constants of EE algorithm
const int M_steps=4000;
const int m2_steps=100;
const double K0=1e-5;
const double c1=0.01;
const double c2=0.001;
const double p2=0.5;
int main(int argc, char *argv[])
{
 printf("A demonstration of the Equilibrium Expectation algorithm\n");
 printf("by estimation of parameters of the Ising model\n");
 init_prng(1);
 int i,j,k, h;
 int Lx, Ly;
 double parameters[Npar];
 double SumChangeStats[Npar];
 double Ki[Npar]; 
 change_stats_func_t* change_stats_funcs[100];
 //Specify model by statitics
 change_stats_funcs[1]=IsingChangeField;
 change_stats_funcs[0]=IsingChangeInteraction;

 FILE * stat_file=fopen("statistics.out","w");
 char hat[1000];
 sprintf(hat, "t  ");
 for (i = 0; i < Npar; i++)
     snprintf(hat+strlen(hat), 1000,"par%i  ", i+1 ); 
 fprintf(stat_file, "%s\n", hat);
 DataType **Data=safe_malloc(Lmax*sizeof(DataType*));
 for(i=0;i<Lmax;i++) Data[i]=safe_malloc(Lmax*sizeof(DataType));

//Load empirial data to study
 Lx=Ly=1000;
 readMatrix(argv[1], Lmax, &Ly, &Lx, Data);
 printf("Image size %i x %i pixels\n", Lx, Ly);
 fflush(stdout);
 for(i=0;i<Lx;i++) for(j=0;j<Ly;j++) 
    {h=Data[i][j]; Data[i][j]=1; if(h) Data[i][j]=-1; }

//Metropolis algorithm to generate simulated data
/*
printf("generate simulated data by Metropolis algorithm\n");
for(i=0;i<Npar;i++) SumChangeStats[i]=0;
parameters[0]=0.461241*1.01;
for(i=0;i<100000;i++) {
 double a= MetropolisSampler(Data, Lx,Ly,  parameters, Npar,
      change_stats_funcs, SumChangeStats, m_steps*10, NULL, FALSE);
  fprintf(stat_file,"%i	%g \n",i,SumChangeStats[0]); 
 }
*/

//Open output files
 FILE * par_file=fopen("parameters.out","w");

//mean_sd.out may be used to check if the approoximate equality (A1) holds
 FILE * mean_sd_file=NULL;
// mean_sd_file=fopen("mean_sd.out","w");
 fprintf(par_file, "%s acc\n", hat);
if(mean_sd_file) fprintf(mean_sd_file, "%s\n", hat);

 printf("Contrastive divergence estimates without averaging:\n");
 fflush(stdout);
 algorithm_S(Data, Lx, Ly, parameters, Npar,  change_stats_funcs,
                 CDsteps, m_steps,  K_CD, Ki, par_file);
 for(i=0;i<Npar;i++) printf("%g  ", parameters[i]);
 printf("\n");


//Set the starting values of Ki
// for (i = 0; i < Npar; i++) Ki[i]=1e-10;
 //A more robust approach is to compute deriviatives, as below
double *Dervi = (double *)safe_malloc(Npar*sizeof(double));
 Deriviative( 10, 1.0, Dervi, Data, Lx, Ly, Npar, 
		change_stats_funcs, parameters, m_steps);
for (i = 0; i < Npar; i++)
  { //printf("%g	", Ki[i]) ;
   Ki[i]=K0/(Dervi[i]*Dervi[i]);
   //printf("%g   \n", Ki[i]);
  }
  free(Dervi);

printf("Computing MLE by EE algorithm. See progress in output files\n");
double Res[Npar];
//Here statistics are difference between actual statistics and statistics in empirical data 
for(i=0;i<Npar;i++) SumChangeStats[i]=0;
EE_algorithm(Data, Lx, Ly, parameters, Npar, change_stats_funcs,
                   m_steps,   c2, M_steps, m2_steps, Ki, p2,  c1,
		   par_file,  stat_file, mean_sd_file, SumChangeStats, Res,  0);
//Average parameters value after burn-in. Statistics should fluctuate around zero
EE_algorithm(Data, Lx, Ly, parameters, Npar,change_stats_funcs,
                   m_steps,   c2, M_steps, m2_steps, Ki, p2,  c1,
		   par_file,  stat_file, mean_sd_file, SumChangeStats, Res, M_steps*m2_steps);

printf("Results of Maximum Likelihood estimation:\n");
for(i=0;i<Npar;i++) printf("%g  ", Res[i]);
printf("\n");
fflush(stdout);
printf("The convergence test was suggested in Byshkin et al,\n"); 
printf("Fast Maximum Likelihood estimation via Equilibrium Expectation for Large Network Data \n");
printf("To check if the estimated parameters are inside the confidence interval\n");
printf("It is possible to compute the expected statistics at limits of this interval\n");
printf("The statistics of observed data should be in-between %\n");

//MLE Test.
//If c2 is zero then paparemter values are not adjusted
//In this case EE algorithm is equivalent to Metropolis algorithm
//The eroor of the obtained estimates is smaller than x % if
//at parameter*(1+x/100) and parameters*(1-x/100) statistics have different sign

/*
double x=1;
printf("A test of accuracy is running\n");
double dp=1+x*0.01;
for(i=0;i<Npar;i++) parameters[i]=Res[i]/dp;
EE_algorithm(Data, Lx, Ly, parameters, Npar, change_stats_funcs,
                   m_steps,  0, M_steps, m2_steps, Ki, p2,  c1,
		   par_file,  stat_file, mean_sd_file, SumChangeStats, Res,2* M_steps*m2_steps);

for(i=0;i<Npar;i++) parameters[i]*=dp*dp;
EE_algorithm(Data, Lx, Ly, parameters, Npar, change_stats_funcs,
                   m_steps,  0, M_steps, m2_steps, Ki, p2,  c1,
		   par_file,  stat_file, mean_sd_file, SumChangeStats, Res,3* M_steps*m2_steps);

*/
  for(i=0;i<Lx;i++) free(Data[i]);
  free(Data);
  fclose(stat_file);
  fclose(par_file);
  if(mean_sd_file) fclose(mean_sd_file);
  return(0);
}
