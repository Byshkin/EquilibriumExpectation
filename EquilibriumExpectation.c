//A demonstration of the Equilibrium Expectation algorithm
//for Monte Carlo Maximum Likelihood estimation
//Copyright (C) <2018>  <Maksym Byshkin>
//This program is free software: you can redistribute it and/or modify it
//under the terms of the GPL-2 license <https://www.gnu.org/licenses/

/*****************************************************************************
 Reference for the Equilibrium Expectation algorithm
 Byshkin et al, Fast Maximum Likelihood estimation via Equilibrium Expectation for Large Network Data. 
 Scientific Reports, under review, 2018 , preprint https://arxiv.org/abs/1802.10311
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

void EE_algorithm(DataType **Data, int Lx, int Ly, double *theta, int Nparameters,
                   change_stats_func_t *change_stats_funcs[],
                   int m_steps,  double c2, int Mouter, int m2,
                   double D0[], double p2, double c1,FILE * theta_outfile, FILE *dzA_outfile,
		   FILE *Kafile, long double *dzA, double *Res, int t0)
{
  int touter, tinner, l;
  int pmn=0;
  double acceptance_rate;
  double theta_mean, theta_sd;
  long double *parmean=safe_malloc(Nparameters*sizeof(long double));
  double *theta_step = (double *)safe_malloc(Nparameters*sizeof(double));
  double **thetamatrix = (double **)safe_malloc(Nparameters*sizeof(double *));
  for (l = 0; l < Nparameters; l++)
    thetamatrix[l] = (double *)safe_malloc(m2*sizeof(double));
   for (l = 0; l < Nparameters; l++) parmean[l]=0;

  for (touter = 0; touter < Mouter; touter++) {
        fprintf(theta_outfile, "%u ", touter*m2+t0);
        fprintf(dzA_outfile, "%u ", touter*m2+t0);
        if(Kafile) fprintf(Kafile, "%u ", touter*m2+t0);
    for (tinner = 0; tinner < m2; tinner++) {
 	acceptance_rate= MetropolisSampler(Data, Lx,Ly,  theta, Nparameters,
          change_stats_funcs, dzA, m_steps, NULL, FALSE);    
      for (l = 0; l < Nparameters; l++) {
        theta_step[l] = (dzA[l] < 0 ? 1 : -1) * D0[l] * dzA[l]*dzA[l];

        if(c2!=0) theta[l] += theta_step[l];
        thetamatrix[l][tinner] = theta[l];
      }
    for (l = 0; l < Nparameters; l++) parmean[l]+=theta[l];
    }

   for (l = 0; l < Nparameters; l++)
    {
    fprintf(dzA_outfile, "%Lg ", dzA[l]); 
    fprintf(theta_outfile, "%g ", theta[l]);
    }
    fprintf(theta_outfile, "%g\n", acceptance_rate);
    fprintf(dzA_outfile, "\n");
  
    // get mean and sd of each theta value over inner loop iterations
    //   and adjust D0 to limit variance of theta 
    for (l = 0; l < Nparameters; l++) {
      theta_mean = mean_and_sd(thetamatrix[l], m2, &theta_sd);
      // force minimum magnitude to stop theta sticking at zero 
      if (fabs(theta_mean) < c1)  theta_mean = c1;
      if(Kafile) fprintf(Kafile, "%g ", theta_sd / fabs(theta_mean));
      D0[l] *= pow(c2 / (theta_sd / fabs(theta_mean)),p2);
    }
      if(Kafile) { fprintf(Kafile, "\n"); fflush(Kafile);  } 
      fflush(dzA_outfile);
      fflush(theta_outfile); 
  }
  for (l = 0; l < Nparameters; l++)   free(thetamatrix[l]);
  for (l = 0; l < Nparameters; l++) Res[l]=parmean[l]/(Mouter*m2);
  free(thetamatrix);
  free(parmean);
  free(theta_step);
}

