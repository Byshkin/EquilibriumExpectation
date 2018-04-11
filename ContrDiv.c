//This software is free, open source and is distributed under the GPL-2 license
//Maksym Byshkin, 2018

#include "headers.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
// 1-step Contrastive Divergence (CD)
/*
References for the CD algorithm
Fellows, I. E. "Why (and When and How) Contrastive Divergence Works." arXiv preprint arXiv:1405.0602 (2014).  
 (Found via Krivitsky, P. N.  Using contrastive divergence to seed Monte Carlo MLE for exponential-family
 random graph models. Computational Statistics & Data Analysis, 107, 149-161 (2017)*/

void algorithm_S(DataType **Data, int Lx, int Ly, double *parameters, int Nparameters,
                 change_stats_func_t *change_stats_funcs[],
                 int M1, int m_steps,  double K0,
                 double Dmean[], FILE * theta_outfile)
{
  int t, l;
  double *Vari = (double *)safe_malloc(Nparameters*sizeof(double));
  double *D0 = (double *)safe_malloc(Nparameters* sizeof(double));
  double *dzA = (double *)safe_malloc(Nparameters*sizeof(double));
  double parameters_change,da, acceptance_rate; 
  /* 1/D0 approximates squared derivatives */  
  for (l = 0; l < Nparameters; l++)    parameters[l] = 0;
 
  for (t = 0; t < M1; t++) {
    fprintf(theta_outfile, "%d ", t-M1);
    acceptance_rate= MetropolisSampler(Data, Lx,Ly,  parameters, Nparameters,
        change_stats_funcs, dzA, m_steps, Vari, TRUE);
    for (l = 0; l < Nparameters; l++) {
      // next 3 lines are for S algorithm
      D0[l] += dzA[l]*dzA[l];
      da = 0;
      if (dzA[l] < 0 || dzA[l] > 0)  da = K0 / (dzA[l]*dzA[l]);

      //I.Fellows algorithm
      double iv=0;
      if(fabs(Vari[l])>1e-10) iv=K0/Vari[l]*0.001;
	  
      //parameters_change = (dzA[l] > 0 ? -1 : 1) * da * dzA[l]*dzA[l];//S algorithm
      parameters_change=-dzA[l]*iv;// I.Fellows algorithm
      parameters[l] += parameters_change;
      fprintf(theta_outfile, "%g ", parameters[l]);
    }
    fprintf(theta_outfile, "%g\n", acceptance_rate);
  }
  
  for (l = 0; l < Nparameters; l++)    Dmean[l] = m_steps / D0[l];
  free(D0);
  free(dzA);
  free(Vari);
}
