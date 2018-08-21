# Appication of the EE algorithm to the Ising model

Compile by any C compiler, e.g.
gcc -o IsingEE  main.c Metropolis.c ContrDiv.c EquilibriumExpectation.c utils.c FinDif.c

Execute by command ./IsingEE Lugano.txt

Here Lugano.txt is a binary image in text format (genarated with ImageJ program)

For this data with 3006 x 4800 points Maximum Likelihood estimation took 2 minutes on MAC PRO

First 1-step Contrastive Divergance is computed. Then MLE.

Output is written to files parameters.out, statistics.out

parameters.out contains values of parameters as a function of step t 

statistic.out contains values of statistics dz(t).

See the paper Byshkin et al (2018) [Fast Maximum Likelihood Estimation via Equilibrium Expectation for Large Network Data](https://www.nature.com/articles/s41598-018-29725-8). *Scientific Reports* 8:11509. https://doi.org/10.1038/s41598-018-29725-8 and the presentation slides at www.estimnet.org 
