# Appication of the Equilibrium Expectation algorithm to the Ising model

Compile by any C compiler, e.g.
gcc -o IsingEE  main.c Metropolis.c ContrDiv.c EquilibriumExpectation.c utils.c FinDif.c

Execute by command ./IsingEE Lugano.txt

Here Lugano.txt is a binary image in text format (genarated with ImageJ program)

For this data with 3006 x 4800 points Maximum Likelihood estimation took 2 minutes on MAC PRO

Output is written to files parameters.out, statistics.out

parameters.out contains values of parameters as a function of step t 

statistic.out contains values of statistics dz(t).

See the paper Byshkin et al Fast Maximum Likelihood estimation via Equilibrium Expectation for Large Network Data 2018 and the presentation slides at www.estimnet.org 
