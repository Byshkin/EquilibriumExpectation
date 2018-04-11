#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "headers.h"
#include <string.h>
#include <limits.h>

//Pseudo random numbers
#define LLONG long long
#define ULLONG unsigned LLONG
static float Twoto_52  = 1.0/1048576/1048576/4096;
static float Pi2to_52  = 2*3.14159265358979323846/1048576/1048576/4096;

#define MULT_64_0      0x6765C793FA10079DL
#define MULT_64_1      0x147A6DA2B7F86475L

static  ULLONG rnd64_0=0x3A8DFA08E48B2827L;
static  ULLONG rnd64_1=0x85A36366EB71F043L;

void  SeedR64(int time) {
        rnd64_0 = time;
        rnd64_0 = rnd64_0<<32;
	rnd64_0 += ((ULLONG)clock()<<1) +1;
	rnd64_1 = rnd64_0;
        rnd64_0 *= MULT_64_0;
        rnd64_1 *= MULT_64_1;
}
static inline ULLONG ulR64Inline(void) {
    rnd64_0 *= MULT_64_0;
    rnd64_1 *= MULT_64_1;
    return  (rnd64_0 + rnd64_1);
}

double  dR64(void) {
    LLONG    U = ulR64Inline()>>12;
    return   U*Twoto_52; /* 1/2^52 */
}
int iR64 (void)  {
    unsigned  u = ulR64Inline()>>32;
    return   u >> 1;
}

// memory allocation
void *safe_malloc(size_t size)
{
  void *p = malloc(size);
  if (!p)  {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }
  return p;
}

// simple stats functions
/*
 * Compute mean and standard deviation of an array of doubles
 *
 * Parameters:
 *    values  - array of values to compute mean and sd of
 *    nvalues - length of values array
 *    sd      - (Out) standard deviation
 *
 * Return value:
 *   mean of the values
 *
 */
double mean_and_sd(double values[], int nvalues, double *sd)
{
  int  i;
  double  mean = 0;
  double  meandiff;

  for (i = 0; i < nvalues; i++)
    mean += values[i];
  mean /= nvalues;
  *sd = 0;
  for (i = 0; i < nvalues; i++) {
    meandiff = values[i] - mean;
    *sd += meandiff * meandiff;
  }
  *sd = sqrt(*sd / nvalues);
  return mean;
}



/* maximum size of an input line : adapt to your needs */
#define SIZE 100000

int readMatrix(char *infile, int Lmax, int *Lx, int *Ly, DataType ** Data ) {
    
    FILE *file=fopen(infile, "r");
    if(file) printf("Reading matrix from file %s\n", infile);
    if(!file) {printf("no file %s", infile); exit(1);}
    int numberOfNumbs=0,value, valsRead;
    float average;
    char  *val;
    char *line=safe_malloc(SIZE*sizeof(char));
    char delims[] = " \t\r\n";
   int i=0;  
   while(fgets(line, SIZE, file) != NULL)
   {
    val = strtok(line, delims);
    valsRead = sscanf(val, "%d",&value);
    int j=1;
    numberOfNumbs=1;
    Data[i][0]=value;
    while(valsRead>0)
    {
        numberOfNumbs++;
        val = strtok(NULL, delims);
        valsRead = (val == NULL) ? 0 : sscanf(val, "%d",&value);
        Data[i][j]=value;
        j++;
        if((i==Lmax)||(j==Lmax)) {printf("matrix too large, Increase Lmax and recompile\n"); exit(1);}
    }
    if(i==0) *Lx=numberOfNumbs-1;
    if((i!=0)&&(numberOfNumbs-1!=0)&&(*Lx!=numberOfNumbs-1)) {printf("%s file is not a matrix", infile);exit(1);}
    i++;
  }
    *Ly=i;
    fclose(file);
    free(line);
    return (0);
}
