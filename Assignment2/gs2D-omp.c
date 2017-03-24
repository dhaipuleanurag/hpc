#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

void copy_array(double **source, double **destination, int n){
  int i;
  int j;
 // #pragma omp parallel for
  for(i = 0; i < (n+2); i++) {
    #pragma omp parallel for
    for(j = 0; j < (n+2); j++) {
		destination[i][j] = source[i][j];
    }
  }
}

/*double compute_residue(double *u, int len){
  double sum = 0;
  int i = 0;
  for(i = 0; i < len; i++){
    sum = sum + u[i]*u[i];  
  }
  double sq = fabs(sqrt(sum));
  return sq;
}

double compute_residue2(double *u, int len, double h) {
  double sum = 0;
  int i = 0;
  for(i = 1; i < len; i++) {
     if(i == len -1) {
      double val =  fabs(u[i] * 2 *(1/h) *(1/h) - u[i-1] * (1/h) * (1/h)- 1);
      sum = sum + val * val;

     } else {
     
      double val = fabs(u[i] * 2 *(1/h) * (1/h)-  u[i-1]*(1/h)*(1/h) - u[i+1] *(1/h)*(1/h)- 1);
      sum = sum + val * val;
     }   
   }
   double sq = fabs(sqrt(sum));
   return sq;
}
*/

void printu(double **u, int n) {
  int i = 0;
  int j = 0;
  printf("uvalue:\n");
  for(i = 0; i < (n+2); i++) {
	//#pragma omp parallel for 
	for(j = 0; j < (n+2); j++){
		//printf("%.10e ", u[i][j]);
		printf("%f ",u[i][j]);
	}
	printf("\n");
  }
  printf("\n");
}

int main (int argc, char **argv)
{
  int  i, j, n;
  //double *u,*prev_u,*f;
  int iterations = 1000;
  int pass;
  
  if (argc != 2) {
    fprintf(stderr, "Function needs vector size and number of passes as input arguments!\n");
    abort();
  }
  
  n = atol(argv[1]);
  double h = 1.0/(n+1);
  double invh = 1.0/h;
  
  
  double **u = (double **)malloc((n+2) * sizeof(double *));
  for (i=0; i < (n+2); i++)
	u[i] = (double *)malloc((n+2) * sizeof(double));
  
  double **f = (double **)malloc((n+2) * sizeof(double *));
  for (i=0; i < (n+2); i++)
	f[i] = (double *)malloc((n+2) * sizeof(double));
  
  double **prev_u = (double **)malloc((n+2) * sizeof(double *));
  for (i=0; i < (n+2); i++)
	prev_u[i] = (double *)malloc((n+2) * sizeof(double));
  
  for (i = 0; i < (n+2); i++) {
	for(j = 0; j <  (n+2); j++) {
		u[i][j] = 0;
		f[i][j] = 1;
		prev_u[i][j] = 0;
	}
  }

  timestamp_type time1, time2;
  get_timestamp(&time1);

  
  for (pass = 0; pass < iterations; pass++) {
		  //printu(u,n);
	//#pragma omp parallel for
	for (i = 1; i < n+1; i++) {
		
		if(i % 2 == 1) {
			#pragma omp parallel for
			for(j = 1; j < n+1; j= j+2) {
				u[i][j] = 0.25 * (h * h + prev_u[i-1][j] + prev_u[i][j-1] + prev_u[i+1][j] + prev_u[i][j+1]);
	  		
			}
			#pragma omp parallel for
			for(j = 2; j < n+1; j= j+2) {
                        	u[i][j] = 0.25 * (h * h + u[i-1][j] + u[i][j-1] + u[i+1][j] + u[i][j+1]);
			
              		}
		
		} else {
                        #pragma omp parallel for
                        for(j = 2; j < n+1; j= j+2) {
                                u[i][j] = 0.25 * (h * h + prev_u[i-1][j] + prev_u[i][j-1] + prev_u[i+1][j] + prev_u[i][j+1]);

                        }
                        #pragma omp parallel for
                        for(j = 1; j < n+1; j= j+2) {
                                u[i][j] = 0.25 * (h * h + u[i-1][j] + u[i][j-1] + u[i+1][j] + u[i][j+1]);

                        }

                } 
	}
    	copy_array(u, prev_u, n);
  }
  

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printu(u, n);
  printf("Time elapsed is %f seconds.\n", elapsed);
  
  for(i = 0; i < n+2; i++) {
	free(u[i]);
	free(f[i]);
	free(prev_u[i]);	
  }
  free(u);
  free(f);
  free(prev_u);
  return 0;
}
