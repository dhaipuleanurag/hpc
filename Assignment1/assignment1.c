/* multiply vector components, write into a vector,
 *  and compute the inner product  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
void copy_array(double *source, double *destination, int len){
  int i;
  for(i = 0; i < len; i++) {
    destination[i] = source[i];
  }
}

double compute_residue(double *u, int len){
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

void printui(double *u, int len) {
  int i = 0;
  printf("uvalue:");
  for(i = 1; i < len; i++) {
     printf("%.10e ", u[i]);
  }
  printf("\n");
}

int main (int argc, char **argv)
{
  int  i, n;
  double *u,*prev_u,*f;
  int iterations = 1000;
  int pass;

  //int **mat = (int **)malloc(rows * sizeof(int*));
  //for(i = 0; i < rows; i++) mat[i] = (int *)malloc(cols * sizeof(int));


  if (argc != 2) {
    fprintf(stderr, "Function needs vector size and number of passes as input arguments!\n");
    abort();
  }
  n = atol(argv[1]);
  double h = 1.0/(n+1);
 // printf("%d", n);
 // printf("%f", h);  
  prev_u = (double *) malloc(sizeof(double) * (n+1));
  u = (double *) malloc(sizeof(double) * (n+1)); 
  f = (double *) malloc(sizeof(double) * (n+1));
  /*Check*/
  /* fill vectors */
  for (i = 0; i < n+1; ++i) {
    u[i] = 0;
    f[i] = 1;
    prev_u[i] = 0;
  }

  timestamp_type time1, time2;
  get_timestamp(&time1);

  for (pass = 0; pass < iterations; ++pass) {
  // printui(u, n+1);
  // printui(prev_u, n+1);
   double residue = compute_residue2(u, n+1, h);


  // printf("iteration %d :residue%f \n", pass+1, residue);
   printf("%f\n", residue);
  // a = (double *) malloc(sizeof(double) * n);
 
   for (i = 1; i < n+1; ++i) {
      if(i == 1){
        u[i] =h * h*0.5*(f[i] +(1/h) * (1/h) * prev_u[i+1]);	
      }
      else if(i == n){
	u[i] = h * h *0.5*(f[i] + (1/h)* (1/h)* prev_u[i-1]);	
      }
      else {
	u[i] = h * h *0.5*(f[i] + (1/h)* (1/h)* prev_u[i-1] + (1/h)* (1/h )* prev_u[i+1]);	
      }
     // printf("%f ", u[i]);
    }
  // printf("\n"); 
  // printf("u cur value %f",u[2]);
  // printf("\n");

  // printf("cur value %f",prev_u[2]);
    copy_array(u, prev_u, n+1);
  // double residue = compute_residue2(u, n+1);
  // printf("residue %f \n", residue);

//  printf("\n");

  // printf("next value %f", prev_u[2]);
 //  printf("\n");
  }

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Time elapsed is %f seconds.\n", elapsed);
  free(u);
  free(f);
  free(prev_u);
  return 0;
}
