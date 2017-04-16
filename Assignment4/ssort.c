/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

int main( int argc, char *argv[])
{
  int rank;
  int i, j, N, p;
  int *vec, *splits;
  MPI_Status status1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
  sscanf(argv[1], "%d", &N);
    
  if(N % p != 0 && N < 12){
    printf("Exiting. N must be a multiple of p and greater than 12\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  
  vec = calloc(N, sizeof(int));
  splits = calloc(p-1, sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  
  /* sort locally */
  qsort(vec, N, sizeof(int), compare);

  /* Get splitters locally */
  for(i = N/p, j = 0; i < N; i= i+N/p, j++) {
     splits[j] = vec[i];
  }
 
  /* Gather all splitters to process 0 */ 
  int *allsplits;
  if(rank == 0){
    //allsplits = calloc((p -1)*p, sizeof(int));
    allsplits = (int *)malloc(p *(p-1) * sizeof(int)); 
  }
  MPI_Gather(splits, p-1, MPI_INT, allsplits, (p-1), MPI_INT, 0, MPI_COMM_WORLD);
  

  /* Process 0 decides final splitters based on which numbers will be sent to different 
  * buckets or processes */
  int *finalsplits = calloc(p-1, sizeof(int));
  if(rank == 0){
    qsort(allsplits, p*(p-1), sizeof(int), compare);
    for(i = p-1,j = 0; i <= p*(p-1) ; i= i+p-1, j++){
      finalsplits[j] = allsplits[i];  
    }
  }
  MPI_Bcast(finalsplits, p-1, MPI_INT, 0, MPI_COMM_WORLD);
  

  
  int **sr = (int **)malloc(p * sizeof(int *));
  for (i=0; i < p; i++)
     sr[i] = (int *)malloc(p * sizeof(int));
  
  
  /* Every process decides how many to values it should expect from it. */
  int *myarray = calloc(2*N, sizeof(int));
  int myiter = 0;
  int count = 0;
  for(i = 0, j = 0; i < N && j < p-1; i++){
     if(vec[i] <= finalsplits[j]){
       count++;
       if(rank == j){
         myarray[myiter++] = vec[i];
       }
     } else {
       if(rank != j){      
         MPI_Send(&count, 1, MPI_INT, j, 100 , MPI_COMM_WORLD);
       }
       j++;
       //printf("\njvalue %d    %d\n",  finalsplits[j-1], j-1);
       count = 1;
       if(rank == j){
         myarray[myiter++] = vec[i];
       }
    }
  }
  count = N-i+1; 
  if(rank != j){
    MPI_Send(&count, 1, MPI_INT, j, 100 , MPI_COMM_WORLD);
  } else {
    for(i = i; i < N; i++){
      myarray[myiter++] = vec[i];
    }
  }
  /* Receive and store the number of numbers to expect from other processes */  
  for(j = 0; j < p; j++) {
    if(rank != j){
      MPI_Recv(&sr[rank][j], 1, MPI_INT, j, 100, MPI_COMM_WORLD, &status1); 
    }
  }  


  /* Start sending and receiving number to different buckets */
  for(i = 0, j = 0; i < N && j < p-1; i++){
    if(vec[i] <= finalsplits[j]) {
      if(rank != j){
        MPI_Send(&vec[i], 1, MPI_INT, j, 100 , MPI_COMM_WORLD);
        
      }
    } else {
      j++;
      if(rank != j){
        MPI_Send(&vec[i], 1, MPI_INT, j, 100 , MPI_COMM_WORLD);          
      } 
    }
  }
  for(i = i; i < N; i++){
   if(rank != j){  
      MPI_Send(&vec[i], 1, MPI_INT, j, 100 , MPI_COMM_WORLD);
    }
  }
  
 
  int totalreceived = 0;
  int receivecount = 0;
  for(i = 0; i < p; i++){
    if(i != rank){
      receivecount = sr[rank][i];
      for(j = 0; j < receivecount; j++){
        MPI_Recv(&myarray[myiter++], 1, MPI_INT, i, 100, MPI_COMM_WORLD, &status1);
        totalreceived++;     
      } 
    }
  }
  

  /* Sort the buckets locally */
  qsort(myarray, myiter, sizeof(int), compare);

  /* every processor writes its result to a file */
  char filename[20];
  sprintf(filename, "sorted-%d.txt", rank);
  FILE* fp = fopen(filename, "w");
  for (i=0; i<myiter; i++) {
     fprintf(fp, "%d ", myarray[i]);
  }  

  /*
  printf("\n Sorted values in process %d: \n", rank);
  for(i = 0; i < myiter; i++){
    printf(" %d", myarray[i]);
  }
  */

  free(vec);
  //free(splits);
  free(myarray);
  //free(sr);
  //free(allsplits);  
  MPI_Finalize();
  return 0;
}
