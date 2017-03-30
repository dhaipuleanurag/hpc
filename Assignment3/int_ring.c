/* Illustrate exchange of array of doubles
 * ping-pong style between even and odd processors
 * every processor writes its result to a file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "util.h"
int main( int argc, char *argv[])
{
  int rank, i;
  MPI_Status status;
  int total_threads;
  int N = atoi(argv[1]);
  int a_size = 512000;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total_threads);
  int *message_out = calloc(a_size, sizeof(int));
  int *message_in = calloc(a_size, sizeof(int));
  for(i = 0; i < a_size; ++i)
  {
    message_out[i] = 0;
    message_in[i] = 0;
  }

  int tag = 99;
  int origin, destination;

  timestamp_type time1, time2;

  destination = rank + 1;
  origin = rank - 1;
  if(rank == 0)
	origin = total_threads - 1;
  if(rank == total_threads - 1)
        destination = 0; 

  if(rank == 0)
  {  
    MPI_Send(message_out, a_size, MPI_INT, destination, tag, MPI_COMM_WORLD);
    printf("Message sent by %d: %d\n", rank, message_out[0]);  
    get_timestamp(&time1);
  }
  for(i = 0; i < N; i++)
  {
      if(rank == 0 && i == N-1) {
        break;
      }
      MPI_Recv(message_in,  a_size, MPI_INT, origin,      tag, MPI_COMM_WORLD, &status);
      if(rank == 0)
      {
         get_timestamp(&time2);
         double elapsed = timestamp_diff_in_seconds(time1,time2);
         printf("Time elapsed is %f seconds.\n", elapsed);
      }
      printf("Message received at %d: %d\n", rank, message_in[0]);
      message_out[0] = message_in[0] + rank; 
      message_out[1] = message_in[1] + rank + rank; 
      MPI_Send(message_out, a_size, MPI_INT, destination, tag, MPI_COMM_WORLD);
      printf("Message sent by %d: %d\n", rank, message_out[0]);
      if(rank == 0)
      {
         get_timestamp(&time1);
      }
  }
  if(rank == 0)
  {
   MPI_Recv(message_in, a_size, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
   printf("Final Message received at %d index 0: %d\n", rank, message_in[0]);
   printf("Final Message received at %d index 1: %d\n", rank, message_in[1]);
  }
  MPI_Finalize();
  return 0;
}
