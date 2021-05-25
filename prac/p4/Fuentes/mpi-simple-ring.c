/*
** Sending simple, point-to-point messages.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h" 

#define MASTER 0

int main(int argc, char* argv[])
{
  int myrank;
  int size;
  int dest;              /* destination rank for message */
  int source;            /* source rank of a message */
  int tag = 0;           /* scope for adding extra information to a message */
  MPI_Status status;     /* struct used by MPI_Recv */
  MPI_Request req;       /* struct used by MPI_Isend */
  char message[BUFSIZ];

  /* MPI_Init returns once it has started up processes */
  MPI_Init( &argc, &argv );

  /* size and rank will become ubiquitous */ 
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  /* 
  ** SPMD - conditionals based upon rank
  ** will also become ubiquitous
  */


  /* create a message to send, in this case a character array */
  sprintf(message, "Hola Mundo! Soy el proceso %d!", myrank);
  /* send to the master process */
  // dest = MASTER;
  /* send to the next process */
  dest = (myrank + 1) % size;
  /* 
  ** Send our first message!
  ** use strlen()+1, so that we include the string terminator, '\0'
  */
  MPI_Isend(message,strlen(message)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD, &req);

  if (myrank!=MASTER){
    /* Not master recieves from previous*/
    source = myrank - 1;
  } else {
    /* Master recieves from the last*/
    source = size-1;
  }
  /* recieving messages.. */
  MPI_Recv(message, BUFSIZ, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  /* and then printing them */
  printf("Mensaje recibido por el proceso %d: %s\n", myrank, message);

  /* don't forget to tidy up when we're done */
  MPI_Finalize();

  /* and exit the program */
  return EXIT_SUCCESS;
}