#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static void proceso(int,int);

int main(int argc, char* argv[]){
 int myRank, cantidadDeProcesos;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&cantidadDeProcesos);

 	proceso(myRank,cantidadDeProcesos);

	MPI_Finalize();
	return(0);
}

static void proceso(int myRank,int cantidadDeProcesos){
 MPI_Status status;
 int mensaje;
 int i;
 int p;
 char hostname[255];

	gethostname(hostname,255);
	printf("Proceso %i >> LISTO EN HOST %s\n",myRank,hostname);

	// Envia un mensaje a cada uno de los otros procesos
	for(p=0;p<cantidadDeProcesos;p++){
		if(p!=myRank){
			MPI_Send(&myRank,1,MPI_INT,p,0,MPI_COMM_WORLD);
			printf("Proceso %i >> ENVIO A %i \n",myRank,p);
		}
	}

	// Recibe el mensaje de cada uno de los otros procesos
	for(p=0;p<cantidadDeProcesos-1;p++){
		if(p!=myRank){
			MPI_Recv(&mensaje,1,MPI_INT,p,0,MPI_COMM_WORLD,&status);
			printf("Proceso %i >> RECIBIO DE %i \n",myRank,mensaje);
		}
	}

}
