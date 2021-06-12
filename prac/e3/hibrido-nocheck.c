#include <stdio.h>
#include <stdlib.h>   /* malloc() */
#include <sys/time.h> /* gettimeofday */
#include <math.h>     /* sin, cos, fabs, ceil */
#include <time.h>     /* srand((unsigned) time(&t)) */
#include <omp.h>      /* hilos */
#include <mpi.h>

/* Para calcular el tiempo */
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

/* Generación de números aleatorios */
double randFP(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

#define PI 3.14159265358979323846

#define COORDINATOR 0

/************** MAIN *************/
int main(int argc, char *argv[])
{ 
  int numProcs, rank, n, bs, bs2, procRows, thRows, provided, numThreads, id, size, workerSize;
  MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Status status;

  	/* Lee parámetros de la línea de comando */
	if ((argc != 4) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((numThreads = atoi(argv[3])) <= 0)) {
    printf("\nUsar: %s size bs\n  size: Dimension de la matriz\n  bs: Dimension del bloque\n", argv[0]);
		exit(1);
	}
  if ((n % bs) != 0)
  {
    printf("\nError en los parámetros. Usage: ./%s size bs (size debe ser múltiplo de bs)\n", argv[0]);
    exit(1);
  }
	if (n % numProcs != 0) {
		printf("El tamaño de la matriz debe ser múltiplo del numero de procesos.\n");
		exit(1);
	}

  /* Para números aleatorios */
  time_t t;
  srand((unsigned) time(&t));

  /* Punteros */
  double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, *ablk, *bblk, *cblk, avgR1, avgR2, lavgR1, lavgR2, sinPhi, cosPhi;

  /* Índices */
  int i, j, k, ii, jj, kk;

  /* Para medir el tiempo */
  double timetick;
  int COMM_SIZE = 6;
  double commTimes[COMM_SIZE], maxCommTimes[COMM_SIZE], minCommTimes[COMM_SIZE], commTime, totalTime;

  /* Calcular porción de cada worker */
  procRows = n / numProcs;
  thRows = procRows / numThreads;
  /**/
  size = n * n;
  workerSize = procRows * n;

	// Reservar memoria
	if (rank == COORDINATOR) {
    printf("n:                 %d\n", n);
    printf("numProcs:          %d\n", numProcs);
    printf("numThreads:        %d\n", numThreads);
    printf("bs:                %d\n", bs);
    printf("Tiras(n/bs):       %d\n", n/bs);
    printf("procRows(n/nP):    %d\n", procRows);
    printf("thRows(sS/nT):     %d\n", thRows);
  }
  // NO puede ser thRows<bs
  // N=512  P=4 bs=64 numThreads=8
  // Strip=128
  // thRows=16
  // ==> bs := 16
  bs2 = bs;
  if (thRows < bs) {
    bs = thRows;
    if (rank == COORDINATOR) {
      printf("New bs:            %d\n", bs);
    }
  }

  if (rank == COORDINATOR) {
    /* Matrices completas, debe inicializarlas el COORDINATOR */
    T = (double*) malloc(sizeof(double)*size);
    M = (double*) malloc(sizeof(double)*size);
    /* Matrices completa, para totalizar */
    C = (double*) malloc(sizeof(double)*size);

  } else  {
    /* Los workers solo necesitan su parte para trabajar */
    T = (double*) malloc(sizeof(double)*workerSize);
    M = (double*) malloc(sizeof(double)*workerSize);
    C = (double*) malloc(sizeof(double)*workerSize);
  }
  /* Las matrices de la izquierda de la multiplicación
  se necesita solo la parte de cada proceso, incluído
  el COORDINADOR */
  R1 = (double*) malloc(sizeof(double)*workerSize);
  R2 = (double*) malloc(sizeof(double)*workerSize);

  // A y B al estar a la derecha de la multiplicación se necesitan completas en todos los procesos
  A = (double*) malloc(sizeof(double)*size);
  B = (double*) malloc(sizeof(double)*size);

  /* El resultado de los productos se necesita solo para hacer la suma, no hay que "totalizarlo",
  se puede calcular cada strip de R y luego se totaliza esta solamente. */
  R1A = (double*) malloc(sizeof(double)*workerSize);
  R2B = (double*) malloc(sizeof(double)*workerSize);


  /* Inicializar datos */
  if (rank == COORDINATOR) {
    // por col
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        A[i+j*n] = 1; 
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        B[i+j*n] = 1;
    // por filas
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        T[i*n+j] = 1;
    // En cero para +=
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        C[i*n+j] = 0;
    /* Rellenar M con valores aleatorios entre 0 y 2Pi */
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        M[i*n+j] = randFP(0, 2*PI);
  }
  /* Estas las tienen que inicializar todos */
  for (i=0; i<procRows ; i++)
    for (j=0; j<n ; j++)
      R1A[i*n+j] = 0;
  for (i=0; i<procRows ; i++)
    for (j=0; j<n ; j++)
      R2B[i*n+j] = 0;
  
  /* Inicializar promedios */
  avgR1 = 0.0;
  avgR2 = 0.0;
  lavgR1 = 0.0;
  lavgR2 = 0.0;

  /* Barera */
  MPI_Barrier(MPI_COMM_WORLD); /* ????????????????????????????????????????????????? */

  /* Distribuir datos*/
  if (rank == COORDINATOR) commTimes[0] =  MPI_Wtime();
  MPI_Scatter(M, workerSize, MPI_DOUBLE, M, workerSize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Scatter(T, workerSize, MPI_DOUBLE, T, workerSize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Bcast(A, size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Bcast(B, size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  if (rank == COORDINATOR) commTimes[1] =  MPI_Wtime();

  /* Comienza cálculo */
  #pragma omp parallel
  {
    numThreads = omp_get_num_threads();

    #pragma omp for private (i,j,k,cosPhi,sinPhi) reduction(+:lavgR1) reduction(+:lavgR2)
    for(i=0;i<procRows;i++){
      for(j=0;j<n;j++){
        k = i*n+j;
        sinPhi = sin(M[k]);
        cosPhi = cos(M[k]);
        R1[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
        lavgR1 += R1[k];
        R2[k] = (1-T[k])*(1-sinPhi)+T[k]*cosPhi;
        lavgR2 += R2[k];
      }
    }
    #pragma omp single
    { 
      /* Un hilo de cada nodo reduce el promedio*/
      if (rank == COORDINATOR) commTimes[2] =  MPI_Wtime();
      MPI_Allreduce(&lavgR1, &avgR1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&lavgR2, &avgR2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (rank == COORDINATOR) commTimes[3] =  MPI_Wtime();
      /* El promedio se comparte entre los hilos del nodo */
      avgR1 = avgR1 / (size);
      avgR2 = avgR2 / (size);
    }
    /* Calcular R1 * A */
    #pragma omp for private (i,j,k,ii,jj,kk,ablk,bblk,cblk) nowait
    for (i = 0; i < procRows; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &R1A[i*n + j];
        for (k = 0; k < n; k+=bs)
        { 
          ablk = &R1[i*n + k];
          bblk = &A[j*n + k];
          for (ii=0; ii < bs; ii++)
          {
            for (jj = 0; jj < bs; jj++)
            {
              for (kk = 0; kk < bs; kk++)
              {
                cblk[ii*n + jj] += ablk[ii*n + kk] * bblk[jj*n + kk];
              }
            }
          }
        }
      }
    }
    /* Calcular R2 * B */
    #pragma omp for private (i,j,k,ii,jj,kk,ablk,bblk,cblk)
    for (i = 0; i < procRows; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &R2B[i*n + j];
        for (k = 0; k < n; k+=bs)
        { 
          ablk = &R2[i*n + k];
          bblk = &B[j*n + k];
          for (ii=0; ii < bs; ii++)
          {
            for (jj = 0; jj < bs; jj++)
            {
              for (kk = 0; kk < bs; kk++)
              {
                cblk[ii*n + jj] += ablk[ii*n + kk] * bblk[jj*n + kk];
              }
            }
          }
        }
      }
    }
    /* Calcular C */
    #pragma omp for private (i,j,k)
    for(i = 0; i < procRows; i++)
    {
      for(j=0;j<n;j++)
      {
        k = i*n+j;
        C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
      }
    }
  }
  if (rank == COORDINATOR) commTimes[4] = MPI_Wtime();
  MPI_Gather(C, workerSize, MPI_DOUBLE, C, workerSize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  if (rank == COORDINATOR) commTimes[5] = MPI_Wtime();


  if (rank == COORDINATOR) {

    totalTime = commTimes[5] - commTimes[0];
    commTime = (commTimes[1] - commTimes[0]) + (commTimes[3] - commTimes[2]) + (commTimes[5] - commTimes[4]);

    printf("totalTime:         %lf\ncommTime:          %lf\n",totalTime,commTime);
  }


  free(A);
  free(B);
  free(M);
  free(R1);
  free(R2);
  free(C);
  free(T);
  free(R1A);
  free(R2B);

  MPI_Finalize();

  return 0;
}
