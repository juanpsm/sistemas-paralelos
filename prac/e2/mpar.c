#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin y cos */
#include<time.h>     /* srand((unsigned) time(&t)) */
#include<pthread.h>  /* hilos */

/* For pthreads */
void * calculate (void * ptr);

/* Time calculation */
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

/* Random number generation */
double randFP(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

#define PI 3.14159265358979323846

/* Shared variables */
double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, avgR1, avgR2;
int n, Th;

pthread_mutex_t mutex_avgR1;
pthread_mutex_t mutex_avgR2;

pthread_barrier_t   barrier_R_ready;
pthread_barrier_t   barrier_averages_ready;

/************** MAIN *************/
int main(int argc, char *argv[])
{

  /* Check command line parameters */
  if (argc < 3){
	  printf("\nFaltan argumentos. Usar %s n Th",argv[0]);
	  exit(1);
  }

  n = atoi(argv[1]);
  Th = atoi(argv[2]);

  /* Random numbers */
  time_t t;
  srand((unsigned) time(&t));

  /* Indexes */
  int i, j;

  /* Time measurement */
  double timetick;

  /* Getting memory */
  A   = (double*)malloc(sizeof(double)*n*n); 
  B   = (double*)malloc(sizeof(double)*n*n); 
  C   = (double*)malloc(sizeof(double)*n*n); 
  R1  = (double*)malloc(sizeof(double)*n*n); 
  R2  = (double*)malloc(sizeof(double)*n*n); 
  T   = (double*)malloc(sizeof(double)*n*n); 
  M   = (double*)malloc(sizeof(double)*n*n); 
  R1A = (double*)malloc(sizeof(double)*n*n); 
  R2B = (double*)malloc(sizeof(double)*n*n); 
  
  printf("Incializando matrices %d x %d\n", n, n);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      /* A and B by column */
      A[i+j*n]=1.0;
      B[i+j*n]=1.0;
      /* The rest by rows */
      T[i*n+j]=1.0;
      /* Zero then += */
      C[i*n+j]=0.0;
      R1A[i*n+j]=0.0;
      R2B[i*n+j]=0.0;
      /* Fill M matrix with random values beetween 0 an 2*Pi */
      M[i*n+j] = randFP(0, 2*PI);
    }
  }
  
  /* Averages initialization */
  avgR1 = 0.0;
  avgR2 = 0.0;

  /* Threads initialization */
  int id, ids[Th];
  pthread_attr_t attr;
  pthread_t threads[Th];
  pthread_attr_init(&attr);

  pthread_mutex_init(&mutex_avgR1, NULL);
  pthread_mutex_init(&mutex_avgR2, NULL);

  pthread_barrier_init (&barrier_R_ready, NULL, Th);
  pthread_barrier_init (&barrier_averages_ready, NULL, Th);

  printf("Calculando ... \n");
  printf("  HILOS:   %d\n", Th);
  printf("  %.2f filas x hilo\n\n", n / (double) Th);

  /* Start time measurement */
  timetick = dwalltime();

  /* Create threads */
  for (id = 0; id < Th; id++) {
    ids[id] = id;
    pthread_create(&threads[id], &attr, calculate, &ids[id]);
  }
  /* Join */
  for (id = 0; id < Th; id++) {
    pthread_join(threads[id], NULL);
  }

  printf(" TIEMPO = %f\n", dwalltime() - timetick);

  free(A);
  free(B);
  free(M);
  free(R1);
  free(R2);
  free(C);
  free(T);
  free(R1A);
  free(R2B);
  
  return 0;
}

/*****************************************************************/

/* Función para calcular el producto de matrices partido por filas.
   Calcula las filas de A*B desde start_row a end_row. */
void * calculate (void * ptr) {
  int id;
  id = *((int *) ptr);
  int i,j,k,start_row,end_row;
  double local_avgR1, local_avgR2, sinPhi, cosPhi;
  
  /* Each thread get some rows to operate with */
  start_row = id*n/Th;
  end_row = (id+1)*n/Th;

  local_avgR1 = 0.0;
  local_avgR2 = 0.0;
  
  // debug info
  // printf("(%d) El hilo %d hará %d filas  ->  for i = %d .. %d (no incl)\n",id, id, end_row-start_row, start_row, end_row);

  /* Calculate R1, R2 and their averages */
  // también se puede hacer cada cosa en for distintos
  for(i=start_row;i<end_row;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      sinPhi = sin(M[k]);
      cosPhi = cos(M[k]);

      R1[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
      local_avgR1 += R1[k];
      
      R2[k] = (1-T[k])*(1-sinPhi)+T[k]*cosPhi;
      local_avgR2 += R2[k];
    }
  }

  /* Calculate R1 * A */
  for(i = start_row; i < end_row; i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R1A[i*n+j] += R1[i*n+k]*A[k+j*n];
      }
    }
  }

  /* Calculate R2 * B */
  for(i = start_row; i < end_row; i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R2B[i*n+j] += R2[i*n+k]*B[k+j*n];
      }
    }
  }

  /* Update shared average */
  pthread_mutex_lock(&mutex_avgR1);
    avgR1 += local_avgR1;
  pthread_mutex_unlock(&mutex_avgR1);

  pthread_mutex_lock(&mutex_avgR2);
    avgR2 += local_avgR2;
  pthread_mutex_unlock(&mutex_avgR2);

  /* Only one divides the accumulated, butwe need a barrier. Every thread must have finished accumulating the average*/
  pthread_barrier_wait (&barrier_R_ready);

  /* Thanks to the previous barrier, you can be sure no one is accessing the
    averages, so no need to lock. */
  if (id = 0) {
    avgR1 = avgR1 / (n*n);
    avgR2 = avgR2 / (n*n);
  }

  /* Now we need another barrier in order to calculate C, the averages must be complete */
  pthread_barrier_wait (&barrier_averages_ready);
  // EN REALIDAD SI HAGO LA DIVISION DEL PROMEDIO JUNTO CON LA SUMA, ME ALCANZA CON UNA SOLA BARRERA, PERO HAGO MAS DIVISIONES, NO SE QUE SERÁ MÁS RÁPIDO

  /* Calculate C */
  for(i = start_row; i < end_row; i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
    }
  }
  pthread_exit(0);
}

/*****************************************************************/