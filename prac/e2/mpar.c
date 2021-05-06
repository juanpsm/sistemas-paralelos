#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin y cos */
#include<time.h>     /* srand((unsigned) time(&t)) */
#include<pthread.h>  /* hilos */

/* For pthreads */
void * calculate (void * ptr);

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

/* Variables compartidas */
double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, avgR1, avgR2;
int n, Th;

pthread_mutex_t mutex_avgR1;
pthread_mutex_t mutex_avgR2;

pthread_barrier_t   barrier_averages_ready;

/************** MAIN *************/
int main(int argc, char *argv[])
{

  /* Verificar parámetros */
  if (argc < 3){
	  printf("\nFaltan argumentos. Usar %s n Th",argv[0]);
	  exit(1);
  }

  n = atoi(argv[1]);
  Th = atoi(argv[2]);

  /* Para números aleatorios */
  time_t t;
  srand((unsigned) time(&t));

  /* Índices */
  int i, j;

  /* Para medir el tiempo */
  double timetick;

  /* Alocar memoria */
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
      /* A y B por columna */
      A[i+j*n]=1.0;
      B[i+j*n]=1.0;
      /* El resto por filas */
      T[i*n+j]=1.0;
      /* Zero then += */
      C[i*n+j]=0.0;
      R1A[i*n+j]=0.0;
      R2B[i*n+j]=0.0;
      /* Rellenar M con valores aleatorios entre 0 y 2Pi */
      M[i*n+j] = randFP(0, 2*PI);
    }
  }
  
  /* Inicializar promedios */
  avgR1 = 0.0;
  avgR2 = 0.0;

  /* Inicializar hilos*/
  int id, ids[Th];
  pthread_attr_t attr;
  pthread_t threads[Th];
  pthread_attr_init(&attr);

  pthread_mutex_init(&mutex_avgR1, NULL);
  pthread_mutex_init(&mutex_avgR2, NULL);

  pthread_barrier_init (&barrier_averages_ready, NULL, Th);

  printf("Calculando ... \n");
  printf("  HILOS:   %d\n", Th);
  printf("  %.2f filas x hilo\n\n", n / (double) Th);

  /* Empieza a medir el tiempo */
  timetick = dwalltime();

  /* Crear hilos */
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

  pthread_mutex_destroy(&mutex_avgR1);
  pthread_mutex_destroy(&mutex_avgR2);
  pthread_barrier_destroy(&barrier_averages_ready);
  
  return 0;
}

/*****************************************************************/

/* Función para los hilos. Cada uno calcula las filas de C desde start_row a end_row. */
void * calculate (void * ptr) {
  int id;
  id = *((int *) ptr);
  int i,j,k,start_row,end_row;
  double local_avgR1, local_avgR2, sinPhi, cosPhi;
  
  /* Cada hilo obtiene ciertas filas sobre las que operará */
  start_row = id*n/Th;
  end_row = (id+1)*n/Th;

  /* Inicializar acumuladores para promedios*/
  local_avgR1 = 0.0;
  local_avgR2 = 0.0;
  
  // debug info
  // printf("(%d) El hilo %d hará %d filas  ->  for i = %d .. %d\n",id, id, end_row-start_row, start_row, end_row);

  /* Calcular R1, R2 y acumular para los promedios */
  for(i = start_row; i < end_row; i++){
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

  /* Calcular R1 * A */
  for(i = start_row; i < end_row; i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R1A[i*n+j] += R1[i*n+k]*A[k+j*n];
      }
    }
  }

  /* Calcular R2 * B */
  for(i = start_row; i < end_row; i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R2B[i*n+j] += R2[i*n+k]*B[k+j*n];
      }
    }
  }

  /* Actualizar promedio compartido */
  pthread_mutex_lock(&mutex_avgR1);
    avgR1 += local_avgR1 / (n*n);
  pthread_mutex_unlock(&mutex_avgR1);

  pthread_mutex_lock(&mutex_avgR2);
    avgR2 += local_avgR2 / (n*n);
  pthread_mutex_unlock(&mutex_avgR2);

  /* Ahora se necesita una barrera, ya que para hacer los calculos de C deben estar completos los promedios */
  pthread_barrier_wait (&barrier_averages_ready);

  /* Calcular C */
  for(i = start_row; i < end_row; i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
    }
  }
  pthread_exit(0);
}

/*****************************************************************/