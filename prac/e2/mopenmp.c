#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin, cos, fabs */
#include<time.h>     /* srand((unsigned) time(&t)) */
#include<omp.h>      /* hilos */

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

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int n, Th;

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

  /* Punteros */
  double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, avgR1, avgR2, sinPhi, cosPhi, *C_CHECK;

  /* Índices */
  int i, j, k;

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
      /* EN cero para luego += */
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

  /* Setear el numero de hilos */
  omp_set_num_threads(Th);

  printf("Calculando con");
  
  #pragma omp parallel
  {
    // Test hilos
    // printf("Nº:   %d\n", omp_get_thread_num());
    #pragma omp single
    {
      printf(" %d hilos\n", omp_get_num_threads());
    }
  }

  /* Empieza a medir el tiempo */
  timetick = dwalltime();

  /* Calcular R1, R2 y acumular para los promedios */
  #pragma omp parallel for private (i,j,k,cosPhi,sinPhi) reduction (+:avgR1) reduction (+:avgR2)
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      sinPhi = sin(M[k]);
      cosPhi = cos(M[k]);
      R1[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
      avgR1 += R1[k];
      R2[k] = (1-T[k])*(1-sinPhi)+T[k]*cosPhi;
      avgR2 += R2[k];
    }
  }
  #pragma omp single
  {
    avgR1 = avgR1 / (n*n);
    avgR2 = avgR2 / (n*n);
  }
  /* Calcular R1 * A */
  #pragma omp parallel for private (i,j,k,cosPhi,sinPhi)
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R1A[i*n+j] += R1[i*n+k]*A[k+j*n];
      }
    }
  }
  /* Calcular R2 * B */
  #pragma omp parallel for private (i,j,k,cosPhi,sinPhi) 
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R2B[i*n+j] += R2[i*n+k]*B[k+j*n];
      }
    }
  }

  /* Calcular C */
  #pragma omp parallel for private (i,j,k)
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
    }
  }
  printf(" TIEMPO = %f\n", dwalltime() - timetick);

    /************* Secuencial ***************/

  C_CHECK   = (double*)malloc(sizeof(double)*n*n);

  /* Resetear matrices*/
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      C_CHECK[i*n+j]=0.0;
      R1A[i*n+j]=0.0;
      R2B[i*n+j]=0.0;
    }
  }

  /* Resetear promedios */
  avgR1 = 0.0;
  avgR2 = 0.0;

  printf("Calculando secuencialmente... \n");

  /* Empieza a medir el tiempo */
  timetick = dwalltime();

  /* Calcular R1, R2 y acumular para los promedios */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      sinPhi = sin(M[k]);
      cosPhi = cos(M[k]);
      R1[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
      avgR1 += R1[k];
      R2[k] = (1-T[k])*(1-sinPhi)+T[k]*cosPhi;
      avgR2 += R2[k];
    }
  }
  avgR1 = avgR1 / (n*n);
  avgR2 = avgR2 / (n*n);

  /* Calcular R1 * A */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R1A[i*n+j] += R1[i*n+k]*A[k+j*n];
      }
    }
  }
  /* Calcular R2 * B */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R2B[i*n+j] += R2[i*n+k]*B[k+j*n];
      }
    }
  }

  /* Calcular C_CHECK */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      C_CHECK[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
    }
  }
  printf(" TIEMPO = %f\n", dwalltime() - timetick);
/*********** END Secuencial ***************/

  /* Comprobación */
  int error = 0;
  for(i=0; i < n*n; i++){
    if (fabs(C[i] - C_CHECK[i]) > 0.000001){
      error = 1;
      break;
    }
  }
  printf("Resultado ");
  if (error) printf("erroneo.\n"); else printf("correcto.\n");

  free(A);
  free(B);
  free(M);
  free(R1);
  free(R2);
  free(C);
  free(C_CHECK);
  free(T);
  free(R1A);
  free(R2B);
  
  return 0;
}

/*****************************************************************/