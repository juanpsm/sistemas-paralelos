/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<math.h> /* sin y cos */
#include<sys/time.h>  /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */

/* Init square matrix with a specific value */
void initvalmat(double *mat, int N, double val, int transpose); 
 
/* Multiply square matrices, blocked version */
// void * matmulblks(void * ptr);

/* Multiply (block)submatrices */
void * blkmul(void * ptr);

// Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

double *A,*B,*AB;
int N, T, bs;

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int i, j, k;

  double timetick;

  /* Check command line parameters */
  if ( (argc != 4) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0) || ((bs = atoi(argv[3])) <= 0) || ((N % bs) != 0))
  {
    printf("\nError en los parámetros. Usage: ./%s N T BS (N debe ser multiplo de BS)\n", argv[0]);
    exit(1);
  }
  
  /* Getting memory */  
  A=(double*)malloc(sizeof(double)*N*N); 
  B=(double*)malloc(sizeof(double)*N*N); 
  AB=(double*)malloc(sizeof(double)*N*N); 

  printf("Incializando matrices %d x %d...\n", N, N);
  initvalmat(A, N, 1.0, 0);
  initvalmat(B, N, 1.0, 1);
  initvalmat(AB, N, 0.0, 0);

  int id, ids[T];
  pthread_attr_t attr;
  pthread_t threads[T] ;
  pthread_attr_init(&attr);

  printf("Calculando A*B usando %d bloques y %d threads... \n", bs, T);
  timetick = dwalltime();

  for (i = 0; i < N; i += bs)
  {
    for (j = 0; j < N; j += bs)
    {
      for  (k = 0; k < N; k += bs)
      {
        /* Crea los hilos */
        for (id = 0; id < T; id++) {
          ids[id] =  id;
          pthread_create(&threads[id], &attr, blkmul, &ids[id]);
        }
        /* Espera a que los hilos terminen */
        for (id = 0; id < T; id++)
          pthread_join(threads[id], NULL);
      }
      // printf("(%d) AB %d %d = %f \n", id, i, j, AB[i*N+j]);
    }
  }

  printf(" TIEMPO = %f\n", dwalltime() - timetick);

  // Check
  int error = 0;
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
          if (AB[i*N+j]!=N && !error){
              // printf(" Error en  %d %d \n", i, j);
              error = 1;
          }
          // Uncomment to print AB
          if ( N <= 16 ) printf(" %f ", AB[i*N+j]);
      }
      if ( N <= 16 ) printf("\n");
  }
  printf("Resultado ");
  if (error) printf("erroneo.\n"); else printf("correcto.\n");
  
  free(A);
  free(B);
  free(AB);

  return 0;
}

/*****************************************************************/

/* Init square matrix with a specific value */
void initvalmat(double *mat, int N, double val, int transpose)
{
  int i, j;      /* Indexes */

	if (transpose == 0) {
	  for (i = 0; i < N; i++)
	  {
      for (j = 0; j < N; j++)
      {
        mat[i*N + j] = val;
      }
	  }
	} else {
	  for (i = 0; i < N; i++)
	  {
      for (j = 0; j < N; j++)
      {
        mat[j*N + i] = val;
      }
	  }
	}
}

/*****************************************************************/

/* Multiply square matrices, blocked version */
void * matmulblks(void * ptr)
{
  // ahora las matrices, N y bs son globales, esto solo recibe el id

  int id;
  id = *((int *) ptr);
  int i,j,k,inicial,final; /* Indices... */
  
  inicial = id*N/T;
  final = (id+1)*N/T;
  
  // printf(" Inicia Hilo %d (%d - %d)\n", id, inicial, final); 

  for (i = 0; i < N; i += bs)
  {
    for (j = 0; j < N; j += bs)
    {
      for  (k = 0; k < N; k += bs)
      {
        //blkmul(k);
      }
      printf("(%d) AB %d %d = %f \n", id, i, j, AB[i*N+j]);
    }
  }
  pthread_exit(0);
}

/*****************************************************************/

/* Multiply (block)submatrices */
void * blkmul(void * ptr)
{
  int id;
  id = *((int *) ptr);
  int i,j,k; /* Indices... */

  for (i = 0; i < id; i++)
  {
    for (j = 0; j < id; j++)
    {
      for  (k = 0; k < id; k++)
      {
        AB[i*N + j] += A[i*N + k] * B[j*N + k];
        // &A[i*N + k], &B[j*N + k], &AB[i*N + j]
        printf("(%d) AB %d %d = %f \n", id, i, j, AB[i*N+j]);
      }
    }
  }
  pthread_exit(0);
}
    
/*****************************************************************/
