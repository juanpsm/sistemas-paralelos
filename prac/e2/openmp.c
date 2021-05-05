/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin, cos, ceil */
#include<time.h>     /* srand((unsigned) time(&t)) */
#include<omp.h>      /* hilos */

/* Init square matrix with a specific value */
void initvalmat(double *mat, int n, double val, int transpose); 
 
/* Multiply square matrices, blocked version, for OpenMP */
void calculate();

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
int n, Th, bs;

/************** MAIN *************/
int main(int argc, char *argv[])
{

  /* Check command line parameters */
  if  ( (argc != 4) ||
        ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((Th = atoi(argv[3])) <= 0)
      )
  {
    printf("\nError en los parámetros. Usage: ./%s n BS Th\n", argv[0]);
    exit(1);
  }
  if ((n % bs) != 0)
    {
    printf("\nError en los parámetros. Usage: ./%s n BS Th (n debe ser multiplo de BS)\n", argv[0]);
    exit(1);
  }

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
  /* A and B by column */
  initvalmat(A,   n, 1.0, 1);
  initvalmat(B,   n, 1.0, 1);
  /* The rest by rows */
  initvalmat(T,   n, 1.0, 0);
  initvalmat(C,   n, 0.0, 0);
  initvalmat(R1A, n, 0.0, 0);
  initvalmat(R2B, n, 0.0, 0);

  /* Fill M matrix with random values beetween 0 an 2*Pi */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      M[i*n+j] = randFP(0, 2*PI);
    }
  }

  /* Averages initialization */
  avgR1 = 0.0;
  avgR2 = 0.0;

  printf("Calculando con bloques de %dx%d\n", bs, bs);
  printf("  Tiras:   %d\n", n/bs);
  printf("  HILOS:   %d\n", Th);
  printf("  %.2f tiras x hilo\n\n", n/bs / (double) Th);

  /* Start time measurement */
  timetick = dwalltime();
  
  /* Calcular */
  calculate();


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

/* Init square matrix with a specific value */
void initvalmat(double *mat, int n, double val, int transpose)
{
  int i, j;      /* Indexes */

	if (transpose == 0) {
	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  mat[i*n + j] = val;
		}
	  }
	} else {
	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  mat[j*n + i] = val;
		}
	  }
	}
}

/*****************************************************************/

/* Multiply square matrices, blocked version */
void calculate()
{
  int i,j,k,ii,jj,kk, start_row, end_row;
  double local_avgR1, local_avgR2, sinPhi, cosPhi;
  double *ablk, *bblk, *cblk;

  #pragma omp parallel private(i,j,k,ii,jj,kk, start_row, end_row)
  {
    int id = omp_get_thread_num();

    int tiras = (int) ceil(n/bs / (double) Th);
    start_row = id * tiras * bs;
    end_row = (id+1) * tiras * bs;
    // Si acotamos el end_row, no entrará a los for los hilos que sobren
    if (end_row > n) end_row = n;

    // debug info
    // printf("(%d) El hilo %d hará %d filas  ->  for i = %d .. %d\n",id, id, end_row-start_row, start_row, end_row);


    /* Calculate R1, R2 and their averages */
    for (i = start_row; i < end_row; i++){
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
    /* Block iteration */
    for (i = start_row; i < end_row; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &R1A[i*n + j];
        for (k = 0; k < n; k+=bs)
        { 
          ablk = &R1[i*n + k];
          bblk = &A[j*n + k];
          /* Inner row itetarions */
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

    /* Calculate R2 * B */
    /* Block iteration */
    for (i = start_row; i < end_row; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &R2B[i*n + j];
        for (k = 0; k < n; k+=bs)
        { 
          ablk = &R2[i*n + k];
          bblk = &B[j*n + k];
          /* Inner row itetarions */
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

    /* Update shared average */

    /* Update shared average */

    /* Barrier */

    /* Calculate C */
    for(i = start_row; i < end_row; i++){
      for(j=0;j<n;j++){
        k = i*n+j;
        C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
      }
    }
  }
}

/*****************************************************************/