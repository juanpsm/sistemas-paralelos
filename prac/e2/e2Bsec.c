/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<math.h> /* sin y cos */
#include<sys/time.h>  /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
/* Init square matrix with a specific value */
void initvalmat(double *mat, int n, double val, int transpose); 
 
/* Multiply square matrices, blocked version */
void matmulblks(double *a, double *b, double *c, int n, int bs);

/* Multiply (block)submatrices */
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);

// Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

double randFP(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

#define PI 3.14159265358979323846

/************** MAIN *************/
int main(int argc, char *argv[])
{
  double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, avgR1, avgR2, sinPhi, cosPhi;
  int n, bs, i, j, k;

  double timetick;

  /* Check command line parameters */
  if ( (argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0))
  {
    printf("\nError en los parámetros. Usage: ./%s n BS (n debe ser multiplo de BS)\n", argv[0]);
    exit(1);
  }
  
  // Para crear números aleatorios
  time_t t;
  srand((unsigned) time(&t));

  /* Getting memory */  
  A=(double*)malloc(sizeof(double)*n*n); 
  B=(double*)malloc(sizeof(double)*n*n); 
  C=(double*)malloc(sizeof(double)*n*n); 
  R1=(double*)malloc(sizeof(double)*n*n); 
  R2=(double*)malloc(sizeof(double)*n*n); 
  T=(double*)malloc(sizeof(double)*n*n); 
  M=(double*)malloc(sizeof(double)*n*n); 
  R1A=(double*)malloc(sizeof(double)*n*n); 
  R2B=(double*)malloc(sizeof(double)*n*n); 
  
  printf("Incializando matrices ...\n");
  initvalmat(A, n, 1.0, 1);
  initvalmat(B, n, 1.0, 1);
  initvalmat(T, n, 1.0, 0);
  initvalmat(C, n, 0.0, 0);
  initvalmat(R1A, n, 0.0, 0);
  initvalmat(R2B, n, 0.0, 0);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      M[i*n+j]=randFP(0, 2*PI);
    }
  }
  
  avgR1 = 0;
  avgR2 = 0;

  printf("Calculando ... \n");

  timetick = dwalltime();

  // Calcular R y promedio
  // tambien se puede hacer cada cosa en for distintos
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

  // Calc R1 * A
  matmulblks(R1, A, R1A, n, bs);

  // Calc R2 * B
  matmulblks(R2, B, R2B, n, bs);

  // printf("Calculando C...\n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
    }
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
void matmulblks(double *a, double *b, double *c, int n, int bs)
{
  int i, j, k;    /* Guess what... */

  /* Init matrix c, just in case */  
  initvalmat(c, n, 0.0, 0);
  
  for (i = 0; i < n; i += bs)
  {
    for (j = 0; j < n; j += bs)
    {
      for  (k = 0; k < n; k += bs)
      {
        blkmul(&a[i*n + k], &b[j*n + k], &c[i*n + j], n, bs);
      }
    }
  }
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs)
{
  int i, j, k;    /* Guess what... again... */

  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < bs; j++)
    {
      for  (k = 0; k < bs; k++)
      {
        cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
      }
    }
  }
}
  