/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<math.h> /* sin y cos */
#include<sys/time.h>  /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
/* Inicializar matriz cuadrada con un valor específico */
void initvalmat(double *mat, int n, double val, int transpose); 
 
/* Multiplicar matrices cuadradas, por bloques */
void matmulblks(double *a, double *b, double *c, int n, int bs);

/* Multiplicar submatrices (bloques) */
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
  double *A,*B,*C,*R,*T,*M,*RA,*RB, avgR;
  int n, bs, i, j, k;

  double timetick;

  /* Verificar parámetros */
  if ( (argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0))
  {
    printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
    exit(1);
  }
  
  // Para crear números aleatorios
  time_t t;
  srand((unsigned) time(&t));

  /* Alocar memoria */  
  A=(double*)malloc(sizeof(double)*n*n); 
  B=(double*)malloc(sizeof(double)*n*n); 
  C=(double*)malloc(sizeof(double)*n*n); 
  R=(double*)malloc(sizeof(double)*n*n); 
  T=(double*)malloc(sizeof(double)*n*n); 
  M=(double*)malloc(sizeof(double)*n*n); 
  RA=(double*)malloc(sizeof(double)*n*n); 
  RB=(double*)malloc(sizeof(double)*n*n); 

  printf("Incializando matrices ...\n");
  initvalmat(A, n, 1.0, 1);
  initvalmat(B, n, 1.0, 1);
  initvalmat(T, n, 1.0, 0);
  initvalmat(C, n, 0.0, 0);
  initvalmat(RA, n, 0.0, 0);
  initvalmat(RB, n, 0.0, 0);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      M[i*n+j]=randFP(0, 2*PI);
    }
  }

  avgR = 0;

  printf("Calculando ... \n");

  timetick = dwalltime();

  // Calcular R y promedio
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      sinPhi = sin(M[k]);
      cosPhi = cos(M[k]);
      R[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
      avgR += R[k];
    }
  }
  avgR = avgR / (n*n);

  // Calc R * A
  matmulblks(R, A, RA, n, bs);

  // Calc R * B
  matmulblks(R, B, RB, n, bs);
  
  printf("Calculando C...\n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      C[i*n+j] = T[i*n+j] + avgR * (RA[i*n+j] + RB[i*n+j]);
    }
  }
  printf(" TIEMPO = %f\n", dwalltime() - timetick);

  free(A);
  free(B);
  free(M);
  free(R);
  free(C);
  free(T);
  free(RA);
  free(RB);

  return 0;
}

/*****************************************************************/

/* Inicializar matriz cuadrada con un valor específico */
void initvalmat(double *mat, int n, double val, int transpose)
{
  int i, j;      /* Índices */

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

/* Multiplicar matrices cuadradas, por bloques */
void matmulblks(double *a, double *b, double *c, int n, int bs)
{
  int i, j, k;

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

/* Multiplicar submatrices (bloques) */
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs)
{
  int i, j, k;

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
    
/*****************************************************************/
