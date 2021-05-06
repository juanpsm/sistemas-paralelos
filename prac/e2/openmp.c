/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin, cos, ceil */
#include<time.h>     /* srand((unsigned) time(&t)) */
#include<omp.h>      /* hilos */

/* Inicializar matriz cuadrada con un valor específico */
void initvalmat(double *mat, int n, double val, int transpose); 
 
/* Multiplicar matrices cuadradas, por bloques, para for OpenMP */
void calculate();

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
int n, Th, bs;

/************** MAIN *************/
int main(int argc, char *argv[])
{

  /* Verificar parámetros */
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
  /* A y B por columna */
  initvalmat(A,   n, 1.0, 1);
  initvalmat(B,   n, 1.0, 1);
  /* El resto por filas */
  initvalmat(T,   n, 1.0, 0);
  initvalmat(C,   n, 0.0, 0);
  initvalmat(R1A, n, 0.0, 0);
  initvalmat(R2B, n, 0.0, 0);

  /* Rellenar M con valores aleatorios entre 0 y 2Pi */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      M[i*n+j] = randFP(0, 2*PI);
    }
  }

  /* Inicializar promedios */
  avgR1 = 0.0;
  avgR2 = 0.0;

  printf("Calculando con bloques de %dx%d\n", bs, bs);
  printf("  Tiras:   %d\n", n/bs);
  printf("  HILOS:   %d\n", Th);
  printf("  %.2f tiras x hilo\n\n", n/bs / (double) Th);

  /* Empieza a medir el tiempo */
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
void calculate()
{
  int i,j,k,ii,jj,kk, start_row, end_row;
  double local_avgR1, local_avgR2, sinPhi, cosPhi;
  double *ablk, *bblk, *cblk;
  int tiras = (int) ceil(n/bs / (double) Th);
  
  #pragma omp parallel private(i,j,k,ii,jj,kk,local_avgR1, local_avgR2, sinPhi, cosPhi, start_row, end_row, ablk, bblk, cblk)
  {
    int id = omp_get_thread_num();
    
    /* Cada hilo obtiene ciertas filas sobre las que operará */
    start_row = id * tiras * bs;
    end_row = (id+1) * tiras * bs;
    // Si acotamos el end_row, no entrará a los for los hilos que sobren
    if (end_row > n) end_row = n;

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
    /* Iteraciones por bloques  */
    for (i = start_row; i < end_row; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &R1A[i*n + j];
        for (k = 0; k < n; k+=bs)
        { 
          ablk = &R1[i*n + k];
          bblk = &A[j*n + k];
          /* Iteraciones dentro de cada  bloque  */
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
    /* Iteraciones por bloques  */
    for (i = start_row; i < end_row; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &R2B[i*n + j];
        for (k = 0; k < n; k+=bs)
        { 
          ablk = &R2[i*n + k];
          bblk = &B[j*n + k];
          /* Iteraciones dentro de cada  bloque  */
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

    /* Actualizar promedio compartido */
    #pragma omp critical
    {
      avgR1 += local_avgR1 / (n*n);
    }
    #pragma omp critical
    {
     avgR2 += local_avgR2 / (n*n);
    }

    /* Barrera */
    #pragma omp barrier

    /* Calcular C */
    for(i = start_row; i < end_row; i++){
      for(j=0;j<n;j++){
        k = i*n+j;
        C[k] = T[k] + avgR1 * avgR2 * (R1A[k] + R2B[k]);
      }
    }
  }
}

/*****************************************************************/
