/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */

/* Inicializar matriz cuadrada con un valor específico */
void initvalmat(int *mat, int n, int val, int transpose); 
 
/* Multiplicar matrices cuadradas, por bloques */
void matmulblks(int *a, int *b, int *c, int n, int bs);

/* Multiplicar submatrices (bloques) */
void blkmul(int *ablk, int *bblk, int *cblk, int n, int bs);

/* Para calcular el tiempo */
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int n, bs;

  /* Verificar parámetros */
  if  ( (argc != 3) ||
        ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) ||
        ((n % bs) != 0)
      )
  {
    printf("\nError en los parámetros. Usage: ./%s n bs (n debe ser multiplo de bs)\n", argv[0]);
    exit(1);
  }

  /* Punteros */
  int *A,*B,*AB;

  /* Índices */
  int i, j, k;

  /* Para medir el tiempo */
  double timetick;

  /* Alocar memoria */  
  A  = (int*)malloc(sizeof(int)*n*n); 
  B  = (int*)malloc(sizeof(int)*n*n); 
  AB = (int*)malloc(sizeof(int)*n*n); 

  printf("Incializando matrices %d x %d\n", n, n);
  // A por filas
  initvalmat(A, n, 0, 0);
  // B por columnas
  initvalmat(B, n, 0, 1);
  // AB por filas
  initvalmat(AB, n, 0, 0);
  
  // Fill with known pattern for later check
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      A[i*n+j] = i+j+1;
    }
  }
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      B[i+n*j] = j+1;
    }
  }

  printf("Calculando A x B con bloques de %dx%d\n", bs, bs);

  /* Empieza a medir el tiempo */
  timetick = dwalltime();

  /* Calculate A x B */
  matmulblks(A, B, AB, n, bs);

  printf(" TIEMPO = %f\n", dwalltime() - timetick);

  /* Check */
  int * sumfila;
  sumfila=(int*)malloc(sizeof(int)*n); 
  int print = (n<=32);
  if (print) printf("  A:\n");
  for(i=0;i<n;i++){
    sumfila[i]=0;
    for(j=0;j<n;j++){
      sumfila[i] += A[i*n+j];
      if (print) printf(" %d ", A[i*n+j]);
    }
    if (print) printf("Suma fila = %d\n", sumfila[i]);
  }
  if (print) printf("\n  B:\n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (print) printf(" %d ", B[i+j*n]);
    }
    if (print) printf("\n");
  }

  if (print) printf("\n  AB:\n");
  int error = 0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (AB[i*n+j]!=sumfila[i]*(j+1) && !error){
        printf(" Error en  AB_%d_%d \n", i, j);
        error = 1;
      }
      if (print)  printf(" %d ", AB[i*n+j]);
    }
    if (print)  printf("\n");
  }
  printf("Resultado ");
  if (error) printf("erroneo.\n"); else printf("correcto.\n");

  free(A);
  free(B);
  free(AB);

  return 0;
}

/*****************************************************************/

/* Inicializar matriz cuadrada con un valor específico */
void initvalmat(int *mat, int n, int val, int transpose)
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
void matmulblks(int *a, int *b, int *c, int n, int bs)
{
  int i, j, k;

  /* Init matrix c, just in case */  
  //initvalmat(c, n, 0.0, 0);

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
void blkmul(int *ablk, int *bblk, int *cblk, int n, int bs)
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