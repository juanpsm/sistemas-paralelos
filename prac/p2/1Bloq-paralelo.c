/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<math.h> /* sin y cos */
#include<sys/time.h>  /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */

/* Init square matrix with a specific value */
void initvalmat(int *mat, int n, int val, int transpose); 
 
/* Multiply square matrices, blocked version */
void matmulblks(void);

/* Multiply (block)submatrices */
void blkmul(int id);
// void blkmul(int *ablk, int *bblk, int *cblk, int n, int bs);

// Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

int *A,*B,*AB;
int n, T, bs;
int ig, jg, kg;

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int i, j, k;

  double timetick;

  /* Check command line parameters */
  if ( (argc != 4) || ((n = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0) || ((bs = atoi(argv[3])) <= 0) || ((n % bs) != 0))
  {
    printf("\nError en los parámetros. Usage: ./%s n T BS (n debe ser multiplo de BS)\n", argv[0]);
    exit(1);
  }

  /* Getting memory */  
  A=(int*)malloc(sizeof(int)*n*n); 
  B=(int*)malloc(sizeof(int)*n*n); 
  AB=(int*)malloc(sizeof(int)*n*n); 

  printf("Incializando matrices %d x %d\n", n, n);
  initvalmat(A, n, 0, 0);
  initvalmat(B, n, 0, 1);
  initvalmat(AB, n, 0, 0);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      A[i*n+j]=i+j+1;
    }
  }
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      B[i+n*j]=j+1;
    }
  }

  int id, ids[T];
  pthread_attr_t attr;
  pthread_t threads[T] ;
  pthread_attr_init(&attr);

  printf("Calculando A x B con %d bloques\n", bs);
  timetick = dwalltime();
  // Calc A * B
  matmulblks();

  printf(" TIEMPO = %f\n", dwalltime() - timetick);
  
  // Check
  if ( n <= 16 ){
    printf("  A:\n");
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        printf(" %d ", A[i*n+j]);
      }
      printf("\n");
    }
    printf("\n  B:\n");
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        printf(" %d ", B[i+j*n]);
      }
      printf("\n");
    }
  }
  printf("\n  AB:\n");
  int error = 0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (AB[i*n+j]!=n && !error){
        // printf(" Error en  %d %d \n", i, j);
        error = 1;
      }
      // Uncomment to print AB
      if ( n <= 16 ) printf(" %d ", AB[i*n+j]);
    }
    if ( n <= 16 ) printf("\n");
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
void initvalmat(int *mat, int n, int val, int transpose)
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
void matmulblks(void)
{
  int id;
//printf("| i j k |\n | i' j' k' |  AB   |   A  |  B   |    AB    |\n");
//| 0  0  0  | AB_00 | A_00 | B_00 | 0.000000 | 1.000000 |

  for (ig = 0; ig < n; ig+=bs)
  {
    for (jg = 0; jg < n; jg+=bs)
    {
      for  (kg = 0; kg < n; kg+=bs)
      {
        printf("\n| %d %d |\n", ig,jg);
        // printf("i=%d j=%d k=%d\n", i,j,k);
        for (id=0; id < T; id++)
          blkmul(id); // id: 0...T-1
        for (i = 0; i < bs; i++)
        {
          for (j = 0; j < bs; j++)
          {
            for  (k = 0; k < bs; k++)
            {
              printf(" | %d%d %d%d |\n", i,k, k,j);
              AB[(i+ig)*n + j+jg] += A[(i+ig)*n + k+kg] * B[(j+jg)*n+k+kg];
              // cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
              // ablk : [      a          ]
              //               ^ik
              // ig=0 jg=0 kg=1
              // i=0 j=0 k=0
              // AB [0] += A[1] * B[1]
            }
          }
        }
      }
    }
  }
}

/*
a b c d  a b c d    | aa+ba+  ab+bb  cc.. dd.. |     | aa+ab+ac+ad
a b c d  a b c d    | aa+ba  ab+bb            |     |
a b c d  a b c d    | 
a b c d  a b c d    | 

una es hacer un bloque por hilo
N=4 bs=2 -> T=4 = 4*4/2*2

T = N^2/bs^2

(N, T) -> bs = sqrt(N^2/T)
no se puede optimizar bs

Otra forma
La cantidad de bloques es N^2/bs^2, darle a cada hilo una cantidad de bloques
cada hilo recibe (N^2/bs^2)/T

lo que falta:
como divir los hilos
mutex para escribir porque los bloques se pisan en la matriz resultado

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmul(int id)
{
  int i, j, k; // hacer los otros ijk globales y estos llamarlos de otra forma

  //int *a = &A[], *b = &B[], *c = &AB[];
  //blkmul(&a[i*n + k], &b[j*n + k], &c[i*n + j], n, bs);

  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < bs; j++)
    {
      for  (k = 0; k < bs; k++)
      {
        printf(" | %d%d %d%d |\n", i,k, k,j);
        AB[(i+ig)*n + j+jg] += A[(i+ig)*n + k+kg] * B[(j+jg)*n+k+kg];
        // cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
        // ablk : [      a          ]
        //               ^ik
      }
    }
  }
}
    
/*****************************************************************/

