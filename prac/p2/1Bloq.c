/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<math.h> /* sin y cos */
#include<sys/time.h>  /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
/* Init square matrix with a specific value */
void initvalmat(int *mat, int n, int val, int transpose); 
 
/* Multiply square matrices, blocked version */
void matmulblks(int *a, int *b, int *c, int n, int bs);

/* Multiply (block)submatrices */
void blkmul(int *ablk, int *bblk, int *cblk, int n, int bs);

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int *A,*B,*AB, avgR;
  int n, bs, i, j, k;

  /* Check command line parameters */
  if ( (argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0))
  {
    printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
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

  printf("Calculando A x B con %d bloques\n", bs);

  // Calc A * B
  matmulblks(A, B, AB, n, bs);

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
void matmulblks(int *a, int *b, int *c, int n, int bs)
{
  int i, j, k;    /* Guess what... */

  /* Init matrix c, just in case */  
  initvalmat(c, n, 0.0, 0);
  
  printf("| i j k |\n | i' j' k' |  AB   |   A  |  B   |    AB    |\n");
//| 0  0  0  | AB_00 | A_00 | B_00 | 0.000000 | 1.000000 |

  for (i = 0; i < n; i += bs)
  {
    for (j = 0; j < n; j += bs)
    {
      for  (k = 0; k < n; k += bs)
      {
        printf("| %d %d %d |\n", i,j,k);
        // printf("i=%d j=%d k=%d\n", i,j,k);
        blkmul(&a[i*n + k], &b[j*n + k], &c[i*n + j], n, bs);
      }
    }
  }
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmul(int *ablk, int *bblk, int *cblk, int n, int bs)
{
  int i, j, k;    /* Guess what... again... */

  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < bs; j++)
    {
      for  (k = 0; k < bs; k++)
      {
        printf(" | %d  %d  %d  |", i,j,k);
        printf(" AB_%d%d | A_%d%d | B_%d%d |", i,j, i,k, k,j);
        cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
        printf(" %d |\n", cblk[i*n+j]);
      }
    }
  }
}
    
/*****************************************************************/

