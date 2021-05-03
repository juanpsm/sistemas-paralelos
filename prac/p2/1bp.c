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
void * matmulblks(void * ptr);

/* Multiply (block)submatrices */
// no hace falta
// void blkmul(int id);
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

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int i, j, k;

  double timetick;

  /* Check command line parameters */
  if  ( (argc != 4) ||
        ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((T = atoi(argv[3])) <= 0) ||
        ((n % bs) != 0)
      )
  {
    printf("\nError en los parámetros. Usage: ./%s n BS T\n", argv[0]);
    exit(1);
  }
  if ((n % bs) != 0)
    {
    printf("\nError en los parámetros. Usage: ./%s n BS T (n debe ser multiplo de BS)\n", argv[0]);
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

  printf("Calculando A x B con bloques de %dx%d\n", bs, bs);
  printf("  Tiras:   %d\n", n/bs);
  printf("  HILOS:   %d\n", T);
  printf("  %.2f tiras x hilo\n\n", n/bs / (double) T);

  timetick = dwalltime();
  // Calc A * B
  for (i = 0; i < T; i++) {
    ids[i] = i;
    pthread_create(&threads[i], NULL, matmulblks, &ids[i]);
  }

  // Unir
  for (i = 0; i < T; i++) {
    pthread_join(threads[i], NULL);
  }


  printf(" TIEMPO = %f\n", dwalltime() - timetick);
  
  // Check
  int * sumfila;
  sumfila=(int*)malloc(sizeof(int)*n); 
  int print = n<=32;
  printf("  A:\n");
  for(i=0;i<n;i++){
    sumfila[i]=0;
    for(j=0;j<n;j++){
      sumfila[i] += A[i*n+j];
      if (print) printf(" %d ", A[i*n+j]);
    }
    if (print) printf("Suma fila = %d\n", sumfila[i]);
  }
  printf("\n  B:\n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (print) printf(" %d ", B[i+j*n]);
    }
    if (print) printf("\n");
  }

  printf("\n  AB:\n");
  int error = 0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (AB[i*n+j]!=sumfila[i]*(j+1) && !error){
        printf(" Error en  %d %d \n", i, j);
        error = 1;
      }
      // Uncomment to print AB
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
void * matmulblks(void * ptr)
{
  int id;
  id = *((int *) ptr);
  int i,j,k,f,c,h, inicial, final;
  int *ablk, *bblk, *cblk;
  double cant;

  // La cantidad de "tiras" de bloques es N/BS. Tengo que repartir entre los procesos estas tiras.
  // Los "saltos" de esta tira es i: [ 0, BS, 2*BS, ... , ((N/BS)-1)*BS ]     (N/BS)*BS = N ya no entra
  //                              id:[ 0, 1, ... , T-1]
  
  // La cantidad de tiras que puedo hacer lo calculo redondeando para arriba
  int tiras = (int) ceil(n/bs / (double) T);
  // El primer bloque sera respecto del ide del hilo
  int bloqi = id * tiras;
  // y finalizara en el siguiente
  int bloqf = (id+1) * tiras;
  // Si alguno se pasa del maximo numero de tiras:
  if (bloqi > n/bs) bloqi = n/bs;
  if (bloqf > n/bs) bloqf = n/bs;
  // Para el numero de fila tengo que multiplicar el nuemro de bloque por el ancho del mismo
  inicial = bloqi * bs;
  final = bloqf * bs;
  // Si acotamos el final, no entrará a los for los hilos que sobren
  if (final > n) final = n;

  // info
  printf("(%d) El hilo %d hará %d tiras  ->",id, id, bloqf-bloqi);
  printf("  bloque inicial: %d", bloqi);
  printf("  bloque final (no incl): %d  ->", bloqf);
  printf("  for i = %d .. %d (no incl)\n", inicial, final);

  // FOR PARA BLOQUES
  for (i = inicial; i < final; i+=bs)
  { 
    for (j = 0; j < n; j+=bs)
    {
      //printf("    (%d) for j = %d; j < %d\n",id, j, n);
      //printf("    (%d) acceso a AB %d %d\n",id, i, j);
      cblk = &AB[i*n + j];
      for (k = 0; k < n; k+=bs)
      { 
        ablk = &A[i*n + k];
        bblk = &B[j*n + k];
        // printf("i=%d j=%d k=%d\n", i,j,k);
        //FOR PARA PROCESAR BLOQUE
        for (f=0; f < bs; f++)
        {
          for (c = 0; c < bs; c++)
          {
            for (h = 0; h < bs; h++)
            {
             
              cblk[f*n + c] += ablk[f*n + h] * bblk[c*n + h];
              // printf("(%d) i:%d j:%d k:%d f:%d c:%d h:%d\n",id,i,j,k,f,c,h);

              // AB[(i+ig)*n + j+jg] += A[(i+ig)*n + k+kg] * B[(j+jg)*n+k+kg];
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
/*****************************************************************/

