/* Using Blocked matrix multiplication example by Fernando G. Tinetti                   */

#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin y cos */
#include<omp.h>      /* hilos */

/* Inicializar matriz cuadrada con un valor específico */
void initvalmat(int *mat, int n, int val, int transpose); 
 
/* Multiplicar matrices cuadradas, por bloques, for OpenMP */
void matmulblks();

/* Para calcular el tiempo */
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

/* Variables compartidas */
int *A,*B,*AB;
int n, T, bs;

/************** MAIN *************/
int main(int argc, char *argv[])
{
  /* Verificar parámetros */
  if  ( (argc != 4) ||
        ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((T = atoi(argv[3])) <= 0)
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
  printf("  Tiras:   %d\n", n/bs);
  printf("  HILOS:   %d\n", T);
  printf("  %.2f tiras x hilo\n\n", n/bs / (double) T);

  /* Empieza a medir el tiempo */
  timetick = dwalltime();
  
  /* Calcular */
  matmulblks();


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
void matmulblks()
{
  int i,j,k,ii,jj,kk, start_row, end_row;
  int *ablk, *bblk, *cblk;

  // La cantidad de "tiras" de bloques es N/BS. Tengo que repartir entre los procesos estas tiras.
  // Los "saltos" de esta tira es i: [ 0, BS, 2*BS, ... , ((N/BS)-1)*BS ]     (N/BS)*BS = N ya no entra
  //                              id:[ 0, 1, ... , T-1]
   #pragma omp parallel private(i,j,k,ii,jj,kk, start_row, end_row)
   {
    int id = omp_get_thread_num();

    // La cantidad de tiras que puedo hacer lo calculo redondeando para arriba
    int tiras = (int) ceil(n/bs / (double) T);
    // El primer bloque sera respecto del id del hilo
    int bloqi = id * tiras;
    // y finalizara en el siguiente
    int bloqf = (id+1) * tiras;
    // Si alguno se pasa del maximo numero de tiras:
    if (bloqi > n/bs) bloqi = n/bs;
    if (bloqf > n/bs) bloqf = n/bs;
    // Para el numero de fila tengo que multiplicar el nuemro de bloque por el ancho del mismo
    start_row = bloqi * bs;
    end_row = bloqf * bs;
    // Si acotamos el end_row, no entrará a los for los hilos que sobren
    if (end_row > n) end_row = n;

    // debug info
    printf("(%d) El hilo %d hará %d tiras  ->  bloque inicial: %d  bloque final: %d (no incl) ->  for i = %d .. %d (no incl)\n",id, id, bloqf-bloqi, bloqi, bloqf, start_row, end_row);

    /* Iteraciones por bloques  */
    for (i = start_row; i < end_row; i+=bs)
    { 
      for (j = 0; j < n; j+=bs)
      {
        cblk = &AB[i*n + j];
        for (k = 0; k < n; k+=bs)
        {
          ablk = &A[i*n + k];
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
  }
}

/*****************************************************************/