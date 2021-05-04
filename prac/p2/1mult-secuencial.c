#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */

/* Multiply matrices */
void multip(int *a, int *b, int *c, int n);

/* Time calculation */
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
  int n;

  /* Check command line parameters */
  if (argc < 2){
	  printf("\nFaltan argumentos. Usar %s n",argv[0]);
	  exit(1);
  }

  n = atoi(argv[1]);

  /* Pointers */
  int *A,*B,*AB;

  /* Indexes */
  int i, j, k;

  /* Time measurement */
  double timetick;

  /* Getting memory */
  A  = (int*)malloc(sizeof(int)*n*n);
  B  = (int*)malloc(sizeof(int)*n*n); 
  AB = (int*)malloc(sizeof(int)*n*n);

  printf("Incializando matrices %d x %d\n", n, n);
  // Fill with known pattern for later check
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
    // A por filas
    A[i*n+j] = i+j+1;
    // B por columnas
    B[i+j*n] = j+1;
    // AB en cero para despues poner +=
    AB[i*n+j] = 0;
    }
  }

  printf("Calculando A x B... \n");

  /* Start time measurement */
  timetick = dwalltime();

  multip(A, B, AB, n);

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

/* Multiply square matrices */
void multip(int *a, int *b, int *c, int n)
{
  int i,j,k;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        c[i*n+j] += a[i*n+k]*b[k+j*n];
      }
    }
  }
}

/*****************************************************************/