#include<stdlib.h>
#include<stdio.h>
#include<math.h> /* sin y cos */
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */

//Para calcular tiempo
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

void multipThread(double *A, double *B, int inicial, int final, int N, double *AB){
  int i,j,k;
  // Calc A * B de las matriz "partida" inicial-final
  // Limito el recorrido de filas solamente
  for(i=inicial;i<final;i++){
    for(j=0;j<N;j++){
      printf(" AB[%d] = %f , \n", i*N+j, AB[i*N+j]);
      for(k=0;k<N;k++){
        AB[i*N+j] += A[i*N+k]*B[k+j*N];
        printf(" AB %d %d += A %d %d * B %d %d = %f * %f = %f = %f, \n", i, j, i, k, k, j, A[i*N+k], B[k+j*N], A[i*N+k]*B[k+j*N], AB[i*N+j]);
      }
      printf(" AB %d %d = %f , \n", i, j, AB[i*N+j]);
    }
    printf("\n");
  }
}

int main(int argc, char* argv[]){

  if (argc < 3) {
	  printf("\nFaltan argumentos. Usar %s N T",argv[0]);
	  exit(1);
  }

  unsigned long long N = atoi(argv[1]);
  unsigned long long T = atoi(argv[2]);

  time_t t;
  srand((unsigned) time(&t));

  double timetick;

  int i,j,k;
  double *A,*B,*AB;

  A=(double*)malloc(sizeof(double)*N*N);
  B=(double*)malloc(sizeof(double)*N*N); 
  AB=(double*)malloc(sizeof(double)*N*N);

  printf("Incializando matrices ...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      // A por filas
      A[i*N+j]=1.0;
      // B por columnas
      B[i+j*N]=1.0;
      // AB en cero para despues poner +=
      AB[i*N+j]=0.0;
    }
  }

  printf("Calculando A*B... \n");
  timetick = dwalltime();

  // Calc A * B
  /*
    Para dividir la matris en partes iguales..
    T = 2
    i=0 -> i*N/T i*N/T+N/T -> 0 N/2
    i=1 -> i*N/T i*N/T+N/T -> N/2 N

    T = 4
    i=0 -> i*N/T i*N/T+N/T -> 0 N/4
    i=1 -> i*N/T i*N/T+N/T -> N/4 N/2
    i=2 -> i*N/T i*N/T+N/T -> N/2 N*3/4
    i=3 -> i*N/T i*N/T+N/T -> N*3/4 N
  */
  for(i=0;i<T;i++){
    // Con cada thread llamar a la funcion que calcula la submatriz
    //multipThread(A, B, 0, N/T, AB)
    printf("%lld - %lld \n", i*N/T,(i+1)*N/T);
    multipThread(A, B, i*N/T, (i+1)*N/T, N, AB);
  }

  printf(" TIEMPO = %f\n", dwalltime() - timetick);
  
  // Check
  for(i=0;i<N;i++){
    printf(" | ");
    for(j=0;j<N;j++){
      printf(" %f ", AB[i*N+j]);
    }
    printf("\n");
  }
  

  free(A);
  free(B);
  free(AB);
  
  return 0;
}