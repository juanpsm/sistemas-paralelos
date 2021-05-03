#include<stdlib.h>
#include<stdio.h>
#include<math.h> /* sin y cos */
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */

#define PI 3.14159265358979323846

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

int main(int argc, char* argv[]){

  if (argc < 2) {
	  printf("\nFalta argumento. %s N",argv[0]);
	  exit(1);
  }

  int N = atoi(argv[1]);

  time_t t;
  srand((unsigned) time(&t));

  double timetick;

  int i,j,k;
  double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, avgR1, avgR2, sinPhi, cosPhi;

  A=(double*)malloc(sizeof(double)*N*N); 
  B=(double*)malloc(sizeof(double)*N*N); 
  C=(double*)malloc(sizeof(double)*N*N); 
  R1=(double*)malloc(sizeof(double)*N*N); 
  R2=(double*)malloc(sizeof(double)*N*N); 
  T=(double*)malloc(sizeof(double)*N*N); 
  M=(double*)malloc(sizeof(double)*N*N); 
  R1A=(double*)malloc(sizeof(double)*N*N); 
  R2B=(double*)malloc(sizeof(double)*N*N); 
  
  printf("Incializando matrices ...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      // a y b por columnas
      A[i+j*N]=1.0;
      B[i+j*N]=1.0;
      // Las demas por filas
      T[i*N+j]=1.0;
      M[i*N+j]=randFP(0, 2*PI);

      // En cero para despues poner +=
      C[i*N+j]=0.0;
      R1A[i*N+j]=0.0;
      R2B[i*N+j]=0.0;
    }
  }
  
  avgR1 = 0;
  avgR2 = 0;

  printf("Calculando ... \n");
  timetick = dwalltime();

  // Calcular R y promedio
  // tambien se puede hacer cada cosa en for distintos
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      k = i*N+j;
      sinPhi = sin(M[k]);
      cosPhi = cos(M[k]);
      R1[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
      avgR1 += R1[k];
      R2[k] = (1-T[k])*(1-sinPhi)+T[k]*cosPhi;
      avgR2 += R2[k];
    }
  }
  avgR1 = avgR1 / (N*N);
  avgR2 = avgR2 / (N*N);
// en el paralelo si divido x filas los promedios quedan "parciales"

  // printf("   t = %f\n", dwalltime() - timetick);

  // printf("Calculando C con multiplicacion...\n");
  // timetick = dwalltime();
  // Calc R1 * A
  // printf("Calculando R * A...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        R1A[i*N+j] += R1[i*N+k]*A[k+j*N];
      }
    }
  }
  // Calc R2 * B
  // printf("Calculando R * B...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        R2B[i*N+j] += R2[i*N+k]*B[k+j*N];
      }
    }
  }

  // aca iria una barrera que debe abrirse cuando esten listos los promedios

  // printf("Calculando C...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      k = i*N+j;
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