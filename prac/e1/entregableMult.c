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

  unsigned long long N = atoi(argv[1]);

  time_t t;
  srand((unsigned) time(&t));

  double timetick;

  int i,j,k;
  double *A,*B,*C,*R,*T,*M,*RA,*RB, avgR;

  A=(double*)malloc(sizeof(double)*N*N); 
  B=(double*)malloc(sizeof(double)*N*N); 
  C=(double*)malloc(sizeof(double)*N*N); 
  R=(double*)malloc(sizeof(double)*N*N); 
  T=(double*)malloc(sizeof(double)*N*N); 
  M=(double*)malloc(sizeof(double)*N*N); 
  RA=(double*)malloc(sizeof(double)*N*N); 
  RB=(double*)malloc(sizeof(double)*N*N); 

  printf("Incializando matrices ...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      // a y b por columnas
      A[i+j*N]=1.0;
      B[i+j*N]=1.0;
      // Las demas por filas
      T[i*N+j]=1.0;
      M[i*N+j]=randFP(0, 2*PI);

      // EN cero para despues poner +=
      C[i*N+j]=0.0;
      RA[i*N+j]=0.0;
      RB[i*N+j]=0.0;
    }
  }
  
  avgR = 0;

  printf("Calculando ... \n");
  timetick = dwalltime();

  // Calcular R y promedio
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        R[i*N+j] = (1-T[i*N+j])*(1-cos(M[i*N+j]))+T[i*N+j]*sin(M[i*N+j]);
      }
      avgR += R[i*N+j] / (N*N);
    }
  }
  
  // printf("   t = %f\n", dwalltime() - timetick);

  // printf("Calculando C con multiplicacion...\n");
  // timetick = dwalltime();
  // Calc R * A
  // printf("Calculando R * A...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        RA[i*N+j] += R[i*N+k]*A[k+j*N];
      }
    }
  }
  // Calc R * B
  // printf("Calculando R * B...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        RB[i*N+j] += R[i*N+k]*B[k+j*N];
      }
    }
  }

  // printf("Calculando C...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      C[i*N+j] = T[i*N+j] + avgR * (RA[i*N+j] + RB[i*N+j]);
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