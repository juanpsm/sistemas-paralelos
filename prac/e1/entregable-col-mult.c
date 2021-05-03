#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#include <sys/time.h>

#define PI 3.14159265358979323846
#define MIN 0
#define MAX 1

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

  int i=0,j=0,k=0;
  double *A,*B,*C,*R,*T,*M,*RA,*RB, avgR;

  A=(double*)malloc(sizeof(double)*N*N); 
  B=(double*)malloc(sizeof(double)*N*N); 
  C=(double*)malloc(sizeof(double)*N*N); 
  R=(double*)malloc(sizeof(double)*N*N); 
  T=(double*)malloc(sizeof(double)*N*N); 
  M=(double*)malloc(sizeof(double)*N*N); 
  RA=(double*)malloc(sizeof(double)*N*N); 
  RB=(double*)malloc(sizeof(double)*N*N); 

  printf("Incializando matrices X col...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A[i+j*N]=1.0;
      B[i+j*N]=1.0;
      T[i+j*N]=1.0;
      M[i+j*N]=randFP(0, 2*PI);
    }
  }

  // Calcular R
  printf("Calculando R...");
  timetick = dwalltime();
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        R[i+j*N] = (1-T[i+j*N])*(1-cos(M[i+j*N]))+T[i+j*N]*sin(M[i+j*N]);
      }
    }
  }
  printf("   t = %f\n", dwalltime() - timetick);

  // Calcular avgR
  printf("Calculando avgR...");
  timetick = dwalltime();
  avgR = 0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      avgR += R[i+j*N];
    }
  }
  avgR = avgR / (N*N);
  printf("t = %f\n", dwalltime() - timetick);

  printf("Calculando C con multiplicacion...\n");
  timetick = dwalltime();
  // Calc R * A
  printf("Calculando R * A...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      RA[i+j*N]=0;
      for(k=0;k<N;k++){
        RA[i+j*N] += R[i+k*N]*A[k+j*N];
      }
    }
  }
  // Calc R * B
  printf("Calculando R * B...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      RB[i+j*N]=0;
      for(k=0;k<N;k++){
        RB[i+j*N] += R[i+k*N]*B[k+j*N];
      }
    }
  }
  printf("Calculando C...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      C[i+j*N] = 0;
      C[i+j*N] = T[i+j*N] + avgR * (RA[i+j*N] + RB[i+j*N]);
    }
  }
  printf(" TIEMPO = %f\n", dwalltime() - timetick);

  free(A);
  free(B);
  free(C);
  free(R);
  free(M);
  free(T);
  free(RA);
  free(RB);
}