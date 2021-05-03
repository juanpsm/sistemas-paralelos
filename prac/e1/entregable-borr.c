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

double elememtoMult(double *matriz1, double *matriz2 ,int i, int j, int n, int orden){
  double result=0;
  int k=0;

  for(k=0;k<n;k++){
    if (orden==1){
      result += matriz1[i*n+k]*matriz2[k*n+j];
    }else{
      result += matriz1[i+k*n]*matriz2[k+j*n];
    }
  }
  return result;
}

void printM(double *m, int n, int orden){
  int i,j;
  for(i=0;i<n;i++){
    printf(" | ");
    for(j=0;j<n;j++){
      if(orden==1){
        printf(" %f ", m[i*n+j]);
      } else {
        printf(" %f ", m[i+j*n]);
      }
    }
    printf("\n");
  }
}

int main(int argc, char* argv[]){

  if (argc < 2) {
	  printf("\nFalta argumento. %s N",argv[0]);
	  exit(1);
  }

  unsigned long long N = atoi(argv[1]);
  int show = atoi(argv[2]);

  time_t t;
  srand((unsigned) time(&t));
  
  //printf("srand %f\n", randFP(MIN,MAX));

  double timetick;

  int i=0,j=0,k=0;
  double *A,*B,*C,*R,*T,*M,*RA,*RB,*Ac,*Bc,*Tc,*Mc, avgR;

  A=(double*)malloc(sizeof(double)*N*N); 
  B=(double*)malloc(sizeof(double)*N*N); 
  C=(double*)malloc(sizeof(double)*N*N); 
  R=(double*)malloc(sizeof(double)*N*N); 
  T=(double*)malloc(sizeof(double)*N*N); 
  M=(double*)malloc(sizeof(double)*N*N); 
  RA=(double*)malloc(sizeof(double)*N*N); 
  RB=(double*)malloc(sizeof(double)*N*N); 
  Ac=(double*)malloc(sizeof(double)*N*N); 
  Bc=(double*)malloc(sizeof(double)*N*N); 
  Tc=(double*)malloc(sizeof(double)*N*N); 
  Mc=(double*)malloc(sizeof(double)*N*N);

  printf("Incializando matrices x filas...\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A[i*N+j]=0.0;
      B[i*N+j]=0.0;
      T[i*N+j]=0.0;
      M[i*N+j]=PI;
      // A[i*N+j]=randFP(MIN,MAX);
      // B[i*N+j]=randFP(MIN,MAX);
      // T[i*N+j]=randFP(MIN,MAX);
      // M[i*N+j]=randFP(0, 2*PI);
    }
    for(j=N/2;j<N;j++){
      A[i*N+j]=1.0;
      B[i*N+j]=2.0;
      T[i*N+j]=3.0;
      M[i*N+j]=PI;
    }
  }

  if (show){
    printf("A:\n");
    printM(A,N,1);
    printf("B:\n");
    printM(B,N,1);
    printf("T:\n");
    printM(T,N,1);
    printf("M:\n");
    printM(M,N,1);
  }

  // Calcular R
  printf("Calculando R... ");
  timetick = dwalltime();
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        R[i*N+j] = (1-T[i*N+j])*(1-cos(M[i*N+j]))+T[i*N+j]*sin(M[i*N+j]);
      }
    }
  }
  printf("   t = %f\n", dwalltime() - timetick);

  if (show){
    printf("R:\n");
    printM(R,N,1);
  }

  // Calcular avgR
  printf("Calculando avgR... ");
  timetick = dwalltime();
  avgR = 0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      avgR += R[i*N+j];
    }
  }
  avgR = avgR / (N*N);
  printf("t = %f\n", dwalltime() - timetick);
  if (show){
    printf("avgR: %f\n", avgR);
  }

  // Calc C
  printf("Calculando C con funcion... \n");
  timetick = dwalltime();
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      C[i*N+j] = 0;
      C[i*N+j] = T[i*N+j] + avgR * ( elememtoMult(R,A,i,j,N,1) + elememtoMult(R,B,i,j,N,1));
    }
  }
  printf(" TIEMPO = %f\n\n", dwalltime() - timetick);

  if (show){
    printf("C:\n");
    printM(C,N,1);
  }

  printf("Calculando C con multiplicacion... \n");
  timetick = dwalltime();
  // Calc R * A
  printf("Calculando R * A... \n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      RA[i*N+j]=0;
      for(k=0;k<N;k++){
        RA[i*N+j] += R[i*N+k]*A[k*N+j];
      }
    }
  }
  // Calc R * B
  printf("Calculando R * B... \n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      RB[i*N+j]=0;
      for(k=0;k<N;k++){
        RB[i*N+j] += R[i*N+k]*B[k*N+j];
      }
    }
  }
  printf("Calculando C... \n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      C[i*N+j] = 0;
      C[i*N+j] = T[i*N+j] + avgR * (RA[i*N+j] + RB[i*N+j]);
    }
  }
  printf(" TIEMPO = %f\n", dwalltime() - timetick);

  if (show){
    printf("RA:\n");
    printM(RA,N,1);
    printf("RB:\n");
    printM(RB,N,1);
    printf("C:\n");
    printM(C,N,1);
  }

  printf("\n\nIncializando matrices X col... \n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      Ac[i+j*N]=A[i*N+j];
      Bc[i+j*N]=B[i*N+j];
      Tc[i+j*N]=T[i*N+j];
      Mc[i+j*N]=M[i*N+j];
    }
  }
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A[i+j*N]=Ac[i+j*N];
      B[i+j*N]=Bc[i+j*N];
      T[i+j*N]=Tc[i+j*N];
      M[i+j*N]=Mc[i+j*N];
    }
  }

  if (show){
    printf("A:\n");
    printM(A,N,0);
    printf("B:\n");
    printM(B,N,0);
    printf("T:\n");
    printM(T,N,0);
    printf("M:\n");
    printM(M,N,0);
  }

  // Calcular R
  printf("Calculando R... ");
  timetick = dwalltime();
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        R[i+j*N] = (1-T[i+j*N])*(1-cos(M[i+j*N]))+T[i+j*N]*sin(M[i+j*N]);
      }
    }
  }
  printf("   t = %f\n", dwalltime() - timetick);

  if (show){
    printf("R:\n");
    printM(R,N,0);
  }

  // Calcular avgR
  printf("Calculando avgR... ");
  timetick = dwalltime();
  avgR = 0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      avgR += R[i+j*N];
    }
  }
  avgR = avgR / (N*N);
  printf("t = %f\n", dwalltime() - timetick);
  if (show){
    printf("avgR: %f\n", avgR);
  }

  // Calc C
  printf("Calculando C con funcion... \n");
  timetick = dwalltime();
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      C[i+j*N] = 0;
      C[i+j*N] = T[i+j*N] + avgR * ( elememtoMult(R,A,i,j,N,0) + elememtoMult(R,B,i,j,N,0));
    }
  }
  printf(" TIEMPO = %f\n\n", dwalltime() - timetick);

  if (show){
    printf("C:\n");
    printM(C,N,0);
  }

  printf("Calculando C con multiplicacion... \n");
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

  if (show){
    printf("RA:\n");
    printM(RA,N,0);
    printf("RB:\n");
    printM(RB,N,0);
    printf("C:\n");
    printM(C,N,0);
  }

  free(A);
  free(B);
  free(C);
  free(R);
  free(M);
  free(T);
  free(RA);
  free(RB);
}