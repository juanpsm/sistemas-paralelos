#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<math.h>     /* sin y cos */
#include<time.h>     /* srand((unsigned) time(&t)) */

/* Time calculation */
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

/* Random number generation */
double randFP(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

#define PI 3.14159265358979323846

/************** MAIN *************/
int main(int argc, char *argv[])
{
  int n;

  /* Check command line parameters */
  if (argc < 2) {
	  printf("\nFalta argumento. Usar %s n",argv[0]);
	  exit(1);
  }

  n = atoi(argv[1]);

  /* Random numbers */
  time_t t;
  srand((unsigned) time(&t));

  /* Pointers */
  double *A,*B,*C,*R1,*R2,*T,*M,*R1A,*R2B, avgR1, avgR2, sinPhi, cosPhi;

  /* Indexes */
  int i, j, k;

  /* Time measurement */
  double timetick;

  /* Getting memory */
  A   = (double*)malloc(sizeof(double)*n*n); 
  B   = (double*)malloc(sizeof(double)*n*n); 
  C   = (double*)malloc(sizeof(double)*n*n); 
  R1  = (double*)malloc(sizeof(double)*n*n); 
  R2  = (double*)malloc(sizeof(double)*n*n); 
  T   = (double*)malloc(sizeof(double)*n*n); 
  M   = (double*)malloc(sizeof(double)*n*n); 
  R1A = (double*)malloc(sizeof(double)*n*n); 
  R2B = (double*)malloc(sizeof(double)*n*n); 
  
  printf("Incializando matrices %d x %d\n", n, n);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      /* A and B by column */
      A[i+j*n]=1.0;
      B[i+j*n]=1.0;
      /* The rest by rows */
      T[i*n+j]=1.0;
      /* Zero then += */
      C[i*n+j]=0.0;
      R1A[i*n+j]=0.0;
      R2B[i*n+j]=0.0;
      /* Fill M matrix with random values beetween 0 an 2*Pi */
      M[i*n+j] = randFP(0, 2*PI);
    }
  }
  
  /* Averages initialization */
  avgR1 = 0;
  avgR2 = 0;

  printf("Calculando ... \n");

  /* Start time measurement */
  timetick = dwalltime();

  /* Calculate R1, R2 and their averages */
  // también se puede hacer cada cosa en for distintos
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
      sinPhi = sin(M[k]);
      cosPhi = cos(M[k]);
      R1[k] = (1-T[k])*(1-cosPhi)+T[k]*sinPhi;
      avgR1 += R1[k];
      R2[k] = (1-T[k])*(1-sinPhi)+T[k]*cosPhi;
      avgR2 += R2[k];
    }
  }
  avgR1 = avgR1 / (n*n);
  avgR2 = avgR2 / (n*n);
// en el paralelo si divido x filas los promedios quedan "parciales"

  /* Calculate R1 * A */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R1A[i*n+j] += R1[i*n+k]*A[k+j*n];
      }
    }
  }
  /* Calculate R2 * B */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        R2B[i*n+j] += R2[i*n+k]*B[k+j*n];
      }
    }
  }

  // aca iria una barrera que debe abrirse cuando esten listos los promedios

  /* Calculate C */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      k = i*n+j;
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

/*****************************************************************/