//Ejercicio 1
#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

int main(int argc,char*argv[]){
 double *A;
 int i,j;
 int check = 1; 
 
 int numThreads = atoi(argv[2]);
 int N=atoi(argv[1]);
 omp_set_num_threads(numThreads);

 //Aloca memoria para la matriz
  A=(double*)malloc(sizeof(double)*N*N);
 
 //Inicializa la matriz. Cada posicion debe quedar con el valor I*J
 // I => fila J=> columna. 
 #pragma omp parallel private(i,j) shared(A) 
 {
	  for(i=0;i<N;i++){
	   printf("I = %d \n", i);
	   printf("(HILO %d) \n",omp_get_thread_num());
	   #pragma omp for 
	   for(j=0;j<N;j++){

			A[i*N+j]=i*j;
			printf("(HILO %d)Ejecuto for A[%d*%d+%d] = %f \n",omp_get_thread_num(),i,N,j, A[i*N+j]);
	   }
	  }   
  }
 //Verifica el resultado
  for(i=0;i<N;i++){
   for(j=0;j<N;j++){
	check=check&&(A[i*N+j]==i*j);
   }
  }   

  if(check){
   printf("Resultado correcto\n");
  }else{
   printf("Resultado erroneo \n");
  }

 free(A);

 return(0);
}

