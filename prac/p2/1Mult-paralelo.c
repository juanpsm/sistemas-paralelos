#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */

void * multipThread (void * ptr);
double *A,*B,*AB;
int N, T;

// Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

// Para números aleatorios
double randFP(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

int main(int argc, char* argv[]){

    if (argc < 3){
	    printf("\nFaltan argumentos. Usar %s N T",argv[0]);
	    exit(1);
    }
    N = atoi(argv[1]);
    T = atoi(argv[2]);

    time_t t;
    srand((unsigned) time(&t));
    double timetick;
    int i,j;

    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N); 
    AB=(double*)malloc(sizeof(double)*N*N);

    printf("Incializando matrices %d x %d...\n", N, N);
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

    int id, ids[T];
    pthread_attr_t attr;
    pthread_t threads[T] ;
    pthread_attr_init(&attr);

    printf("Calculando A*B con %d threads... \n", T);
    timetick = dwalltime();

    /* Crea los hilos */
    for (id = 0; id < T; id++) {
        ids[id] =  id;
        pthread_create(&threads[id], &attr, multipThread, &ids[id]);
    }
    /* Espera a que los hilos terminen */
    for (id = 0; id < T; id++)
        pthread_join(threads[id], NULL);

    printf(" TIEMPO = %f\n", dwalltime() - timetick);    

    // Check
    int error = 0;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            if (AB[i*N+j]!=N && !error){
                // printf(" Error en  %d %d \n", i, j);
                error = 1;
            }
            // Uncomment to print AB
            if ( N <= 16 ) printf(" %f ", AB[i*N+j]);
        }
        if ( N <= 16 ) printf("\n");
    }
    printf("Resultado ");
    if (error) printf("erroneo.\n"); else printf("correcto.\n");

    free(A);
    free(B);
    free(AB);
    return 0;
}

// Función para calcular el producto de matrices partido por filas.
// Calcula las filas de A*B desde inicial a final.
void * multipThread (void * ptr) {
    int id;
    id = *((int *) ptr);
    int i,j,k,inicial,final;
    
    inicial = id*N/T;
    final = (id+1)*N/T;
    
    // printf(" Inicia Hilo %d (%d - %d)\n", id, inicial, final);

    for(i=inicial;i<final;i++){
        for(j=0;j<N;j++){
            for(k=0;k<N;k++){
                AB[i*N+j] += A[i*N+k]*B[k+j*N];
            }
            // printf("(%d) AB %d %d = %f \n", id, i, j, AB[i*N+j]);
        }
    }
    pthread_exit(0);
}