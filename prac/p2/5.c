#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */

void * intersect (void * ptr);
int *A, *B, *AB;
int N, T, shared_index;

pthread_mutex_t lock;

// Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

// Para números aleatorios
int randFP(int min, int max) {
  int range = (max - min);
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

    A = (int*)malloc(sizeof(int)*N);
    B = (int*)malloc(sizeof(int)*N);
    // AB en realidad no se cuanto ocupara pero a lo sumo N
    AB = (int*)malloc(sizeof(int)*N);

    int i,j,k, len;
    // Para ver el tamaño de un arreglo (ahora esta vacío)
    // igual no anda, lo chequeo luego recorriendolo... Aguante c!
    len = sizeof AB / (sizeof(int)*N);
    // printf("len(AB) = %d\n", len);
    // printf("sizeof(int) = %ld\n", sizeof(int));

    printf("Incializando arreglos de %d elementos\n", N);
    for(i=0;i<N;i++){
        A[i]=randFP(0, 10);
        B[i]=randFP(0, 10);
    }

    int id, ids[T];
    pthread_attr_t attr;
    pthread_t threads[T] ;
    pthread_attr_init(&attr);

    // este se va a usar para ir agregando los elementos comunes a A B
    shared_index=0;

    printf("Buscando elementos comunes en arreglos\n");
    timetick = dwalltime();

    /* Crea los hilos */
    for (id = 0; id < T; id++) {
        ids[id] =  id;
        pthread_create(&threads[id], &attr, intersect, &ids[id]);
    }
    /* Espera a que los hilos terminen */
    for (id = 0; id < T; id++)
        pthread_join(threads[id], NULL);

    printf(" TIEMPO paralelo : %f\n", dwalltime() - timetick);

    // Check seq
    k=0;
    timetick = dwalltime();
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            if (A[i] == B[j]){
                AB[k] = A[i];
                k++;
            }
        }
    }
    printf(" TIEMPO seq : %f\n", dwalltime() - timetick);

    // Imprimir resultados, si son chicos
    
    if (N<64){
        printf("A = [");
        for (i=0; i<N; i++) printf(" %d ", A[i]);
        printf("]\nB = [");
        for (j=0; j<N; j++) printf(" %d ", B[j]);
        printf("]\nAB = [");
        for (i=0; i<k; i++) printf(" %d ", AB[i]);
    }
    printf("]\nResultado ");
    if (shared_index != k) printf("erroneo.\n"); else printf("correcto.\n");
    printf("par: %d seq: %d \n", shared_index, k);

    free(A);
    free(B);
    free(AB);
    return 0;
}

// Función buscar elemento en el arreglo partido por la cantidad de hilos.

// Función para encontrar la interseccion
void * intersect (void * ptr) {
    int id;
    id = *((int *) ptr);
    int i,j, inicial, final;

    inicial = id*N/T;
    final = (id+1)*N/T;
    
    // printf(" Inicia Hilo %d (%d - %d)\n", id, inicial, final);

    for(i=inicial;i<final;i++){
        // cada hilo accede exclusivamente a una parte de A, pero todos leen B y escriben en AB
        for(j=0; j<N; j++){
            if (A[i] == B[j]) {
                // fijarse si el elem ya esta en AB y no agregarlo o eliminar despues los repetidos?
                // tambien se puede ahorrar comprobando si el elem actual de a ya "Lo hice"

                pthread_mutex_lock(&lock);
                // agregar_sin_repeticion(AB,shared_index,A[i])
                AB[shared_index] = A[i];
                shared_index++;
                pthread_mutex_unlock(&lock);
                // Para ver cuando encuentra un comun. Ojo con N grandes!
                // printf("(%d) A[%d]=%d  B[%d]=%d\n", id, i,A[i],j,B[j]);
            }
        }
    }
    pthread_exit(0);
}