#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */

void * contar (void * ptr);
int *array;
int *array2;
int N, T, X, ocurrencias;

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

    array = (int*)malloc(sizeof(int)*N);
    array2 = (int*)malloc(sizeof(int)*N);

    X = 5;
    int i, res, elem;

    printf("Incializando arreglo de %d elementos\n", N);
    for(i=0;i<N;i++){
        //array[i]=randFP(0, 10);
        elem = randFP(0, 10);
        array[i]=elem;
        array2[i]=elem;
    }

    int id, ids[T];
    pthread_attr_t attr;
    pthread_t threads[T] ;
    pthread_attr_init(&attr);
    pthread_mutex_init(&lock, NULL);

    printf("Buscando elemento %d en arreglo\n", X);
    timetick = dwalltime();

    /* Crea los hilos */
    for (id = 0; id < T; id++) {
        ids[id] =  id;
        pthread_create(&threads[id], &attr, contar, &ids[id]);
    }
    /* Espera a que los hilos terminen */
    for (id = 0; id < T; id++)
        pthread_join(threads[id], NULL);

    printf(" TIEMPO paralelo : %f\n", dwalltime() - timetick);    

    // Check seq
    res=0;
    timetick = dwalltime();
    for (i=0; i<N; i++){
        if (array2[i] == X) res++;
        // printf(" %d ", array[i]);
    }
    printf(" TIEMPO seq      : %f\n", dwalltime() - timetick);
    
    printf("Resultado ");
    if (ocurrencias != res) printf("erroneo.\n"); else printf("correcto.\n");
    printf("ocu: %d res: %d\n", ocurrencias, res);
    free(array);
    free(array2);
    return 0;
}

// Función buscar elemento en el arreglo partido por la cantidad de hilos.
void * contar (void * ptr) {
    int id;
    id = *((int *) ptr);
    int i, inicial, final;
    int ocurrencias_local;

    inicial = id*N/T;
    final = (id+1)*N/T;
    
    // printf(" Inicia Hilo %d (%d - %d)\n", id, inicial, final);

    for(i=inicial;i<final;i++){
        if (array[i] == X) {
            ocurrencias_local++;
            // printf("(%d) en %d : %d\n", id, i, ocurrencias);
        }
    }
    // sumar localmente y pedir mutex para escribir
    pthread_mutex_lock(&lock);
    ocurrencias += ocurrencias_local;
    pthread_mutex_unlock(&lock);

    pthread_exit(0);
}