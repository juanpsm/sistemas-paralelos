#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */
#include <semaphore.h>

void * maxmin (void * ptr);
int *array;
int N, T, MAX, MIN;

sem_t mutex_max;
sem_t mutex_min;

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

    int i, res;

    printf("Incializando arreglo de %d elementos\n", N);
    for(i=0;i<N;i++){
        array[i]=randFP(-5770, 10015);
    }

    int id, ids[T];
    pthread_attr_t attr;
    pthread_t threads[T] ;
    pthread_attr_init(&attr);

    sem_init(&mutex_max, 0, 1);
    sem_init(&mutex_min, 0, 1);

    printf("Buscando elemento MAX y MIN en arreglo\n");
    timetick = dwalltime();

    /* Crea los hilos */
    for (id = 0; id < T; id++) {
        ids[id] =  id;
        pthread_create(&threads[id], &attr, maxmin, &ids[id]);
    }
    /* Espera a que los hilos terminen */
    for (id = 0; id < T; id++)
        pthread_join(threads[id], NULL);

    printf(" TIEMPO paralelo : %f\n", dwalltime() - timetick);    

    // Check seq
    int max=0, min=0;
    timetick = dwalltime();
    for (i=0; i<N; i++){
        if (array[i] < min) min=array[i];
        if (array[i] > max) max=array[i];
    }
    printf(" TIEMPO seq      : %f\n", dwalltime() - timetick);
    
    printf("Resultado ");
    if (MAX != max) printf("erroneo.\n"); else printf("correcto.\n");
    printf("MAX: %d max: %d\n", MAX, max);
    printf("MIN: %d min: %d\n", MIN, min);

    free(array);
    return 0;
}

// Función buscar elemento en el arreglo partido por la cantidad de hilos.
void * maxmin (void * ptr) {
    int id;
    id = *((int *) ptr);
    int i, inicial, final;

    inicial = id*N/T;
    final = (id+1)*N/T;
    
    int min_loc, max_loc;

    // printf(" Inicia Hilo %d (%d - %d)\n", id, inicial, final);

    for(i=inicial;i<final;i++){
        if (array[i] > max_loc) max_loc = array[i];
        if (array[i] < min_loc) min_loc = array[i];
    }

    sem_wait(&mutex_max);
    if (max_loc > MAX) MAX = max_loc;
    sem_post(&mutex_max);
    sem_wait(&mutex_min);
    if (min_loc < MIN) MIN = min_loc;
    sem_post(&mutex_min);

    pthread_exit(0);
}