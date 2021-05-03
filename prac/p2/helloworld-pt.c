#include <stdio.h>
#include <pthread.h>
#define NUM_THREADS 10

void * hello_world (void * ptr);
int main(){
    int i, ids[NUM_THREADS];
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS] ;
    pthread_attr_init(&attr);
    /* Crea los hilos */
    for (i = 0; i < NUM_THREADS; i++) {
        ids[i] =  i;
        pthread_create(&threads[i], &attr, hello_world, &ids[i]);
    }
    /* Espera a que los hilos terminen */
    for (i = 0; i < NUM_THREADS; i++)
        pthread_join(threads[i], NULL);
    return 0;
}

void * hello_world (void * ptr) {
    int id;
    id = *((int *) ptr);
    printf("¡Hola Mundo! Soy el hilo %d \n",id);

    pthread_exit(0);
}