#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h> /* gettimeofday */
#include<time.h> /* srand((unsigned) time(&t)) */
#include<pthread.h> /* hilos */
#include<math.h>

void * calc_promedio (void * ptr);
double *array;
int N, T;

double promedio;

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

    array = (double*)malloc(sizeof(double)*N);

    int i;

    printf("Incializando arreglo de %d elementos\n", N);
    for(i=0;i<N;i++){
        array[i]=randFP(0.0, 100.0);
        //array[i]= 4.8759;
    }

    int id, ids[T];
    pthread_attr_t attr;
    pthread_t threads[T] ;
    pthread_attr_init(&attr);

    printf("Calculando promedio\n");
    timetick = dwalltime();

    /* Crea los hilos */
    for (id = 0; id < T; id++) {
        ids[id] =  id;
        pthread_create(&threads[id], &attr, calc_promedio, &ids[id]);
    }
    /* Espera a que los hilos terminen */
    for (id = 0; id < T; id++)
        pthread_join(threads[id], NULL);

    printf(" TIEMPO paralelo : %f\n", dwalltime() - timetick);    

    // Check seq
    double res = 0.0;
    timetick = dwalltime();
    for (i=0; i<N; i++){
        res += array[i];
    }
    res = res / N;
    printf(" TIEMPO seq      : %f\n", dwalltime() - timetick);
    
    // Si se quiere ver el arreglo:
    // printf(" array = [");
    // for (i=0; i<N; i++){
    //     printf(" %f ", array[double]);
    // }
    // printf("]\n");


    printf("Resultado ");
    // defino un margen de error, porque es dificil comparar doubles.
    if (fabs(promedio - res) > 0.000001) printf("erroneo.\n"); else printf("correcto.\n");
    printf("par: %f seq: %f\n", promedio, res);

    free(array);
    return 0;
}

// Función buscar elemento en el arreglo partido por la cantidad de hilos.
void * calc_promedio (void * ptr) {
    int id;
    id = *((int *) ptr);
    int i, inicial, final;

    inicial = id*N/T;
    final = (id+1)*N/T;
    
    // printf(" Inicia Hilo %d (%d - %d)\n", id, inicial, final);
    double parcial = 0.0;

    for(i=inicial;i<final;i++){
        // printf("id %d inic %d fin %d i %d\n",id,inicial,final,i);
        parcial += array[i];
    }
    parcial = parcial / N;
    
    pthread_mutex_lock(&lock);
    promedio += parcial;
    pthread_mutex_unlock(&lock);

    pthread_exit(0);
}