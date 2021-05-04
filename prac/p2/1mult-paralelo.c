#include<stdio.h>
#include<stdlib.h>   /* malloc() */
#include<sys/time.h> /* gettimeofday */
#include<pthread.h>  /* hilos */

/* Multiply matrices, for pthreads */
void * multipThread (void * ptr);

/* Time calculation */
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

// Shared variables
int *A,*B,*AB;
int n, T;

/************** MAIN *************/
int main(int argc, char *argv[]){

    /* Check command line parameters */
    if (argc < 3){
	    printf("\nFaltan argumentos. Usar %s n T",argv[0]);
	    exit(1);
    }
    n = atoi(argv[1]);
    T = atoi(argv[2]);

    /* Indexes */
    int i, j;

    double timetick;

    /* Getting memory */
    A=(int*)malloc(sizeof(int)*n*n);
    B=(int*)malloc(sizeof(int)*n*n); 
    AB=(int*)malloc(sizeof(int)*n*n);

    printf("Incializando matrices %d x %d\n", n, n);
    // Fill with known pattern for later check
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
        // A por filas
        A[i*n+j] = i+j+1;
        // B por columnas
        B[i+j*n] = j+1;
        // AB en cero para despues poner +=
        AB[i*n+j] = 0;
        }
    }

    int id, ids[T];
    pthread_attr_t attr;
    pthread_t threads[T] ;
    pthread_attr_init(&attr);

    printf("Calculando A x B con %d threads... \n", T);
    
    /* Start time mesurement */
    timetick = dwalltime();

    /* Create threads */
    for (id = 0; id < T; id++) {
        ids[id] = id;
        pthread_create(&threads[id], &attr, multipThread, &ids[id]);
    }
    /* Join */
    for (id = 0; id < T; id++) {
        pthread_join(threads[id], NULL);
    }

    printf(" TIEMPO = %f\n", dwalltime() - timetick);    

    /* Check */
    int * sumfila;
    sumfila=(int*)malloc(sizeof(int)*n); 
    int print = (n<=32);
    if (print) printf("  A:\n");
    for(i=0;i<n;i++){
        sumfila[i]=0;
        for(j=0;j<n;j++){
        sumfila[i] += A[i*n+j];
        if (print) printf(" %d ", A[i*n+j]);
        }
        if (print) printf("Suma fila = %d\n", sumfila[i]);
    }
    if (print) printf("\n  B:\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
        if (print) printf(" %d ", B[i+j*n]);
        }
        if (print) printf("\n");
    }

    if (print) printf("\n  AB:\n");
    int error = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
        if (AB[i*n+j]!=sumfila[i]*(j+1) && !error){
            printf(" Error en  AB_%d_%d \n", i, j);
            error = 1;
        }
        if (print)  printf(" %d ", AB[i*n+j]);
        }
        if (print)  printf("\n");
    }
    printf("Resultado ");
    if (error) printf("erroneo.\n"); else printf("correcto.\n");

    free(A);
    free(B);
    free(AB);

    return 0;
}

/*****************************************************************/

/* Función para calcular el producto de matrices partido por filas.
   Calcula las filas de A*B desde start_row a end_row. */
void * multipThread (void * ptr) {
    int id;
    id = *((int *) ptr);
    int i,j,k,start_row,end_row;
    
    /* Each thread get some rows to operate with */
    start_row = id*n/T;
    end_row = (id+1)*n/T;
    
    // debug info
    printf("(%d) El hilo %d hará %d filas  ->",id, id, end_row-start_row);
    printf("  for i = %d .. %d (no incl)\n", start_row, end_row);

    for(i = start_row; i < end_row; i++){
        for(j=0;j<n;j++){
            for(k=0;k<n;k++){
                AB[i*n+j] += A[i*n+k]*B[k+j*n];
            }
        }
    }
    pthread_exit(0);
}
/*****************************************************************/