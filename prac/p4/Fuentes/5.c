#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define COORDINATOR 0

/*R = AB + CD + EF, donde A, B, C, D , E y F son matrices cuadradas de NxN. Ejecute para N = {512, 1024, 2048}.
*/
int main(int argc, char* argv[]){
	int i, j, k, numProcs, rank, n, stripSize, check=1;
	double *a, *b, *c, *d, *e, *f, *r, *ab, *cd, *ef;
	MPI_Status status;
	double commTimes[4], maxCommTimes[4], minCommTimes[4], commTime, totalTime;

	/* Lee parámetros de la línea de comando */
	if ((argc != 2) || ((n = atoi(argv[1])) <= 0) ) {
	    printf("\nUsar: %s size \n  size: Dimension de la matriz y el vector\n", argv[0]);
		exit(1);
	}

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (n % numProcs != 0) {
		printf("El tamaño de la matriz debe ser multiplo del numero de procesos.\n");
		exit(1);
	}

	// calcular porcion de cada worker
	stripSize = n / numProcs;

	// Reservar memoria
	if (rank == COORDINATOR) {
		printf("numprocs: %d\n", numProcs);
		a  = (double*) malloc(sizeof(double)*n*n);
		c  = (double*) malloc(sizeof(double)*n*n);
    e  = (double*) malloc(sizeof(double)*n*n);
    r  = (double*) malloc(sizeof(double)*n*n);
	}
	else  {
    a  = (double*) malloc(sizeof(double)*n*stripSize);
		c  = (double*) malloc(sizeof(double)*n*stripSize);
    e  = (double*) malloc(sizeof(double)*n*stripSize);
    r  = (double*) malloc(sizeof(double)*n*stripSize);
	}
	
	b = (double*) malloc(sizeof(double)*n*n);
  d = (double*) malloc(sizeof(double)*n*n);
  f = (double*) malloc(sizeof(double)*n*n);

  /* El resultado de los productos se necesita solo para hacer la suma, no hay que "totalizarlo",
  se puede calcular cada strip de R y luego se totaliza esta solamente. */
  ab = (double*) malloc(sizeof(double)*n*stripSize);
	cd = (double*) malloc(sizeof(double)*n*stripSize);
  ef = (double*) malloc(sizeof(double)*n*stripSize);


	// inicializar datos
	if (rank == COORDINATOR) {
		for (i=0; i<n ; i++)
			for (j=0; j<n ; j++)
				a[i*n+j] = 1;
		for (i=0; i<n ; i++)
			for (j=0; j<n ; j++)
				b[i*n+j] = 1;
    for (i=0; i<n ; i++)
			for (j=0; j<n ; j++)
				c[i*n+j] = 1;
    for (i=0; i<n ; i++)
			for (j=0; j<n ; j++)
				d[i*n+j] = 1;
    for (i=0; i<n ; i++)
			for (j=0; j<n ; j++)
				e[i*n+j] = 1;
    for (i=0; i<n ; i++)
			for (j=0; j<n ; j++)
				f[i*n+j] = 1;

	}

	MPI_Barrier(MPI_COMM_WORLD);

	commTimes[0] = MPI_Wtime();

	/* distribuir datos*/
	MPI_Scatter(a, stripSize*n, MPI_DOUBLE, a, stripSize*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Scatter(c, stripSize*n, MPI_DOUBLE, c, stripSize*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Scatter(e, stripSize*n, MPI_DOUBLE, e, stripSize*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
	MPI_Bcast(b, n*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Bcast(d, n*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Bcast(f, n*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

	commTimes[1] = MPI_Wtime();

	/* computar multiplicaciones parciales */
	for (i=0; i<stripSize; i++) {
		for (j=0; j<n ;j++ ) {
			ab[i*n+j]=0;
			for (k=0; k<n ;k++ ) { 
				ab[i*n+j] += (a[i*n+k]*b[j*n+k]); 
			}
		}
	}

  for (i=0; i<stripSize; i++) {
		for (j=0; j<n ;j++ ) {
			cd[i*n+j]=0;
			for (k=0; k<n ;k++ ) { 
				cd[i*n+j] += (c[i*n+k]*d[j*n+k]); 
			}
		}
	}

  for (i=0; i<stripSize; i++) {
		for (j=0; j<n ;j++ ) {
			ef[i*n+j]=0;
			for (k=0; k<n ;k++ ) { 
				ef[i*n+j] += (e[i*n+k]*f[j*n+k]); 
			}
		}
	}

  for (i=0; i<stripSize; i++) {
		for (j=0; j<n ;j++ ) {
			r[i*n+j] = (ab[i*n+j]+cd[i*n+j]+ef[i*n+j]); 
		}
	}

	commTimes[2] = MPI_Wtime();

	// recolectar resultados parciales

	MPI_Gather(r, stripSize*n, MPI_DOUBLE, r, stripSize*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

	commTimes[3] = MPI_Wtime();

	MPI_Reduce(commTimes, minCommTimes, 4, MPI_DOUBLE, MPI_MIN, COORDINATOR, MPI_COMM_WORLD);
	MPI_Reduce(commTimes, maxCommTimes, 4, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);

	MPI_Finalize();

	if (rank==COORDINATOR) {

		// Check results
		for(i=0;i<n;i++)
			for(j=0;j<n;j++)
				check=check&&(c[i*n+j]==n);

		if(check){
			printf("Multiplicación de matrices resultado correcto\n");
		}else{
			printf("Multiplicación de matrices resultado erroneo\n");
		}
		
		totalTime = maxCommTimes[3] - minCommTimes[0];
		commTime = (maxCommTimes[1] - minCommTimes[0]) + (maxCommTimes[3] - minCommTimes[2]);

		printf("Multiplicacion de matrices (N=%d)\tTiempo total=%lf\tTiempo comunicacion=%lf\n",n,totalTime,commTime);
	}

  int error = 0;
	if (rank==COORDINATOR){
    for (i=0; i<n; i++) {
      for (j=0; j<n ;j++ ) {
        if (r[i*n+j]!=3*n) {
          error = 1;
        }
        // printf(" %f ", a[i*n+k]); 
      }
  }
  printf("error = %d", error);

  free(a);
  free(b);
  free(c);
  free(d);
  free(e);
  free(f);
  free(r);
  free(ab);
  free(cd);
  free(ef);

  return 0;
}