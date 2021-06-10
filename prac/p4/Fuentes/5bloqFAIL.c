#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h> /* ceil */

#define COORDINATOR 0

/*R = AB + CD + EF, donde A, B, C, D , E y F son matrices cuadradas de NxN. Ejecute para N = {512, 1024, 2048}.
*/
int main(int argc, char* argv[]){
	int i, j, k, ii, jj, kk, numProcs, rank, n, bs, stripSize, check=1;
	double *a, *b, *c, *d, *e, *f, *r, *ab, *cd, *ef, *ablk, *bblk, *cblk;
	MPI_Status status;
	double commTimes[4], maxCommTimes[4], minCommTimes[4], commTime, totalTime;

	/* Lee parámetros de la línea de comando */
	if ((argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) ) {
    printf("\nUsar: %s size bs\n  size: Dimension de la matriz\n  bs: Dimension del bloque\n", argv[0]);
		exit(1);
	}
  if ((n % bs) != 0)
  {
    printf("\nError en los parámetros. Usage: ./%s size bs (size debe ser múltiplo de bs)\n", argv[0]);
    exit(1);
  }

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  printf("(%d) hello \n", rank);

	if (n % numProcs != 0) {
		printf("El tamaño de la matriz debe ser múltiplo del numero de procesos.\n");
		exit(1);
	}

	// calcular porción de cada worker
	stripSize = n / numProcs;
  int strip = (int) ceil(n/bs / (double) numProcs);

  /* Preparar bloques */
  int start_row = rank * strip * bs;
  int end_row = (rank+1) * strip * bs;
  if (end_row > n) end_row = n;

	// Reservar memoria
	if (rank == COORDINATOR) {
		printf("  numprocs:       %d\n", numProcs);
    printf("  Bloques:        %dx%d\n", bs, bs);
    printf("  Tiras:          %d\n", n/bs);
    printf("  Strips:         %d\n", strip);
    printf("  StripSize:      %d\n", stripSize);
    /* La matriz se dividirá en tiras que dependen del ancho del bloque 
      hay un total de n/bs tiras, en cada hilo tendré:*/
    printf("  %.2f tiras x hilo\n\n", n/bs / (double) numProcs);
		a  = (double*) malloc(sizeof(double)*n*n);
		c  = (double*) malloc(sizeof(double)*n*n);
    e  = (double*) malloc(sizeof(double)*n*n);
    r  = (double*) malloc(sizeof(double)*n*n);
	}
	else  {
    a  = (double*) malloc(sizeof(double)*n*strip);
		c  = (double*) malloc(sizeof(double)*n*strip);
    e  = (double*) malloc(sizeof(double)*n*strip);
    r  = (double*) malloc(sizeof(double)*n*strip);
	}
	
	b = (double*) malloc(sizeof(double)*n*n);
  d = (double*) malloc(sizeof(double)*n*n);
  f = (double*) malloc(sizeof(double)*n*n);

  /* El resultado de los productos se necesita solo para hacer la suma, no hay que "totalizarlo",
  se puede calcular cada strip de R y luego se totaliza esta solamente. */
  ab = (double*) malloc(sizeof(double)*n*strip);
	cd = (double*) malloc(sizeof(double)*n*strip);
  ef = (double*) malloc(sizeof(double)*n*strip);

  printf("(%d) start: %d  end: %d \n", rank, start_row, end_row);

  // inicializar datos
  if (rank == COORDINATOR) {
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        a[i*n+j] = 1; // por filas
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        b[i+j*n] = 1; // por col
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        c[i*n+j] = 1; // por filas
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        d[i+j*n] = 1; // por col
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        e[i*n+j] = 1; // por filas
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        f[i+j*n] = 1; // por col
    // r en cero
    for (i=0; i<n ; i++)
      for (j=0; j<n ; j++)
        r[i*n+j] = 0; // por filas
  }
  // Estas las tienen que inicializar todos
    for (i=0; i<strip ; i++)
			for (j=0; j<n ; j++)
				ab[i*n+j] = 0;
    for (i=0; i<strip ; i++)
			for (j=0; j<n ; j++)
				cd[i*n+j] = 0;
    for (i=0; i<strip ; i++)
			for (j=0; j<n ; j++)
				ef[i*n+j] = 0;

  printf("(%d) barr \n", rank);
	MPI_Barrier(MPI_COMM_WORLD);

	commTimes[0] = MPI_Wtime();

	/* distribuir datos*/
	MPI_Scatter(a, strip*n, MPI_DOUBLE, a, strip*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Scatter(c, strip*n, MPI_DOUBLE, c, strip*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Scatter(e, strip*n, MPI_DOUBLE, e, strip*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
	MPI_Bcast(b, n*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Bcast(d, n*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
  MPI_Bcast(f, n*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

  printf("(%d) comm \n", rank);
	commTimes[1] = MPI_Wtime();

	/* computar multiplicaciones parciales */

  for (i = 0; i < strip; i+=bs)
  {
    for (j = 0; j < n; j+=bs)
    {
      cblk = &ab[i*n + j];
      for (k = 0; k < n; k+=bs)
      { 
        ablk = &a[i*n + k];
        bblk = &b[j*n + k];
        /* Iteraciones dentro de cada  bloque  */
        for (ii=0; ii < bs; ii++)
        {
          for (jj = 0; jj < bs; jj++)
          {
            for (kk = 0; kk < bs; kk++)
            {
              cblk[ii*n + jj] += ablk[ii*n + kk] * bblk[jj*n + kk];
            }
          }
        }
      }
    }
  }
  printf("(%d) termina ab \n", rank);
  for (i = 0; i < strip; i+=bs)
  { 
    for (j = 0; j < n; j+=bs)
    {
      cblk = &cd[i*n + j];
      for (k = 0; k < n; k+=bs)
      { 
        ablk = &c[i*n + k];
        bblk = &d[j*n + k];
        /* Iteraciones dentro de cada  bloque  */
        for (ii=0; ii < bs; ii++)
        {
          for (jj = 0; jj < bs; jj++)
          {
            for (kk = 0; kk < bs; kk++)
            {
              cblk[ii*n + jj] += ablk[ii*n + kk] * bblk[jj*n + kk];
            }
          }
        }
      }
    }
  }
  printf("(%d) termina cd \n", rank);
    for (i = 0; i < strip; i+=bs)
  { 
    for (j = 0; j < n; j+=bs)
    {
      cblk = &ef[i*n + j];
      for (k = 0; k < n; k+=bs)
      { 
        ablk = &e[i*n + k];
        bblk = &f[j*n + k];
        /* Iteraciones dentro de cada  bloque  */
        for (ii=0; ii < bs; ii++)
        {
          for (jj = 0; jj < bs; jj++)
          {
            for (kk = 0; kk < bs; kk++)
            {
              cblk[ii*n + jj] += ablk[ii*n + kk] * bblk[jj*n + kk];
            }
          }
        }
      }
    }
  }
  printf("(%d) termina ef \n", rank);
  for (i=0; i<strip; i++) {
		for (j=0; j<n ;j++ ) {
			r[i*n+j] = (ab[i*n+j]+cd[i*n+j]+ef[i*n+j]); 
		}
	}
  printf("(%d) termina R \n", rank);
	commTimes[2] = MPI_Wtime();

	// recolectar resultados parciales

	MPI_Gather(r, strip*n, MPI_DOUBLE, r, strip*n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

	commTimes[3] = MPI_Wtime();

	MPI_Reduce(commTimes, minCommTimes, 4, MPI_DOUBLE, MPI_MIN, COORDINATOR, MPI_COMM_WORLD);
	MPI_Reduce(commTimes, maxCommTimes, 4, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);

  printf("(%d) comm2 \n", rank);
	MPI_Finalize();

	if (rank==COORDINATOR) {

		// Check results
		for(i=0;i<n;i++)
			for(j=0;j<n;j++)
				check=check&&(r[i*n+j]==3*n);

		if(check){
			printf("Multiplicación de matrices resultado correcto\n");
		}else{
			printf("Multiplicación de matrices resultado erroneo\n");
		}
		
		totalTime = maxCommTimes[3] - minCommTimes[0];
		commTime = (maxCommTimes[1] - minCommTimes[0]) + (maxCommTimes[3] - minCommTimes[2]);

		printf("Multiplicacion de matrices (N=%d)\tTiempo total=%lf\tTiempo comunicacion=%lf\n",n,totalTime,commTime);
	}

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