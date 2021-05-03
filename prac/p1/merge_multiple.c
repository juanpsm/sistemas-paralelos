/*
Dado como entrada un Array de Listas el Algoritmo compara 3 algoritmos de mezcla distintos:

 1) mergeMultiple: mezcla todas las Listas a la vez.
 2) mergeMultipleIncremental: inicialmente mezcla las primeras dos Listas. Luego, mezcla la siguiente Lista con la Lista obtenida en la mezcla anterior.
 3) mergeMultiplePares: va mezclando de a pares de listas. Luego, sigue la misma estrategia una y otra vez con las listas resultantes.

La estructura de datos es un Array de Listas. Cada Lista está a su vez implementada como un Array.
Para simplificar el algoritmo, la longitud del Array de Listas debe ser siempre una potencia de 2.
Por defecto el algoritmo no imprime las listas. Si se quieren imprimir se debe compilar con la opcion -DPRINT. 
 
Para ejecutar:
	./programa <cantidad de listas> <tipo de inicializacion> <Longitud de lista mínima> <Longitud de lista máxima> 
	
	Parámetros:
		<cantidad de listas>: la longitud del Array de Listas. Debe ser potencia de 2 (2,4,8,16,32,64 ...etc).
		
		<tipo de inicializacion>: como serán inicializados los datos en cada lista. Hay 7 formas de inicializarlos:
		 
			SECUENCIA: en el Array de Listas, cada lista puede no tener el mismo tamaño. Cada posición se completa de la siguiente forma:
						posicion 0 tiene una lista con los números 1,2,3,4
						posicion 1 tiene una lista con los números 5,6,7,8,9
						posición 2 tiene una lista con los números 10,11
						etc. 
			PRIMOS: en el Array de Listas, cada lista puede no tener el mismo tamaño. Cada posición se completa de la siguiente forma:
						posicion 0 tiene una lista con los múltiplos del número primo 2 
						posicion 1 tiene una lista con los múltiplos del número primo 3
						posicion 2 tiene una lista con los múltiplos del número primo 5
						etc.
			PARES:  en el Array de Listas, cada lista puede no tener el mismo tamaño. Cada posición se completa de la siguiente forma:
						posicion 0 tiene una lista con los numeros 2,4,6,8 
						posicion 1 tiene una lista con los numeros 10,12,14,16,18
						posicion 2 tiene una lista con los numeros 20,22
						etc.
			IMPARES: en el Array de Listas, cada lista puede no tener el mismo tamaño. Cada posición se completa de la siguiente forma:
						posicion 0 tiene una lista con los numeros 1,3,5,7 
						posicion 1 tiene una lista con los numeros 9,11,13,15,17
						posicion 2 tiene una lista con los numeros 19,21,23
						etc.
			UNIFORME_CRECIENTE: en el Array de Listas, cada lista tiene el mismo tamaño. 
					Suponiendo un Array de 4 listas de longitud 3, cada posición se completa de la siguiente forma:
						posicion 0: 1 5 9
						posicion 1: 2 6 10
						posicion 2: 3 7 11
						posicion 3: 4 8 12
			UNIFORME_DECRECIENTE: en el Array de Listas, cada lista tiene el mismo tamaño. 
					Suponiendo un Array de 4 listas de longitud 3, cada posición se completa de la siguiente forma:
						posicion 0: 12 8 4
						posicion 1: 11 7 3
						posicion 2: 10 6 2
						posicion 3: 9 5 1
			ALEATORIA: en el Array de Listas, Cada lista puede no tener el mismo tamaño.
					Cada lista se genera aleatoriamente tanto en longitud como en contenido.
					Para generar una Lista, el primer número se genera aleatoriamente entre 1 y 10.
					El siguiente número se generar con el número anterior más un desplazamiento aleatorio entre 1 y 10.
					
		<Longitud de lista mínima>: La mínima cantidad de elementos que tendrá una lista.
		
		<Longitud de lista máxima>: La máxima cantidad de elementos que tendrá una lista.
			
*/
#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#define TRUE 1
#define FALSE 0

/*** BEGIN HACK ***/

#define CREATEENUM(name, first, ...) \
    typedef enum { first = 0, __VA_ARGS__ } name; \
    char name##_s[] = #first ", " #__VA_ARGS__;
#define TOSTR(x) #x
#define TOENUM(name, x) ((name) _toenum(name##_s, x))

long _toenum(char *enum_s, const char *x) {
    long i = 0;
    size_t len = strlen(enum_s);

    char *copy = (char*) malloc(sizeof(char) * (len + 1));
    strncpy(copy, enum_s, len);
    copy[len] = '\0';

    char *saveptr = NULL;
    char *s = strtok_r(copy, ", ", &saveptr);
    do {
        if (strcmp(s, x) == 0) {
            free(copy);
            return i;
        }
        i++;
    } while((s = strtok_r(NULL, ", ", &saveptr)) != NULL);

    free(copy);
    return -1;
}

/*** END HACK ***/

//Crea un enumerativo con nombre Tinicializadores
// Equivalente a crear:
//	enum{
//		SECUENCIA,
//		PRIMOS,
//		PARES,
//		IMPARES,
//		SIZE_INICIALIZADORES
//	}Tinicializadores;
CREATEENUM(Tinicializadores,
	SECUENCIA,
	PRIMOS,
	PARES,
	IMPARES,
	UNIFORME_CRECIENTE,
	UNIFORME_DECRECIENTE,
	ALEATORIA,
	SIZE_INICIALIZADORES
)

typedef struct lista{
 unsigned int N;
 unsigned int indiceActual;
 unsigned int* lista;
} Tlista;

/*FUNCIONES AUXILIARES*/
double dwalltime(){
 double sec;
 struct timeval tv;

	gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

void imprimirVector(unsigned int *a, unsigned int N, char* s){
unsigned int i;

 printf("%s ",s);
 for(i=0;i<N;i++){
 	printf("%d ",a[i]); 
 }
 printf("\n");

}

//Retorna si un numero es potencia de 2
unsigned int esPotenciaDe2(unsigned int v){
    return v && !(v & (v - 1));
}

//Retorna los primeros Na numeros primos en un array
unsigned int* primos(unsigned int* a, unsigned int Na){
unsigned int i=0; // contador de numeros encontrados 
unsigned int n=2; // numeros primos 
unsigned int q=1; // numeros por los q se divide n 
		
	while(i<Na){ // cantidad de numeros a imprimir
		if(n>=q){
			q++;
			
			if(((n%q)==0) && (n!=q)){ // si no es primo
				n++;
				q=1;
			}
			
			if((n%q)==0 && (n==q)){ // si es primo 
				a[i]=n;
				n++;
				i++;
				q=1;
			}	

		}else n++;	
	}
	
 return a;	
}

/*END FUNCIONES AUXILIARES*/

/*FUNCIONES DE INICIALIZACION Y VALIDACION*/
//Inicializacion donde cada posicion del arreglo lo almacena en secuencia de menor a mayor
//La longitud de cada lista se obtiene aleatoriamente entre min y max
//EJ: si la posicion 0 tiene una lista con los numeros 1,2,3,4 y la posicion 1 tiene 5,6,7,8,9 la posicion 2 10,11,... etc.
Tlista* inicializarSecuencia(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;
 unsigned int valorActual=1;

 srand(time(NULL)); 
 unsigned int r;
 
  //Inicializa Listas
 for(i=0;i<nroListas;i++){
	r = min + (rand() % (max-min+1));		
	container[i].N = r;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	for(j=0;j<container[i].N;j++){
		container[i].lista[j]=valorActual;
		valorActual++;
	}	
 }
 
 return container;
}

int validarSecuencia(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if( a[i]+1 != a[i+1] ){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }	
 return 0;
}

//Incializacion donde cada posicion del arreglo lo almacena con multiplos de un numero primo diferente
//EJ: la posicion 0 tiene una lista con los multiplos de 2, la 1 con los de 3, la 2 con los de 5 etc.
Tlista* inicializarPrimos(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;
 unsigned int* arrayPrimos = (unsigned int*)malloc(sizeof(unsigned int)*nroListas);
 arrayPrimos = primos(arrayPrimos,nroListas);

 srand(time(NULL)); 
 unsigned int r;
 
 //Inicializa Listas 
 for(i=0;i<nroListas;i++){
	r = min + (rand() % (max-min+1)); 
	container[i].N = r;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	for(j=0;j<container[i].N;j++){
		container[i].lista[j]=arrayPrimos[i]*(j+1);
	}	
 }
 
 free(arrayPrimos);
 return container;	
}

int validarPrimos(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if(a[i]>a[i+1]){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }
 return 0;
}

//Incializacion donde cada posicion del arreglo lo almacena en secuencia con numeros pares de menor a mayor
//EJ: si la posicion 0 tiene una lista con los numeros 2,4,6,8 y la posicion 1 tiene 10,12,14,16,18 la posicion 2 20,22,... etc.
Tlista* inicializarPares(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;
 unsigned int valorActual=1;
 
 srand(time(NULL)); 
 unsigned int r;
 
  //Inicializa Listas 
 for(i=0;i<nroListas;i++){
	r = min + (rand() % (max-min+1));	
	container[i].N = r;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	for(j=0;j<container[i].N;j++){
		container[i].lista[j]=2*valorActual;
		valorActual++;
	}	
 }

 return container;
}

int validarPares(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if( a[i]+2 != a[i+1] ){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }		
 return 0;
}

//Incializacion donde cada posicion del arreglo lo almacena en secuencia con numeros impares de menor a mayor
//EJ: si la posicion 0 tiene una lista con los numeros 1,3,5 y la posicion 1 tiene 7,9 la posicion 2 11,... etc.
Tlista* inicializarImpares(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;
 unsigned int valorActual=1;
 
 srand(time(NULL)); 
 unsigned int r;
 
  //Inicializa Listas
 for(i=0;i<nroListas;i++){
	r = min + (rand() % (max-min+1)); 		
	container[i].N = r;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	for(j=0;j<container[i].N;j++){
		container[i].lista[j]=2*valorActual - 1;
		valorActual++;
	}	
 }

 return container;
}

int validarImpares(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if( a[i]+2 != a[i+1] ){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }	 
 return 0;
}

//Incializacion donde todas las listas tienen el mismo tamaño, dado por max.
// Luego, el primer una lista es mayor que el primer elemento de la lista anterior y menor que el primer elemento de la lista siguiente
//EJ: si el primer elemento de la lista en la posicion 1 del arreglo es 2 entonces:
// -El primer elemento de la lista en la posicion 0 del arreglo es 1
// -Y el primer elemento de la lista en la posicion 2 del arreglo es 3
Tlista* inicializarUniformeCreciente(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;
 unsigned int valorActual=1;

  //Inicializa Listas
 for(i=0;i<nroListas;i++){
	container[i].N = max;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	for(j=0;j<container[i].N;j++)
		container[i].lista[j]=valorActual + j*nroListas;
	valorActual++;
 }
 
 return container;
}

int validarUniformeCreciente(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if( a[i]+1 != a[i+1] ){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }	
 return 0;
}

//Incializacion donde todas las listas tienen el mismo tamaño, dado por max.
// Luego, el primer una lista es menor que el primer elemento de la lista anterior y mayor que el primer elemento de la lista siguiente
//EJ: si el primer elemento de la lista en la posicion 1 del arreglo es 4 entonces:
// -El primer elemento de la lista en la posicion 0 del arreglo es 5
// -Y el primer elemento de la lista en la posicion 2 del arreglo es 3
Tlista* inicializarUniformeDecreciente(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;
 unsigned int valorActual=nroListas;

  //Inicializa Listas
 for(i=0;i<nroListas;i++){
	container[i].N = max;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	for(j=0;j<container[i].N;j++)
		container[i].lista[j]=valorActual + j*nroListas;
	valorActual--;
 }
 
 return container;
}

int validarUniformeDecreciente(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if( a[i]+1 != a[i+1] ){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }	
 return 0;
}

//Inicializacion aleatoria de los elementos de cada lista. El tamaño de cada lista es aleatorio en el rango de min y max
//El primer elemento de cada lista se genera aleatoriamente entre 1 y 10, el siguiente es el valor anterior mas un numero aleatorio entre 1 y 10
Tlista* inicializarAleatoria(Tlista* container,unsigned int nroListas, unsigned int min, unsigned int max){
 unsigned int i,j;

 srand(time(NULL)); 
 unsigned int r;
 
  //Inicializa Listas
 for(i=0;i<nroListas;i++){
	r = min + (rand() % (max-min+1));	
	container[i].N = r;		
	container[i].indiceActual = 0;
	container[i].lista = (unsigned int*)malloc(sizeof(unsigned int)*container[i].N);
	container[i].lista[0]= 1 + (rand() % 10);
	for(j=1;j<container[i].N;j++)
		container[i].lista[j]= container[i].lista[j-1] + (1 + (rand() % 10));
 }
 
 return container;
}

int validarAleatoria(unsigned int* a, unsigned int N){
unsigned int i;
 for(i=0;i<N-1;i++){
		if( a[i] > a[i+1] ){
			 printf("ERROR: creo que esto está mal mezclado en las posiciones: %d y %d\n",i,i+1);
			 return -1;
		}
 }	
 return 0;
}

/*END FUNCIONES DE INICIALIZACION Y VALIDACION*/


unsigned int* merge(unsigned int* a,unsigned int Na, unsigned int*b, unsigned int Nb, unsigned int* mezcla){
  unsigned int i, j, k;
 
  j = k = 0;
 
  for (i = 0; i < Na + Nb;) {
    if (j < Na && k < Nb) {
      if (a[j] < b[k]) {
        mezcla[i] = a[j];
        j++;
      }
      else {
        mezcla[i] = b[k];
        k++;
      }
      i++;
    }
    else if (j == Na) {
      for (; i < Na + Nb;) {
        mezcla[i] = b[k];
        k++;
        i++;
      }
    }
    else {
      for (; i < Na + Nb;) {
        mezcla[i] = a[j];
        j++;
        i++;
      }
    }
  }

  return mezcla;
}

Tlista* reiniciarIndiceActual(Tlista* container,unsigned int nroListas){
 unsigned int i;

 for(i=0;i<nroListas;i++)
	container[i].indiceActual = 0;

 return container;	
}

void imprimirMultiplesArrays(Tlista* container,unsigned int nroListas){
 unsigned int i,j;

 for(i=0;i<nroListas;i++){
	printf("Lista:%d N:%d i:%d ->",i,container[i].N, container[i].indiceActual);
	for(j=0;j<container[i].N;j++)
		printf("%d ",container[i].lista[j]);
	printf("\n");
 }
}

void freeMultiplesArrays(Tlista* container,unsigned int nroListas){
 unsigned int i;

 for(i=0;i<nroListas;i++)
	free(container[i].lista);

 free(container);
}

//Recibe una arreglo de listas y retorna la mezcla de todas las listas
//La estrategia mezcla todas las listas al mismo tiempo
unsigned int* mergeMultiple(unsigned int* mezcla,Tlista* container,unsigned int nroListas){
unsigned int min;
unsigned int hayMinimo=FALSE;
unsigned int indiceMinimo;
unsigned int i;
unsigned int indiceActualMezcla=0;
unsigned int N=0;
unsigned int it;

 //Obtiene la longitud del array mezclado
 for(i=0;i<nroListas;i++)
	N+=container[i].N;
 //Estoy seguro que la cantidad de iteraciones a realizar es la suma del tamaño de las listas
 for(it=0;it<N;it++){ 
	
	for(i=0;i<nroListas;i++){
		if(container[i].indiceActual < container[i].N){
			if ( (!hayMinimo) || (container[i].lista[container[i].indiceActual]<min) ){
				min=container[i].lista[container[i].indiceActual];
				indiceMinimo=i;
				hayMinimo=TRUE;
			}
		}
	}

	mezcla[indiceActualMezcla] = container[indiceMinimo].lista[container[indiceMinimo].indiceActual];
	container[indiceMinimo].indiceActual++;
	indiceActualMezcla++;
	hayMinimo=FALSE;
}	 
 
 return mezcla;	
}

//Recibe una arreglo de listas y retorna la mezcla de todas las listas
//La estrategia mezcla primero un par de arreglos y luego va mezclando cada arreglo con el que obtuvo antes
unsigned int* mergeMultipleIncremental(unsigned int* buffer2,Tlista* container,unsigned int nroListas){
  unsigned int i;
  unsigned int NParte=0; //Acumula el tamaño de cada parte mezclada
  unsigned int* buffer1;
  unsigned int N=0;
  unsigned int* pbaseOUT; //Puntero a esa posicion del offset
  unsigned int* pbaseIN;
  unsigned int bufferDeSalida=2; //Indica cual es el buffer de salida
  
  //Obtiene la longitud del array mezclado
  for(i=0;i<nroListas;i++)
	N+=container[i].N;
  buffer1 = (unsigned int*)malloc(sizeof(unsigned int)*N);
  pbaseOUT = buffer1;
  
  //Primera mezcla: mezcla un par de listas, las dos primeras, desde el container a un arreglo (buffer1)
  pbaseOUT = merge(container[0].lista,container[0].N, container[1].lista, container[1].N, pbaseOUT);
  NParte = container[0].N + container[1].N;
  
  //Segunda mezcla: mezcla lo que quedo en el arreglo (buffer1)
  pbaseOUT = buffer2;
  pbaseIN=buffer1;
  for(i=2;i<nroListas;i++){
	pbaseOUT = merge(pbaseIN,NParte,container[i].lista,container[i].N,pbaseOUT);
	NParte+=container[i].N;	
	
	//Intercambia el buffer de entrada por el de salida
	if (bufferDeSalida==1){
		bufferDeSalida=2;
		pbaseOUT = buffer2;
		pbaseIN = buffer1;
	}else{
		bufferDeSalida=1;
		pbaseOUT = buffer1;
		pbaseIN = buffer2;
	}
  }
 
 if(bufferDeSalida==1) {
	free(buffer1);
	return buffer2;
 }else{
	free(buffer2);
	return buffer1;
 }

}

//Recibe una arreglo de listas y retorna la mezcla de todas las listas
//La estrategia mezcla de a pares de arreglos
unsigned int* mergeMultiplePares(unsigned int* buffer2,Tlista* container,unsigned int nroListas){
  unsigned int i;
  unsigned int pares = nroListas/2; 
  unsigned int NParte[pares]; //Acumula el tamaño de cada parte mezclada
  unsigned int* buffer1;
  unsigned int N=0;
  unsigned int* pbaseOUT; //Puntero a esa posicion del offset
  unsigned int* pbaseIN1;
  unsigned int* pbaseIN2;
  unsigned int bufferDeSalida=2; //Indica cual es el buffer de salida
  unsigned int p;
  
  //Obtiene la longitud del array mezclado
  for(i=0;i<nroListas;i++)
	N+=container[i].N;
  buffer1 = (unsigned int*)malloc(sizeof(unsigned int)*N);
  pbaseOUT = buffer1;
  
  //Primera mezcla: de a pares de listas y va dejando una mezcla despues de otra en un arreglo (buffer1)
  for(i=0;i<nroListas;i+=2){
	pbaseOUT = merge(container[i].lista,container[i].N, container[i+1].lista, container[i+1].N, pbaseOUT);
	pbaseOUT+= container[i].N + container[i+1].N;
  }
  //Actualiza los tamaños de cada par
  for(i=0;i<pares;i++){
	NParte[i] = container[2*i].N + container[2*i+1].N;
  }

  //Segunda mezcla: mezcla lo que quedo en el arreglo (buffer1)
  pbaseOUT = buffer2;
  pbaseIN1=buffer1;

  for(p=pares;p>1;p=p/2){
	for(i=0;i<pares;i+=2){
		pbaseIN2 = pbaseIN1 + NParte[i];
		pbaseOUT = merge(pbaseIN1, NParte[i], pbaseIN2, NParte[i+1], pbaseOUT);
		pbaseIN1+=NParte[i] + NParte[i+1];
		pbaseOUT+= NParte[i] + NParte[i+1];
	}
 
	//Actualiza el vector de partes
	pares = pares/2;
	for(i=0;i<pares;i++){
			NParte[i] = NParte[2*i] + NParte[2*i+1];
	}
	
	//Intercambia el buffer de entrada por el de salida
	if (bufferDeSalida==1){
		bufferDeSalida=2;
		pbaseOUT = buffer2;
		pbaseIN1 = buffer1;
	}else{
		bufferDeSalida=1;
		pbaseOUT = buffer1;
		pbaseIN1 = buffer2;
	}

  }	

 if(bufferDeSalida==1) {
	free(buffer1);
	return buffer2;
 }else{
	free(buffer2);
	return buffer1;
 }
 
}

int main(int argc, char* argv[]){
Tlista* (*funcionDeInicializacion[SIZE_INICIALIZADORES]) (Tlista*, unsigned int, unsigned int, unsigned int);
int (*funcionDeValidacion[SIZE_INICIALIZADORES]) (unsigned int*, unsigned int);
Tlista* container;
unsigned int* mezcla;
unsigned int N=0;
double timetick; 
unsigned int i;
unsigned int min, max;
 
 if(argc < 5){
 	printf("Faltan argumentos, correr: ./programa <cantidad de listas> <tipo de inicializacion> <Longitud de lista mínima> <Longitud de lista máxima>\n");
	return -1;
 }	 

 unsigned int nroListas=atol(argv[1]);
 if (!esPotenciaDe2(nroListas)){
 	printf("La cantidad de arrays debe ser un número potencia de 2\n");
 	return -1;
 }

 Tinicializadores inicializador = TOENUM(Tinicializadores, argv[2]);
 if( (inicializador == -1) || (inicializador>=SIZE_INICIALIZADORES) ){
	printf("Error en la opcion de inicializacion, debe ser: SECUENCIA, PRIMOS, PARES, IMPARES, UNIFORME_CRECIENTE, UNIFORME_DECRECIENTE o ALEATORIA\n"); 
	return -1;
 }
 
 if( (inicializador==UNIFORME_CRECIENTE) || (inicializador==UNIFORME_DECRECIENTE)){
	printf("Aclaración: para la inicializacion UNIFORME_CRECIENTE o UNIFORME_DECRECIENTE solo se tomara en cuenta la longitud de lista maxima\n");	
 }
 
 min=atol(argv[3]);
 max=atol(argv[4]);
 if ( (min==0) || (max==0) || (min > max)){
	printf("Error en los parámetros de longitudes mínima y máxima: minima debe ser <> 0 y min < max\n"); 
	return -1;
 }
 
 funcionDeInicializacion[SECUENCIA]=inicializarSecuencia;
 funcionDeInicializacion[PRIMOS]=inicializarPrimos;
 funcionDeInicializacion[PARES]=inicializarPares;
 funcionDeInicializacion[IMPARES]=inicializarImpares;
 funcionDeInicializacion[UNIFORME_CRECIENTE]=inicializarUniformeCreciente;
 funcionDeInicializacion[UNIFORME_DECRECIENTE]=inicializarUniformeDecreciente;
 funcionDeInicializacion[ALEATORIA]=inicializarAleatoria;
 
 funcionDeValidacion[SECUENCIA]=validarSecuencia;
 funcionDeValidacion[PRIMOS]=validarPrimos;
 funcionDeValidacion[PARES]=validarPares;
 funcionDeValidacion[IMPARES]=validarImpares; 
 funcionDeValidacion[UNIFORME_CRECIENTE]=validarUniformeCreciente;
 funcionDeValidacion[UNIFORME_DECRECIENTE]=validarUniformeDecreciente;
 funcionDeValidacion[ALEATORIA]=validarAleatoria;
 
 //Crea e inicializa el arreglo de listas
 container = (Tlista*)malloc(sizeof(Tlista)*nroListas);
 container = funcionDeInicializacion[inicializador](container,nroListas,min,max);
 #ifdef PRINT
	imprimirMultiplesArrays(container,nroListas);
 #endif
 
 //Obtiene la longitud del arreglo mezclado y reserva memoria para este
 for(i=0;i<nroListas;i++)
	N+=container[i].N;
 mezcla = (unsigned int*)malloc(sizeof(unsigned int)*N);

 //Obtiene la mezcla de todas las listas al mismo tiempo
 timetick = dwalltime();
 mezcla = mergeMultiple(mezcla,container,nroListas);
 printf("Tiempo Mezcla multiple %f \n", dwalltime() - timetick);
 #ifdef PRINT
	imprimirVector(mezcla, N, "Multiple: ");
 #endif	
 if (funcionDeValidacion[inicializador](mezcla,N)==0) printf("Mezcla OK\n");
 
 //Reinicializa los indices para iterar las listas y pone en cero el arreglo de la mezcla
 container = reiniciarIndiceActual(container, nroListas);
 memset(mezcla, 0, sizeof(unsigned int)*N);
 
 //Obtiene la mezcla haciendolo incrementalmente
 timetick = dwalltime();
 mezcla = mergeMultipleIncremental(mezcla,container,nroListas);
 printf("Tiempo Mezcla de incremental %f \n", dwalltime() - timetick);
 #ifdef PRINT
	imprimirVector(mezcla,N, "Pares: ");
 #endif	
 if (funcionDeValidacion[inicializador](mezcla,N)==0) printf("Mezcla OK\n");

 
//Reinicializa los indices para iterar las listas y pone en cero el arreglo de la mezcla
 container = reiniciarIndiceActual(container, nroListas);
 memset(mezcla, 0, sizeof(unsigned int)*N);

 //Obtiene la mezcla haciendolo de a pares de listas
 timetick = dwalltime();
 mezcla = mergeMultiplePares(mezcla,container,nroListas);
 printf("Tiempo Mezcla de a pares %f \n", dwalltime() - timetick);
 #ifdef PRINT
  imprimirVector(mezcla,N, "Pares: ");
 #endif
 if (funcionDeValidacion[inicializador](mezcla,N)==0) printf("Mezcla OK\n");

 //Libera recursos
 free(mezcla);
 freeMultiplesArrays(container,nroListas);

 return 0;	
}

