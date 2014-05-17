/***************************************************************************************************
 * MAC0438 Programação Concorrente
 * EP2 22/05/2014
 *
 * Daniel Augusto Cortez 2960291
 * 
 * Arquivo: ep2.c
 * Compilação: gcc -o ep2 -Wall -pedantic ep2.c -lpthread -lgmp -lm
 * Última atualização: 14/05/2014
 **************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>  
#include <time.h> 
#include <unistd.h>
#include <sys/unistd.h>
#include <gmp.h>

/*
#define _XOPEN_SOURCE 600
*/

/* Utilizado para inicialização das threads */
#define SHARED 1 

/* Valor utilizado no critério de parada */
long double eps;

/* Critério de parada selecionado */
char stop_criteria;

/* Opção de debug para o programa */
int debug = 0;

/* Número de threads que calculam termos */
int num_threads;

/* Vetor das threads que calculam os termos */
pthread_t* threads;

/* Thread coordenadora para a barreira de sincronização */
pthread_t coordinator;

/* Vetores de semaforos utilizados para implementar a barreira de sincronização */
sem_t* arrive;
sem_t* go_on;

/* Vetor onde os termos serão depositados pelas threads */
double *terms;

/* Número da iteração atual */
long iteraction = 0;

/* Número total de termos calculados até a iteração atual */
long total_terms = 0;

/* Sinaliza que o critério de parada foi satisfeito */
int stop = 0;

/* Número de Euler */
long double e = 0.0;

/* Ponteiro para função que faz o cálculo de um termo da série */
long double (*get_term)(long k);

/* Barreira de sincronização pronta */
/* pthread_barrier_t barrier; */

/* Mutex compartilhado */
/* pthread_mutex_t mutex; */

/**************************************************************************************************/

void print_usage();

void sequential();

void parallel();
void create_terms_vector();	
void create_barrier_semaphores();
void create_coordinator_thread();
void* barrier_sync(void* arg);
void create_worker_threads();
void* evaluate(void* arg);
void join_threads();
void clean_up();	

long double taylor(long k);

/**************************************************************************************************/

int main(int argc, char** argv)
{
	if (argc <= 3 || argc >= 6) {
		print_usage();
		return EXIT_FAILURE;
	}

	/* Está linha define a função que irá calcular um termo da série */ 
	get_term = taylor;

	num_threads = atoi(argv[1]);
	if (num_threads == 0) 
		num_threads = sysconf(_SC_NPROCESSORS_ONLN);

	stop_criteria = argv[2][0];

	eps = atof(argv[3]);

	if (argc == 5 && argv[4][0] == 'd') 
		debug = 1; 

	printf("num_threads = %d\n", num_threads);
	printf("       stop = %c\n", stop_criteria);
	printf("        eps = %Lf\n", eps);
	printf("      debug = %d\n", debug);

	if (argc == 5 && argv[4][0] == 's') {
		printf("[sequencial]\n");	
		sequential();
	}
	else { 
		printf("[paralelo]\n");		
		parallel();
	}

	printf("e = %20.19Lf\n", e);

	return EXIT_SUCCESS;
}

/**************************************************************************************************/

void print_usage()
{
	printf("Uso: ep2 <num_threads> <f|m> <eps> [d|s]\n");
	printf("  num_threads \t número de threads (0 indica número de threads = número de núcles)\n");
	printf("            f \t para com diferença menor do que eps\n");	
	printf("            m \t para com último termo menor do que eps\n");
	printf("          eps \t valor que define a parada\n");
	printf("            d \t informações de debug\n");
	printf("            s \t execução sequencial\n");
}

/**************************************************************************************************/

long double taylor(long k)
{
	long i, fat = 1;

	for (i = 2; i <= k; ++i)
		fat *= i;

	return (long double) 1.0 / fat;
}

/**************************************************************************************************/

void sequential()
{
	long double term;
	long k = 0;
	
	while (1) {
			e += (term = get_term(k++));
			if (term < eps) return;
	}
}

/**************************************************************************************************/

void parallel()
{
	create_terms_vector();	
	create_barrier_semaphores();
	create_coordinator_thread();
	create_worker_threads();
	join_threads();
	clean_up();	
}

/**************************************************************************************************/

void create_terms_vector()
{
	if (!(terms = malloc(num_threads * sizeof(double)))) {
		fprintf(stderr, "Não foi possível alocar o vetor de termos.\n");
    exit(EXIT_FAILURE);	
	}
}

/**************************************************************************************************/

void create_barrier_semaphores()
{
	int i;

	if (!(arrive = malloc(num_threads * sizeof(sem_t)))) {
		fprintf(stderr, "Não foi possível alocar semáforo.\n");
    exit(EXIT_FAILURE);	
	}
	if (!(go_on = malloc(num_threads * sizeof(sem_t)))) {
		fprintf(stderr, "Não foi possível alocar semáforo.\n");
    exit(EXIT_FAILURE);	
	}
	for (i = 0; i < num_threads; ++i) {
		if (sem_init(&arrive[i], SHARED, 0)) {
			fprintf(stderr, "Não foi possível inicializar semáforo.\n");
    	exit(EXIT_FAILURE);	
		}
		if (sem_init(&go_on[i], SHARED, 0)) {
			fprintf(stderr, "Não foi possível inicializar semáforo.\n");
    	exit(EXIT_FAILURE);
		}
	}
}

/**************************************************************************************************/

void create_coordinator_thread()
{
	pthread_attr_t attr;
	
	pthread_attr_init(&attr);
	pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

	if (pthread_create(&coordinator, &attr, barrier_sync, NULL)) {
		fprintf(stderr, "Não foi possível criar thread coordenadora.\n");
		exit(EXIT_FAILURE); 
	}
}

/**************************************************************************************************/

void* barrier_sync(void* arg)
{
	int i;
	long double previous_e;

	while (!stop) {
		/* Espera todos terminarem a iteração */
		for (i = 0; i < num_threads; ++i)
			sem_wait(&arrive[i]);				

		/* Rendezvous */
		total_terms = (++iteraction) * num_threads;
		if (stop_criteria == 'f') {	
			previous_e = e;
			for (i = 0; i < num_threads; i++)
				e += terms[i];
			if (e - previous_e < eps) 
				stop = 1;			
		}
		else if (stop_criteria == 'm') {
			for (i = 0; i < num_threads; i++) {
				e += terms[i];			
				if (terms[i] < eps)
					stop = 1;
			}
		}

		/* Libera para continuar a próxima iteração */
		for (i = 0; i < num_threads; ++i)
			sem_post(&go_on[i]);	
	}

	return NULL;
}

/**************************************************************************************************/

void create_worker_threads()
{
	pthread_attr_t attr;
	int* ip;	
	int i;
	
	pthread_attr_init(&attr);
	pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
	if (!(threads = malloc(num_threads * sizeof(pthread_t)))) {
		fprintf(stderr, "Não foi possível alocar threads.\n");
    exit(EXIT_FAILURE);
	}
	for (i = 0; i < num_threads; ++i) {
		ip = malloc(sizeof(int)); *ip = i; 
		if (pthread_create(&threads[i], &attr, evaluate, (void*) ip)) {
			fprintf(stderr, "Não foi possível criar thread %d.\n", i);
      exit(EXIT_FAILURE); 
		}
	}
}

/**************************************************************************************************/

void* evaluate(void* arg)
{
	int i = *((int*) arg);
	
	if (arg) free(arg);

	while (!stop) {	
		terms[i] = get_term(total_terms + i);

		printf("[thread %d calculou na iteracao %ld]\n", i, iteraction);	

		sem_post(&arrive[i]);
		sem_wait(&go_on[i]);
	}

	return NULL;
}

/**************************************************************************************************/

void join_threads()
{
	int i;

	if (pthread_join(coordinator, NULL)) {
		fprintf(stderr, "Não é possível juntar a thread coordenadora.\n");
	  exit(EXIT_FAILURE);
	}
	for (i = 0; i < num_threads; ++i) {
		if (pthread_join(threads[i], NULL)) {
	  	fprintf(stderr, "Não é possível juntar a thread %d.\n", i);
	  	exit(EXIT_FAILURE);
		}
	}
}
/**************************************************************************************************/

void clean_up()
{
	int i;
	
	if (terms) 
		free(terms);

	if (arrive) {
		for (i = 0; i < num_threads; ++i)
			sem_destroy(&arrive[i]);
		free(arrive);
	}
	
	if (go_on) {
		for (i = 0; i < num_threads; ++i)
			sem_destroy(&go_on[i]);
		free(go_on);
	}
}

/**************************************************************************************************/


/**************************************************************************************************

if (pthread_barrier_init(&barrier, NULL, THREADS)) {
	fprintf(stderr, "Não é possível criar a barreira.\n");
	return EXIT_FAILURE;
}
pthread_barrier_destroy(&barrier);

if (pthread_mutex_init(&mutex, NULL)) {
	fprintf(stderr, "Não é possível inicializar o mutex.\n");
	return EXIT_FAILURE;
}
pthread_mutex_destroy(&mutex);
		
pthread_mutex_lock(&mutex); 
pthread_mutex_unlock(&mutex); 
rc = pthread_barrier_wait(&barrier);
if (rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
	fprintf(stderr, "Could not wait on barrier\n");
	exit(EXIT_FAILURE);
}

***************************************************************************************************/
