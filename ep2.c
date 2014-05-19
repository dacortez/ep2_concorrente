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
#include <math.h>
#include <pthread.h>
#include <semaphore.h> 
#include <string.h>   
#include <time.h> 
#include <unistd.h>
#include <sys/unistd.h>
#include <gmp.h>

#define _XOPEN_SOURCE 600

/* Utilizado para inicialização das threads */
#define SHARED 1 

/* Número de repetições de um experimento */
#define TRIALS 10

/* Número de Euler */
mpf_t e;

/* Valor utilizado no critério de parada */
mpf_t eps;

/* Número de casas decimais de eps */
int decimals;

/* Critério de parada selecionado */
char stop_criteria;

/* Mode de utilização do programa */
char mode;

/* Número de threads que calculam termos */
int num_threads;

/* Vetor das threads que calculam os termos */
pthread_t* threads;

/* Barrier variable */
pthread_barrier_t barr;

/* Vetor onde os termos serão depositados pelas threads */
mpf_t *terms;

/* Número da iteração atual */
int iteraction = 0;

/* Número total de termos calculados até a iteração atual */
int total_terms = 0;

/* Sinaliza que o critério de parada foi satisfeito */
int stop = 0;

/* Ponteiro para função que faz o cálculo de um termo da série */
void (*set_term)(mpf_t term, const int k);

/**************************************************************************************************/

void print_usage();
void set_decimals(const char* eps_str);
void sequential_experiment();
void parallel_experiment();
double average(double* x, int size);
double sdv(double* x, int size);
int sequential();
void print_e(int* k);
void parallel();
void create_terms_vector();	
void create_barrier();
void create_worker_threads();
void* evaluate(void* arg);
void join_threads();
void clean_up();	

void taylor(mpf_t term, const int k);

/**************************************************************************************************/

int main(int argc, char** argv)
{
	int k;

	if (argc <= 3 || argc >= 6) {
		print_usage();
		return EXIT_FAILURE;
	}

	/* Define a função que irá calcular um termo da série */ 
	set_term = taylor;

	num_threads = atoi(argv[1]);
	if (num_threads == 0) 
		num_threads = sysconf(_SC_NPROCESSORS_ONLN);

	stop_criteria = argv[2][0];

	set_decimals(argv[3]);

	mpf_set_default_prec(4 * decimals);

	if (mpf_init_set_str(eps, argv[3], 10) == -1) {
		fprintf(stderr, "Precisão inváida.\n");
    exit(EXIT_FAILURE);	
	}	

	if (argc == 5) 
		mode = argv[4][0]; 

	mpf_init(e);

	if (mode == 'x') {
		sequential_experiment();
	}
	else if (mode == 'y') {
		parallel_experiment();
	}
	else if (mode == 's') {	
		k = sequential();
		print_e(NULL);	
		printf("termos = %d\n", k);
	}
	else { 	
		parallel();		
		print_e(NULL);
		printf("iterações = %d\n", iteraction);
	}

	mpf_clear(e);
	mpf_clear(eps);

	return EXIT_SUCCESS;
}

/**************************************************************************************************/

void print_usage()
{
	printf("Uso: ep2 <threads> <f|m> <eps> [d|s|x|y]\n");
	printf("  threads \t número de threads (0 utiliza o número de núcles)\n");
	printf("        f \t para com diferença menor do que eps\n");	
	printf("        m \t para com último termo menor do que eps\n");
	printf("      eps \t valor que define a parada\n");
	printf("        d \t informações de debug\n");
	printf("        s \t execução sequencial\n");
	printf("        x \t experimentos sequencial\n");
	printf("        y \t experimentos paralelo\n");
}

/**************************************************************************************************/

void set_decimals(const char* eps_str)
{
	char* p;
	
	if ((p = strrchr(eps_str, '.'))) {
		decimals = strlen(p + 1);
	}
	else if ((p = strrchr(eps_str, 'e'))) {
		decimals = atoi(p + 1);
		if (decimals < 0) decimals = -decimals;
	}	 
	else {
		fprintf(stderr, "Precisão inváida.\n");
  	exit(EXIT_FAILURE);
	}
}

/**************************************************************************************************/

void sequential_experiment()
{
	int i;
	clock_t begin, end;
	double time_spent[TRIALS];

	for (i = 0; i < TRIALS; ++i) {
		stop = 0;		
		mpf_set_d(e, 0.0);
		begin = clock();
		sequential();
		end = clock();
		time_spent[i] = (double)(end - begin) / CLOCKS_PER_SEC;
	}
	printf(" média = %.5fs\n", average(time_spent, TRIALS));
	printf("desvio = %.5fs\n", sdv(time_spent, TRIALS));
}

/**************************************************************************************************/

void parallel_experiment()
{
	int i;
	clock_t begin, end;
	double time_spent[TRIALS];

	for (i = 0; i < TRIALS; ++i) {
		stop = total_terms = iteraction = 0;
		mpf_set_d(e, 0.0);
		create_terms_vector();	
		create_barrier();
		begin = clock();
		create_worker_threads();
		join_threads();
		end = clock();
		clean_up();
		time_spent[i] = (double)(end - begin) / CLOCKS_PER_SEC;
	}
	printf(" média = %.5fs\n", average(time_spent, TRIALS));
	printf("desvio = %.5fs\n", sdv(time_spent, TRIALS));
}

/**************************************************************************************************/

double average(double* x, int size)
{
	int i = 0;
	double sum = 0.0;

	for (i = 0; i < TRIALS; i++)
		sum += x[i];
	return sum / size;
}

/**************************************************************************************************/

double sdv(double* x, int size)
{
	int i = 0;
	double sum = 0.0;
	double avg = average(x, size);

	for (i = 0; i < TRIALS; i++)
		sum += (x[i] - avg) * (x[i] - avg);
	return sqrt(sum / size);
}

/**************************************************************************************************/

int sequential()
{
	mpf_t term;
	int k = 0;
	
	mpf_init(term);
	while (!stop) {
		set_term(term, k++);
		mpf_add(e, e, term);
		if (mpf_cmp(term, eps) < 0) 
			stop = 1;
		if (mode != 'x')
			print_e(&k);	
	}
	mpf_clear(term); 		

	return k;
}

/**************************************************************************************************/

void print_e(int* k)
{
	if (k) 
		gmp_printf("e[%5lu] = %.*Ff\n", *k, decimals, e);
	else 
		gmp_printf("e[final] = %.*Ff\n", decimals, e);
}

/**************************************************************************************************/

void parallel()
{
	create_terms_vector();	
	create_barrier();
	create_worker_threads();
	join_threads();
	clean_up();
}

/**************************************************************************************************/

void create_terms_vector()
{
	int i;

	if (!(terms = malloc(num_threads * sizeof(mpf_t)))) {
		fprintf(stderr, "Não foi possível alocar o vetor de termos.\n");
    exit(EXIT_FAILURE);	
	}
	for (i = 0; i < num_threads; ++i)
		mpf_init(terms[i]);
}

/**************************************************************************************************/

void create_barrier()
{
	if (pthread_barrier_init(&barr, NULL, num_threads)) {
		printf("Não foi possĩvel criar a barreira.\n");
		exit(EXIT_FAILURE);
	}
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
	mpf_t sum;
	int j, i = *((int*) arg);
	
	if (arg) free(arg);

	if (i == 0) mpf_init(sum);
	
	while (!stop) {			
		set_term(terms[i], total_terms + i);
		
		if (mode == 'd');
			printf("[Thread %d chegou na barreira na iteração %d]\n", i, iteraction + 1);	
		
		/* Rendezvous 1 */
		pthread_barrier_wait(&barr);

		/* A primeira thread é responsável pela atualização de e */
		if (i == 0) {
			iteraction++; 
			total_terms += num_threads;
			if (stop_criteria == 'f') {	
				mpf_set_d(sum, 0.0);
			for (j = 0; j < num_threads; ++j)
				mpf_add(sum, sum, terms[j]);
			mpf_add(e, e, sum);
			if (mpf_cmp(sum, eps) < 0)
				stop = 1;			
			}
			else if (stop_criteria == 'm') {
				for (j = 0; j < num_threads; ++j) {
					mpf_add(e, e, terms[j]);			
					if (mpf_cmp(terms[j], eps) < 0)
						stop = 1;
				}
			}
			if (mode == 'd'); 
				print_e(&iteraction);
		}

		/* Rendezvous 2 */
		pthread_barrier_wait(&barr);
	}

	if (i == 0) mpf_clear(sum);	

	return NULL;
}

/**************************************************************************************************/

void join_threads()
{
	int i;

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
	
	if (terms) {
		for (i = 0; i < num_threads; ++i)
			mpf_clear(terms[i]);	
		free(terms);
	}

	pthread_barrier_destroy(&barr);

	if (threads) 
		free(threads); 	
}

/**************************************************************************************************/

void taylor(mpf_t term, const int k)
{
	int i;
 	mpf_t fatorial;

	mpf_init_set_ui(fatorial, 1);

	for (i = 2; i <= k; ++i)
		mpf_mul_ui(fatorial, fatorial, i);

	mpf_ui_div(term, 1, fatorial);
}

/**************************************************************************************************/


