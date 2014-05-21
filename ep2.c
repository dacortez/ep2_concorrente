/***************************************************************************************************
 * MAC0438 Programação Concorrente
 * EP2 22/05/2014
 *
 * Daniel Augusto Cortez 2960291
 * 
 * Arquivo: ep2.c
 * Compilação: gcc -o ep2 -Wall -pedantic -O3 ep2.c -lpthread -lgmp -lm
 * Última atualização: 20/05/2014
 **************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h> 
#include <string.h>
#include <sys/time.h>   
#include <time.h> 
#include <unistd.h>
#include <sys/unistd.h>
#include <gmp.h>

/* Utilizado para inicialização da barreira */
#define _XOPEN_SOURCE 600

/* Utilizado para inicialização das threads */
#define SHARED 1 

/* Número de repetições de um experimento */
#define TRIALS 10

/* Número de Euler */
mpf_t e;

/* Precisão utilizada no critério de parada */
mpf_t eps;

/* Número de casas decimais de eps */
int decimal_places;

/* Critério de parada selecionado */
char stop_criteria;

/* Mode de utilização do programa */
char mode;

/* Número de threads trabalhadoras (uma para cada termo) */
int num_threads;

/* Vetor das threads que calculam os termos */
pthread_t* threads;

/* Barreira de sincronizção */
pthread_barrier_t barr;

/* Vetor onde os termos serão depositados pelas threads */
mpf_t *terms;

/* Número da iteração atual */
int iteraction = 0;

/* Número total de termos calculados até a iteração atual */
int total_terms = 0;

/* Sinaliza que o critério de parada foi satisfeito */
int stop = 0;

/* Parâmetro para expansão de Taylor somando de p em p termos */ 
int p;

/* Armazena maior fatoria calculado até a presente iteração */
mpf_t max_fat; 

/* Maior fatorial calculado pela última thread */
mpf_t last_fat;

/* Ponteiro para função que verifica o critério de parada e atualiza e */
void (*test_and_update_e)();

/* Ponteiro para função que faz o cálculo de um termo da série */
void (*set_term)(mpf_t term, int k, int id);

/**************************************************************************************************/

void print_usage();
void set_decimal_places(char* eps_str);
void do_sequential_experiment();
void do_parallel_experiment();
double my_clock();
double average(double* x, int size);
double sdv(double* x, int size);
void sequential();
void print_e(int k);
void print_final();
void parallel();
void create_terms();		
void create_barrier();
void create_threads();
void* work(void* arg);
void testf_and_update_e();
void testm_and_update_e();
void join_threads();
void destroy_terms();	
void destroy_barrier();
void destroy_threads();	
	
/* Três funções para o cálculo do termo */
void taylor_1(mpf_t term, int k, int id);
void taylor_3(mpf_t term, int k, int id);
void taylor_p(mpf_t term, int k, int id);

/**************************************************************************************************/

int main(int argc, char** argv)
{
	if (argc <= 3 || argc >= 7) {
		print_usage();
		return EXIT_FAILURE;
	}

	num_threads = atoi(argv[1]);
	if (num_threads == 0) 
		num_threads = sysconf(_SC_NPROCESSORS_ONLN);

	test_and_update_e = (argv[2][0] == 'f') ? testf_and_update_e : testm_and_update_e;

	set_decimal_places(argv[3]);

	/* precison = - log_2 (eps) => precision = decimal_places * log_2 (10) */ 
	mpf_set_default_prec(4 * decimal_places);

	if (mpf_init_set_str(eps, argv[3], 10) == -1) {
		fprintf(stderr, "Precisão inváida.\n");
    exit(EXIT_FAILURE);	
	}	

	if (argc >= 5) 
		mode = argv[4][0]; 

	/* Define a função que irá calcular um termo da série */ 
	if (argc == 6) {
		p = atoi(argv[5]);	
		set_term = taylor_p;
	}
	else {
		p = decimal_places / 100; 
		set_term = (p <= 3) ? taylor_3 : taylor_p;  
	}	

	mpf_init(e);

	if (mode == 'x') {
		do_sequential_experiment();
	}
	else if (mode == 'y') {
		do_parallel_experiment();
	}
	else if (mode == 's') {
		sequential();
		print_final();	
		printf("termos = %d\n", total_terms);
	}
	else { 	
		parallel();		
		print_final();
		printf("iterações = %d\n", iteraction);
	}
	
	mpf_clear(e);
	mpf_clear(eps);

	return EXIT_SUCCESS;
}

/**************************************************************************************************/

void print_usage()
{
	printf("Uso: ep2 <threads> <f|m> <eps> [d|s|x|y] [p]\n");
	printf("  threads \t número de threads (0 utiliza o número de núcles)\n");
	printf("        f \t para com diferença menor do que eps\n");	
	printf("        m \t para com termo menor do que eps\n");
	printf("      eps \t valor que define a precisão (ex: 1e-100)\n");
	printf("        d \t informações de debug\n");
	printf("        s \t execução sequencial\n");
	printf("        x \t experimentos sequencial\n");
	printf("        y \t experimentos paralelo\n");
	printf("        p \t parâmetro de compressão da série de Taylor\n");
}

/**************************************************************************************************/

void set_decimal_places(char* eps_str)
{
	char* p;
	
	if ((p = strrchr(eps_str, '.')))
		decimal_places = strlen(p + 1);
	else if ((p = strrchr(eps_str, 'e'))) {
		decimal_places = atoi(p + 1);
		if (decimal_places < 0) 
			decimal_places = -decimal_places;
	}	 
	else {
		fprintf(stderr, "Precisão inváida.\n");
  	exit(EXIT_FAILURE);
	}
}

/**************************************************************************************************/

void do_sequential_experiment()
{
	int i;
	double begin, end, time_spent[TRIALS];

	for (i = 0; i < TRIALS; ++i) {
		stop = 0;		
		begin = my_clock();
		sequential();
		end = my_clock(); 
		time_spent[i] = end - begin;
		/* print_final(); */
	}
	printf(" média[s] = %.6fs\n", average(time_spent, TRIALS));
	printf("desvio[s] = %.6fs\n", sdv(time_spent, TRIALS));
}

/**************************************************************************************************/

void do_parallel_experiment()
{
	int i;
	double begin, end, time_spent[TRIALS];

	for (i = 0; i < TRIALS; ++i) {
		stop = total_terms = iteraction = 0;
		begin = my_clock();
		parallel();		
		end = my_clock();
		time_spent[i] = end - begin;
		/* print_final(); */ 
	}
	printf(" média[p] = %.6fs\n", average(time_spent, TRIALS));
	printf("desvio[p] = %.6fs\n", sdv(time_spent, TRIALS));
}

/**************************************************************************************************/

inline double my_clock() 
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (1.0e-6 * t.tv_usec + t.tv_sec);
}

/**************************************************************************************************/

double average(double* x, int size)
{
	int i;
	double sum = 0.0;

	for (i = 0; i < size; i++)
		sum += x[i];
	return sum / size;
}

/**************************************************************************************************/

double sdv(double* x, int size)
{
	int i;
	double sum = 0.0;
	double avg = average(x, size);

	for (i = 0; i < size; i++)
		sum += (x[i] - avg) * (x[i] - avg);
	return sqrt(sum / size);
}

/**************************************************************************************************/

void sequential()
{
	mpf_t term;
	
	num_threads = 1; /* deve ser 1 para poder utilizar max_fat em taylor_p() */
	total_terms = 0;	

	mpf_init(term);
	mpf_init(last_fat);
	mpf_init_set_ui(max_fat, 1);
	
	mpf_set_d(e, 0.0);
	while (!stop) {
		set_term(term, total_terms, 0); /* id deve ser 0 para poder utilizar last_fat em taylor_p() */
		mpf_add(e, e, term);
		mpf_set(max_fat, last_fat);
		if (mpf_cmp(term, eps) < 0) 
			stop = 1;
		if (mode == 's')
			print_e(total_terms);
		total_terms++;	
	}

	mpf_clear(term);
	mpf_clear(last_fat);
	mpf_clear(max_fat);	
}

/**************************************************************************************************/

void print_e(int k)
{
	gmp_printf("e[%5lu] = %.*Ff\n", k, decimal_places, e);
}

/**************************************************************************************************/

void print_final()
{
	gmp_printf("e[final] = %.*Ff\n", decimal_places, e);
}

/**************************************************************************************************/

void parallel()
{
	mpf_init(last_fat);
	mpf_init_set_ui(max_fat, 1);
	
	mpf_set_d(e, 0.0);

	create_terms();	
	create_barrier();
	create_threads();
	join_threads();
	destroy_terms();
	destroy_barrier();
	destroy_threads();

	mpf_clear(last_fat);	
	mpf_clear(max_fat);	
}

/**************************************************************************************************/

void create_terms()
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

void create_threads()
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
		if (pthread_create(&threads[i], &attr, work, (void*) ip)) {
			fprintf(stderr, "Não foi possível criar thread %d.\n", i);
      exit(EXIT_FAILURE); 
		}
	}
}

/**************************************************************************************************/

void* work(void* arg)
{
	int i = *((int*) arg);
	
	if (arg) 
		free(arg);

	while (!stop) {			
		set_term(terms[i], total_terms + i, i);
		if (mode == 'd')
			printf("[Thread %d chegou na barreira na iteração %d]\n", i, iteraction + 1);	
		
		/* Barreira 1 */
		pthread_barrier_wait(&barr);

		/* A primeira thread é responsável pela atualização da iteração */
		if (i == 0) {
			iteraction++; 
			total_terms += num_threads;
			mpf_set(max_fat, last_fat);
			test_and_update_e();
			if (mode == 'd')
				print_e(iteraction);
		}

		/* Barreira 2 */
		pthread_barrier_wait(&barr);
	}

	return NULL;
}

/**************************************************************************************************/

void testf_and_update_e()
{
	int i;
	mpf_t sum;

	mpf_init(sum);
	
	for (i = 0; i < num_threads; ++i)
		mpf_add(sum, sum, terms[i]);
	if (mpf_cmp(sum, eps) < 0)
		stop = 1;
	mpf_add(e, e, sum);
	
	mpf_clear(sum);
}

/**************************************************************************************************/

void testm_and_update_e()
{
	int i;

	for (i = 0; i < num_threads; ++i) {
		mpf_add(e, e, terms[i]);			
		if (mpf_cmp(terms[i], eps) < 0)
			stop = 1;
	}
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

void destroy_terms()
{
	int i;
	
	if (terms) {
		for (i = 0; i < num_threads; ++i)
			mpf_clear(terms[i]);	
		free(terms);
	}
}

/**************************************************************************************************/

void destroy_barrier()
{
	pthread_barrier_destroy(&barr);
}

/**************************************************************************************************/

void destroy_threads()
{
	if (threads) 
		free(threads); 	
}

/**************************************************************************************************/

void taylor_1(mpf_t term, int k, int id)
{
	int i, factors;
	mpf_t aux_fat;

  /* k! = k (k-1) ... max_fat */
	factors = total_terms == 0 ? k : 1 + k - total_terms;	
	mpf_init_set(aux_fat, max_fat);	
	for (i = 0; i < factors; i++)
		mpf_mul_ui(aux_fat, aux_fat, k - i);	

	/* 1 / k! */
	mpf_ui_div(term, 1, aux_fat);

	/* última thread cálcula o maior fatorial (último termo) */
	if (id == num_threads - 1)
		mpf_set(last_fat, aux_fat);

	mpf_clear(aux_fat);
}

/**************************************************************************************************/

void taylor_3(mpf_t term, int k, int id)
{
	int i, factors, three_k = 3 * k;
	mpf_t aux_fat;
	
	/* (3k)! = 3k (3k-1) ... max_fat */
	factors = total_terms == 0 ? three_k : 3 + three_k - 3 * total_terms;	
	mpf_init_set(aux_fat, max_fat);	
	for (i = 0; i < factors; i++)
		mpf_mul_ui(aux_fat, aux_fat, three_k - i);	

	/* (1 + (3k)^2) / (3k)!  */
	mpf_ui_div(term, three_k * three_k + 1, aux_fat);

	/* última thread cálcula o maior fatorial (último termo) */
	if (id == num_threads - 1)
		mpf_set(last_fat, aux_fat);

	mpf_clear(aux_fat);
}

/**************************************************************************************************/

void taylor_p(mpf_t term, int k, int id)
{
	int i, factors, p_k = p * k;
	mpf_t aux_fat, prod, sum;

	/* denominador: (pk)! = pk (pk-1) ... max_fat*/
	factors = total_terms == 0 ? p_k : p + p_k - p * total_terms;	
	mpf_init_set(aux_fat, max_fat);	
	for (i = 0; i < factors; i++)
		mpf_mul_ui(aux_fat, aux_fat, p_k - i);	

	/* numerador: 1 + pk + pk (pk-1) + ... + pk (pk-1) ... (pk-p+2) */
	mpf_init_set_ui(prod, p_k); 
	mpf_init_set_ui(sum, 1 + p_k);
	for (i = 1; i <= p - 2; ++i) {
		mpf_mul_ui(prod, prod, p_k - i);
		mpf_add(sum, sum, prod);		
	}

	/* numerador / denominador */
	mpf_div(term, sum, aux_fat);	

	/* última thread cálcula o maior fatorial (último termo) */
	if (id == num_threads - 1)
		mpf_set(last_fat, aux_fat);	

	mpf_clear(aux_fat);
	mpf_clear(prod);
	mpf_clear(sum);
}

/**************************************************************************************************/

