
#ifndef GLOBALS_H
#define GLOBALS_H

#include<iostream>
#include <vector>
#include <iostream>
#include<ctime>
#include <stdio.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include<math.h>
#include <map>
#include<cmath>
#include<chrono>
#include<random>
#include<string>
#include <openssl/rand.h>
#include <cstdint>
#include <array>
#include <cstdint>
#include <vector>
#include <unistd.h>
#include<cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <cblas.h>


using namespace std;

//global variables
extern int m,K,ii,*snptype,*snppos,known;
extern long long int V;
extern char buff[256],*newWorkingDirectory;
extern FILE *latfile;
extern std::random_device rd;
extern map<int, int**> variant_dict;
extern float *alpha;
extern float *beta;
extern short ***n_m_t,*Docs_1,*idf;
extern unsigned int **n_m_z, **n_z_t;

//global functions
void set();
void set_variables();
void gsl_ran_multinomial_new (const gsl_rng * r, const size_t K, const unsigned int N, const double p[], short n[]);
void initialize(std::uint64_t seed );
void epochs(std::uint64_t e_seed,int iter);
void construct_dict();
void test_initialize(std::uint64_t seed);
void test_set();
#endif // GLOBALS_H
