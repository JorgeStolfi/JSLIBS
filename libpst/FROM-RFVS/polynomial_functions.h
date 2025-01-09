/*Callback functions for the polinomial interpolation !*/
#ifndef __POLY_MODEL_H_
#define __POLY_MODEL_H_

#include <least_squares_nd.h>
#include <approx_system.h>

typedef struct poly_function_options_t{
  int dimentions;
  int degree;
  bool_t homogeneous;
} poly_function_options_t;

approx_model_t* create_ls_polynomial_model(void);

poly_function_options_t* poly_model_parse(argparser_t* pp);
poly_function_t* poly_init_components(int dimensions, int degree, bool_t homogeneous);
poly_function_t* poly_model_init_components(poly_function_options_t o);

double poly_phi(int r, double* x,void* l_data);
int poly_get_number_components(void* l_data);
void poly_set_alphas(double* C, int n,void* l_data);
void poly_get_alphas(void* l_data,double* C, int n);
double poly_evaluate(double* x,void* l_data);
void poly_write_parameters(FILE* arq,void* l_data);
void* poly_read_parameters(FILE* arq);
void* poly_copy_data(void* l_data);
void poly_release_data(void* l_data);

#endif