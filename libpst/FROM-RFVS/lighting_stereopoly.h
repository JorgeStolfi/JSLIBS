#ifndef __LIGHTING_STEREOPOLY_H__
#define __LIGHTING_STEREOPOLY_H__
#include <r2.h>
#include <r3.h>
#include <argparser.h>
#include <approx_system.h>


/* Stereographic polinomial Model Functions*/

typedef struct lighting_stereopoly_options_t{
  /*harmonic*/
  int degree;

} lighting_stereopoly_options_t;

typedef struct lighting_stereopoly_data_t{
  int num_comp; /*Number of terms in the polynomial*/
  int** tuples; /*Degree pairs of each term*/
  double* coeficients /*Coefficient of each term*/;
  r3_t view_dir;
} lighting_stereopoly_data_t;

lighting_stereopoly_options_t* lighting_stereopoly_parse(argparser_t* pp);
approx_model_t* lighting_stereopoly_create_approx_lighting_model(void);

lighting_stereopoly_data_t* lighting_stereopoly_init_components(lighting_stereopoly_options_t* o,r3_t view_dir);
double lighting_stereopoly_phi(int r, double* x, void* l_data);
int lighting_stereopoly_get_number_components(void* l_data);
void lighting_stereopoly_set_alphas(double* C,int n,void* l_data);
void lighting_stereopoly_get_alphas(void* l_data,double* C,int n);
void* lighting_stereopoly_copy_lighting_data(void* l_data);
void lighting_stereopoly_release_lighting_data(void* l_data);
double lighting_stereopoly_shading(double* x,void* l_data);
void lighting_stereopoly_write_parameters(FILE* arq,void* l_data);
void* lighting_stereopoly_read_parameters(FILE* arq);



#endif
