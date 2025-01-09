#ifndef __LIGHTING_HARMONIC_H__
#define __LIGHTING_HARMONIC_H__
#include <r2.h>
#include <r3.h>
#include <argparser.h>
#include <approx_system.h>


typedef struct lighting_harmonic_options_t{
  /*harmonic*/
  int degree;
} lighting_harmonic_options_t;


typedef struct lighting_harmonic_data_t{
  int num_comp; /*Number of triplests in the harmonic function*/
  int** triplets; /*Set of triplets of the harmonic function*/
  double* coeficients /*Value of each coeficient of the harmonic function*/;
  r3_t view_dir;
} lighting_harmonic_data_t ;

lighting_harmonic_options_t* lighting_harmonic_parse(argparser_t* pp);
approx_model_t* lighting_harmonic_create_approx_lighting_model(void);

/*approx model functions */
lighting_harmonic_data_t* lighting_harmonic_init_components(lighting_harmonic_options_t* o,r3_t view_dir);
double lighting_harmonic_phi(int r, double* x, void* l_data);
int lighting_harmonic_get_number_components(void* l_data);
void lighting_harmonic_set_alphas(double* C,int n,void* l_data);
void lighting_harmonic_get_alphas(void* l_data,double* C,int n);
void* lighting_harmonic_copy_lighting_data(void* l_data);
void lighting_harmonic_release_lighting_data(void* l_data);
double lighting_harmonic_shading(double* x,void* l_data);
void lighting_harmonic_write_parameters(FILE* arq,void* l_data);
void* lighting_harmonic_read_parameters(FILE* arq);



#endif
