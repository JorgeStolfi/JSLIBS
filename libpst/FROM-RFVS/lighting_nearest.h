#ifndef __LIGHTING_NEAREST_H__
#define __LIGHTING_NEAREST_H__
#include <r2.h>
#include <r3.h>
#include <argparser.h>
#include <tabela.h>
#include <approx_system.h>


typedef struct lighting_nearest_options_t{
  gauge_data_t* gin;
  int channel;
}lighting_nearest_options_t;


typedef struct lighting_nearest_data_t{
  double* val;
  r3_t* normal;
  int n;
  double alpha; /*Coeficient of linear model*/
} lighting_nearest_data_t ;

lighting_nearest_options_t* lighting_nearest_parse(argparser_t* pp);
approx_model_t* lighting_nearest_create_approx_lighting_model(void);

/*approx model functions */
lighting_nearest_data_t* lighting_nearest_init_components(lighting_nearest_options_t* o,r3_t view_dir);
double lighting_nearest_phi(int r, double* x, void* l_data);
int lighting_nearest_get_number_components(void* l_data);
void lighting_nearest_set_alphas(double* C,int n,void* l_data);
void lighting_nearest_get_alphas(void* l_data,double* C,int n);
void* lighting_nearest_copy_lighting_data(void* l_data);
void lighting_nearest_release_lighting_data(void* l_data);
double lighting_nearest_shading(double* x,void* l_data);
void lighting_nearest_write_parameters(FILE* arq,void* l_data);
void* lighting_nearest_read_parameters(FILE* arq);



#endif
