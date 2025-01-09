#ifndef __LIGHTING_RADIAL_H__
#define __LIGHTING_RADIAL_H__
#include <r2.h>
#include <r3.h>
#include <argparser.h>
#include <approx_system.h>

/*Radial Basis model*/
typedef struct lighting_radial_options_t{
  /*radial*/
  int resolution;
  double span;
  
} lighting_radial_options_t ;


typedef struct lighting_radial_data_t{
  int num_comp; /*Number of radial bases*/
  r3_t* bcenter; /*centers of each radial basis coeficient*/
  double* bradius;/*radius of each radial basis coeficient*/
  double* coeficients; /*ditto*/
  r3_t view_dir;
} lighting_radial_data_t;

lighting_radial_options_t* lighting_radial_parse(argparser_t* pp);
approx_model_t* lighting_radial_create_approx_lighting_model(void);

/*Approx model functions*/
lighting_radial_data_t* lighting_radial_init_components(lighting_radial_options_t* o,r2_t center,double raio, r2_t estica,r3_t pole,double thetaMax);
double lighting_radial_phi(int r, double* x, void* l_data);
int lighting_radial_get_number_components(void* l_data);
void lighting_radial_set_alphas(double* C,int n, void* l_data);
void lighting_radial_get_alphas(void* l_data,double* C,int n);
void* lighting_radial_copy_lighting_data(void* l_data);
void lighting_radial_release_lighting_data(void* l_data);
double lighting_radial_shading(double* x,void* l_data);
void lighting_radial_write_parameters(FILE* arq,void* l_data);
void* lighting_radial_read_parameters(FILE* arq);

#endif
