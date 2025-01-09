#ifndef __MODELS_H__
#define __MODELS_H__
#include <r2.h>
#include <r3.h>
#include <argparser.h>
#include <approx_system.h>




/*Generic Lighting model functions*/

approx_model_t* create_approx_lighting_model(model_type_t type);

double lighting_error_parse(argparser_t* pp,double emax);







/*
Given an arparser {pp} read the type of the model and returns a lighting_data_model_t structure woth its initializers
except view_dir.
The initializer parameters are readen only if {determineParams} as TRUE.
Otherwise it will initialize only the field lightingModelType

  the input format is the following:
  POINTLIKE_MODEL 
  pointlike [backplane|plain] [highlightK {K}] [terminatorRho {RHO}]
  HARMONIC_MODEL 
  harmonic degree {DEGREE}
  HARMONICSG_MODEL
  harmonicSG degree {DEGREE}
  RADIALBASIS_MODEL
  radial resolution {RESOLUTION} span {SPAN}
  GLOSSYLIKE_MODEL
  glossy glossiness {G0} {G1}
  COMPACT_MODEL
  compact dirlight {Ux Uy Uz} dirbackplane {Xix Xiy Xiz } rho {RHO} K {K}
 
*/

void lighting_model_generate_plot(FILE* arq, approx_model_t* lm_test,void* l_datatest,long int num_samples, r3_t view_dir, r3_t ref_dir);
void lighting_model_generate_SG_plot(char* prefix, approx_model_t* lm_test,void* l_datatest, long int num_samples_axis,r3_t view_dir);

void generate_compare_plot(FILE* arq, approx_model_t* lm_test,void* l_datatest, approx_model_t* lm_approx, void* l_datapprox, long int num_samples, r3_t view_dir, r3_t ref_dir);

void   generate_SG_plot(char* prefix, approx_model_t* lm_test,void* l_datatest, approx_model_t* lm_approx, void* l_datapprox, long int num_samples_axis,r3_t view_dir);


#endif
