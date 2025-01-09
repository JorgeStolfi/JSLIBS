#ifndef __LIGHTING_COMPACT_H__
#define __LIGHTING_COMPACT_H__
#include <r2.h>
#include <r3.h>
#include <argparser.h>
#include <tabela.h>
#include <approx_system.h>


typedef struct lighting_compact_options_t{
  /*Compact*/
  r3_t dirlight;
  double err_dirlight;
  double radlight;
  double err_radlight;
  double shine;
  double err_shine;
  bool_t has_isotropic;
  bool_t has_backlight;
  r3_t backnormal;
  double err_backnormal;
} lighting_compact_options_t;


typedef struct lighting_compact_data_t{
  r3_t view_dir; /*Direction of the point of view (FIXED)*/

/*Non linear coefficients*/
  r3_t dirlight; /*Light source direction*/
  double radlight; /*Angular radius of the light source*/
  double shine; /*Intensity of the phong factor */
  r3_t backnormal; /*Background*/
  
  /*Uncertainties of the non linear coefficients*/
  double err_dirlight;
  double err_radlight;
  double err_shine;
  double err_backnormal;
  
  /*Linear coefficients*/
  double Ec; /*Intensity of light source*/
  double El; /*Intensity of the phong shading factor*/
  double Ei; /*Isotropic light source intensity*/
  double Ef; /*Backplane light source intensity*/
  
  /*Which basis elements to include*/ 
  bool_t Hc; /*Compact light (dirlight,radlight)*/
  bool_t Hl; /*Highlight (dirlight,shine)*/
  bool_t Hi; /*isotropic light*/
  bool_t Hf; /*backplane lighting (backnormal)*/
    
} lighting_compact_data_t;

lighting_compact_options_t* lighting_compact_parse(argparser_t* pp);
approx_model_t* lighting_compact_create_approx_lighting_model(void);

lighting_compact_data_t* lighting_compact_init_components(lighting_compact_options_t* o, r3_t view_dir);

/* Approx Model functions */
void lighting_compact_determine_basis_indices(void* l_data,int* ne, int* num_compac, int* num_phong, int* num_isotropic, int* num_backplane);
double lighting_compact_phi(int r, double* x,void* l_data);
int lighting_compact_get_number_components(void* l_data);
void lighting_compact_set_alphas(double* C,int n, void* l_data);
void lighting_compact_get_alphas(void* l_data,double* C,int n);
void* lighting_compact_copy_lighting_data(void* l_data);
void lighting_compact_release_lighting_data(void* l_data);
double lighting_compact_shading(double* x,void* l_data);
void lighting_compact_write_parameters(FILE* arq,void* l_data);
void* lighting_compact_read_parameters(FILE* arq);
// void lighting_compact_update_weights(void* l_data,double** X, double* F,double* weights, int n,double* sigma);
double  light_compact_simple_compute_weight(lighting_compact_data_t* cl,r3_t normal);
double lighting_compact_compare_lightings(void* l_data1,void* l_data2);
int lighting_compact_num_nl_parameters(void* l_data);
void lighting_compact_pack_nl_parameters(void* l_data,double* parameters,int n);
void lighting_compact_get_nl_centers_and_deltas(void* l_data,double* center,double* delta,int n);
void lighting_compact_unpack_nl_parameters(double* parameters,int n,void* l_data);
void lighting_compact_update_params_and_errors(double* parameters,int n,void* l_data, double alpha, double beta, double gamma);
void lighting_compact_generate_simplex(void* l_data, int n, double* S);

/*Those are to estimate the light source direction*/
double lighting_compact_simple_phi(int r, double* x,void* l_data);
int lighting_compact_simple_get_number_components(void* l_data);
void lighting_compact_simple_set_alphas(double* C,int n, void* l_data);
void lighting_compact_simple_get_alphas(void* l_data,double* C,int n);
approx_model_t* lighting_compact_simple_create_approx_lighting_model(void);

void light_compact_estimate_lightdir(double** X, double* F,double* wpos, int n,void* l_data,int update_steps);
/*Ad-Hoc light source direction estimation  with {update_steps} iterations- Given the set points X,F(X) and its fixed weights {wpos}
computes a estimative for lighdir,Ec and Ef OR Ei by Least Squares */
void light_compact_simple_update_weights(void* l_data,double** X,double* weights, int n);

#endif
