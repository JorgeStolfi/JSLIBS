#ifndef __ANALYINV_H__
#define __ANALYINV_H__

#include <r3.h>
#include <approx_system.h>
#include <float_image.h>

typedef struct analytic_inversion_t{
  int type_gauge;
  int m; /*number of illumination fields*/
  /*Direct models */  
  approx_model_t* gauge_model; /*analytic model of a gauge under each illumination*/
  void** gauge_data; /*gauge model parameters - must be of a compatible type specified within {gauge_model}*/
  /*Inverse model*/
  double* u; /*Main direction of signature set*/
  double* v; /*Secondary direction of the signature set*/
  double *b; /*Baricenter of the signature set*/
  double bu,bv; /*{U,V} coordinates of the barycenter*/
  double radius; /*Radius of the plane*/
  int type_inv;
  approx_model_t* inv_model; /*Maps a projected signature to a normal coordinate*/
  void* inv_data[3]; /*parameters of the 2D approximation model for each normal coordinate- must be compatible with {inv_model}*/
  void* mag_data; /*Parameters of the 2D approximation model for each magntude*/
} analytic_inversion_t;

analytic_inversion_t* analytic_inversion_create_model
( 
 int m,
 int type_gauge,
 approx_model_t* gauge_model,
 void** gauge_data,
 r3_t view_dir,
 double thetaMax,
 int resolution,
 int type_inv,
 approx_model_t* inv_model,
 void** inv_data,
 void* mag_data,
 char* plot_pov_prefix,
 char* plot_txt_prefix
);
/* Returns a inverse model of a set of ${m} gauges adjusted to the hemisphere with center {view_dir} and angular radius ${thetaMax},
using {resolution} rings of samples */

void analytic_inversion_compute_normal(double* S, analytic_inversion_t* ai,r3_t* normal, double* albedo, double G[], double* dist);
/*
  Computes the normal {normal} and albedo {albedo} of a observation vector {S} from the analytic model {ai}.
  Also, returns in {G} the best matching gauge observation vector. Returns in {dist} the distance between the scene and gauge signatures.
*/

float_image_t* analytic_inversion_plot_map(analytic_inversion_t* ai,int imSize);

void analytic_inversion_compute_gauge_OV(analytic_inversion_t* ai, r3_t normal, double G[]);
/*Given a analytic inversion model and a normal, computes the gauve OV and stores into {G}*/

void analytic_inversion_write(FILE* arq,analytic_inversion_t* ai);
/*Writes a self-contained file with the contents of {ai} data structure into {arq}.*/
analytic_inversion_t* analytic_inversion_read(FILE* arq);
/*Reads a self-contained file with the contents of analytic_inversion_t data structure from {arq}. The file should be
self-coinatined (that is, written by {analytic_inversion_write}) with all information relevant to do not require any user intervention*/
void   analytic_inversion_generate_SG_plot(char* prefix, analytic_inversion_t* ai, long int num_samples_axis,r3_t view_dir);
void analytic_inversion_generate_UV_plot(char* prefix, analytic_inversion_t* ai, long int num_samples_axis,r3_t view_dir);

#endif
