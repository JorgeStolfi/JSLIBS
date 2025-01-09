#ifndef __ANALYALPHA_H__
#define __ANALYALPHA_H__

#include <analytic_inversion.h>
#include <r3.h>
#include <float_image.h>

typedef double analytic_alpha_log_prob_function_t(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 );

analytic_alpha_log_prob_function_t* escolhe_analytic_alpha_log_prob_function_t(int nLogProb);

double analytic_alpha_log_prob_withalbedo(double alb, double S[],double G[],int m,double sigma, double omg0, double omg1);
/*
  Computes the negated logarithm of the probability density function for the observation vector {S[0..m-1]}, assuming 
  that the gauge observation vector of the same normal is {G[0..m-1]}. Assumes that an inlier observation S[i]
  has gaussian distribution  with mean {alb*G[i]} and deviation {sigma}, where {alb} is the scene albedo. Also assumes 
  S[i] maybe an outlier lower than {alb*G[i]} with probability {omg0} and an outlier higher than {alb*G[i]} 
  with probability {omg1}
 
 */

double analytic_alpha_log_prob_noalbedo07(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 );
double analytic_alpha_log_prob_noalbedo08(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 );
double analytic_alpha_log_prob_noalbedo09(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 );
/*
  Sames as {analytic_alpha_log_prob_withalbedo}, but assumes an unknow albedo with uniform distribution
  between 0 and 1.  Uses Est_LogProb0X, where X is the identifier of the function
 */



void analytic_alpha_get_nl_centers_and_deltas(
  r3_t s, double err_s,
  double a,double err_a,
  double* center,double* delta,
  int num_parms);
  
  
void analytic_alpha_generate_simplex( double center[],double delta[], int n, double* S);

void analytic_alpha_compute_normal(
  double S[],int m,analytic_inversion_t* ai, bool_t optimize_albedo,/*Input*/
  r3_t* normal, double* albedo,double* dist,double G[], /*Output*/
  int nLogProb,double sigma,double omg0, double omg1, int update_steps,double tol_albedo,double tol_normal, r3_t view_dir,char* debug_prefix
);

void analytic_alpha_optmize_parameters(
  double S[],int m, r3_t normal_guess, double tol_normal,bool_t optimize_albedo, double albedo_guess, double tol_albedo,
  approx_model_t * gauge_model, void** gauge_data,
  double sigma, double omg0, double omg1, int update_steps,
  r3_t* normal, double* albedo, double* dist,
  r3_t view_dir,char* debug_prefix
);

double analytic_alpha_goal_function(analytic_alpha_log_prob_function_t* est_logprob ,r3_t normal,bool_t known_albedo,double albedo, double S[] ,int m, approx_model_t* lm, void** l_data, double sigma, double omg0, double omg1 );
void analytic_alpha_compute_normal_aux(
  double S[],analytic_inversion_t* ai,int m,r3_t normal_guess,double rad_normal,  bool_t optimize_albedo,double albedo_guess,double rad_albedo,/*Input*/
  r3_t* normal, double* albedo,double* dist,double G[], /*Output*/
  int nLogProb,double sigma,double omg0, double omg1, int update_steps,double tol_albedo,double tol_normal, r3_t view_dir,char* debug_prefix
);

float_image_t* analytic_alpha_plot_goal_2D(analytic_alpha_log_prob_function_t* est_logprob, r3_t correct_normal, bool_t known_albedo, double albedo,double S[], int m , approx_model_t* gauge_model, void** gauge_data,double sigma,double omg0, double omg1,r3_t view_dir, int imSize);
// float_image_t* analytic_alpha_plot_goal(r3_t correct_normal, bool_t known_albedo, double albedo,double S[], int m , approx_model_t* gauge_model, void** gauge_data,double sigma,double omg0, double omg1,r3_t view_dir, int imSize);

#endif