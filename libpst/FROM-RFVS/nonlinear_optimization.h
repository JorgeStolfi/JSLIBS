#ifndef __NONLINEAR_H__
#define __NONLINEAR_H__

typedef double goal_function_t(double parameters[], int num_parms );
typedef void clip_function_t(double center[], double delta[], int num_parms );

void nonlinear_optimize(
  double center_original[], double delta_original[], int num_parms, /*Input*/
  double r, double rmin, double alpha, int max_iters,double epsilon,int which, /*Adjust parameters*/
  goal_function_t* goal_function, clip_function_t* clip_function, /*Callback functions*/
  double optmized_parameters[], /*Output*/
  char* debug_prefix
);

void nonlinear_optimize_z(
  double initial_z[], int num_parms, /*Input*/
  double r,double theta, double alpha, int max_iters,double epsilon,int which, /*Adjust parameters*/
  goal_function_t* goal_function, /*Callback functions*/
  double opt_z[], /*Output*/
  char* debug_prefix
);
void plot_nonlinear_goalfunction(char* prefix,double center[],double delta[], int num_parms, goal_function_t goal_function);
#endif