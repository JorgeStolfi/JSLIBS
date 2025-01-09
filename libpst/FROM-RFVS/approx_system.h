#ifndef __LSSYSTEM_H__
#define __LSSYSTEM_H__

#include <stdio.h>
#include <bool.h>
#include <argparser.h>

/*This library provides a LS fitting function including support to weights.
It rellies in callback functions defined bellow, it is up to the library which
defines the desired model to initialize ALL the attributes of the ls_model_t structure.
Please reffer to lighting_models.h and lighting_models.c for a example */


typedef enum{ HARMONIC_MODEL, RADIALBASIS_MODEL, STEREOPOLY_MODEL , COMPACT_MODEL,NEAREST_MODEL,GENERICPOLY_MODEL} model_type_t;

typedef struct model_options_t{
  model_type_t modelType;
  void* options;
} model_options_t;

model_options_t model_parse(argparser_t* pp, bool_t parseValues);

/*Types of each function*/
typedef double phi_function(int i, double* x,void* l_data);
/*Given the coeficient number{i}, the argument vector {x} and given lighting model {l_data} computes the value of Phi{i}*/
typedef int number_components_function(void* l_data);
/*Give the number of the coeficients for the LS system*/
typedef void  set_alphas_function(double* C, int n, void* l_data);
/*Stores a array of linear coeficients {C} with ${n} basis into {L_data} structure*/
typedef void get_alphas_function(void* l_data, double* C, int n);
/*Stores internal linear coeficients of {l_data} into the {C} array */
typedef double evaluate_function(double* x,void* l_data);
/*Given a normal and the lighting model, returns its light intensity*/
typedef void write_function(FILE* arq,void* l_data);
/*writes {l_data} parameters into {arq} the parameters in humam-readable fashion*/
typedef void* read_function(FILE* arq);
/*reads {l_data} parameters from {arq} the parameters in humam-readable fashion written by write_function*/
typedef void* copy_data_function(void* l_data);
/*Given a lighting_data structure creates a copy of it with same parameters*/
typedef double compare_function(void* l_data1,void* l_data2);
/*Given two lighting data structures , compute its difference*/
typedef void release_data_function(void* l_data);
/*Releases a given lighting data structure*/
typedef int num_nl_parameters_function(void* l_data);
/*This function returns the number of uncertain non-linear parameters. It has to be the size of the array returned by {unpack_nl_parameters_function} */
typedef void pack_nl_parameters_function(void* l_data,double* parameters, int n);
/*This function type converts the uncertain non linear parameters to an array of parameters values */
typedef void get_nl_centers_and_deltas_function(void* l_data,double* center, double* delta, int n);
/*Convert the uncertain non linear parameters and their errors to arrays of center-delta intervals*/
typedef void unpack_nl_parameters_function(double* parameters,int n, void* l_data);
/*This function type sets the uncertain non linear  parameters from array of parameter values.*/
typedef void update_nl_parameters_and_errors_function(double* parameters,int n, void* l_data,double alpha, double beta, double gamma);
/*This function performs as unpack_parameters but also updates the uncertainties during the non linear optimization*/
typedef void generate_nl_parameter_simplex_function(void* l_data, int n, double* S);
/*Generates a simplex of packed parameter vectors, within current uncertainty bounds*/


/*LS generic model fitting data structure, contain all needed callback functions*/
typedef struct approx_model_t{
  /*Mandatory members -  without then the fitting procedures will crash ! */
  
  int type; /*This is meant to be used for debugging and enum types*/
//  void* ls_data; /*this contain the data structure which defines the modeled function F, stores initially the initial  guess and in the end the fitted function*/
  phi_function* phi; /*Linear basis function*/ 
  number_components_function* get_num_components; /*retrieve the number of linear basis functions*/
  set_alphas_function* set_alphas; /*Sets the linear coeficientes from the model from a vector*/
  get_alphas_function* get_alphas; /*Gets the linear coeficients from the model and stores into a vector*/
  copy_data_function* copy_data; /*Make a copy of ls_data, it is needed to "save" the results of an previous iteration*/
  release_data_function* release_data; /*Release a ls_data pointer from memory, mandatory if you dont wish to have your memory taken by multiple iterations !*/
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  evaluate_function* evaluate; /* compute the value of the fitting function, it is optional*/
  write_function*  write_parameters; /*write the params in human readable manner into a file*/
  read_function*  read_parameters; /*read the params from a file in human readable manner*/
  /*for non linear models*/
  compare_function* compare; /*compare the non linear parameters of two data records*/
  num_nl_parameters_function* get_num_nl_parameters; /*Number of non-linear parameters*/
  pack_nl_parameters_function* pack_nl_parameters; /*Stores non-linear parameters from the data structure into a vector */
  get_nl_centers_and_deltas_function* get_nl_centers_and_deltas;
  unpack_nl_parameters_function* unpack_nl_parameters; /* Sets non linear parameters from a vector into a data structure */
  update_nl_parameters_and_errors_function* update_nl_parameters_and_errors; /*Udate the parameters and errors after a non linear iteration optimization*/
  generate_nl_parameter_simplex_function* generate_nl_parameter_simplex; /*Generates simplex in the non-linear parameters space within the error bounds*/
  
  
} approx_model_t;


void computeLeastSquaresMatrix(double* A,phi_function* phi, int basis_size, double** X,int n, double w[], double wpos[],void* l_data);
/*Computes the Least Square matrix A using the weights,n-dimensional X values and current values of the function.*/

void computeLeastSquaresRHSVector(double* b,phi_function* phi, int basis_size,double** X,double* F,int n , double w[], double wpos[],void* l_data);
/*Same above, but for RHS vector*/

void computeLSTerms(double* A,
	double* b ,
	double* c,
	phi_function* phi,
	int basis_size,
	double** X,
	double* F,
	int n,
	double* w,
	double* wpos,
	void* l_data
	);
	
/*Compute the LS terms for the given matrix {A}, RHS vector {b} and stores it in {c}, weights and current l*/

void approx_model_fit_linear_model(double** X, double* F,double* w,double* wpos, int n, approx_model_t* am, void* l_data);
/* 
Adjusts the linear coeficients of {lm} by Least Squares method, using the current values 
of the non linear parameters. If {update_steps} and {outliers_fraction} are positive, iterates the adjustment at most {update_steps},
recomputing the a posteriori weights {Wo} according to the outlier model given fraction {outliers_fraction} of outliers.
*/

double approx_model_compute_S_star(approx_model_t* am, double** X, double* F,double* w, double* wpos, int n, double* parameters,void* l_data_original);

void plot_non_linear_goal_function(char* prefix,double** X, double* F,double* w,double* wpos, int n, approx_model_t* am, void* l_data);

void approx_model_fit_non_linear_model_simplex(double** X, double* F, double* w,double* wpos, int n, approx_model_t* am,  void* l_data,  int update_steps,double alpha, double beta, double gamma);
/*
Adjusts the linear and non-linear parameters of {lm} by iterated quadratic minimization. The non linear parameters are
adjusted by fitting a quadratic function to the mean quadratic error and obtaining its minimum, repeated at maximum {num_iterations} iterations.
At each iteration adjusts the linear coeficients of {lm} by Least Squares method.
If {update_steps} and {outliers_fraction} are positive, iterates the adjustment at most {update_steps},
recomputing the a posteriori weights {Wo} according to the outlier model given fraction {outliers_fraction} of outliers
*/

#endif     