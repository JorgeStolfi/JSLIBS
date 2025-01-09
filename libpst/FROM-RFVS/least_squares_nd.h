
#ifndef ls_squares_nd_H
#define ls_squares_nd_H

#include <stdio.h>
#include <rn.h>

struct poly_function_t {
  /*Describes an polynomial function P(X) = w0Xe0 + w1Xe1...+ wmXm
  P(X): Rd=>R
  where m is the number of possible coeficients
  Xei is X[0]^coef[i][0] + *X[1]^coef[i][1] + .. + X[d]^coef[i][d]
  where d is the number of dimensions */
  int num_coefs; //how much expoent coeficients we have
  double** coefs; //store the expoent of the dimentions, it is {num_coefs}x{dimensions}
  int dimensions; //dimentionality of P(X)
  double* weights; //the weight of each coeficient
};

typedef struct poly_function_t poly_function_t;

double** CreateCoeficientsVector(int degree,int dimensions,int* num_coefs);
/* Given the degree of the polynon, it generates a list of ${num_coefs} 
with ${degree} positions of possible coeficients for the polynon in R-{dimensions}.
EG:
degree = 2 dimensions=3
it returns 
(  0,  0,  0) (  0,  0,  1) (  0,  0,  2) (  0,  1,  0) (  0,  1,  1)
(  0,  2,  0) (  1,  0,  0) (  1,  0,  1) (  1,  1,  0) (  2,  0,  0) 

degree = 2 dimensions=2
it returns (0,0) (0,1) (0,2) (1,0) (1,1) (2,0) */

double** CreateCoeficientsVectorHomogeneous(int degree,int dimensions,int* num_coefs);

double EvaluateXE(int dimensions, double* x, double* coefs);
/*Given an array {x} and {coef} with dimension {dimensions},
evaluates \sum_{i=0 to dimensions} x_i^coef_i   */

double EvaluatePX(int dimensions,double* x,int num_coefs,double** coefs,double* weights);
/* Given an array {x}, an coeficent array {coefs}, and a weight array {weights} with dimension {dimensions} 
Evaluates P(x) = \sum{i=0 to num_coefs} weights_i*Xe*/


//Remember : A(yT) = Bt !!!

double* EvaluateB(int dimensions, int N,int m, double* F, double** X,double** coefs);
/*Evaluates the B 1xm array from the system */

double* EvaluateA(int dimensions,int N, int m, double** X, double** coefs);
/*Evaluates the A mxm* matrix of the system. Remember that the result is a linearized matrix !*/

double* SolveSystem(double* A,double*B,int m);
/*Given {A} (linearized mxm) and B (mx1) , where AX = B. Computes X (mx1)*/


double EvaluatePvalue(poly_function_t* p,double* x);
/*evaluates the polynomial function of {x}. {x} MUST have same dimension as {p}->dimensions*/

poly_function_t* LS_PolyFitting(int dimensions, int N, int degree,double** X,double* F,int isHomogeneous);
/*Fits and polynomial weighted function and returns its parameters in data structure.
The input data must be:
{N} which is the number of samples
{dimensions} which is how many dimensions {X} have
{degree} which is the degree of the polynomial function to be fitted
an matrix {X} (Nxdimensions) where each row corresponds to a sample coordinate
and array {F} where F[i] = f(X[i])
*/

void PrintfPolyFunction(FILE* arq ,poly_function_t* p);
/* Prints in an FILE the structure contents*/
poly_function_t* ReadPolyFunction(FILE* arq);

#endif