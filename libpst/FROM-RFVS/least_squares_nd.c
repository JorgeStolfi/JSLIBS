#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <rn.h>
#include <rmxn.h>
#include <math.h>
#include <least_squares_nd.h>

/*
This library implements a least squares fitting of a polynomial function f:Rn->R given N data points
The method is described in the article "Multi-Dimensional Least-Squares Polynomial Curve Fitting
(Communications of the ACM, Volume 2 ,  Issue 9  (September 1959) Pages: 29 - 30  Year of Publication: 1959)
*/

/*Internal Stack implementation  for use in the generationf of the coeficients*/

struct double_stack_t{
    double* array;
    struct double_stack_t* next;
};

typedef struct double_stack_t double_stack_t;

int stack_count(double_stack_t* s);
/*returns how many elements we have on stack*/

void stack_push(double_stack_t** sp,double* a, int n);
/*Given an pointer to a stack {sp}, an array {a} with {n} elements
push an new element on top of stack, updating ${sp}. The content of *sp
must be NULL for an empty stack.*/

double* stack_pop(double_stack_t** sp);
/*Given an pointer to a stack ${sp} it pops out an element of the stack,
updating it. An empty stack will update *sp to NULL.*/

/* Internal functions of the library*/

void write_elem_rec(double* coef, int col, int total, int n,double_stack_t** sp);
/*This is a recursive function to generate all the possible coeficients of an
polynomial with degree {n}. Each computed coeficient vector is pushed into the stack {sp}
{total} indicates the value of the coeficiens already writen in to the buffer array {coef}.
{col} indicates the current position on {coef} where the function evaluates the coeficient.
It works in the following manner:
1)Stop  Conditions
a) coef is full (col == n) . Just pushes {coef} into the stack
b) the computed coeficients fulfill the condition of an polinomial of degree {n}.
   {coef} is filled with zeros in the remaining positions and pushed into teh stack
2)Recursion
  the function is called recursively {total} times, from i = {total} to 0 and for each call the {col} position
  in {coef} have i value assigned. For each call i is subtracted from total and col is incremented.
*/


/* Internal functions*/

void EvaluateBi(int dimensions, double*x ,int num_coefs, double** coefs,double* B);
/*Computes Bi in {B}, which is an {num_coefs} dimensional vector 
Bi is defined as (Xe1, Xe2,...,Xem)*/

int stack_count(double_stack_t* s){
  if(s != NULL){
    return 1 + stack_count(s->next);
  }else{
    return 0;
  }
}

void stack_push(double_stack_t** sp,double* a, int n){
  double_stack_t* novo = (double_stack_t*)malloc(sizeof(double_stack_t));
  novo->next = *sp;
  novo->array =  (double*)malloc(sizeof(double)*n);
  int i;
  for(i = 0; i < n;i++){	
    novo->array[i] = a[i];
  }
 // rn_gen_print(stdout,n,a,"%3.0lf","(", ",", ")");
//  rn_gen_print(stdout,n,novo->array,"%3.0lf","(", ",", ")\n");
  *sp=novo;
}

double* stack_pop(double_stack_t** sp){
  double_stack_t* old = *sp;
  if(old == NULL) return NULL;
  double* a = old->array;
  *sp = old->next;
  free(old);
  return a;
}

/* end of internal stack implementation*/

/* Internal functions */
void write_elem_rec(double* coef, int col, int total, int n,double_stack_t** sp){
  if (total == 0){
    int i;
    for(i = col ; i < n; i++){
      coef[i] = 0;
    }
   // rn_gen_print(stdout,n,coef,"%3.0lf","(", ",", ")");
    stack_push(sp,coef,n);
 //   fprintf(stdout,"\n");
    return;
  }

  if(col == n){
 //   rn_gen_print(stdout,n,coef,"%3.0lf","(", ",", ")");
//    fprintf(stdout,"\n");
    stack_push(sp,coef,n);
    return;
  }

  int i;
  for(i = total; i >= 0; i--){
	coef[col] = i;
	write_elem_rec(coef,col + 1, total-i, n,sp);
  }
  return ;
}

void write_elem_rec_homogeneous(double* coef, int col, int total, int n,double_stack_t** sp);
void write_elem_rec_homogeneous(double* coef, int col, int total, int n,double_stack_t** sp){

  if(col == n){
//     rn_gen_print(stdout,n,coef,"%3.0lf","(", ",", ")");
   // fprintf(stdout,"\n");
    stack_push(sp,coef,n);
    return;
  }

  int i;
  for(i = total ; i >= 0; i--){
	coef[col] = i;
	write_elem_rec_homogeneous(coef,col + 1, total, n,sp);
  }
  return ;
}

void EvaluateBi(int dimensions, double*x ,int num_coefs, double** coefs,double* B){
  int i;
  for(i = 0; i < num_coefs; i++){
    B[i] = EvaluateXE(dimensions,x,coefs[i]);
  }
}

/* External functions */


double** CreateCoeficientsVector(int degree,int dimensions,int* num_coefs){
  double** ep;
  double vec[dimensions];
 // rn_zero(degree,vec);
  double_stack_t* s = NULL;
  write_elem_rec(vec, 0, degree, dimensions,&s);
  int n = stack_count(s);
  //fprintf(stderr,"%d\n",n);
  *num_coefs = n;
  ep = (double**)malloc(sizeof(double*)*n);
  int i = 0;
  while(s!= NULL){
     double* coef = stack_pop(&s);
     ep[i] = coef;
     assert(i < n);
     i++;
    
  }
  
  return ep;
}




double** CreateCoeficientsVectorHomogeneous(int degree,int dimensions,int* num_coefs){
  double** ep;
  double vec[dimensions];
 // rn_zero(degree,vec);
  double_stack_t* s = NULL;
  write_elem_rec_homogeneous(vec, 0, degree, dimensions,&s);
  int n = stack_count(s);
  //fprintf(stderr,"%d\n",n);
  *num_coefs = n;
  ep = (double**)malloc(sizeof(double*)*n);
  int i = 0;
  while(s!= NULL){
     double* coef = stack_pop(&s);
     ep[i] = coef;
     assert(i < n);
     i++;
    
  }
  
  return ep;
}


double EvaluateXE(int dimensions, double* x, double* coefs){
  double val = 1;
  int i;
  for(i = 0; i < dimensions; i++){
    val*= pow(x[i],coefs[i]);
  }
  return val;
}

double EvaluatePX(int dimensions,double* x,int num_coefs,double** coefs,double* weights){
  double val = 0;
  int i;
  for(i = 0; i < num_coefs; i++){
    double* coef = coefs[i];
    val+= EvaluateXE(dimensions,x,coef)*weights[i];
  }
  return val;
}

double EvaluatePvalue(poly_function_t* p,double* x){
  double val = 0;
  int i;
  for(i = 0; i < p->num_coefs; i++){
    double* coef = p->coefs[i];
    val+= EvaluateXE(p->dimensions,x,coef)*(p->weights[i]);
  }
  return val;
}

double* EvaluateB(int dimensions, int N,int m, double* F, double** X,double** coefs){
  double* b = rn_alloc(m);
  int i;
  rn_zero(m,b);
  for(i = 0; i < N; i++){
    double Bi[m];
    EvaluateBi(dimensions,X[i] ,m,coefs,Bi);
    double FiBi[m];
    rn_scale(m,F[i],Bi,FiBi);
    rn_add(m,b,FiBi,b);
  }
  return b;
}

double* EvaluateA(int dimensions,int N, int m, double** X, double** coefs){
  double* A = rmxn_alloc(m, m);
  double* BtB = rmxn_alloc(m, m);
  int i;
  for(i = 0; i < N; i++){
    double Bi[m];
    EvaluateBi(dimensions,X[i] ,m,coefs,Bi);
    //Compute Bt * B, just consider A as Bt (mx1) and B as (1xm)
    rmxn_mul(m, 1,m, Bi, Bi, BtB);
    rn_add(m*m,A,BtB,A);
  }
  free(BtB);
  return A;
}

double* SolveSystem(double* A,double*B,int m){
  double* inverse = rmxn_alloc(m,m);
  double det = rmxn_inv_full(m, A, inverse);
  assert(det != 0);
  double* solution = rn_alloc(m);
  rmxn_map_col(m,m,inverse,B,solution);
  free(inverse);
  return solution;
  
}

poly_function_t* LS_PolyFitting(int dimensions, int N, int degree,double** X,double* F,int isHomogeneous){
  
  int m;
  
  double** coefs = NULL;
  if(!isHomogeneous){
    coefs = CreateCoeficientsVector(degree,dimensions,&m);
  }else{
    coefs = CreateCoeficientsVectorHomogeneous(degree,dimensions,&m);
  }
  double* B = EvaluateB(dimensions, N,m,F,X,coefs);
  double* A = EvaluateA(dimensions,N,m, X,coefs);
  double* weights = SolveSystem(A,B,m);
  poly_function_t* p =  (poly_function_t*)malloc(sizeof(poly_function_t));
  p->dimensions = dimensions;
  /* counts how much valid coeficients are valid (AKA no zero value)*/
  int i = 0;
  int p_m = 0;
  for(i = 0; i< m; i++){
//     if(fabs(weights[i]) >  1.0e-10)
      p_m++;
  }
  p->coefs = (double**)malloc(sizeof(double*)*p_m);
  p->weights = (double*)malloc(sizeof(double)*p_m);
  
  int count_m= 0;
  for(i = 0; i< m; i++){
   // if(fabs(weights[i]) >  1.0e-10){
      p->weights[count_m] = weights[i];
      p->coefs[count_m] = (double*)malloc(sizeof(double)*dimensions);
      int j;
      for(j = 0; j < dimensions; j++) p->coefs[count_m][j] = coefs[i][j];
      count_m++;
      
    //}
    free(coefs[i]);
  }

  p->num_coefs = p_m;
  free(coefs);
  free(weights);
  free(A);
  free(B);
  
  return p;
}

void PrintfPolyFunction(FILE* arq,poly_function_t* p){
//   fprintf(arq,"poly_struct_t\n");
//   fprintf(arq,"Dimentionality: %d Num Coefs: %d\n",p->dimensions,p->num_coefs);
  fprintf(arq,"%d %d\n",p->dimensions,p->num_coefs);
  int i;
  for(i = 0; i < p->num_coefs; i++){
//     rn_gen_print(arq,p->dimensions,p->coefs[i],"%3.4lf","(", ",", ")");
    int j;
    for(j = 0; j < p->dimensions; j++){
//     fprintf(arq," x %9.6e\n",p->weights[i]);
      fprintf(arq,"%9.6lf ",p->coefs[i][j]);
    }
    fprintf(arq,"   %9.6lf\n",p->weights[i]);
  }
    
}


poly_function_t* ReadPolyFunction(FILE* arq){
  poly_function_t* p = (poly_function_t*)malloc(sizeof(poly_function_t));
  /*Ignore the first header*/
  fscanf(arq,"%d %d",&(p->dimensions),&(p->num_coefs));
  p->coefs = (double**)malloc(sizeof(double*)*(p->num_coefs));
  p->weights = (double*)malloc(sizeof(double)*(p->num_coefs));
  int i,j;
  for(i = 0; i < p->num_coefs; i++){
    p->coefs[i] = (double*)malloc(sizeof(double)*(p->dimensions));
    for(j = 0; j< p->dimensions ; j++){
      double c;
      fscanf(arq,"%lf" ,&c);
      p->coefs[i][j] = c;
    }
    double w;
    fscanf(arq,"%lf",&w);
    p->weights[i] = w;
  }
  return p;
}
