#define PROG_NAME "polinomial_interpolator"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <argparser.h>
#include <least_squares_nd.h>
#include <jsfile.h>
#include <rn.h>
#include <rmxn.h>
#include <gausol_solve.h>
#include <least_squares_nd.h>
#include <math.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "	[-in {INPUT IMAGE}] \\\n" \
  "	[-out {INPUT IMAGE}] \\\n" \
  "	-degree {N} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  Etc. etc..\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "OPTIONS" \
  "  Etc. etc.."


struct options_t{
   char* infile;
   char* outfile;
   int degree;
   bool_t invert;
   bool_t evaluate;
   char* evaluateinfile;
   char* evaluateoutfile;
};

typedef struct options_t options_t;

options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  o->infile = "-";
  if(argparser_keyword_present(pp, "-in")){
    o->infile = argparser_get_next(pp);
  }
  
  o->outfile = "-";
  if(argparser_keyword_present(pp, "-out")){
    o->outfile = argparser_get_next(pp);
  }
  
  argparser_get_keyword(pp, "-degree");
  o->degree = argparser_get_next_int(pp, 1, 100000);
  
  o->invert = argparser_keyword_present(pp, "-invert");
  
  o->evaluate = argparser_keyword_present(pp, "-evaluate");
  if(o->evaluate){
     o->evaluateinfile = "-";
    if(argparser_keyword_present_next(pp,"infile")){
      o->evaluateinfile = argparser_get_next(pp);
    }
    o->evaluateoutfile = "-";
    if(argparser_keyword_present_next(pp,"outfile")){
      o->evaluateoutfile = argparser_get_next(pp);
    }
  }
  
  argparser_finish(pp);
  return o;
}

void WritePolyData(FILE* arq,poly_function_t* polym){
 	
  fprintf(arq,"%d ",(polym->dimensions));
  fprintf(arq,"%d\n",(polym->num_coefs));
  int icf;
  for(icf = 0; icf < (polym->num_coefs);icf++){
    int d;
    for(d = 0; d < (polym->dimensions); d++){
      fprintf(arq,"%3.4lf ",(polym->coefs[icf][d]));
    }
    fprintf(arq,"%9.6e\n",(polym->weights[icf]));
  }
}

void ReadData(FILE* arq,double*** Xp, double** Fp, int* Np, int* dp,  bool_t invert){
  
    int N;
    int d;
    double** X;
    double* F;
    
    fscanf(arq,"%d %d",&d,&N);
    X = (double**)malloc(sizeof(double*)*N);
    F = (double*)malloc(sizeof(double)*N);
    int i;
    for(i = 0; i < N; i++){
      X[i] = (double*)malloc(sizeof(double)*d);
      int j;
      for(j = 0; j < d; j++){
	fscanf(arq,"%lf",&(X[i][j]));
      }
      fscanf(arq,"%lf",&(F[i]));
      if(invert) F[i] =  1/F[i];
    }
    
   *dp = d;
   *Xp = X;
   *Fp = F;
   *Np = N;
}



void ReadEvaluateInterval(FILE* arq,double* sta_interv, double* end_interv, double* spacing,int dimensions){
  int i;
  for(i = 0; i < dimensions;i++){
    fscanf(arq,"%lf %lf %lf",&(sta_interv[i]),&(end_interv[i]),&(spacing[i]));
    assert(spacing[i] > 0);
  }  
}

void evaluate_recursive(FILE* arq,
			double* index,
			double* sta_interv,
			double* end_interv,
			double* spacing,
			int pos,
			poly_function_t* pol){
  
  if(pos == pol->dimensions){
    int i;
    for(i = 0; i < pol->dimensions; i++){
      fprintf(arq,"%4.0lf ",index[i]);
    }
    double val = EvaluatePvalue(pol,index);
    fprintf(arq,"%lf\n",val);
    return ;
  }else{
    double ind = sta_interv[pos];
    
    while(ind <= end_interv[pos]){
      index[pos] = ind;
      evaluate_recursive(arq,index,sta_interv,end_interv,spacing,pos+1,pol);
      ind+= spacing[pos];
    }

  }
}

void evaluate_recursive_FNI(FILE* arq,
			double* index,
			double* sta_interv,
			double* end_interv,
			double* spacing,
			int pos,
			poly_function_t* pol){
  
  if(pos < 0){
    int i;
    for(i = 0; i < pol->dimensions; i++){
      fprintf(arq,"%4.0lf ",index[i]);
    }
    double val = EvaluatePvalue(pol,index);
    fprintf(arq,"%lf\n",val);
    return;
  }else{
    double ind = sta_interv[pos];
    while(ind <= end_interv[pos]){
      index[pos] = ind;
      evaluate_recursive_FNI(arq,index,sta_interv,end_interv,spacing,pos-1,pol);
      ind+=spacing[pos];
    }
    if(pos != pol->dimensions -1) fprintf(arq,"\n");
  }
  
  
}

void WriteEvaluateData(FILE* arq,double* sta_interv, double* end_interv, double* spacing,poly_function_t* pol,bool_t FNI_FORMAT){
   double index[pol->dimensions];
   if(!FNI_FORMAT){
    evaluate_recursive(arq,index,sta_interv,end_interv,spacing,0,pol);
   }else{
     evaluate_recursive_FNI(arq,index,sta_interv,end_interv,spacing,pol->dimensions-1,pol);
   }
}


typedef double basis_function(int i, double* X) ;

void computeLeastSquaresMatrix(double* A,basis_function* phi, int basis_size, double** X, int N);
/*
Given and basis function phi, compute and writes inside {A}, the LS matrix using the {X} set
*/

void computeLeastSquaresMatrix(double* A,basis_function* phi, int basis_size, double** X, int N){
	int r,s;
	
	int count = 0;
	//fprintf(stderr,"\n\n");
	for( r = 0; r < basis_size; r++){
		for( s = 0; s < basis_size; s++){
			int iA = (r*basis_size) + s;
			double valueA = 0;
			int i;
			for(i = 0; i < N; i++){
				double* Xi = X[i];
				double phiR = phi(r,Xi);
				double phiS = phi(s,Xi);
				valueA+= (phiR*phiS);
			}
			
			A[iA] = valueA;
		}
		//fprintf(stderr,"\033[1A");
		//fprintf(stderr,"Processed [%04d of %04d] - %4.3f%%\n",count,total,count*100.0/(float)total);
		count++;
	}
	
}

void computeLeastSquaresRHSVector(double* b,basis_function* phi, int basis_size,double** X,double* F, int N );
/*Computes the RHS vector for a given function phi and the F(X)*/

void computeLeastSquaresRHSVector(double* b,basis_function* phi, int basis_size,double** X,double* F, int N ){
	int r;
	int count = 0;
	//fprintf(stderr,"\n\n");
	for( r = 0; r < basis_size; r++){
		double valueB = 0;
		int i;
		for(i = 0; i < N; i++){
			double* Xi = X[i];
			double Di = F[i];
			double phiR = phi(r,Xi);
			valueB+=(phiR*Di);
		}
		b[r] = valueB;
	//	fprintf(stderr,"\033[1A");
	//	fprintf(stderr,"Processed [%04d of %04d] - %4.3f%%\n",count,total,count*100.0/(float)total);
		count++;
	}
	

}



poly_function_t*  computePolinomialApproximation(double** X,double* F, int N,int dimensions,int degree){
    
  poly_function_t* p =  (poly_function_t*)malloc(sizeof(poly_function_t));
  
  
  auto double phiPol(int r,double* x);
  
  double phiPol(int r,double* x){
    int i;
    double val = 0;
    for(i = 0; i < p->dimensions; i++){
      val+= pow(x[i],p->coefs[r][i])*(p->weights[r]);
    }
   return val;
  }
  
  int basis_size;
  
  double** coefs = CreateCoeficientsVector(degree,dimensions,&basis_size);
  double* weights = (double*)malloc(sizeof(double)*basis_size);
  p->dimensions = dimensions;
  p->weights = weights;
  p->coefs = coefs;
  p->num_coefs = basis_size;
  int i;
  for(i = 0; i < basis_size; i++){
    weights[i] = 1;
  }
  
  double* A = rmxn_alloc(basis_size,basis_size);
  double* B = rmxn_alloc(basis_size,1);
  double* C = rmxn_alloc(basis_size,1);
  
  double epsilon = 10e-6;
  double diff = 100*epsilon;
  
  int iterations = 0;
  int MAX_ITERATIONS = 10;
  
  double*  last_weights = (double*)malloc(sizeof(double)*basis_size);
  
  fprintf(stderr,"Starting \n");
  do{
    
    rn_copy(basis_size,weights,last_weights);
    PrintfPolyFunction(stderr,p);
    computeLeastSquaresMatrix(A,phiPol,basis_size, X, N);
    fprintf(stderr,"A\n");
    rmxn_gen_print( stderr, basis_size,basis_size,A,"%3.4lf",NULL,NULL,NULL,NULL,NULL,NULL); 

    computeLeastSquaresRHSVector(B,phiPol,basis_size,X,F,N);
    fprintf(stderr,"B\n");
    rmxn_gen_print( stderr, basis_size,1,B,"%3.4lf",NULL,NULL,NULL,NULL,NULL,NULL); 
    // AC = B;
   // double* tecs =  SolveSystemBeta(A,B,basis_size);
    gausol_solve(basis_size, basis_size, A , 1, B, C, TRUE,TRUE, 1.0e-14, NULL,NULL);
 //   rn_copy(basis_size,tecs,C);
   // free(tecs);
    rn_copy(basis_size,C,weights);
    diff = rn_dist(basis_size,weights,last_weights);
    iterations++;
    
    fprintf(stderr,"[%d]: %lf\n",iterations,diff);
  }while((diff > epsilon) && (iterations < MAX_ITERATIONS));
  
  free(last_weights);
  free(A);
  free(B);
  free(C);
  
  return p;
}



int main(int argc, char** argv){

  options_t* o = parse_args(argc,argv);
  FILE* infile = open_read(o->infile,TRUE);
  FILE* outfile = open_write(o->outfile,TRUE);
  double** X;
  double* F;
  int N;
  int d;
  ReadData(infile,&X,&F,&N,&d,o->invert);
  poly_function_t* p = LS_PolyFitting(d,N,o->degree,X,F);
 // poly_function_t* p = computePolinomialApproximation( X,F,N,d,o->degree);
  WritePolyData(outfile,p);
  if(o->evaluate){
    FILE* evalinfile = open_read(o->evaluateinfile,TRUE);
    FILE* evaloutfile = open_write(o->evaluateoutfile,TRUE);
    double sta_interv[d];
    double end_interv[d];
    double spacing[d];
    ReadEvaluateInterval(evalinfile,sta_interv,end_interv,spacing,d);
    WriteEvaluateData(evaloutfile,sta_interv,end_interv,spacing,p,TRUE);
    if(evalinfile != stdin) fclose(evalinfile);
    if(evaloutfile != stdout) fclose(evaloutfile);
  }
  fclose(infile);
  fclose(outfile);
  
  return 0;
}

