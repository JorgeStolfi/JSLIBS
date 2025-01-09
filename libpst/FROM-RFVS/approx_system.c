#define _GNU_SOURCE
#include <approx_system.h>
#include <stdlib.h>
#include <assert.h>
#include <gauss_elim.h>
#include <rmxn.h>
#include <sve_minn.h>
#include <tabela.h>
#include <jsfile.h>
#include <lighting_models.h>
#include <lighting_harmonic.h>
#include <lighting_radial.h>
#include <lighting_stereopoly.h>
#include <lighting_compact.h>
#include <lighting_nearest.h>
#include <polynomial_functions.h>

model_options_t model_parse(argparser_t* pp, bool_t parseValues){
  
  model_options_t o;
  o.options = NULL;
  if (argparser_keyword_present_next(pp, "harmonic")){
    o.modelType = HARMONIC_MODEL;
    if(parseValues)  { o.options = lighting_harmonic_parse(pp); }
  }else if (argparser_keyword_present_next(pp, "radial")){
    o.modelType = RADIALBASIS_MODEL;
    if(parseValues) { o.options = lighting_radial_parse(pp); }
  }else if (argparser_keyword_present_next(pp, "stereopoly")){
    o.modelType = STEREOPOLY_MODEL;
    if(parseValues) { o.options = lighting_stereopoly_parse(pp); }
  }else if (argparser_keyword_present_next(pp, "compact")){
     o.modelType = COMPACT_MODEL;
     if(parseValues) { o.options = lighting_compact_parse(pp); }
  }else if (argparser_keyword_present_next(pp, "nearest")){
     o.modelType = NEAREST_MODEL;
     if(parseValues) { o.options = lighting_nearest_parse(pp); } 
  }else if (argparser_keyword_present_next(pp, "genericpoly")){
     o.modelType = GENERICPOLY_MODEL;
     if(parseValues) { o.options = poly_model_parse(pp); } 
  }else{
    fprintf(stderr,"Keyword was %s\n",argparser_get_next(pp));
    demand(FALSE, "unknown model");
  }
  
  return o;
}


void computeLeastSquaresMatrix(double* A,phi_function* phi, int basis_size, double** X,int n, double w[], double wpos[],void* l_data){
	int r,s;
	int count = 0;
// 	int total = basis_size;
	//fprintf(stderr,"\n\n");
	for( r = 0; r < basis_size; r++){
		for( s = 0; s < basis_size; s++){
			int iA = (r*basis_size) + s;
			double valueA = 0;
			int i;
			for(i = 0; i < n; i++){
			
				//r3_t normal = X[i];
				double* x = X[i];
                                // Change to view-depedent coordinates:
				double w_i = (w != NULL ? w[i] : 1);
				double w_pos = (wpos != NULL ? wpos[i] : 1);
// 				double wi = w[i]*wpos[i];
				double wi = w_i*w_pos;
// 				double phiR = phi(r,&normal,l_data);
// 				double phiS = phi(s,&normal,l_data);
				double phiR = phi(r,x,l_data);
 				double phiS = phi(s,x,l_data);
				valueA+= (phiR*phiS)*wi;
			}
			
			A[iA] = valueA;
		}
		
		//fprintf(stderr,"\033[1A");
		//fprintf(stderr,"Processed [%04d of %04d] - %4.3f%%\n",count,total,count*100.0/(float)total);
		count++;
	}
	
}

void computeLeastSquaresRHSVector(double* b,phi_function* phi, int basis_size,double** X,double* F,int n , double w[], double wpos[],void* l_data){
	int r;
	int count = 0;
// 	int total = basis_size;
	//fprintf(stderr,"\n\n");
	for( r = 0; r < basis_size; r++){
		
		double valueB = 0;
		int i;
		for(i = 0; i < n; i++){
			//r3_t normal = X[i];
			double* x = X[i];
                        // Change to view-depedent coordinates
			double w_i = (w != NULL ? w[i] : 1);
			double w_pos = (wpos != NULL ? wpos[i] : 1);
// 			double wi = w[i]*wpos[i];
			double wi = w_i*w_pos;
			double Di = F[i];
			
// 			double phiR = phi(r,&normal,l_data);
			double phiR = phi(r,x,l_data);
			valueB+=(phiR*Di)*wi;
		}
		b[r] = valueB;
		//fprintf(stderr,"\033[1A");
		//fprintf(stderr,"Processed [%04d of %04d] - %4.3f%%\n",count,total,count*100.0/(float)total);
		count++;
	}
	

}



void approx_model_fit_linear_model(double** X, double* F,double* w,double* wpos, int n, approx_model_t* am, void* l_data){
  int basis_size = am->get_num_components(l_data);
  double* A = rmxn_alloc(basis_size,basis_size);
  double* b = rn_alloc(basis_size);
  double* c = rn_alloc(basis_size);
  computeLeastSquaresMatrix(A,am->phi, basis_size, X,n,w, wpos, l_data);
  computeLeastSquaresRHSVector(b,am->phi,basis_size,X,F,n,w,wpos,l_data);
  double tiny = 1.0e-10;
  gsel_solve(basis_size,basis_size, A, 1, b, c, tiny);
  am->set_alphas(c,basis_size,l_data);
  free(A);
  free(b);
  free(c);
  
}

double approx_model_compute_S_star(
    approx_model_t* am,
    double** X,
    double* F,
    double* w,
    double* wpos,
    int n,
    double* parameters,
    void* l_data_original
){
  void* l_data = am->copy_data(l_data_original);
  int num_parms = am->get_num_nl_parameters(l_data);
  am->unpack_nl_parameters(parameters,num_parms,l_data);
  approx_model_fit_linear_model(X, F,w,wpos,n, am, l_data);
  /*Now we have to compute the square error*/
  double sum = 0;
  double sumw = 10e-5;
  int i;
  for(i = 0; i < n; i++){
    double w_i = (w != NULL ? w[i] : 1);
    double w_pos = (wpos != NULL ? wpos[i] : 1);
    double wi = w_i*w_pos;
//     double wi = w[i]*wpos[i];
    double f = am->evaluate(X[i],l_data);
    sum+= (F[i] - f)*(F[i] - f)*wi;
    sumw+= wi;
  }
  double q_error = sum/sumw;
  
  fprintf(stderr," E: %9.6lf   ",q_error);
  am->write_parameters(stderr,l_data);
  //free(weights);
  am->release_data(l_data);
  
  
  return q_error;
  
}


void plot_non_linear_goal_function(char* prefix,double** X, double* F,double* w,double* wpos, int n, approx_model_t* am, void* l_data)
{
  int num_parms = am->get_num_nl_parameters(l_data);
  double parameter[num_parms];
  double center[num_parms];
  double delta[num_parms];
  int nsteps = 20;
  
  
  am->pack_nl_parameters(l_data,parameter,num_parms);
  
  int i_par;
  for(i_par = 0; i_par < num_parms; i_par++){
    /*plots variation with respect to parameter i_par*/
    char* filename = NULL;
    char *filename = jsprintf("%s-P%02d.txt",prefix,i_par);
    FILE* arq = open_write(filename,TRUE);
    
    
    
    am->get_nl_centers_and_deltas(l_data,center,delta,num_parms);
    double v_old = parameter[i_par];
    void* l_data_test = am->copy_data(l_data);
    int j;
    for( j = -nsteps; j<= +nsteps; j++){
      parameter[i_par] = center[i_par]+delta[i_par]*(j/(double)nsteps);
      double f = approx_model_compute_S_star(am,X,F,w,wpos,n,parameter,l_data_test);
      fprintf(arq,"%15.9lf %15.9lf\n",parameter[i_par],f);
    }
    parameter[i_par] = v_old;
    am->release_data(l_data_test);
    fclose(arq);
    free(filename);
  }
  
}

void approx_model_fit_non_linear_model_simplex(double** X, double* F,double* w,double* wpos, int n,
				  approx_model_t* am,
				  void* l_data,
				  int update_steps,
				  double alpha, double beta, double gamma)
{
  
  assert((alpha > 1.0) );
  assert((beta < 1.0) && (beta > 0.0) );
  assert((gamma < 1.0) && (gamma > 0.0) );
  
  int num_parms = am->get_num_nl_parameters(l_data);
  int num_vertices = num_parms+1;
  double* S = rmxn_alloc(num_vertices,num_parms);
  int num_samples = (num_vertices*(num_vertices+1))/2;
  double* Fv = rn_alloc(num_samples);
  double parameters[num_parms];
  
  
  int num_iters = 0;
  do{
     num_parms = am->get_num_nl_parameters(l_data);
     am->generate_nl_parameter_simplex(l_data, num_parms, S);
     int i,j,k;
     k = 0;
     for(i = 0; i < num_vertices; i++){
       double* Si = &(S[i*num_parms]);
       for(j = 0; j <= i; j++){
	 double* Sj = &(S[j*num_parms]);
	 int r;
	 for(r = 0; r < num_parms; r++){
	   parameters[r] = (i == j ? Si[r] : (Si[r] + Sj[r])/2.0);
	 }
	 fprintf(stderr,"[%d][%d] -",i,j);
	 Fv[k] = approx_model_compute_S_star(am,X,F,w,wpos,n,parameters,l_data);
	 k++;
       }
     }
     
     double baricentric[num_vertices];

     sve_minn_step(num_parms,Fv, baricentric);
     rmxn_map_row (num_vertices, num_parms, baricentric,S, parameters);
     am->update_nl_parameters_and_errors(parameters,num_parms, l_data,alpha,beta,gamma);
     fprintf(stderr,"Baricenter\n");
     (void) approx_model_compute_S_star(am,X,F,w,wpos,n,parameters,l_data);
//      am->write_parameters(stderr,l_data);
     num_iters++;
  }while( (num_iters < update_steps) && (num_parms > 0));
  if(num_parms == 0){
    fprintf(stderr,"!! ALL NL parameters converged !!\n");
  }
}

	
