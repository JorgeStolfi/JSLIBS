/*Callback functions for the polinomial interpolation !*/
#include <polynomial_functions.h>
#include <least_squares_nd.h>
#include <affirm.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


int NoValidation(double* c, int basis_size, bool_t* validity_array);
// int NoValidation(double* c, int basis_size, bool_t* validity_array){
//   return basis_size;
// }

poly_function_options_t* poly_model_parse(argparser_t* pp){
  
  poly_function_options_t* o = (poly_function_options_t*)malloc(sizeof(poly_function_options_t));
  fprintf(stderr,"  GENERICPOLY_MODEL\n");
  
  argparser_get_keyword_next(pp, "dimentions");
  o->dimentions = (int)argparser_get_next_int(pp, 1, 10000);
  fprintf(stderr,"  dimentions %d\n",o->dimentions);
 
  argparser_get_keyword_next(pp, "degree");
  o->degree = (int)argparser_get_next_int(pp, 1, 10000);
  fprintf(stderr,"  degree %d\n",o->degree);
  
  o->homogeneous = argparser_keyword_present_next(pp,"homogeneous");
  
  return o;
}


approx_model_t* create_ls_polynomial_model(void){
  
  approx_model_t* lm = (approx_model_t*) malloc(sizeof(approx_model_t));
  
  lm->type = GENERICPOLY_MODEL;

  lm->phi = poly_phi;
  lm->get_num_components = poly_get_number_components;
  lm->set_alphas = poly_set_alphas;
  lm->get_alphas = poly_get_alphas;
  lm->copy_data = poly_copy_data;
  lm->release_data =  poly_release_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  lm->evaluate = poly_evaluate; 
  lm->write_parameters = poly_write_parameters;
  lm->read_parameters = poly_read_parameters; 
  /*for non linear models - NOT USED*/
  lm->compare = NULL; 
  lm->get_num_nl_parameters = NULL;
  lm->pack_nl_parameters =  NULL; 
  lm->unpack_nl_parameters = NULL;
  lm->update_nl_parameters_and_errors = NULL;
  lm->generate_nl_parameter_simplex = NULL; 
  
  return lm;
  
}

poly_function_t* poly_model_init_components(poly_function_options_t o){
  poly_function_t* pl = poly_init_components(o.dimentions,o.degree,o.homogeneous);
  return pl;
}


poly_function_t* poly_init_components(int dimensions,int degree,bool_t homogeneous){
  poly_function_t* pl = (poly_function_t*)malloc(sizeof(poly_function_t));
  pl->dimensions = dimensions;
  
  if(!homogeneous){
    pl->coefs = CreateCoeficientsVector(degree,pl->dimensions,&(pl->num_coefs));
  }else{
    pl->coefs = CreateCoeficientsVectorHomogeneous(degree,pl->dimensions,&(pl->num_coefs));
  }
  pl->weights = (double*)malloc(sizeof(double)*(pl->num_coefs));
  int i;
  for(i = 0; i < pl->num_coefs; i++){
     pl->weights[i] = 1;
  }
  return pl;
}

double poly_phi(int r, double* x,void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  demand(r < pl->num_coefs, "poly_phi: Invalid coeficient index !");
  double* coefs = pl->coefs[r];
  return EvaluateXE(pl->dimensions, x, coefs);
}

int poly_get_number_components(void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  return pl->num_coefs;
}


void poly_get_alphas(void* l_data,double* C, int n){
  poly_function_t* pl = (poly_function_t*)l_data;
  int i ;
  demand(n <= pl->num_coefs,"poly_retrieve_components: Invalid number of coeficients !");
  for(i = 0; i < pl->num_coefs; i++){
    if(i < n){
      C[i] = pl->weights[i] ;
    }else{
      C[i] = 0;
    }
  }
}


void poly_set_alphas(double* C, int n,void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  int i ;
  demand(n <= pl->num_coefs,"poly_retrieve_components: Invalid number of coeficients !");
  for(i = 0; i < pl->num_coefs; i++){
    if(i < n){
      pl->weights[i] = C[i];
    }else{
      pl->weights[i] = 0;
    }
  }
}

double poly_evaluate(double* x,void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  int dimensions = pl->dimensions;
  
  int degree = 0;
  int i;
  for(i = 0; i < pl->num_coefs; i++){
    int j;
    for(j  = 0; j < pl->dimensions; j++){
      if(degree < pl->coefs[i][j]) { degree = (int)pl->coefs[i][j]; }
    }
  }
  
  double degree_val[(degree+1)*(dimensions)];

    

  for(i = 0; i <= degree; i++){
    int j;
    for(j  = 0; j < pl->dimensions; j++){
      degree_val[i*dimensions + j] = (i == 0 ? 1.0 : degree_val[((i-1)*dimensions) + j]*x[j]);
    }
  }
  
  double val = 0;
  for(i = 0; i < pl->num_coefs; i++){
    double val_alpha = 1.0;
    int j;
    for(j  = 0; j < pl->dimensions; j++){
      double ind = pl->coefs[i][j];
      val_alpha*= degree_val[ind*dimensions+ j];
    }
     double* coef = pl->coefs[i];
//      double val_alpha2 = EvaluateXE( dimensions, x,coef);
//     if( abs(val_alpha - val_alpha2) > 10e-10){
//       fprintf(stderr,"VALUE MISMATCH\nVal1 = %e\nVal2 = %e\n",val_alpha,val_alpha2);
//     }
    val+= pl->weights[i]*val_alpha;
  }
    
//   return EvaluatePvalue(pl,x);
//   double val2 = EvaluatePvalue(pl,x);
//   if( abs(val2 - val) > 10e-100){
//     for(i = 0; i < degree; i++){
//       int j;
//       fprintf(stderr,"D[%02d] ",i);
//       for(j  = 0; j < pl->dimensions; j++){
// 	fprintf(stderr,"%e vs %e ",degree_val[i*dimensions+ j],pow(x[j],i));
//       }
//       fprintf(stderr,"\n");
//     }
//     fprintf(stderr,"VALUE MISMATCH\nVal1 = %12.7lf\nVal2 = %12.7lf\n",val,val2);
//     assert(FALSE);
//   }
  return val;

}

void* poly_read_parameters(FILE* arq){
  poly_function_t* pl = ReadPolyFunction(arq);
  return pl;
}


void poly_write_parameters(FILE* arq,void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  PrintfPolyFunction(arq ,pl);
}



void* poly_copy_data(void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  poly_function_t* pl_new = (poly_function_t*)malloc(sizeof(poly_function_t));
  pl_new->dimensions = pl->dimensions;
  pl_new->num_coefs = pl->num_coefs;
  pl_new->weights = (double*)malloc(sizeof(double)*(pl->num_coefs));
  int i;
  pl_new->coefs = (double**)malloc(sizeof(double*)*(pl->num_coefs));
  for(i = 0; i < pl->num_coefs;i++){
    pl_new->coefs[i] = (double*)malloc(sizeof(double)*(pl->dimensions));
    int j;
    for(j = 0;j < pl->dimensions; j++){
      pl_new->coefs[i][j] = pl->coefs[i][j];
    }
    pl_new->weights[i] = pl->weights[i];
  }
  
  return pl_new;
}

void poly_release_data(void* l_data){
  poly_function_t* pl = (poly_function_t*)l_data;
  int i;
  for(i = 0; pl->num_coefs; i++){
    free(pl->coefs[i]);
  }
  free(pl->weights);
  free(pl->coefs);
  free(pl);
}
