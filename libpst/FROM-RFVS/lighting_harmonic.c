#include <lighting_harmonic.h>
#include <lighting_models.h>
#include <assert.h>
/*Polynomial model*/

lighting_harmonic_options_t* lighting_harmonic_parse(argparser_t* pp){
  
  lighting_harmonic_options_t* o = (lighting_harmonic_options_t*)malloc(sizeof(lighting_harmonic_options_t));
  fprintf(stderr,"  HARMONIC_MODEL\n");
  argparser_get_keyword_next(pp, "degree");
  o->degree = (int)argparser_get_next_int(pp, 1, 10000);
  fprintf(stderr,"  degree %d\n",o->degree);
  return o;
}

approx_model_t* lighting_harmonic_create_approx_lighting_model(void){
  approx_model_t* am = (approx_model_t*)malloc(sizeof(approx_model_t));
  am->type = HARMONIC_MODEL;

  am->phi = lighting_harmonic_phi;
  am->get_num_components = lighting_harmonic_get_number_components;
  am->set_alphas = lighting_harmonic_set_alphas;
  am->get_alphas = lighting_harmonic_get_alphas;
  am->copy_data = lighting_harmonic_copy_lighting_data;
  am->release_data =  lighting_harmonic_release_lighting_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  am->evaluate = lighting_harmonic_shading; 
  am->write_parameters = lighting_harmonic_write_parameters;
  am->read_parameters = lighting_harmonic_read_parameters; 
  /*for non linear models - NOT USED*/
  am->compare = NULL; 
  am->get_num_nl_parameters = NULL;
  am->pack_nl_parameters =  NULL;
  am->get_nl_centers_and_deltas = NULL;
  am->unpack_nl_parameters = NULL;
  am->update_nl_parameters_and_errors = NULL;
  am->generate_nl_parameter_simplex = NULL; 
  return am;
}


lighting_harmonic_data_t* lighting_harmonic_init_components(lighting_harmonic_options_t* o,r3_t view_dir){
  
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*)malloc(sizeof(lighting_harmonic_data_t));
  int g = o->degree;
  int num_comp = (g+1)*(g+1);
  pl->num_comp = num_comp;
  pl->view_dir = view_dir;
  pl->coeficients = (double*)malloc(sizeof(double)*num_comp);
  pl->triplets = (int**)malloc(sizeof(int*)*num_comp);
  
  int i,j,k;
  for(i = 0; i < num_comp; i++){
    pl->coeficients[i] = 1.0;
    pl->triplets[i] = (int*)malloc(sizeof(int)*3);
  }
  
  int count = 0;
  for(i = g; i >= 0; i--){
    for( j = g-i; j >= 0; j--){
      k =  g - i - j;
      pl->triplets[count][0] = i;
      pl->triplets[count][1] = j;
      pl->triplets[count][2] = k;
      count++;
    }
  }
  g = g -1;
  for(i = g; i >= 0; i--){
    for( j = g-i; j >= 0; j--){
      k =  g - i - j;
      pl->triplets[count][0] = i;
      pl->triplets[count][1] = j;
      pl->triplets[count][2] = k;
      count++;
    }
  }
  pl->num_comp = count;
  
  
  return pl;
  
}

double lighting_harmonic_phi(int r, double* x, void* l_data){
  r3_t normal = (r3_t){{x[0],x[1],x[2]}};
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  int ri = pl->triplets[r][0];
  int rj = pl->triplets[r][1];
  int rk = pl->triplets[r][2];
  double phiR = (pow(normal.c[0],ri))*(pow(normal.c[1],rj))*(pow(normal.c[2],rk));
  //phiR = phiR*pl->weights[r];
  return phiR;
}

int lighting_harmonic_get_number_components(void* l_data){
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  return pl->num_comp;
}

void lighting_harmonic_set_alphas(double* C,int n, void* l_data){
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  
  int i;
  for(i = 0; i < pl->num_comp; i++){
    pl->coeficients[i] = C[i];
  }
}


void lighting_harmonic_get_alphas(void* l_data,double* C,int n){
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  assert(n == pl->num_comp);
  int i;
  for(i = 0; i < pl->num_comp; i++){
    C[i] = pl->coeficients[i];
  }
}

void* lighting_harmonic_copy_lighting_data(void* l_data){
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  lighting_harmonic_data_t* npl = (lighting_harmonic_data_t*) malloc(sizeof(lighting_harmonic_data_t));
  npl->num_comp = pl->num_comp;
  npl->triplets = (int**)malloc(sizeof(int*)*(npl->num_comp));
  npl->coeficients  = (double*)malloc(sizeof(double)*(npl->num_comp));
  npl->view_dir = pl->view_dir;
  int i;
  for(i  = 0; i < npl->num_comp;i++){
    npl->coeficients[i] = pl->coeficients[i];
    npl->triplets[i] = (int*)malloc(sizeof(int)*3);
    int j;
    for(j = 0; j < 3;j++){
      npl->triplets[i][j] = pl->triplets[i][j]; 
    }
  }
  return npl;
}


void lighting_harmonic_release_lighting_data(void* l_data){
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  int i;
  for(i  = 0; i < pl->num_comp;i++){
    free(pl->triplets[i]);
  }
  free(pl->triplets);
  free(pl->coeficients);
  free(l_data);
}

double lighting_harmonic_shading(double* x,void* l_data){
  
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  int r;
  double val = 0;
  for(r = 0; r < pl->num_comp; r++){
    val+= (pl->coeficients[r]*lighting_harmonic_phi(r,x,l_data));
  }
  
  return val;
}


void lighting_harmonic_write_parameters(FILE* arq,void* l_data){
  
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) l_data;
  fprintf(arq," %d\n",pl->num_comp);
  int i;
  for(i = 0; i < pl->num_comp;i++){
    fprintf(arq, "%8.6lf %5d %5d %5d\n",pl->coeficients[i],pl->triplets[i][0],pl->triplets[i][1],pl->triplets[i][2] );
  }
  
}

void* lighting_harmonic_read_parameters(FILE* arq){
  
  lighting_harmonic_data_t* pl = (lighting_harmonic_data_t*) malloc(sizeof(lighting_harmonic_data_t));
  fscanf(arq,"%d",&(pl->num_comp));
  int i;
  int num_comp = pl->num_comp;
  pl->coeficients = (double*)malloc(sizeof(double)*num_comp);
  pl->triplets = (int**)malloc(sizeof(int*)*num_comp);
  for(i = 0; i < pl->num_comp;i++){
    pl->triplets[i] = (int*)malloc(sizeof(int)*3);
    fscanf(arq, "%lf %d %d %d",&(pl->coeficients[i]),&(pl->triplets[i][0]),&(pl->triplets[i][1]),&(pl->triplets[i][2]) );
  }
  
  return pl;
  
}

