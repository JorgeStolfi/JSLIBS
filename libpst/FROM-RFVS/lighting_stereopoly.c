#include <lighting_stereopoly.h>
#include <lighting_models.h>
#include <assert.h>

lighting_stereopoly_options_t* lighting_stereopoly_parse(argparser_t* pp){
  lighting_stereopoly_options_t* o = (lighting_stereopoly_options_t*) malloc(sizeof(lighting_stereopoly_options_t));
  fprintf(stderr,"  STEREOPOLY_MODEL\n");
  argparser_get_keyword_next(pp, "degree");
  o->degree = (int)argparser_get_next_int(pp, 1, 10000);
  fprintf(stderr,"  degree %d\n",o->degree);
  return o;
}

approx_model_t* lighting_stereopoly_create_approx_lighting_model(void){
  approx_model_t* am = (approx_model_t*)malloc(sizeof(approx_model_t));
  am->type = STEREOPOLY_MODEL;

  am->phi = lighting_stereopoly_phi;
  am->get_num_components = lighting_stereopoly_get_number_components;
  am->set_alphas = lighting_stereopoly_set_alphas;
  am->get_alphas = lighting_stereopoly_get_alphas;
  am->copy_data = lighting_stereopoly_copy_lighting_data;
  am->release_data =  lighting_stereopoly_release_lighting_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  am->evaluate = lighting_stereopoly_shading; 
  am->write_parameters = lighting_stereopoly_write_parameters;
  am->read_parameters = lighting_stereopoly_read_parameters; 
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

/*Stereographic projection of normals based on polynoms*/

lighting_stereopoly_data_t* lighting_stereopoly_init_components(lighting_stereopoly_options_t* o,r3_t view_dir){
  
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*)malloc(sizeof(lighting_stereopoly_data_t));
  int g = o->degree;
  int num_comp = (g+1)*(g+2)/2;
  pl->num_comp = num_comp;
  pl->view_dir = view_dir;
  pl->coeficients = (double*)malloc(sizeof(double)*num_comp);
  pl->tuples = (int**)malloc(sizeof(int*)*num_comp);
  
  int k,i,j;
  for(i = 0; i < num_comp; i++){
    pl->coeficients[i] = 1.0;
    pl->tuples[i] = (int*)malloc(sizeof(int)*2);
  }
  
  int count = 0;
  
  
  for(k = 0; k <= g; k++){
    for(i = k; i >= 0; i--){
      j = k - i ;
      pl->tuples[count][0] = i;
      pl->tuples[count][1] = j;
      count++;
    }
  }

  pl->num_comp = count;
  
  
  return pl;
}

double lighting_stereopoly_phi(int r, double* x, void* l_data){
  /*Stereographic coordinates*/
  double X,Y;
  X = x[0]/(x[2] + 1.0);
  Y = x[1]/(x[2] + 1.0);
  
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  int ri = pl->tuples[r][0];
  int rj = pl->tuples[r][1];
  double phiR = (pow(X,ri))*(pow(Y,rj));
  //phiR = phiR*pl->weights[r];
  return phiR;
}

int lighting_stereopoly_get_number_components(void* l_data){
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  return pl->num_comp;
}


void lighting_stereopoly_set_alphas(double* C, int n,void* l_data){
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  int i;
  for(i = 0; i < pl->num_comp; i++){
    pl->coeficients[i] = C[i];
  }
  
}

void lighting_stereopoly_get_alphas(void* l_data,double* C,int n){
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  assert(pl->num_comp == n);
  int i;
  for(i = 0; i < pl->num_comp; i++){
    C[i] = pl->coeficients[i];
  }
  
}

void* lighting_stereopoly_copy_lighting_data(void* l_data){
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  lighting_stereopoly_data_t* npl = (lighting_stereopoly_data_t*) malloc(sizeof(lighting_stereopoly_data_t));
  npl->num_comp = pl->num_comp;
  npl->tuples = (int**)malloc(sizeof(int*)*(npl->num_comp));
  npl->coeficients  = (double*)malloc(sizeof(double)*(npl->num_comp));
  int i;
  for(i  = 0; i < npl->num_comp;i++){
    npl->coeficients[i] = pl->coeficients[i];
    npl->tuples[i] = (int*)malloc(sizeof(int)*2);
    int j;
    for(j = 0; j < 2;j++){
      npl->tuples[i][j] = pl->tuples[i][j]; 
    }
  }
  return npl;
}

void lighting_stereopoly_release_lighting_data(void* l_data){
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  int i;
  for(i  = 0; i < pl->num_comp;i++){
    free(pl->tuples[i]);
  }
  free(pl->tuples);
  free(pl->coeficients);
  free(l_data);
}



double lighting_stereopoly_shading(double* x,void* l_data){
  
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  int r;
  double val = 0;
  for(r = 0; r < pl->num_comp; r++){
    val+= (pl->coeficients[r]*lighting_stereopoly_phi(r,x,l_data));
  }
  
  return val;
}

void lighting_stereopoly_write_parameters(FILE* arq,void* l_data){
  
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) l_data;
  fprintf(arq," %d\n",pl->num_comp);
  int i;
  for(i = 0; i < pl->num_comp;i++){
    fprintf(arq, "%8.6lf %5d %5d\n",pl->coeficients[i],pl->tuples[i][0],pl->tuples[i][1] );
  }
  
}

void* lighting_stereopoly_read_parameters(FILE* arq){
  
  lighting_stereopoly_data_t* pl = (lighting_stereopoly_data_t*) malloc(sizeof(lighting_stereopoly_data_t));
  fscanf(arq,"%d",&(pl->num_comp));
  int i;
  int num_comp = pl->num_comp;
  pl->coeficients = (double*)malloc(sizeof(double)*num_comp);
  pl->tuples = (int**)malloc(sizeof(int*)*num_comp);
  for(i = 0; i < pl->num_comp;i++){
    pl->tuples[i] = (int*)malloc(sizeof(int)*2);
    fscanf(arq, "%lf %d %d",&(pl->coeficients[i]),&(pl->tuples[i][0]),&(pl->tuples[i][1]));
  }
  
  return pl;
  
}


