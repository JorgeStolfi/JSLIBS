#include <lighting_nearest.h>
#include <lighting_models.h>
#include <assert.h>
#include <rn.h>
/*Polynomial model*/

lighting_nearest_options_t* lighting_nearest_parse(argparser_t* pp){
  lighting_nearest_options_t* o = (lighting_nearest_options_t*)malloc(sizeof(lighting_nearest_options_t));
  fprintf(stderr,"  NEAREST_MODEL\n");
  o->gin = parse_gauge_args(pp,TRUE);
  o->channel = 0;
  if(argparser_keyword_present_next(pp, "channel")){
    o->channel = (int)argparser_get_next_int(pp,0,2);
  }
  return o;
 
}

approx_model_t* lighting_nearest_create_approx_lighting_model(void){
  approx_model_t* am = (approx_model_t*)malloc(sizeof(approx_model_t));
  am->type = NEAREST_MODEL;

  am->phi = lighting_nearest_phi;
  am->get_num_components = lighting_nearest_get_number_components;
  am->set_alphas = lighting_nearest_set_alphas;
  am->get_alphas = lighting_nearest_get_alphas;
  am->copy_data = lighting_nearest_copy_lighting_data;
  am->release_data =  lighting_nearest_release_lighting_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  am->evaluate = lighting_nearest_shading; 
  am->write_parameters = lighting_nearest_write_parameters;
  am->read_parameters = lighting_nearest_read_parameters; 
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


lighting_nearest_data_t* lighting_nearest_init_components(lighting_nearest_options_t* o,r3_t view_dir){
  
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*)malloc(sizeof(lighting_nearest_data_t));
  gauge_data_t* gin = o->gin;
  assert(! isnan(gin->E.rad));
  assert(! isnan(gin->E.ctr.c[0]));
  
  /* Get the input gauge image {img_in} and its extension {ext}: */
  float_image_t *img_in = read_gauge_image(gin->image, gin->gamma);
  
  int NC = (int)img_in->sz[0];
  assert(NC <= MAX_NC);
  assert(o->channel < NC);
  int NXin = (int)img_in->sz[1];
  int NYin = (int)img_in->sz[2];
  
  
  float_image_t *img_mask = float_image_new(1, NXin, NYin);
  double** F;
   /* Extract the normal and intensity data from the input image: */
  extract_data_from_gauge_image(img_in, img_mask,NULL, &(gin->E),gin->trim, &(gin->view), &(nl->normal), &F, &(nl->n));
  nl->val = (double*)malloc(sizeof(double)*(nl->n));
  rn_copy(nl->n,F[o->channel],nl->val);
  float_image_free(img_mask);
  float_image_free(img_in);
  int i;
  for(i = 0; i < NC; i++){
    free(F[i]);
  }
  free(F);
  return nl;
  
}

double lighting_nearest_phi(int r, double* x, void* l_data){
  assert(r == 0);
  
  r3_t normal = (r3_t){{x[0],x[1],x[2]}};
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
  int i;
  double dbest = +INF;
  double val = NAN;
  for(i = 0; i < nl->n; i++){
    double dist  = r3_dist(&normal,&(nl->normal[i]));
    if(dist < dbest){
      dbest = dist;
      val = nl->val[i];
    }
  }
  assert(!isnan(val));

  return val;
}

int lighting_nearest_get_number_components(void* l_data){
  return 1;
}

void lighting_nearest_set_alphas(double* C,int n, void* l_data){
  assert(n == 1);
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
  nl->alpha = C[0];
}


void lighting_nearest_get_alphas(void* l_data,double* C,int n){
  assert(n == 1);
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
  C[0] = nl->alpha;
}

void* lighting_nearest_copy_lighting_data(void* l_data){
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
  lighting_nearest_data_t* nnl = (lighting_nearest_data_t*) malloc(sizeof(lighting_nearest_data_t));
  nnl->alpha = nl->alpha;
  nnl->n = nl->n;
  int n = nl->n;
  nnl->val = (double*)malloc(sizeof(double)*n);
  nnl->normal = (r3_t*)malloc(sizeof(r3_t)*n);
  int i;
  for(i = 0; i < nl->n; i++){
    nnl->val[i] = nl->val[i];
    nnl->normal[i] = nl->normal[i];
  }
  
  return nnl;
}


void lighting_nearest_release_lighting_data(void* l_data){
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
  free(nl->val);
  free(nl->normal);
  free(l_data);
}

double lighting_nearest_shading(double* x,void* l_data){
  
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
     
  return lighting_nearest_phi(0,x,l_data)*nl->alpha;
}


void lighting_nearest_write_parameters(FILE* arq,void* l_data){
  
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) l_data;
  fprintf(arq,"%d\n",nl->n);
  fprintf(arq,"%8.6lf\n",nl->alpha);
  int i;
  for(i = 0; i < nl->n;i++){
    fprintf(arq, "%8.6lf %8.6lf %8.6lf %8.6lf\n",nl->normal[i].c[0],nl->normal[i].c[1],nl->normal[i].c[2],nl->val[i] );
  }
}

void* lighting_nearest_read_parameters(FILE* arq){
  
  lighting_nearest_data_t* nl = (lighting_nearest_data_t*) malloc(sizeof(lighting_nearest_data_t));
  fscanf(arq,"%d",&(nl->n));
  fscanf(arq,"%lf",&(nl->alpha));
  int n = nl->n;
  int i;
  nl->normal = (r3_t*)malloc(sizeof(r3_t)*n);
  nl->val = (double*)malloc(sizeof(double)*n);
  for(i = 0; i < n;i++){
    fscanf(arq, "%lf %lf %lf %lf",&(nl->normal[i].c[0]),&(nl->normal[i].c[1]),&(nl->normal[i].c[2]),&(nl->val[i]) );
  }
  return nl;
  
}

