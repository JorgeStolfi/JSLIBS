#include <lighting_radial.h>
#include <lighting_models.h>
#include <tabela.h>
#include <assert.h>
/*Radial Basis Model*/

lighting_radial_options_t* lighting_radial_parse(argparser_t* pp){
  lighting_radial_options_t* o = (lighting_radial_options_t*)malloc(sizeof(lighting_radial_options_t));
  fprintf(stderr,"  RADIALBASIS_MODEL\n");
  argparser_get_keyword_next(pp, "resolution");
  o->resolution = argparser_get_next_int(pp, 1, 10000);
  argparser_get_keyword_next(pp, "span");
  o->span = argparser_get_next_double(pp, 0.0, 10000.0);
  fprintf(stderr,"  resolution %d\n  span %4.4lf\n",o->resolution,o->span);
  return o;
}

approx_model_t* lighting_radial_create_approx_lighting_model(void){
  approx_model_t* am = (approx_model_t*)malloc(sizeof(approx_model_t));
  am->type = RADIALBASIS_MODEL;

  am->phi = lighting_radial_phi;
  am->get_num_components = lighting_radial_get_number_components;
  am->set_alphas = lighting_radial_set_alphas;
  am->get_alphas = lighting_radial_get_alphas;
  am->copy_data = lighting_radial_copy_lighting_data;
  am->release_data =  lighting_radial_release_lighting_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  am->evaluate = lighting_radial_shading; 
  am->write_parameters = lighting_radial_write_parameters;
  am->read_parameters = lighting_radial_read_parameters; 
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

lighting_radial_data_t* lighting_radial_init_components(lighting_radial_options_t* o,r2_t center,double raio, r2_t estica,r3_t pole,double thetaMax){
  	r2_t* pontos;
	r3_t* normals;
	int num_comp;
	
	gera_pontos_no_gabarito_eliptico(o->resolution,thetaMax,center,raio,estica,&pontos, &normals,&num_comp,pole);
//	gera_normais_no_gabarito(resolution,center, num_comp,pole);
	lighting_radial_data_t* rl = (lighting_radial_data_t*)malloc(sizeof(lighting_radial_data_t));
	rl->num_comp = num_comp;
	rl->bradius = (double*)malloc(sizeof(double)*num_comp);
	rl->coeficients = (double*)malloc(sizeof(double)*num_comp);
	rl->bcenter = normals;
	rl->view_dir = pole;
	
	int i;
	for(i = 0; i < num_comp; i++){
		r3_t v = normals[i];
		double d = r3_dot(&pole,&v);
		rl->bradius[i] = 3.0*o->span*thetaMax/(o->resolution*(2.0 + d));
		rl->coeficients[i] = 1.0;
		// fprintf(stderr, "  %3d ", i);
		//r3_print(stderr, &((*bcenter)[i]));
		//fprintf(stderr, "  %9.6f\n", br[i]);
	}
	
	free(pontos);
	
	
		
	return rl;
}

double lighting_radial_phi(int r, double* x, void* l_data){
  r3_t normal = (r3_t) {{x[0],x[1],x[2]}};
  lighting_radial_data_t* rl = (lighting_radial_data_t*) l_data;
  r3_t bcenter = rl->bcenter[r];
  double bradius = rl->bradius[r];
  double dot = r3_dot(&normal,&bcenter);
  double d2 = (1.0 - dot)/(bradius*bradius);
  if(d2 >= 1.0) return 0.0;
  return (1 - d2)*(1 - d2);
  
}

int lighting_radial_get_number_components(void* l_data){
  lighting_radial_data_t* rl = (lighting_radial_data_t*) l_data;
  return rl->num_comp;
}


void lighting_radial_set_alphas(double* C,int n, void* l_data){
  
  lighting_radial_data_t* rl = (lighting_radial_data_t*) l_data;
  int i;
  for(i = 0; i < rl->num_comp; i++){
    rl->coeficients[i] = C[i];
  }
}


void lighting_radial_get_alphas(void* l_data,double* C,int n){
  lighting_radial_data_t* rl = (lighting_radial_data_t*) l_data;
  assert(rl->num_comp == n);
  int i;
  for(i = 0; i < rl->num_comp; i++){
    C[i] = rl->coeficients[i];
  }
}

void* lighting_radial_copy_lighting_data(void* l_data){
  lighting_radial_data_t* pl = (lighting_radial_data_t*) l_data;
  lighting_radial_data_t* npl = (lighting_radial_data_t*) malloc(sizeof(lighting_radial_data_t));
  npl->num_comp = pl->num_comp;
  npl->bcenter= (r3_t*)malloc(sizeof(r3_t)*(npl->num_comp));
  npl->coeficients  = (double*)malloc(sizeof(double)*(npl->num_comp));
  npl->bradius  = (double*)malloc(sizeof(double)*(npl->num_comp));
  int i;
  for(i  = 0; i < npl->num_comp;i++){
    npl->bradius[i] = pl->bradius[i];
    npl->bcenter[i] = pl->bcenter[i];
    npl->coeficients[i] = pl->coeficients[i];
  }
  return npl;
}

void lighting_radial_release_lighting_data(void* l_data){
  lighting_radial_data_t* pl = (lighting_radial_data_t*) l_data;
  free(pl->bradius);
  free(pl->bcenter);
  free(pl->coeficients);
  free(l_data);
}

double lighting_radial_shading(double* x,void* l_data){

  lighting_radial_data_t* rl = (lighting_radial_data_t*) l_data;
  int r;
  double val = 0;
  for(r = 0; r < rl->num_comp; r++){
    val+= (rl->coeficients[r]*lighting_radial_phi(r,x,l_data));
  }
  
  return val;
}


void lighting_radial_write_parameters(FILE* arq,void* l_data){
  
  lighting_radial_data_t* rl = (lighting_radial_data_t*) l_data;
  fprintf(arq," %d\n",rl->num_comp);
  int i;
  for(i = 0; i < rl->num_comp;i++){
    fprintf(arq,"%8.6lf %8.6lf",rl->coeficients[i],rl->bradius[i]); 
    r3_gen_print(arq, &(rl->bcenter[i]), "%9.6lf", "  ", " ", "");
    fprintf(arq,"\n");
  }
  
}

void* lighting_radial_read_parameters(FILE* arq){
  
  lighting_radial_data_t* rl = (lighting_radial_data_t*)malloc(sizeof(lighting_radial_data_t));
  fscanf(arq,"%d",&(rl->num_comp));
  int i;
  int num_comp = rl->num_comp;
  rl->bradius = (double*)malloc(sizeof(double)*num_comp);
  rl->coeficients = (double*)malloc(sizeof(double)*num_comp);
  rl->bcenter = (r3_t*)malloc(sizeof(r3_t)*num_comp);
  
  for(i = 0; i < rl->num_comp;i++){
    fscanf(arq,"%lf %lf",&(rl->coeficients[i]),&(rl->bradius[i])); 
    //r3_gen_print(arq, &(rl->bcenter[i]), "%9.6lf", "  ", " ", "");
    fscanf(arq,"%lf %lf %lf",&(rl->bcenter[i].c[0]),&(rl->bcenter[i].c[1]),&(rl->bcenter[i].c[2])); 
//     fprintf(arq,"\n");
  }
  return rl;
  
}

