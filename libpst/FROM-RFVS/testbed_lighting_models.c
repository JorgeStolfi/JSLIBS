#define PROG_NAME "testbed_lighting_models"

#include <stdio.h>
#include <stdlib.h>
#include <jsfile.h>
#include <math.h>
#include <normais.h>
#include <r3.h>
#include <rn.h>
#include <rmxn.h>
#include <ellipse_crs.h>
#include <assert.h>
#include <affirm.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <string.h>
#include <argparser.h>
#include <lighting_models.h>
#include <lighting_compact.h>
#include <lighting_harmonic.h>
#include <lighting_radial.h>
#include <lighting_stereopoly.h>
#include <lighting_nearest.h>
#include <jsrandom.h>
#include <tabela.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "      -inprefix {FILENAME} \\\n" \
  "           {pointlike,harmonic,radial,glossy,harmonicSG} \\\n" \
  "      -outprefix {FILENAME} \\\n" \
  "            [ [backplane] ] \\\n" \
  "            [harmonic degree {DEGREE} ] \\\n" \
  "            [radial  resolution {RES} span {SPAN} ] \\\n" \
  "            [glossy  glossiness {G)} {G1} ] \\\n" \
  "            [harmonicSG  degree {DEGREE}  ] \\\n" \
  "      -plotprefix {FILENAME} \\\n" \
  "      [ -noise {NOISE} ] \\\n" \
  "      [ -viewDir {VX} {VY} VZ} ]\\\n" \
  "      [ -nSteps {STEPS} ] \\\n" \
  "      [ -ignorePosWeight ] \\\n" \
  "      [ -resolution {RES} ] \\\n" \
  "      [ -plotSamples {PLOTSAMPLES} ] \\\n" \
  "      [ -generateSGplot {SGPLOTSAMPLES} ] \\\n" \
  "      [ -plotImage RX] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO ""

#define NOWEIGHT 0
#define COSWEIGHT 1

struct test_data_t{
  double** points;
  double* value;
  double* weights;
  double* wpos;
  int dimension;
  long int num_points;
};

typedef struct test_data_t test_data_t;

struct options_t{
  r3_t view_dir;
  r3_t ref_dir;
  
//   int nSteps;
 int nStepsU;
 int nStepsNL;
 bool_t estimateU;
 bool_t useNonLinear;
  
  double noise;
  
  char* inprefix;
  model_options_t in_lm_data;
   
  char* outprefix;
  model_options_t out_lm_data;
  
  char* plotprefix;
  
  int resolution; /*determine sampling resolution of the generated test set*/
  int plotSamples; /*Number of samples to be used in the plot*/
  bool_t generateSGplot;
  int plotSamplesSG;
  
  bool_t plotImage;
  int radiusImage;
  
  bool_t useCosPosWeight;
  
};

typedef struct options_t options_t;

test_data_t* generate_test_data(approx_model_t* lm_test,void* l_datatest, int resolution, r3_t view_dir,double noise, int type_weight);
test_data_t* generate_test_data(approx_model_t* lm_test,void* l_datatest, int resolution, r3_t view_dir,double noise, int type_weight){
  
  r2_t* points;
  r3_t* normals;
  int num_lines;
  
  r2_t gauge_center = (r2_t){{0,0}};
  r2_t gauge_stretch = (r2_t){{0,0}};
  double radius = 100;

  gera_pontos_no_gabarito_eliptico(resolution,M_PI/2.0, gauge_center,radius,gauge_stretch,&points,&normals,&num_lines,view_dir);
  
  test_data_t* t = (test_data_t*)malloc(sizeof(test_data_t));
  t->points = (double**)malloc(sizeof(double*)*num_lines);
  t->value  = (double*)malloc(sizeof(double)*num_lines);
  t->weights = (double*)malloc(sizeof(double)*num_lines);
  t->wpos = (double*)malloc(sizeof(double)*num_lines);
  t->dimension = 3;
  t->num_points = num_lines;
  
  long int i;
  for(i = 0; i < num_lines; i++){
    t->points[i] = (double*)malloc(sizeof(double)*(t->dimension));
    t->points[i][0] = normals[i].c[0];
    t->points[i][1] = normals[i].c[1];
    t->points[i][2] = normals[i].c[2];
    t->value[i] = lm_test->evaluate(normals[i].c,l_datatest);
    t->value[i]+= noise * dgaussrand();
  }
  
  for(i = 0; i < num_lines; i++){
    t->weights[i] = 1;
    if(type_weight == COSWEIGHT){
      double peso_cos = r3_dot(&view_dir,&(normals[i]));
      t->wpos[i] = (peso_cos < 0 ? 0: peso_cos);
    }else{
      t->wpos[i] = 1;
    }
  }
    
    
  free(points);
  free(normals);
  
  return t;
}


void* approximate_model(test_data_t* t, approx_model_t* lm_approx, model_options_t lmodel, r3_t view_dir, bool_t estimateU, int nStepsU, bool_t useNonLinear, int nStepsNL );
void* approximate_model(test_data_t* t, approx_model_t* lm_approx, model_options_t lmodel, r3_t view_dir, bool_t estimateU, int nStepsU, bool_t useNonLinear, int nStepsNL ){

  void* approx_data;
  
   if(lmodel.modelType == COMPACT_MODEL){
    lighting_compact_options_t* o = lmodel.options;
    approx_data = lighting_compact_init_components(o,view_dir);
    lm_approx->write_parameters(stderr,approx_data);
  }else if(lmodel.modelType == HARMONIC_MODEL){
    lighting_harmonic_options_t* o = lmodel.options;
    approx_data = lighting_harmonic_init_components(o,view_dir);
  }else if(lmodel.modelType == RADIALBASIS_MODEL){
    lighting_radial_options_t* o = lmodel.options;
//     approx_data = radialbasis_init_components(lmodel.radialResolution,lmodel.radialSpan,(r2_t){{0,0}},100,(r2_t){{0,0}},view_dir);
    approx_data = lighting_radial_init_components(o,(r2_t){{0,0}},100,(r2_t){{0,0}},view_dir,M_PI/2.0);    
  }else if(lmodel.modelType == STEREOPOLY_MODEL){
    lighting_stereopoly_options_t* o = lmodel.options;
    approx_data = lighting_stereopoly_init_components(o,view_dir);
  }
//   fitModelToFunction(t->points,t->value, t->num_points,lm_approx,approx_data,nSteps);
  if(lmodel.modelType == COMPACT_MODEL && estimateU){
    fprintf(stderr,"COMPACT MODEL: Estimating Pointlike Light source direction");
    light_compact_estimate_lightdir(t->points,t->value,t->wpos, t->num_points,approx_data,nStepsU);
  }
  if(useNonLinear){
    double alpha = 1.5;
    double beta = 0.7;
    double gamma = 0.2;
    fprintf(stderr,"Estimating Non-Linear Parameters\n");
    approx_model_fit_non_linear_model_simplex(t->points,t->value, t->weights,t->wpos, t->num_points, lm_approx, approx_data, nStepsNL,alpha, beta,gamma);
  }
  approx_model_fit_linear_model(t->points, t->value,t->weights,t->wpos, t->num_points, lm_approx, approx_data);
  
 
  return approx_data;
}







options_t* parse_args(int argc, char** argv);
options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  argparser_get_keyword(pp, "-inprefix");
  o->inprefix = argparser_get_next(pp);
  fprintf(stderr,"INPREFIX %s\n",o->inprefix);
  o->in_lm_data = model_parse(pp,FALSE);
  
  
  argparser_get_keyword(pp, "-outprefix");
  o->outprefix = argparser_get_next(pp);
  fprintf(stderr,"OUTPREFIX %s\n",o->outprefix);
  o->out_lm_data = model_parse(pp,TRUE);
  o->useNonLinear = FALSE;
  o->estimateU = FALSE;
  if( o->out_lm_data.modelType == COMPACT_MODEL){
      o->estimateU = argparser_keyword_present(pp, "-estimateU");
      if(o->estimateU && argparser_keyword_present_next(pp,"steps")){
	o->nStepsU = argparser_get_next_int(pp,1,100000);
      }
      else{
	o->nStepsU = 3;
      }
      o->useNonLinear = argparser_keyword_present(pp, "-useNonLinear");
      if(o->useNonLinear && argparser_keyword_present_next(pp,"steps")){
	o->nStepsNL = argparser_get_next_int(pp,1,100000);
      }
      else{
	o->nStepsNL = 3;
      }
   }
  
  argparser_get_keyword(pp, "-plotprefix");
  o->plotprefix = argparser_get_next(pp);
  fprintf(stderr,"PLOTPREFIX %s\n",o->plotprefix);
  
  

  
  
  
  o->view_dir = (r3_t){{0,0,1}};
  fprintf(stderr,"VIEW DIR ");
  if (argparser_keyword_present(pp, "-viewDir")){
    o->view_dir.c[0] = argparser_get_next_double(pp, -100000.0, 100000.0);
    o->view_dir.c[1] = argparser_get_next_double(pp, -100000.0, 100000.0);
    o->view_dir.c[2] = argparser_get_next_double(pp, -100000.0, 100000.0);
    double len = r3_dir(&(o->view_dir),&(o->view_dir));   
    if(len != 1.0){
      fprintf(stderr,"- Não normalizado: %+3.4lf * ",len);
    }
  }
  fprintf(stderr,"( %3.4f %3.4f %3.4f )\n",o->view_dir.c[0],o->view_dir.c[1],o->view_dir.c[2]);
  
  o->ref_dir = (r3_t){{1,0,0}};
  fprintf(stderr,"REF DIR ");
  if (argparser_keyword_present(pp, "-refDir")){
    o->ref_dir.c[0] = argparser_get_next_double(pp, -100000.0, 100000.0);
    o->ref_dir.c[1] = argparser_get_next_double(pp, -100000.0, 100000.0);
    o->ref_dir.c[2] = argparser_get_next_double(pp, -100000.0, 100000.0);
    double len = r3_dir(&(o->ref_dir),&(o->ref_dir));   
    if(len != 1.0){
      fprintf(stderr,"- Não normalizado: %+3.4lf * ",len);
    }
  }
  fprintf(stderr,"( %3.4f %3.4f %3.4f )\n",o->ref_dir.c[0],o->ref_dir.c[1],o->ref_dir.c[2]);
  
  
  o->resolution= 75;
  if (argparser_keyword_present(pp, "-resolution")){
    o->resolution = argparser_get_next_int(pp, 1.0, 100000.0);
  }
  fprintf(stderr,"RESOLUTION %d\n",o->resolution);
  
  o->plotSamples = 2000;
  if (argparser_keyword_present(pp, "-plotSamples")){
    o->plotSamples = argparser_get_next_int(pp, 1.0, 100000.0);
  }
  
  fprintf(stderr,"PLOTSAMPLES %d\n",o->plotSamples);
  
  o->noise = 0;
  if (argparser_keyword_present(pp, "-noise")){
    o->noise = argparser_get_next_double(pp, 0.0, 100000.0);
  }
  fprintf(stderr,"NOISE %lf\n",o->noise);
  
  o->generateSGplot = argparser_keyword_present(pp,"-generateSGplot");
  if(o->generateSGplot){
    o->plotSamplesSG = argparser_get_next_int(pp, 1.0, 100000.0);
    fprintf(stderr,"USING SGPLOT\n");
  }
  
  
  o->useCosPosWeight = argparser_keyword_present(pp,"-useCosPosWeight");
  if(o->useCosPosWeight){fprintf(stderr,"Cosine weighting function\n");}
  
  o->plotImage = argparser_keyword_present(pp,"-plotImage");
  if(o->plotImage){
    o->radiusImage = argparser_get_next_int(pp, 1.0, 100000.0);
    fprintf(stderr,"PLOTING IMAGE WITH RADIUS %d\n",o->radiusImage);
  }

  argparser_finish(pp);
  return o;

}

int main(int argc, char** argv){
  
  options_t* o = parse_args(argc,argv);
  
  /*First we init the light model to be tested*/
  approx_model_t* lm_test =  create_approx_lighting_model(o->in_lm_data.modelType);
  char* arq_light_param_name = NULL;
  char *arq_light_param_name = jsprintf("%s",o->inprefix);
  FILE* arq_light_param = open_read(arq_light_param_name,TRUE);
  void* l_datatest = lm_test->read_parameters(arq_light_param);
  fclose(arq_light_param);

  /*Generate the testdata*/
  int type_weight = NOWEIGHT;
  if(o->useCosPosWeight) type_weight = COSWEIGHT;
  test_data_t* t = generate_test_data(lm_test,l_datatest, o->resolution, o->view_dir,o->noise,type_weight);
  
  /*We will approximate the model */
  approx_model_t* lm_approx = create_approx_lighting_model(o->out_lm_data.modelType);
  void* l_datapprox = approximate_model(t,lm_approx, o->out_lm_data, o->view_dir,o->estimateU,o->nStepsU, o->useNonLinear,o->nStepsNL);
  if(o->useNonLinear){
    fprintf(stderr,"Plotting goals functions\n");
    plot_non_linear_goal_function(o->plotprefix,t->points,t->value,t->weights,t->wpos, t->num_points, lm_approx, l_datapprox);
  }
  fprintf(stderr,"Modelo aproximado, criando plots\n");
  
  char* arq_light_out_name = NULL;
  char *arq_light_out_name = jsprintf("%s",o->outprefix);
  FILE* arq_out_param = open_write(arq_light_out_name,TRUE);
  lm_approx->write_parameters(arq_out_param,l_datapprox);
  fclose(arq_out_param);
  
  
  /*Plot evrything nicely*/
  char* arq_plot_name = NULL;
  char *arq_plot_name = jsprintf("%s",o->plotprefix);
  FILE* arq_plot = open_write(arq_plot_name,TRUE);
  
  generate_compare_plot(arq_plot,lm_test,l_datatest,lm_approx, l_datapprox, o->plotSamples,o->view_dir,o->ref_dir);
  fclose(arq_plot);
  /*Return zero*/
  if(o->generateSGplot){
    generate_SG_plot(o->plotprefix,lm_test,l_datatest,lm_approx, l_datapprox, o->plotSamplesSG,o->view_dir);
  }
  
  if(o->plotImage){
    int NXplt =  o->radiusImage*2;
    int NYplt =  o->radiusImage*2;
    ellipse_crs_t ellipse_data;
    ellipse_data.rad = o->radiusImage;
    ellipse_data.ctr = (r2_t) {{ NXplt/2.0, NYplt/2.0}};
    ellipse_data.str = (r2_t) {{ 0,0 }};
    
    float_image_t *img_plt = float_image_new(1, NXplt, NYplt);
    fprintf(stderr,"Plotting Synthetic Image of Input Model\n");
    paint_synthetic_gauge_image(img_plt, &(ellipse_data), &(o->view_dir),&l_datatest,lm_test->evaluate);
    char* arq_img_in_name = NULL;
    char *arq_img_in_name = jsprintf("%s-input.pgm",o->outprefix);
    float_pnm_image_write(arq_img_in_name, img_plt,FALSE, 1.0, 0.0,TRUE,TRUE,FALSE);
    
    fprintf(stderr,"Plotting Synthetic Image of Approximated Model\n");
    paint_synthetic_gauge_image(img_plt, &(ellipse_data), &(o->view_dir),&l_datapprox,lm_approx->evaluate);
    char* arq_img_ou_name = NULL;
    char *arq_img_ou_name = jsprintf("%s-approx.pgm",o->outprefix);
    float_pnm_image_write(arq_img_ou_name, img_plt,FALSE, 1.0, 0.0,TRUE,TRUE,FALSE);
    
  }
  
  return 0;
  
}
