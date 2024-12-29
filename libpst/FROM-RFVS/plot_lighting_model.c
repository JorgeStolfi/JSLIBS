#define PROG_NAME "compare_models"
#include <assert.h>
#include <affirm.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <string.h>
#include <argparser.h>
#include <jsfile.h>
#include <lighting_models.h>
#include <lighting_compact.h>
#include <lighting_harmonic.h>
#include <lighting_radial.h>
#include <lighting_stereopoly.h>
#include <lighting_nearest.h>
#include <jsrandom.h>


#define PROG_HELP \
  PROG_NAME " \\\n" \
  "      -model {FILENAME} {MODEL}\\\n" \
  "      -plotprefix {FILENAME} \\\n" \
  "      [ -viewDir {VX} {VY} VZ} ]\\\n" \
  "      [ -refDir {VX} {VY} VZ} ]\\\n" \
  "      [ -plotSamples {PLOTSAMPLES} ] \\\n" \
  "      [ -generateSGplot {SGPLOTSAMPLES} ] \\\n" \
  "      [ -plotImage {RADIUS}] \\\n" \
  "      [ -renderNormalMap {NORMAL MAP}] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO ""

struct options_t{
  r3_t view_dir;
  r3_t ref_dir;
  
//   int nSteps;
  char* modelA;
  model_options_t lm_opt_A;
   
  char* plotprefix;
  
  int plotSamples; /*Number of samples to be used in the plot*/
  bool_t generateSGplot;
  int plotSamplesSG;
  
  bool_t plotImage;
  int radiusImage;
  
  bool_t renderNormalMap;
  char* normalMap;
  
};

typedef struct options_t options_t;

options_t* parse_args(int argc, char** argv);
options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  argparser_get_keyword(pp, "-model");
  o->modelA = argparser_get_next(pp);
  fprintf(stderr,"ModelA %s\n",o->modelA);
  o->lm_opt_A = model_parse(pp,FALSE);
  
  
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
      fprintf(stderr,"!! Não normalizado: %+3.4lf * ",len);
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
      fprintf(stderr,"!! Não normalizado: %+3.4lf * ",len);
    }
  }
  fprintf(stderr,"( %3.4f %3.4f %3.4f )\n",o->ref_dir.c[0],o->ref_dir.c[1],o->ref_dir.c[2]);
  
  o->plotSamples = 2000;
  if (argparser_keyword_present(pp, "-plotSamples")){
    o->plotSamples = argparser_get_next_int(pp, 1, 100000);
  }
  
  fprintf(stderr,"PLOTSAMPLES %d\n",o->plotSamples);
  
  o->generateSGplot = argparser_keyword_present(pp,"-generateSGplot");
  if(o->generateSGplot){
    o->plotSamplesSG = argparser_get_next_int(pp, 1, 100000);
    fprintf(stderr,"USING SGPLOT\n");
  }
  
  
  
  o->plotImage = argparser_keyword_present(pp,"-plotImage");
  if(o->plotImage){
    o->radiusImage = argparser_get_next_int(pp, 1.0, 100000.0);
    fprintf(stderr,"PLOTING IMAGE WITH RADIUS %d\n",o->radiusImage);
  }

  o->renderNormalMap = argparser_keyword_present(pp,"-renderNormalMap");
  if(o->renderNormalMap){
    o->normalMap = argparser_get_next(pp);
  }
  
  argparser_finish(pp);
  return o;

}

int main(int argc, char** argv){
  
  options_t* o = parse_args(argc,argv);
  
  /*First we init the light model to be tested*/
  approx_model_t* lm_test_A =  create_approx_lighting_model(o->lm_opt_A.modelType);
  char* arq_light_param_name_A = NULL;
  char *arq_light_param_name_A = jsprintf("%s",o->modelA);
  FILE* arq_light_param_A = open_read(arq_light_param_name_A,TRUE);
  void* l_data_A = lm_test_A->read_parameters(arq_light_param_A);
  lm_test_A->write_parameters(stderr,l_data_A);
  fclose(arq_light_param_A);
  
  
   
  /*Plot evrything nicely*/
  char* arq_plot_name = NULL;
  char *arq_plot_name = jsprintf("%s",o->plotprefix);
  FILE* arq_plot = open_write(arq_plot_name,TRUE);
  lighting_model_generate_plot(arq_plot, lm_test_A,l_data_A,o->plotSamples, o->view_dir, o->ref_dir);
  fclose(arq_plot);
  /*Return zero*/
  if(o->generateSGplot){
    lighting_model_generate_SG_plot(o->plotprefix, lm_test_A,l_data_A,o->plotSamplesSG,o->view_dir);
  }
  
  if(o->plotImage){
    int NXplt =  o->radiusImage*2;
    int NYplt =  o->radiusImage*2;
    ellipse_crs_t ellipse_data;
    ellipse_data.rad = o->radiusImage;
    ellipse_data.ctr = (r2_t) {{ NXplt/2.0, NYplt/2.0}};
    ellipse_data.str = (r2_t) {{ 0,0 }};
    
    float_image_t *img_plt = float_image_new(1, NXplt, NYplt);
    fprintf(stderr,"Plotting Image of Input Model\n");
    paint_synthetic_gauge_image(img_plt, &(ellipse_data), &(o->view_dir),&l_data_A,lm_test_A->evaluate);
    char* arq_img_in_name = NULL;
    char *arq_img_in_name = jsprintf("%s-modelA.pgm",o->plotprefix);
    float_pnm_image_write(arq_img_in_name, img_plt,FALSE, 1.0, 0.0,TRUE,TRUE,FALSE);
    
  }
  
  if(o->renderNormalMap){
    
    fprintf(stderr,"Loading Normal map Image \n");
    FILE* arq_norm_map = open_read(o->normalMap,TRUE);
    float_image_t* normal_map = float_image_read(arq_norm_map);
    fclose(arq_norm_map);
    int NX,NY,NC;
    NC = normal_map->sz[0];
    NX = normal_map->sz[1];
    NY = normal_map->sz[2];
    assert(NC == 3);
    int x,y;
    float_image_t* render_map = float_image_new(1,NX,NY);
    fprintf(stderr,"Rendering Normal Mapl\n");
    for(x = 0; x < NX; x++){
      for(y = 0; y < NY; y++){
	double dx = float_image_get_sample(normal_map,0,x,y);
	double dy = float_image_get_sample(normal_map,1,x,y);
	double dz = float_image_get_sample(normal_map,2,x,y);
	r3_t normal = (r3_t){{ dx,dy,dz}};
	r3_dir(&normal,&normal);
	double val = lm_test_A->evaluate(normal.c,l_data_A);
	float_image_set_sample(render_map,0,x,y,val);
      }
    }
    char* arq_render_name = NULL;
    char *arq_render_name = jsprintf("%s-renderA.pgm",o->plotprefix);
    float_pnm_image_write(arq_render_name, render_map,FALSE, 1.0, 0.0,TRUE,TRUE,FALSE);
    
  }
  
  return 0;
  
}
