#define PROG_NAME "generate_analytic_inversion"
#define _GNU_SOURCE


#include <jsfile.h>
#include <lighting_models.h>
#include <approx_system.h>
#include <polynomial_functions.h>
#include <analytic_inversion.h>
#include <float_image.h>


#define PROG_HELP \
  PROG_NAME " \\\n" \
  "      -nLights {NLIGHTS} \\\n" \
  "      -outprefix {FILENAME} \\\n" \
  "      -gaugeModel \\\n" \
  "           {compact,radial,stereopoly,harmonic} \\\n" \
  "      -invModel {MODELPARMS}  \\\n" \
  "      -gaugeNames {GAUGE NAMES} \\\n" \
  "      -resolution {RESOLUTION} \\\n" \
  "      [-thetaMax {thetaMax} ] \\\n" \
  "      [-mapSize {MAPSIZE} ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO ""

typedef struct options_t{
  int nLights;
  char** gaugeNames;
  model_options_t lmodel;
  model_options_t imodel;
  r3_t view_dir;
  double thetaMax;
  int mapSize;
  int resolution;
  char* outPrefix;
  
} options_t ;


options_t* parse_args(int argc, char** argv);
options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  argparser_get_keyword(pp, "-nLights");
  o->nLights = argparser_get_next_int(pp, 3, 100000);
  argparser_get_keyword(pp, "-gaugeModel");
  o->lmodel = model_parse(pp,FALSE);
  argparser_get_keyword(pp, "-gaugeNames");
  int i;
  o->gaugeNames = (char**)malloc(sizeof(char*)*(o->nLights));
  for(i = 0; i < o->nLights; i++){
    o->gaugeNames[i] = argparser_get_next(pp);
  }
  argparser_get_keyword(pp, "-invModel");
  o->imodel = model_parse(pp,TRUE);
  
  argparser_get_keyword(pp, "-resolution");
  o->resolution = argparser_get_next_int(pp, 1, 100000);
  
  o->thetaMax = M_PI/2.0;
  if( argparser_keyword_present_next(pp,"-thetaMax") ){
    o->thetaMax = argparser_get_next_double(pp, 0, M_PI/2.0);
  }
  
  o->mapSize = 100;
  if( argparser_keyword_present_next(pp,"-mapSize") ){
    o->mapSize = argparser_get_next_int(pp, 1, 100000);
  }
  
  argparser_get_keyword(pp, "-outPrefix");
  o->outPrefix = argparser_get_next(pp);
  
  argparser_finish(pp);
  return o;
}

int main(int argc, char** argv){
  options_t* o = parse_args(argc,argv);
  
  int i;
  approx_model_t* gauge_amodel = create_approx_lighting_model(o->lmodel.modelType);
  
  void** gauge_ldata = (void**)malloc(sizeof(void*)*(o->nLights));
  for(i = 0; i < o->nLights; i++){
    FILE* arq_parms = open_read(o->gaugeNames[i],TRUE);
    gauge_ldata[i] = gauge_amodel->read_parameters(arq_parms);
  }
  
  approx_model_t* inv_model = NULL;
  void* inv_data[3];
  if( o->imodel.modelType == GENERICPOLY_MODEL){
      inv_model = create_ls_polynomial_model();
      int c;
      for(c = 0 ; c < 3; c++){
	poly_function_options_t* po = o->imodel.options;
	inv_data[c] = poly_model_init_components(*po);
      }
  }
  void* mag_data = (void**)malloc(sizeof(void*));
  if( o->imodel.modelType == GENERICPOLY_MODEL){
    poly_function_options_t* po = o->imodel.options;
    mag_data = poly_model_init_components(*po);
  }
  
  
  
  analytic_inversion_t* ai = analytic_inversion_create_model
  (
    o->nLights,
    o->lmodel.modelType,
    gauge_amodel,
    gauge_ldata,
    o->view_dir,
    o->thetaMax,
    o->resolution,
    o->imodel.modelType,
    inv_model,
    inv_data,
    mag_data,
    NULL,
    o->outPrefix
  );
 
 char* inv_filename_model =  NULL;
 char *inv_filename_model = jsprintf("%s_Imodel.txt",o->outPrefix);
 FILE* inv_arq_model = open_write(inv_filename_model,TRUE);
 analytic_inversion_write(inv_arq_model,ai);
 fclose(inv_arq_model);
 
 
 char* inv_filename_map =  NULL;
 char *inv_filename_map = jsprintf("%s_Imapa.fni",o->outPrefix);
 FILE* inv_arq_map = open_write(inv_filename_map,TRUE);
 float_image_t* im_map = analytic_inversion_plot_map(ai,o->mapSize);
 float_image_write(inv_arq_map,im_map);
 fclose(inv_arq_map);
  
//  analytic_inversion_generate_SG_plot(o->outPrefix, ai, o->mapSize ,o->view_dir);
 
 return 0;
  
}