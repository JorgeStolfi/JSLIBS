#define PROG_NAME "interpolate_normals"

#define _GNU_SOURCE
#include <stdio.h>
#include <normal_interpolation.h>
#include <float_image.h>
#include <ellipse_crs.h>
#include <assert.h>
#include <affirm.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <string.h>
#include <jsfile.h>
#include <argparser.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -nGauges {NUM} \\\n" \
  "    -inPrefix {INPUT PREFIX} \\\n" \
  "    -gaugeTags {gauge0 gauge1... gaugeN} \\\n" \
  "    -outPrefix {OUTPUT PREFIX} \\\n" \
  "    -interpolation {bestLogProb|avgLogProb}  \\\n" \
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
    char* inPrefix;
    char* outPrefix;
    int nGauges;
    char** gaugeTags;
    int polyDegree;
    r2_t* gaugeCenters;
    normal_interp_t interpType;
    double high12K;
    r3_t* highClusterDir;
};

typedef struct options_t options_t;

options_t *parse_args(int argc, char** argv);

options_t *parse_args(int argc, char** argv){
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);      

  options_t *o = (options_t *)malloc(sizeof(options_t));
  
  argparser_get_keyword(pp, "-inPrefix");
  o->inPrefix = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-outPrefix");
  o->outPrefix = argparser_get_next(pp);
    
  argparser_get_keyword(pp, "-nGauges");
  o->nGauges = argparser_get_next_int(pp, 1, 10000);
  
  argparser_get_keyword(pp, "-gaugeTags");
  o->gaugeTags = (char**)malloc(sizeof(char*)*(o->nGauges));
  int i;
  for(i = 0; i < o->nGauges; i++){
    o->gaugeTags[i] = argparser_get_next(pp);
  }
  
  argparser_get_keyword(pp, "-interpolation");
  char* interp  = argparser_get_next(pp);
  o->highClusterDir = NULL;
  
  if(strcmp(interp,"bestLogProb") == 0){
    o->interpType = PROB_BEST;
  }else if(strcmp(interp,"avgLogProb") == 0){
    o->interpType = PROB_AVERAGE;
  }
  else if (strcmp(interp,"high12PS") == 0){
    o->highClusterDir = (r3_t*)malloc(sizeof(r3_t)*(o->nGauges));
    o->interpType = HIGHLIGHT_12PS;
    argparser_get_keyword(pp, "-kFactor");
    o->high12K = argparser_get_next_double(pp, 10e-10, 10000);
    argparser_get_keyword(pp, "-clusterDir");
    int i;
    for(i = 0; i < o->nGauges; i++){
      o->highClusterDir[i].c[0] = argparser_get_next_double(pp, -1.0, 1.0);
      o->highClusterDir[i].c[1] = argparser_get_next_double(pp, -1.0, 1.0);
      o->highClusterDir[i].c[2] = argparser_get_next_double(pp, 0.0, 1.0);
    }
  }
  else if(strcmp(interp,"polyPos") == 0){
    o->interpType = POS_POLY;
    argparser_get_keyword(pp, "-polyDegree");
    o->polyDegree = argparser_get_next_int(pp, 1, 10000);
    double scale = 1.0;
    if(argparser_keyword_present(pp, "-scale")){
      scale = scale*argparser_get_next_double(pp, 0, 10000000);
    }
    argparser_get_keyword(pp, "-gaugeCenters");
    o->gaugeCenters = (r2_t*)malloc(sizeof(r2_t)*(o->nGauges));
    int i;
    for(i = 0; i < o->nGauges; i++){
      o->gaugeCenters[i].c[0] = scale*argparser_get_next_double(pp, 0, 10000000);
      o->gaugeCenters[i].c[1] = scale*argparser_get_next_double(pp, 0, 10000000);
    }
  }else{
    argparser_error(pp, "Invalid interpolation option");
  }
 // o->gaugeCenters = NULL;
  
  argparser_finish(pp);
  return o;
}

float_image_t* FNIRead(char* filename);
float_image_t* FNIRead(char* filename){
  FILE* arq = open_read(filename,TRUE);
  float_image_t* im =  float_image_read(arq);
  fclose(arq);
  return im;
}

void FNIWrite(char* filename,float_image_t* im);
void FNIWrite(char* filename,float_image_t* im){
  FILE* arq = open_write(filename,TRUE);
  float_image_write(arq,im);
  fclose(arq);
}

int main(int argc,char** argv){
  options_t *o = parse_args(argc, argv);
  
  /*First read the normal maps*/
  float_image_t** normal_maps = (float_image_t**)malloc(sizeof(float_image_t*)*(o->nGauges));
  int i;
  fprintf(stderr,"Reading Normal maps\n");
  for(i = 0; i < o->nGauges;i++){
    char* nm_name = NULL;
    char *nm_name = jsprintf("%s_G%s_normals.fni",o->inPrefix,o->gaugeTags[i]);
    normal_maps[i] = FNIRead(nm_name);
    free(nm_name);
  }
  
  float_image_t* interp_normal_map = NULL;
  if((o->interpType == PROB_AVERAGE) || (o->interpType == PROB_BEST) || (o->interpType == HIGHLIGHT_12PS)){
    
    fprintf(stderr,"Reading LogProb maps, required for Probability-based interpolation\n");
    float_image_t** prob_maps = (float_image_t**)malloc(sizeof(float_image_t*)*(o->nGauges));
    
    for(i = 0; i < o->nGauges;i++){
      char* nm_name = NULL;
      char *nm_name = jsprintf("%s_G%s_logPrSG.fni",o->inPrefix,o->gaugeTags[i]);
      prob_maps[i] = FNIRead(nm_name);
      free(nm_name);
    }
    
    fprintf(stderr,"Interpolating maps...");
    if(o->interpType  != HIGHLIGHT_12PS){
      float_image_t* select_map = float_image_new(1,normal_maps[0]->sz[1],normal_maps[0]->sz[2]);
      interp_normal_map = normal_interpolate_prob(normal_maps,prob_maps, o->nGauges, o->interpType,select_map);
      fprintf(stderr,"OK\n");
      
      char* selectfilename = NULL;
      char *selectfilename = jsprintf("%s_gabselect.ppm",o->outPrefix);
      float_pnm_image_write(selectfilename, select_map,FALSE, 1.0, 0.0,TRUE,TRUE,FALSE);
      char *selectfilename = jsprintf("%s_gabselect.fni",o->outPrefix);
      FNIWrite(selectfilename,select_map);
    }else{
      interp_normal_map = normal_interpolate_hightlight12(normal_maps,prob_maps,o->nGauges,o->highClusterDir,o->high12K);
    }
    
  }else if(o->interpType == POS_POLY){
    interp_normal_map = normal_interpolate_pos(normal_maps,o->gaugeCenters, o->nGauges, o->polyDegree);
  }else{
    fprintf(stderr,"Sorry, method not implemented !\n");
    return 0;
  }
  
  fprintf(stderr,"Saving interpolated normal map\n");
  char* outfilename = NULL;
  char *outfilename = jsprintf("%s_normals.fni",o->outPrefix);
  FNIWrite(outfilename,interp_normal_map);
  
  return 0;
}