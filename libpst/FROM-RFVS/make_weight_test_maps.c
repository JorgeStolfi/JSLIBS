#define _GNU_SOURCE
#define PROG_NAME "make_weights_test_maps"
#define PROG_VERS "1.0"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "   -prefix {PREFIX} \\\n" \
  "    [-threshold {THRESHOLD}] "

#define PROG_INFO \
  "NAME\n" \
  "  Etc. etc..\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "OPTIONS" \
  "  Etc. etc.."


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <jsfile.h>
#include <argparser.h>
#include <float_image.h>

typedef struct options_t
  {
    double threshold;
    char* prefix;
  } options_t;
  

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */
  
  

float_image_t* readFNI(char* filename);

float_image_t* readFNI(char* filename){
	FILE* arq = open_read(filename,TRUE);
	float_image_t* img = float_image_read(arq);
	fclose(arq);
	return img;
}

void writeFNI(char* filename, float_image_t* img);

void writeFNI(char* filename, float_image_t* img){
	FILE* arq = open_write(filename,TRUE);
	float_image_write(arq,img);
	fclose(arq);
}

void add_border(float_image_t* im);

void add_border(float_image_t* im){
  int NX,NY,NC;
  NC = im->sz[0];
  NX = im->sz[1];
  NY = im->sz[2];
  
  int x,y;
  
  x = 0;
  for(y = 0; y < NY; y++){
    float_image_set_sample(im,0,x,y,0);
  }
  
  x = NX - 1;
  for(y = 0; y < NY; y++){
    float_image_set_sample(im,0,x,y,0);
  }
  
  y = 0;
  for(x = 0; x < NX; x++){
    float_image_set_sample(im,0,x,y,0);
  }
  
  y = NY - 1;
  for(x = 0; x < NX; x++){
    float_image_set_sample(im,0,x,y,0);
  }
  
}

float_image_t*  create_binary_weights(float_image_t* im,double threshold );
float_image_t*  create_binary_weights(float_image_t* im,double threshold ){
  int NX,NY,NC;
  NC = im->sz[0];
  NX = im->sz[1];
  NY = im->sz[2];
  
  float_image_t* im_b = float_image_new(NC,NX,NY);
  
  int x,y;
  for(x = 0; x < NX; x++){
    for(y = 0; y < NY; y++){
	float val = float_image_get_sample(im,0,x,y);
	if(val >= threshold){
	  float_image_set_sample(im_b,0,x,y,1);
	}else{
	  float_image_set_sample(im_b,0,x,y,0);
	}
	
    }
  }
  
  return im_b;
  
}


float_image_t* create_confidence_map(float_image_t* bw);
float_image_t* create_confidence_map(float_image_t* bw){
  int NX,NY,NC;
  NC = bw->sz[0];
  NX = bw->sz[1];
  NY = bw->sz[2];
  
  float_image_t* im_cm = float_image_new(NC,NX,NY);
  
  int x,y;
  for(x = 0; x < NX; x++){
    for(y = 0; y < NY; y++){
      double bw_val = float_image_get_sample(bw,0,x,y);
      double val = 0.0;
      if(bw_val != 0.0){
	//verify 8 neighbors
	val = 1.0;
	int ix,iy;
	
	int startx = ( (x - 1) < 0 ? 0 : x - 1 );
	int starty = ( (y - 1) < 0 ? 0 : y -1 );
	int endx = ( (x + 1) >= NX  ? NX -1 : x + 1 );
	int endy = ( (y + 1) >= NY  ? NY -1 : y + 1 );
	
	for(ix = startx; ix <= endx; ix++){
	  for(iy = starty; iy <= endy; iy++){
	    double n_val = float_image_get_sample(bw,0,ix,iy);
	    if(n_val == 0) val = 0;
	  }
	}
      }
      float_image_set_sample(im_cm,0,x,y,val);
    }
  }
  
  add_border(im_cm);
  
  return im_cm;
}


options_t *parse_options(int argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    argparser_get_keyword(pp, "-prefix");
    o->prefix = argparser_get_next(pp);
    
    o->threshold = 0.5;
    if (argparser_keyword_present(pp, "-threshold")){
      o->threshold = argparser_get_next_double(pp, -1000000.0, +1000000.0);
    }
 
    argparser_finish(pp);
    
    return o;
  }
  
int main(int argc,char** argv){
  
  options_t* o = parse_options(argc,argv);
  
  char* in_filename = NULL;
  char *in_filename = jsprintf("%s-W.fni",o->prefix);
  
  char* bi_filename = NULL;
  char *bi_filename = jsprintf("%s-BW.fni",o->prefix);
  
  char* cm_filename = NULL;
  char *cm_filename = jsprintf("%s-CM.fni",o->prefix);

  float_image_t* weight_map = readFNI(in_filename);
  
  float_image_t* binary_weight_map = create_binary_weights(weight_map,o->threshold);
  
  float_image_t* confidence_map = create_confidence_map(binary_weight_map);
  
  writeFNI(bi_filename,binary_weight_map);
  writeFNI(cm_filename,confidence_map);
  
  return 0;
}

