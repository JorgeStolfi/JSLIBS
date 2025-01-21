#define PROG_NAME "fni_to_graph"
#define PROG_DESC "Generates a graph file from a gradient or weight map"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright � 2005 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-01-10 07:54:21 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -addDiags ] \\\n" \
  "    [ -shrinkSteps ${NUM STEPS} ] \\\n" \
  "    -slopes {const| random| ramp |IG_FNI_NAME} \\\n" \
  "    -weights {const | random | ramp |IW_FNI_NAME} \\\n" \
  "    [-size {NX} {NY} ] \\\n" \
  "    -outPrefix {PREFIX} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  PROG_INFO_FILE "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  fni_to_pnm(1), pnm_to_fni(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2005-08-15 by Jorge Stolfi, UNICAMP.\n" \
  "MODIFICATION HISTORY\n" \
  "  2008-11-08 by J. Stolfi, IC-UNICAMP: added the weight map option.\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

  

#define PROG_INFO_DESC \
    " "
  
#define PROG_INFO_FILE \
  " ."

#define PROG_INFO_OPTS \
  " "


#define _GNU_SOURCE 

#include <stdio.h>
#include <stdlib.h>
#include <float_image.h>
#include <pst_img_graph.h>
#include <pswr.h>
#include <jsfile.h>
#include <assert.h>
#include <math.h>
#include <argparser.h>
#include <teste_func.h>
#include <sys_stats.h>

#define GRAD_FUNC_NONE 0
#define GRAD_FUNC_CONST 1
#define GRAD_FUNC_RAMP 2
#define GRAD_FUNC_RANDOM 3

#define WGHT_FUNC_NONE 0
#define WGHT_FUNC_CONST 1
#define WGHT_FUNC_RAMP 2
#define WGHT_FUNC_RANDOM 3


#define X_AXIS 0
#define Y_AXIS 1

typedef struct options_t{
  char* inPrefix;
  char* outPrefix;
  char* weightMap;
  char* gradientMap;
  int32_t shrinksSteps;
  bool_t addDiags;
  int32_t weightFunction;
  int32_t gradientFunction;
  int32_t NX,NY;
} options_t;

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


void printGRA(char* filename, pst_img_graph_t* g);

void printGRA(char* filename, pst_img_graph_t* g){
  FILE* arq = open_write(filename,TRUE);
  pst_img_graph_write(arq,g);
  fclose(arq);
}

void writeGRA(char* filename, pst_img_graph_t* g);

void writeGRA(char* filename, pst_img_graph_t* g){
  FILE* arq = open_write(filename,TRUE);
  pst_img_graph_write(arq,g);
  fclose(arq);
}


void writeGraph(pst_img_graph_t* g,char* outPrefix,char* tag, int32_t level);
void writeGraph(pst_img_graph_t* g,char* outPrefix,char* tag, int32_t level){
  char* filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.grf",outPrefix,level,tag);
  writeGRA(filename,g);
  filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.txt",outPrefix,level,tag);
  printGRA(filename,g);
  free(filename);
}

options_t *parse_options(int32_t argc, char **argv);
options_t *parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
  
    o->addDiags = argparser_keyword_present(pp, "-addDiags");
    
    if (argparser_keyword_present(pp, "-shrinkSteps"))
      { o->shrinksSteps = argparser_get_next_int32_t(pp, 0, INT64_MAX); }
    else
      { o->shrinksSteps = 0; }


    if (argparser_keyword_present(pp, "-slopes"))
      { 
	o->gradientMap = NULL;
	o->gradientFunction = GRAD_FUNC_NONE;
	if(argparser_keyword_present_next(pp,"random") )
	{
	  o->gradientFunction = GRAD_FUNC_RANDOM;
	}else if (argparser_keyword_present_next(pp,"const") )
	{
	  o->gradientFunction = GRAD_FUNC_CONST;
	}else if (argparser_keyword_present_next(pp,"ramp") ){
	  o->gradientFunction = GRAD_FUNC_RAMP;
	}else{	o->gradientMap = argparser_get_next_non_keyword(pp);}
    }else { o->gradientMap = NULL ;}
    
    if (argparser_keyword_present(pp, "-weights"))
      {
	o->weightMap = NULL;
	o->weightFunction = WGHT_FUNC_NONE;
	if(argparser_keyword_present_next(pp,"random") )
	{
	  o->weightFunction = WGHT_FUNC_RANDOM;
	}else if (argparser_keyword_present_next(pp,"const") )
	{
	  o->weightFunction = WGHT_FUNC_CONST;
	}else{ o->weightMap = argparser_get_next(pp); }
      }else{ o->weightMap = NULL; }
    
    
    
    if( (o->weightMap == NULL) &&  (o->gradientMap == NULL)){
      argparser_get_keyword(pp, "-size");
      o->NX = argparser_get_next_int32_t(pp, 0, INT64_MAX);
      o->NY = argparser_get_next_int32_t(pp, 0, INT64_MAX);
    }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }
  
  
void pst_graph_interpolate_two_samples
  (  float_image_t* I, float_image_t* W,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *v, double* w
   )
   {
     int32_t NX = I->sz[1]; 
     if( W != NULL ){ assert(W->sz[1] == NX);}
     int32_t NY = I->sz[2]; 
     if( W != NULL ){ assert(W->sz[2] == NY);}

     
     double v0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : float_image_get_sample(I,c,x0,y0));
     double v1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : float_image_get_sample(I,c,x1,y1));
     double w0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W,0,x0,y0)));
     double w1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W,0,x1,y1)));


    
   
     /* First we get the interpolation of one diagonal*/
     *v = (w0*v0+w1*v1)/2;
     *w = (w0+w1)/2;
     
   }

void  pst_img_graph_get_axial_edge_data_from_maps(
    float_image_t* IG,
    float_image_t* IW,
    int32_t x, int32_t y,
    int32_t axis, int32_t dir,
    double *d, double *w
 );


int32_t main(int32_t argc, char** argv){
  
  options_t* o = parse_options(argc,argv);
  
  
  
  float_image_t* IG = (o->gradientMap == NULL ? NULL: readFNI(o->gradientMap));
  
  float_image_t* IW = (o->weightMap == NULL ? NULL : readFNI(o->weightMap) );
  
  fprintf(stderr,"START!\n");
  pst_img_graph_t* g = pst_img_graph_from_gradient_weights(IG, IW,o);
  
  assert(g != NULL);
//   fprintf(stderr,"Generated graph with %ld vertices and %ld edges\n",g->n, g->m);
  
  if(o->addDiags == FALSE) fprintf(stderr,"NO DIAGS\n");
  
  int32_t level = 0;
  do {
    g = pst_img_graph_copy(g);
    fprintf(stderr,"Level[%02ld] - %ld vertexes %ld edges \n",level,g->n_valid,g->m_valid);
    pst_img_graph_mark_vertex_removal(g,NULL);
     writeGraph(g,o->outPrefix,"G",level);
    
    double time_before = user_cpu_time_usec();
    pst_img_graph_shrink(g,NULL,1.0);
    double time_after = user_cpu_time_usec(); 
    double diff_time = time_after - time_before;
    fprintf(stderr,"\nTime level [%ld] : %ld %ld %lf\n",level,g->n_valid, g->m_valid,diff_time/1000000.0);
    
    level++;
  }while( (g->n_valid > 2) && (level < o->shrinksSteps) );
  
  return 0;  
}


