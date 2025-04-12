#define PROG_NAME "graph_to-height"
#define PROG_DESC "Integrates  gradient represented by a graph file"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright � 2011 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-01-10 07:53:58 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -graph {GRAPH_FILE} \\\n" \
  "    -refZ {REF_Z_FILE} \\\n" \
  "    -size {NX} {NY} \\\n" \
  "[    -maxIter {MAXITER} ] \\\n" \
  "[    -convTol {CONVTOL} ] \\\n" \
  "[    -wmag {WMAG} ] \\\n" \
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
  "  Created 2011-02-15 by Jorge Stolfi, UNICAMP.\n" \
  "MODIFICATION HISTORY\n" \
  "  2011-02-18 by R. Saracchini, IC-UNICAMP: changed to robust graph reduction.\n" \
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
#include <argparser.h>
#include <pst_img_graph.h>
#include <jsfile.h>
#include <string.h>
#include <float.h>
#include <sys_stats.h>


#define MAXITER_DEFAULT 200
#define CONVTOL_DEFAULT 0.000005
#define WMAG_DEFAULT 1.0

typedef struct options_t{
  char* refFile;
  char* graphFile;
  char* outPrefix;
  int32_t NX,NY;
  double convTol;
  int32_t maxIter;
  bool_t debug;
  double wmag;
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




options_t *parse_options(int32_t argc, char **argv);
options_t *parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 

    if (argparser_keyword_present(pp, "-refZ"))
      { o->refFile = argparser_get_next(pp); }
    else
      { o->refFile = NULL; }

    argparser_get_keyword(pp, "-graph");
    o->graphFile = argparser_get_next(pp);    

    argparser_get_keyword(pp, "-size");
    o->NX = (int32_t)argparser_get_next_int(pp, 0, INT64_MAX);
    o->NY = (int32_t)argparser_get_next_int(pp, 0, INT64_MAX);
    
    o->debug = argparser_keyword_present(pp,"-debug");
    
    o->maxIter = MAXITER_DEFAULT;
    if(argparser_keyword_present(pp,"-maxIter")){
      o->maxIter = argparser_get_next_int(pp,0,INT64_MAX);
    }
    o->convTol = CONVTOL_DEFAULT;
    if(argparser_keyword_present(pp,"-convTol")){
      o->convTol = argparser_get_next_double(pp,0,DBL_MAX);
    }
    o->wmag = WMAG_DEFAULT;
    if(argparser_keyword_present(pp,"-wmag")){
      o->wmag = argparser_get_next_double(pp,0,DBL_MAX);
    }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }


int main(int argc, char** argv){
  
  options_t* o = parse_options(argc,argv);
  
    
  float_image_t* RZ = NULL;
  if(o->refFile != NULL) RZ = readFNI(o->refFile);
  
  int32_t maxIter = o->maxIter;
  double convTol = o->convTol;
  int para = 0;
  int szero = 1;
  bool_t verbose = FALSE;
  
  FILE* arq = open_read(o->graphFile,TRUE);
  pst_img_graph_t* gr = pst_img_graph_read(arq);
  fclose(arq);
  fprintf(stderr,"Generated graph with %ld vertexes and %ld edges\n",gr->n, gr->m);
  int32_t NX = o->NX;
  int32_t NY = o->NY;
  int32_t  NX_Z = NX+1;
  int32_t  NY_Z = NY+1;
  int32_t  NXY_Z = NX_Z*NY_Z;
  float_image_t* OZ = float_image_new(1,NX_Z,NY_Z);
  double* iZ = (double*)malloc(sizeof(double)*NXY_Z);
  double* iW = (double*)malloc(sizeof(double)*NXY_Z);
  
  process_stats_t stats_before = get_process_status();
  double time_before = user_cpu_time_usec();
  
  pst_img_graph_integration_recursive(gr,iZ,iW,o->wmag,maxIter,convTol,para,szero,verbose,0,OZ,RZ,o->outPrefix,o->debug);
  
  double time_after = user_cpu_time_usec(); 
  process_stats_t stats_after = get_process_status();
  FILE* arq_stats = open_write("process_stats.txt",TRUE);
  double diff_time = time_after - time_before;
  double diff_vmen = stats_after.vsize - stats_before.vsize;
  fprintf(arq_stats,"%ld %lf %lf",NX*NY,diff_time/1000000.0,diff_vmen);
  fclose(arq_stats);
  
  char* hgth_filename = NULL;
  char *hgth_filename = jsprintf("%s-oZ.fni",o->outPrefix);
  writeFNI(hgth_filename,OZ);
  
  
  return 0;  
}
