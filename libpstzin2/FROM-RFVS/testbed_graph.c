#define PROG_NAME "testbed_graph"
#define PROG_DESC "Compares the solution of a graph with the solution of itself without one vertex"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright (c)2011 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -graphPrefix {GRAPH_PREFIX} \\\n" \
  "    -outPrefix {OUT_PREFIX} \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    -removeVertex {XR} {YR} \\\n" \
  "    -weights {NWEIGHTS} {W0} {W1}...{WN-1} \\\n" \
  "    -errFile {ERRFILE} \\\n" \
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

#include <stdio.h>
#include <stdlib.h>
#include <argparser.h>
#include <float_image.h>
#include <pst_img_graph.h>
#include <pswr.h>
#include <fget.h>
#include <jsfile.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <affirm.h>
#include <float.h>


typedef struct options_t{
  char* graphPrefix;
  char* errFile;
  long int NX,NY;
  char* outPrefix;
  r2_t removeVertex;
  long int nWeights;
  double* weights;
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


void printGRA(char* filename, pst_img_graph_t* g);

void printGRA(char* filename, pst_img_graph_t* g){
  FILE* arq = open_write(filename,TRUE);
  pst_img_graph_print(arq,g);
  fclose(arq);
}

void writeGRA(char* filename, pst_img_graph_t* g);

void writeGRA(char* filename, pst_img_graph_t* g){
  FILE* arq = open_write(filename,TRUE);
  pst_img_graph_write(arq,g);
  fclose(arq);
}

pst_img_graph_t* readGRA(char* filename);

pst_img_graph_t* readGRA(char* filename){
  FILE* arq = open_read(filename,TRUE);
  pst_img_graph_t* g = pst_img_graph_read(arq);
  fclose(arq);
  return g;
}

pst_img_graph_t* readGraph(char* inPrefix,char* tag,int level);
pst_img_graph_t* readGraph(char* inPrefix,char* tag,int level){
  char* filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.grf",inPrefix,level,tag);
  pst_img_graph_t* g = readGRA(filename);
  free(filename);
  return g;
  
}



void writeGraph(pst_img_graph_t* g,char* outPrefix,char* tag, int level);
void writeGraph(pst_img_graph_t* g,char* outPrefix,char* tag, int level){
  char* filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.txt",outPrefix,level,tag);
  printGRA(filename,g);
  free(filename);
  filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.grf",outPrefix,level,tag);
  writeGRA(filename,g);
  free(filename);
}

void compute_solution_weights(pst_img_graph_t *ig,double *iW,r2_t rem_coord,double *eW );
void compute_solution_weights(pst_img_graph_t *ig,double *iW,r2_t rem_coord,double *eW){

  long int i;
  
  for(i = 0; i < ig->n; i++){
    pst_vertex_data_t* iv = &(ig->vertex[i]);

    if(iv->id == -1) continue;
    double d2 = r2_dist_sqr(&rem_coord,&(iv->coords));
    eW[i] = iW[i]*d2;
  }
  
}

void normalize_solution_dist2(pst_img_graph_t* ig, double* iZ,double *eW);
void normalize_solution_dist2(pst_img_graph_t* ig, double* iZ,double *eW){
  
  double sW = 0;
  double sWZ = 0;
  long int i;
  
  for(i = 0; i < ig->n; i++){
    pst_vertex_data_t* iv = &(ig->vertex[i]);

    if(iv->id == -1) continue;
    sW+= eW[i];
    sWZ+= (eW[i]*iZ[i]);
  }
  
  if(sW == 0) return;
  
  double avg = sWZ/sW;
  for(i = 0; i < ig->n; i++){
    pst_vertex_data_t* iv = &(ig->vertex[i]);
    if(iv->id == -1) iZ[i] = 0;
    else{ 
      iZ[i]-=avg;
    }
  }
  
  
}

options_t* parse_options(int argc, char **argv);


int main(int argc, char** argv){

  options_t* o = parse_options(argc,argv);
  
 
  long int NX_Z = o->NX +1;
  long int NY_Z = o->NY +1;
  
  char* ig_prefix = NULL;
  char *ig_prefix = jsprintf("%s-I",o->outPrefix);
  char* jg_prefix = NULL;
  char *jg_prefix = jsprintf("%s-J",o->outPrefix);
  
  fprintf(stderr,"START!\n");
  pst_img_graph_t*  ig = readGraph(o->graphPrefix,"G", 0);
  writeGraph(ig,ig_prefix,"G", 0);
  
  assert(ig != NULL);
  fprintf(stderr,"Generated graph with %ld vertices and %ld edges\n",ig->n, ig->m);

  
  
  
  float_image_t* OZ = float_image_new(1,NX_Z,NY_Z);
  float_image_t* SZ = float_image_new(1,NX_Z,NY_Z);
  float_image_t* CZ = float_image_new(1,NX_Z,NY_Z);
  double* iZ = (double*)malloc(sizeof(double)*(ig->n));
  double* jZ = (double*)malloc(sizeof(double)*(ig->n));
  double* iW = (double*)malloc(sizeof(double)*(ig->n));
  double* jW = (double*)malloc(sizeof(double)*(ig->n));
  double* sZ = (double*)malloc(sizeof(double)*(ig->n));
  
  long int maxIter = 50000;
  double convTol = 0.00000005;
  int para = 0;
  int szero = 1;
  bool_t verbose = FALSE;
  pst_img_graph_integration(ig,iZ,iW,maxIter,convTol,para,szero,verbose,0,OZ,NULL,ig_prefix);
  
  
  pst_img_graph_t* jg = pst_img_graph_copy(ig);
    
  long int ix = pst_img_graph_find_nearest_vertex(jg,o->removeVertex);
  demand(ix != -1,"Invalid vertex index !");
  pst_vertex_data_t* v = &(jg->vertex[ix]);
  fprintf(stderr,"Removing vertex [%ld] = %ld = (%lf,%lf) with query (%lf %lf)\n",
	  ix,v->id,v->coords.c[0],v->coords.c[1],o->removeVertex.c[0],o->removeVertex.c[1]
  );
  long int n_edges = pst_img_graph_vertex_count_neighbours(jg,ix);
  if(o->weights != NULL ){
    demand(n_edges == o->nWeights,"Number of weights differ from number of edges.");
  }
  pst_img_graph_vertex_remove(jg, ix,o->weights,o->wmag,TRUE);
  pst_img_graph_remove_paralel_edges(jg);
   
  pst_img_graph_integration(jg,jZ,jW,maxIter,convTol,para,szero,verbose,0,SZ,NULL,jg_prefix);
  pst_img_graph_copy_solution_from_shrunk(jg,jZ,ig,sZ);
  double *eW = (double*)malloc(sizeof(double)*(ig->n));
  
  compute_solution_weights(ig,iW,o->removeVertex,eW);
  
  normalize_solution_dist2(ig,iZ,eW);
  normalize_solution_dist2(ig,sZ,eW);
  
  char* filename_error = NULL;
  char *filename_error = jsprintf("%s-ERR.txt",o->outPrefix);
  FILE* arq_err = open_write(filename_error,TRUE);
  double sumW = 1.0e-300; /*safe*/
  double sumWE2 = 0;
  float_image_fill_channel(CZ,0,0);
  long int i;
  for(i = 0; i < ig->n; i++){
    pst_vertex_data_t* iv = &(ig->vertex[i]);
    if(iv->id != -1){
      if(iv->mark != MARK_VERTEX_REMOVED){
	long int ix,iy;
	pst_img_graph_get_vertex_image_indices(&(iv->coords),NX_Z,NY_Z,&ix,&iy);
	double err = iZ[i] - sZ[i];
	double d2 = r2_dist_sqr(&(o->removeVertex),&(iv->coords));
	double w = iW[i]/(1+d2);
	sumW+=w;
	sumWE2+=w*err*err;
	fprintf(arq_err,"%9.3lf %9.3lf %12.6lf %e\n",iv->coords.c[0],iv->coords.c[1],sqrt(d2),err);
	float_image_set_sample(CZ,0,ix,iy,err);
      }
    }
  }
  
  fclose(arq_err);
  free(filename_error);
  double avgErr = sqrt(sumWE2/sumW);
  FILE* arq_err_avg = fopen(o->errFile,"at");
  demand(arq_err_avg != NULL,"Could not open errFile !");
  fprintf(arq_err_avg,"%12.6lf %e\n",o->wmag,avgErr);
  fprintf(stderr,"WMAG = %12.6lf AVGERR = %e\n",o->wmag,avgErr);
  fclose(arq_err_avg);
  
  
	
  char* diff_filename = NULL;
  char *diff_filename = jsprintf("%s-D.fni",o->outPrefix);
  writeFNI(diff_filename,CZ);
  free(diff_filename);
  free(ig_prefix);
  free(jg_prefix);
  free(iZ);
  free(jZ);
  free(sZ);
  free(iW);
  free(jW);
  
  
  
  return 0;  
}


options_t* parse_options(int argc, char **argv){
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
   
  options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
  
  argparser_get_keyword(pp, "-graphPrefix");
  o->graphPrefix = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-removeVertex");
  o->removeVertex.c[0] = argparser_get_next_double(pp, -DBL_MAX, DBL_MAX);
  o->removeVertex.c[1] = argparser_get_next_double(pp, -DBL_MAX, DBL_MAX);
  
  argparser_get_keyword(pp, "-size");
  o->NX = argparser_get_next_int(pp, 1, INT64_MAX );
  o->NY = argparser_get_next_int(pp, 1, INT64_MAX );
  
   o->wmag = 1.0;
    if(argparser_keyword_present(pp,"-wmag")){
      o->wmag = argparser_get_next_double(pp,0,DBL_MAX);
    }
  
  o->nWeights = 0;
  o->weights = NULL;
  if(argparser_keyword_present(pp, "-weights")){
    o->nWeights = argparser_get_next_int(pp, 1, INT64_MAX );
    o->weights = (double*)malloc(sizeof(double)*(o->nWeights));
    long int i;
    for(i = 0; i < o->nWeights;i++){
      o->weights[i] = argparser_get_next_double(pp, 0, DBL_MAX);
    }
  }
  
  
  argparser_get_keyword(pp, "-outPrefix");
  o->outPrefix = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-errFile");
  o->errFile = argparser_get_next(pp);

  argparser_skip_parsed(pp);
    
  argparser_finish(pp);
    
    
    
  return o;
}
