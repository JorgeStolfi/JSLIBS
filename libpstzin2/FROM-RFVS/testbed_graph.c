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
  int32_t NX,NY;
  char* outPrefix;
  r2_t removeVertex;
  int32_t nWeights;
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


void printGRA(char* filename, pst_img_graph_t* gr);

void printGRA(char* filename, pst_img_graph_t* gr){
  FILE* arq = open_write(filename,TRUE);
  pst_img_graph_write(arq,gr);
  fclose(arq);
}

void writeGRA(char* filename, pst_img_graph_t* gr);

void writeGRA(char* filename, pst_img_graph_t* gr){
  FILE* arq = open_write(filename,TRUE);
  pst_img_graph_write(arq,gr);
  fclose(arq);
}

pst_img_graph_t* readGRA(char* filename);

pst_img_graph_t* readGRA(char* filename){
  FILE* arq = open_read(filename,TRUE);
  pst_img_graph_t* gr = pst_img_graph_read(arq);
  fclose(arq);
  return gr;
}

pst_img_graph_t* readGraph(char* inPrefix,char* tag,int32_t level);
pst_img_graph_t* readGraph(char* inPrefix,char* tag,int32_t level){
  char* filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.grf",inPrefix,level,tag);
  pst_img_graph_t* gr = readGRA(filename);
  free(filename);
  return gr;
  
}



void writeGraph(pst_img_graph_t* gr,char* outPrefix,char* tag, int32_t level);
void writeGraph(pst_img_graph_t* gr,char* outPrefix,char* tag, int32_t level){
  char* filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.txt",outPrefix,level,tag);
  printGRA(filename,gr);
  free(filename);
  filename = NULL;
  char *filename = jsprintf("%s-%02d-%s.grf",outPrefix,level,tag);
  writeGRA(filename,gr);
  free(filename);
}

void compute_solution_weights(pst_img_graph_t *gri,double *iW,r2_t rem_coord,double *eW );
void compute_solution_weights(pst_img_graph_t *gri,double *iW,r2_t rem_coord,double *eW){

  int32_t i;
  
  for(i = 0; i < gri->n; i++){
    pst_vertex_data_t* iv = &(gri->vertex[i]);

    if(iv->id == -1) continue;
    double d2 = r2_dist_sqr(&rem_coord,&(iv->coords));
    eW[i] = iW[i]*d2;
  }
  
}

void normalize_solution_dist2(pst_img_graph_t* gri, double* iZ,double *eW);
void normalize_solution_dist2(pst_img_graph_t* gri, double* iZ,double *eW){
  
  double sW = 0;
  double sWZ = 0;
  int32_t i;
  
  for(i = 0; i < gri->n; i++){
    pst_vertex_data_t* iv = &(gri->vertex[i]);

    if(iv->id == -1) continue;
    sW+= eW[i];
    sWZ+= (eW[i]*iZ[i]);
  }
  
  if(sW == 0) return;
  
  double avg = sWZ/sW;
  for(i = 0; i < gri->n; i++){
    pst_vertex_data_t* iv = &(gri->vertex[i]);
    if(iv->id == -1) iZ[i] = 0;
    else{ 
      iZ[i]-=avg;
    }
  }
  
  
}

options_t* parse_options(int32_t argc, char **argv);


int32_t main(int32_t argc, char** argv){

  options_t* o = parse_options(argc,argv);
  
 
  int32_t NX_Z = o->NX +1;
  int32_t NY_Z = o->NY +1;
  
  char* gri_prefix = NULL;
  char *gri_prefix = jsprintf("%s-I",o->outPrefix);
  char* grj_prefix = NULL;
  char *grj_prefix = jsprintf("%s-J",o->outPrefix);
  
  fprintf(stderr,"START!\n");
  pst_img_graph_t*  gri = readGraph(o->graphPrefix,"G", 0);
  writeGraph(gri,gri_prefix,"G", 0);
  
  assert(gri != NULL);
  fprintf(stderr,"Generated graph with %ld vertices and %ld edges\n",gri->n, gri->m);
  
  float_image_t* OZ = float_image_new(1,NX_Z,NY_Z);
  float_image_t* SZ = float_image_new(1,NX_Z,NY_Z);
  float_image_t* CZ = float_image_new(1,NX_Z,NY_Z);
  double* iZ = (double*)malloc(sizeof(double)*(gri->n));
  double* jZ = (double*)malloc(sizeof(double)*(gri->n));
  double* iW = (double*)malloc(sizeof(double)*(gri->n));
  double* jW = (double*)malloc(sizeof(double)*(gri->n));
  double* sZ = (double*)malloc(sizeof(double)*(gri->n));
  
  int32_t maxIter = 50000;
  double convTol = 0.00000005;
  int32_t para = 0;
  int32_t szero = 1;
  bool_t verbose = FALSE;
  pst_img_graph_integration(gri,iZ,iW,maxIter,convTol,para,szero,verbose,0,OZ,NULL,gri_prefix);
  
  
  pst_img_graph_t* grj = pst_img_graph_copy(gri);
    
  int32_t ix = pst_img_graph_find_nearest_vertex(grj,o->removeVertex);
  demand(ix != -1,"Invalid vertex index !");
  pst_vertex_data_t* v = &(grj->vertex[ix]);
  fprintf(stderr,"Removing vertex [%ld] = %ld = (%lf,%lf) with query (%lf %lf)\n",
	  ix,v->id,v->coords.c[0],v->coords.c[1],o->removeVertex.c[0],o->removeVertex.c[1]
  );
  int32_t n_edges = pst_img_graph_vertex_out_degree(grj,ix);
  if(o->weights != NULL ){
    demand(n_edges == o->nWeights,"Number of weights differ from number of edges.");
  }
  pst_img_graph_vertex_remove(grj, ix,o->weights,o->wmag,TRUE);
  pst_img_graph_remove_paralel_edges(grj);
   
  pst_img_graph_integration(grj,jZ,jW,maxIter,convTol,para,szero,verbose,0,SZ,NULL,grj_prefix);
  pst_img_graph_copy_solution_from_shrunk(grj,jZ,gri,sZ);
  double *eW = (double*)malloc(sizeof(double)*(gri->n));
  
  compute_solution_weights(gri,iW,o->removeVertex,eW);
  
  normalize_solution_dist2(gri,iZ,eW);
  normalize_solution_dist2(gri,sZ,eW);
  
  char* filename_error = NULL;
  char *filename_error = jsprintf("%s-ERR.txt",o->outPrefix);
  FILE* arq_err = open_write(filename_error,TRUE);
  double sumW = 1.0e-300; /*safe*/
  double sumWE2 = 0;
  float_image_fill_channel(CZ,0,0);
  int32_t i;
  for(i = 0; i < gri->n; i++){
    pst_vertex_data_t* iv = &(gri->vertex[i]);
    if(iv->id != -1){
      if(iv->mark != MARK_VERTEX_REMOVED){
	int32_t ix,iy;
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
  free(gri_prefix);
  free(grj_prefix);
  free(iZ);
  free(jZ);
  free(sZ);
  free(iW);
  free(jW);
  
  
  
  return 0;  
}


options_t* parse_options(int32_t argc, char **argv){
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
  o->NX = (int32_t)argparser_get_next_int(pp, 1, INT64_MAX );
  o->NY = (int32_t)argparser_get_next_int(pp, 1, INT64_MAX );
  
   o->wmag = 1.0;
    if(argparser_keyword_present(pp,"-wmag")){
      o->wmag = argparser_get_next_double(pp,0,DBL_MAX);
    }
  
  o->nWeights = 0;
  o->weights = NULL;
  if(argparser_keyword_present(pp, "-weights")){
    o->nWeights = argparser_get_next_int(pp, 1, INT64_MAX );
    o->weights = (double*)malloc(sizeof(double)*(o->nWeights));
    int32_t i;
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
