
#define PROG_NAME "plot_graph_ps"
#define PROG_DESC "Plots a EPS  from a graph file "
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright � 2005 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-01-10 07:53:26 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -in {GRAPH_FILENAME} \\\n" \
  "    -out {EPS_FILENAME} \\\n" \
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


typedef struct options_t{
  char* graphFile;
  char* epsFile;
  double ballSize;
  bool_t useLabels;
  double NX,NY;
} options_t;

void img_graph_write_ps(PSStream *ps, pst_img_graph_t* g,bool_t useLabels,double ballSize, int32_t NX, int32_t NY);

void img_graph_write_ps(PSStream *ps, pst_img_graph_t* g,bool_t useLabels,double ballSize, int32_t NX, int32_t NY)
{
  
  auto void write_ps_vertex(int32_t index,bool_t label);
  
  void write_ps_vertex(int32_t index,bool_t label)
  {
    pst_vertex_data_t* v = &(g->vertex[index]);
    r2_t coords = v->coords;
    if(label){
      char* id_text = NULL;
      char *id_text = jsprintf("%ld",v->id);
      pswr_label(ps,id_text,coords.c[0],coords.c[1],0.5,0.5);
      free(id_text);
    }
    if(!label){
      if(v->mark == 0){pswr_set_fill_color(ps,1,1,1);}
      if(v->mark == 1){pswr_set_fill_color(ps,1,0,0);}
      if(v->mark == 2){pswr_set_fill_color(ps,0,1,0);}
      pswr_circle(ps,coords.c[0],coords.c[1],ballSize,TRUE,TRUE);
    }
  }
  
  auto void write_ps_edge(oct_arc_t e);
  void write_ps_edge(oct_arc_t e){
    
    int32_t o = pst_img_graph_get_arc_origin(g,e);
    int32_t d = pst_img_graph_get_arc_origin(g,oct_sym(e));
    
    write_ps_vertex(o,FALSE);
    write_ps_vertex(d,FALSE);

    
    r2_t org_coords  = g->vertex[o].coords;
    r2_t dst_coords  = g->vertex[d].coords;
    
    pst_path_t path = pst_img_graph_get_edge_path(g,e);
//     path = pst_path_create_empty();
    int32_t i;
    r2_t  p_medio = (r2_t){{ 0,0 }};
    for(i = 0; i <= path.n; i++){
	
	r2_t p = (i == 0 ? org_coords :pst_path_get_vertex(path,i-1));
	r2_t q = (i == path.n ? dst_coords : pst_path_get_vertex(path,i) );
	if( i == (path.n/2)){
	  if((path.n%2) == 0){
	    r2_mix(0.5,&p,0.5,&q,&p_medio);
	  }else{
	    p_medio = q;
	  }
	}
	pswr_segment(ps,p.c[0],p.c[1],q.c[0],q.c[1]);
    }
    
    if( useLabels) {
      char* id_text;
      
      char *id_text = jsprintf("%ld",pst_img_graph_get_edge_num(e));
      pswr_label(ps,id_text,p_medio.c[0],p_medio.c[1],0.5,0.5);
      free(id_text);
    }
    
  }
  
  /* First draw the segments in blue (non tree)*/
  pswr_set_pen(ps,0,0,1,0.1, 0, 0);
  int32_t i;
  for(i = 0; i < g->m; i++){
    if(g->???[i].aout != oct_NULL){
      write_ps_edge(g->{hedge|dedge}@@[i].aout);
    }
  }
  
  /* Now draw the vertices in black*/
   pswr_set_pen(ps,0,0,0,0.1, 0, 0);
   for(i = 0; i < g->n; i++){
     if(g->vertex[i].id != -1 ){
// 	write_ps_vertex(i,g->n <= 256);
      write_ps_vertex(i,useLabels);
     }
   }
  

}



void writePSA(char* outPrefix,pst_img_graph_t* g,bool_t useLabels,double ballSize, int32_t NX, int32_t NY);
void writePSA(char* outPrefix,pst_img_graph_t* g,bool_t useLabels,double ballSize, int32_t NX, int32_t NY){
//   
    double xMin = -0.20*NX;
    double xMax = +1.20*NX;
    double yMin = -0.20*NY;
    double yMax = +1.20*NY;
    
    
    
    double xSize = xMax - xMin;
    double ySize = yMax - yMin;
    double r = sqrt(ySize/xSize);
    double hSize = 8+320/r;
    double vSize = 8+320*r;
    
    char* filename = NULL;
    char *filename = jsprintf("%s.eps",outPrefix);
    FILE* arq = open_write(filename,TRUE);
    
    fprintf(stderr,"\nScaleX %lf ScaleY %lf\n",(hSize-8)/(xMax - xMin),(vSize -8)/(yMax - yMin));
    
    
    PSStream *ps = pswr_new_stream(NULL, arq, TRUE, "fig", NULL, hSize,vSize);
    pswr_new_canvas(ps, NULL);
    pswr_set_window
        ( ps, 
          xMin,xMax,yMin,yMax,
	  4, hSize -4,4, vSize - 4
   );
   pswr_set_fill_color(ps,1.0,1.0,0.95);
   pswr_rectangle(ps,xMin,xMax,yMin,yMax,TRUE,FALSE);
   img_graph_write_ps(ps,g,useLabels,ballSize, NX,NY);
   pswr_close_stream(ps);
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
    o->useLabels = argparser_keyword_present(pp, "-useLabels");
    
    argparser_get_keyword(pp, "-in");
    o->graphFile = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-out");
    o->epsFile = argparser_get_next_non_keyword(pp);

    o->ballSize = 0.3;
    if(argparser_keyword_present(pp,"-ballSize")){
      o->ballSize = argparser_get_next_double(pp,0,100000);
    }
    
    o->NX = -1;
    o->NY = -1;
    if(argparser_keyword_present(pp,"-size")){
      o->NX = argparser_get_next_double(pp,0,100000);
      o->NY = argparser_get_next_double(pp,0,100000);
    }
    
    
    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }
  
  
int32_t main(int32_t argc, char** argv){
  
  options_t* o = parse_options(argc,argv);
  
  FILE* graph_arq = open_read(o->graphFile,TRUE);
  pst_img_graph_t* g = pst_img_graph_read(graph_arq);
  fclose(graph_arq);
  /*We have to determnie NX,NY looking at the vertices*/
  int32_t NX, NY;
  NX = 0; NY = 0;
  int32_t i;
  if(o->NX == -1){
    for(i = 0; i < g->n; i++){
      pst_vertex_data_t* v = &(g->vertex[i]);
      if( v->coords.c[0] > NX) NX = v->coords.c[0];
      if( v->coords.c[1] > NY) NY = v->coords.c[1];
    }
  }else{
    NX = o->NX;
    NY = o->NY;
  }
  fprintf(stderr,"%ld %ld\n",NX,NY);
  writePSA(o->epsFile,g,o->useLabels,o->ballSize,NX,NY);
  
  return 0;
}
