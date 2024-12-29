/* See {pst_img_graph_from_maps.h} */
/* Last edited on 2024-12-25 08:14:18 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <haf.h>
#include <float_image.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <rn.h>

#include <pst_imgsys.h>
#include <pst_img_graph.h>

#include <pst_img_graph_from_maps.h>

void pst_img_graph_get_diagonal_edge_data_from_maps(
  float_image_t* IG,
  float_image_t* IW,
  int32_t x, int32_t y,
  int32_t dirx,int32_t diry,
  double *d, double*w
  );


void get_edge_parameters(
  int32_t NX, int32_t NY,
  float_image_t* IG,
  int32_t gradientFunction,
  float_image_t *IW,
  int32_t weightFunction,
  r2_t u, r2_t v,
  double *d, double *w
  );

/* IMPLEMENTATIONS */

int32_t pst_img_graph_get_vertex_index_from_image_indices(int32_t ix, int32_t iy, int32_t NX, int32_t NY)
  { if ((iy < 0) || (iy > NY)) { return -1; }
    if ((ix < 0) || (ix > NX)) { return -1; }
    return iy*NX + ix;
  }

void pst_img_graph_get_vertex_image_indices(r2_t *p,int32_t NX, int32_t NY, int32_t *ix, int32_t *iy)
  { (*ix) = (int32_t)floor(p->c[0] + 0.5);
    (*iy) = (int32_t)floor(p->c[1] + 0.5);
    if ((*ix) >= NX ) { (*ix) = NX -1; }
    if ((*iy) >= NY ) { (*iy) = NY -1; }
    if ((*ix) < 0 ) { (*ix) = 0; }
    if ((*iy) < 0 ) { (*iy) = 0; }
  }

pst_img_graph_t* pst_img_graph_create_from_gradient_and_weight_maps(float_image_t* IG, float_image_t* IW, bool_t add_diags){
  
  int32_t NX ;
  int32_t NY ;
  if( (IG != NULL) || (IW != NULL) ){
   NX = (IG == NULL ? IW->sz[1]: IG->sz[1]);
   NY = (IG == NULL ? IW->sz[2]: IG->sz[2]);
   if( IG != NULL ) assert(IG->sz[0] == 2);
   if( IW != NULL ) assert(IW->sz[0] == 1);
   if( (IG != NULL) && (IW != NULL) ){
     assert((IG->sz[1] == IW->sz[1]) && (IG->sz[2] == IW->sz[2]));
   }
  }else{
    NX = o->NX;
    NY = o->NY;
  }

  int32_t NX_Z = NX+1;
  int32_t NY_Z = NY+1;
  
  auto oct_arc_t insert_edge(pst_img_graph_t* ig, int32_t orgx, int32_t orgy, int32_t dstx, int32_t dsty);
  
  oct_arc_t insert_edge(pst_img_graph_t* ig, int32_t orgx,int32_t orgy, int32_t dstx, int32_t dsty){
    int32_t org = pst_img_graph_get_vertex_index_from_image_indices(orgx,orgy,NX_Z,NY_Z);
    int32_t dst = pst_img_graph_get_vertex_index_from_image_indices(dstx,dsty,NX_Z,NY_Z);
    double we,de;
    char* label = NULL;
//     char *label = jsprintf("From %ld (%ld,%ld) to %ld(%ld,%ld)", org,orgx,orgy, dst,dstx,dsty);
    r2_t rorg = (r2_t) {{ orgx,orgy }};
    r2_t rdst = (r2_t) {{ dstx,dsty }};
    get_edge_parameters(NX,NY,IG,o->gradientFunction,IW,o->weightFunction,rorg,rdst,&de,&we);
    return  pst_img_graph_edge_insert(ig,org,dst,de,we,label,pst_path_create_empty());
    
  }
  
  

  if( (NX == 0) && (NY == 0) ){ return NULL; }
  

  pst_img_graph_t* g = pst_img_graph_create(((NX+1)*(NY)) + ((NX)*(NY+1)), NX_Z*NY_Z);
  
  
  int32_t x,y;
   for(y = 0; y < NY_Z; y++){
    for(x = 0; x < NX_Z; x++){
      int32_t id = pst_img_graph_get_vertex_index_from_image_indices(x,y,NX_Z,NY_Z);
      r2_t coords = (r2_t){{x,y}};
      pst_img_graph_vertex_add(g,id,oct_NULL,coords);
    }
  }

//   oct_arc_t a = pst_img_graph_insert_edge_from_gradient_weights(g,IG,IW,0,0,X_AXIS,+1);
  oct_arc_t a = insert_edge(g,0,0,1,0);
  oct_arc_t edge_y = a;
    
  for( x= 1; x < NX; x++ ){
//     oct_arc_t e = pst_img_graph_insert_edge_from_gradient_weights(g,IG,IW,x,0,X_AXIS,+1);
    oct_arc_t e = insert_edge(g,x,0,x+1,0);
    oct_splice(oct_sym(edge_y),e);
    edge_y = e;
  }
  
  for( y = 0; y < NY ; y++){
//     oct_arc_t e = pst_img_graph_insert_edge_from_gradient_weights(g,IG,IW,0,y,Y_AXIS,+1);
    oct_arc_t e = insert_edge(g,0,y,0,y+1);
    oct_splice(a,e);
    oct_arc_t b = e;
    oct_arc_t c = a;
    for( x  = 0; x < NX; x++){
//       oct_arc_t ev = pst_img_graph_insert_edge_from_gradient_weights(g,IG,IW,x+1,y+0,Y_AXIS,+1);
      oct_arc_t ev = insert_edge(g,x+1,y,x+1,y+1);
//       oct_arc_t eh = pst_img_graph_insert_edge_from_gradient_weights(g,IG,IW,x+1,y+1,X_AXIS,-1);
      oct_arc_t eh = insert_edge(g,x+1,y+1,x,y+1);
      oct_arc_t t = oct_lnext(c);
      oct_splice(t,ev);
      oct_splice(oct_sym(ev),eh);
      oct_splice(oct_sym(b), oct_sym(eh));
      if(o->addDiags){
// 	oct_arc_t ed = pst_img_graph_insert_diagonal_edge_from_gradient_weights(g,IG,IW,x+1,y,-1,+1);
	oct_arc_t ed = insert_edge(g,x+1,y,x,y+1);
	oct_splice(ed,ev);
	oct_splice(oct_sym(ed),oct_sym(b));
	
      }
      
      c = t;
      b = ev;
    }
      if(!o->addDiags){
	a = oct_sym(oct_lprev(oct_lprev(a)));
      }else{
	a = oct_onext(oct_sym(oct_lnext(a)));
      }
  }
  /*now we have to clean the table from the null edges*/
  int32_t valid_edges = 0;  
  int32_t i;
  for(	i = 0; i < g->m; i++){
    if(g->edge[i].data->weight > 0){ valid_edges++;}
  }
  
  if( valid_edges != g->m){
    fprintf(stderr,"Fixing...\n");
    int32_t count_valid = 0;
    for(i = 0; i < g->m; i++){
      if(g->edge[i].data->weight > 0){
	count_valid++;
      }else{
	pst_img_graph_edge_remove(g,g->edge[i].edge);
      }
    }
  }

  return g;
}

void get_edge_parameters(
  int32_t NX, int32_t NY,
  float_image_t* IG,
  int32_t gradientFunction,
  float_image_t *IW,
  int32_t weightFunction,
  r2_t u, r2_t v,
  double *d, double *w
  ){
  
  auto void get_nearest_grid_corner(r2_t pt, int32_t NX, int32_t NY, int32_t* x, int32_t* y);
  void get_nearest_grid_corner(r2_t pt, int32_t NX, int32_t NY, int32_t* x, int32_t* y){
    assert( (pt.c[0] >= 0) && (pt.c[0] <= NX ));
    assert( (pt.c[1] >= 0) && (pt.c[1] <= NY ));
    *x = pt.c[0];
    *y = pt.c[1];
  }
  
  if( (IG != NULL ) || (IW != NULL) )
  {
      NX =  (IG == NULL ?  IW->sz[1] : IG->sz[1]);
      NY =  (IG == NULL ?  IW->sz[2] : IG->sz[2]);
     if( (IG != NULL) && (IW != NULL) ){ assert( (IG->sz[1] == IW->sz[1]) && (IG->sz[2] == IW->sz[2]));}
     
      int32_t xu,yu,xv,yv;
      get_nearest_grid_corner(u,NX,NY,&xu,&yu);
      get_nearest_grid_corner(v,NX,NY,&xv,&yv);
      int32_t dx = xv - xu;
      int32_t dy = yv - yu;
      if( dx == 0) { pst_img_graph_get_axial_edge_data_from_maps(IG,IW,xu,yu,Y_AXIS,dy,d,w); }
      else if( dy == 0) { pst_img_graph_get_axial_edge_data_from_maps(IG,IW,xu,yu,X_AXIS,dx,d,w); }
      else{ pst_img_graph_get_diagonal_edge_data_from_maps(IG,IW,xu,yu,dx,dy,d,w); }
	
  }
  else{ *d = 0; *w =  1/r2_dist_sqr(&u,&v);  }
  r2_t p;
  r2_mix(0.5,&u,0.5,&v,&p);
  double ww;
  
  switch(weightFunction){
    case WGHT_FUNC_NONE :
      break;
    case WGHT_FUNC_CONST :
      *w = 2.5;
      break;
    case WGHT_FUNC_RAMP :
      *w = (0.25/NX)*p.c[0] + (0.75/NY)*p.c[1] ;
      break;
    case WGHT_FUNC_RANDOM :
      ww = (sin(17*p.c[0]*p.c[0]) + cos(29*p.c[1]*p.c[1]))*47;
      (*w) = ww - floor(ww);
      break;
    default:
      demand(FALSE,"Invalid Weight function");
  };
  
  r2_t uv;
  r2_sub(&v,&u,&uv);
  double dx,dy,dd;
  switch(gradientFunction){
    case GRAD_FUNC_NONE :
      break;
    case GRAD_FUNC_CONST :
      *d = 2.5;
      break;
    case GRAD_FUNC_RAMP :
      dx = 2*((0.25/NX)*p.c[0] + (0.75/NY)*p.c[1]) -1 ;
      dy = 2*((0.75/NX)*p.c[0] + (0.25/NY)*p.c[1]) -1 ;
      *d = dx*uv.c[0] + dy*uv.c[1];
      break;
    case GRAD_FUNC_RANDOM :
      dd = (sin(17*p.c[0]*p.c[0]) + cos(29*p.c[1]*p.c[1]))*47;
      dd = 2*(dd - floor(dd)) -1;
      dd *= r2_dist(&u,&v);
      (*d)= dd; 
      break;
    default:
      demand(FALSE,"Invalid Grad function");
  };
}


void pst_img_graph_get_diagonal_edge_data_from_maps(
  float_image_t* IG,
  float_image_t* IW,
  int32_t x, int32_t y,
  int32_t dirx,int32_t diry,
  double *d, double*w
  ){
  
  if( (IG == NULL) && (IW == NULL) ){
    *d = 0;
    *w = 0.5;
    return;
  }
  
  assert( (abs(dirx) == 1) && (abs(diry) == 1) );
  int32_t NX =  (IG == NULL ?  IW->sz[1] : IG->sz[1]);
  int32_t NY =  (IG == NULL ?  IW->sz[2] : IG->sz[2]);
  
  int32_t ix = (dirx > 0 ? x:x-1);
  int32_t iy = (diry > 0 ? y:y-1);
  
  assert( (ix >=0) && (ix < NX));
  assert( (iy >=0) && (iy < NY));
  
  double dx = dirx*(IG == NULL ? 0 : float_image_get_sample(IG,0,ix,iy));
  double dy = diry*(IG == NULL ? 0 : float_image_get_sample(IG,1,ix,iy));
  *w = 0.5*(IW == NULL ? 1.0 : float_image_get_sample(IW,0,ix,iy));
 
  *d = dx+dy;
  if( *w == 0){ *d = 0; }
}
