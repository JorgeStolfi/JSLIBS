/* See {delaunay_plot_POV.h}. */
/* Last edited on 2013-03-18 00:49:59 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <jsstring.h>
#include <quad.h>
#include <bool.h>

#include <delaunay.h>
#include <delaunay_plot_POV.h>

void delaunay_plot_POV_triangles(FILE *wr, quad_arc_t e, double height[], char* texture)
  {
    auto void draw_delaunay_triangle(quad_arc_t e);
    auto void draw_delaunay_triangle_pair(quad_arc_t e);
    
    quad_arc_vec_t er = quad_arc_vec_new(1);
    er.e[0] = e; 
    quad_enum(&er, draw_delaunay_triangle_pair);
    
    void draw_delaunay_triangle_pair(quad_arc_t e)
      { draw_delaunay_triangle(e);
        draw_delaunay_triangle(quad_sym(e));
      }

    void draw_delaunay_triangle(quad_arc_t e)
      { delaunay_site_t *a = ORG(e);
        delaunay_site_t *b = DST(e);
	delaunay_site_t *c = DST(quad_lnext(e));
	if( (a->index > b->index) || (a->index > c->index ) ) { return ; }
	double az,bz,cz;
	if( delaunay_orient(a,b,c) <= 0 )
          { /* Draw triangle on the reference plane: */
	    az = bz = cz = 0;
	  } 
	else
          { /* Draw triangle on surface: */
            az = height[a->index];
            bz = height[b->index];
            cz = height[c->index];	
	  }
        r2_t ap = delaunay_r2_from_hi2(&(a->pt));
        r2_t bp = delaunay_r2_from_hi2(&(b->pt));
        r2_t cp = delaunay_r2_from_hi2(&(c->pt));
        
	fprintf(wr,"    triangle{ ");
	fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", ap.c[0], ap.c[1], az);
	fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", bp.c[0], bp.c[1], bz);
 	fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>",  cp.c[0], cp.c[1], cz);
	fprintf(wr,"  texture{ %s }", texture);
	fprintf(wr," }\n");
	return ;
      }
  }

void delaunay_plot_POV_skirt(FILE *wr, quad_arc_t e,double height[],char* texture)
  {
    auto bool_t triangle_is_CCW(quad_arc_t e);
    auto void draw_silouette_edge(quad_arc_t e);
    auto void draw_delaunay_triangle_pair(quad_arc_t e);
    
    quad_arc_vec_t er = quad_arc_vec_new(1);
    er.e[0] = e; 
    quad_enum(&er, draw_delaunay_triangle_pair);
    
    void draw_delaunay_triangle_pair(quad_arc_t e)
      {	
        bool_t L = triangle_is_CCW(e);
	bool_t R = triangle_is_CCW(quad_sym(e));
	if(L != R) draw_silouette_edge(e);
      }

    bool_t triangle_is_CCW(quad_arc_t e)
      { 
        delaunay_site_t *a = ORG(e);
        delaunay_site_t *b = DST(e);
	delaunay_site_t *c = DST(quad_lnext(e));
	return delaunay_orient(a,b,c) > 0;
      }

    void draw_silouette_edge( quad_arc_t e)
      {
        delaunay_site_t *a = ORG(e); 
        r2_t ap = delaunay_r2_from_hi2(&(a->pt));
        double az = height[a->index];

        delaunay_site_t *b = DST(e); 
        r2_t bp = delaunay_r2_from_hi2(&(b->pt));
        double bz = height[b->index];

        if (bz != 0)
          { fprintf(wr,"    triangle{ ");
            fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", ap.c[0],ap.c[1], az);
            fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", bp.c[0],bp.c[1], bz);
            fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>",  bp.c[0],bp.c[1], 0.0);
            fprintf(wr,"  texture{ %s } ",texture);
            fprintf(wr," }\n");
          }
        if (az != 0)
          {
            fprintf(wr,"    triangle{ ");
            fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", ap.c[0],ap.c[1], az);
            fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", bp.c[0],bp.c[1], 0.0);
            fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>",  ap.c[0],ap.c[1], 0.0);
            fprintf(wr,"  texture{ %s } ",texture);
            fprintf(wr," }\n");
          }
      }
  }

