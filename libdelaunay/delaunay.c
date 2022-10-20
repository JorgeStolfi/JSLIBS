/* Delaunay triangulation by straightline divide-and-conquer. */
/* Last edited on 2011-12-25 02:13:19 by stolfi */ 

/* 
** Written by J. Stolfi on april 1993, based on an original
** implementation by Jim Roth (DEC CADM Advanced Group, May 1986).  
** See the copyright notice at the end of this file.
*/

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <quad.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>
#include <sign_get.h>
#include <i3.h>
#include <hi2.h>
#include <r2.h>

#include <delaunay.h>
#include <delaunay_debug.h>

#define MC (delaunay_MAX_COORD)

#define DEBUG (DEBUGGING_DELAUNAY || DEBUGGING_QUAD)

/* Internal prototypes: */

sign_t delaunay_cmp_sites(delaunay_site_t *a, delaunay_site_t *b);
  /* Returns {-1}, {0}, or {+1} if the Cartesian X coordinate of
    {a} is less than,equal to, or greater than that of {b},
    respectively. */

void sort_sites(delaunay_site_t sites[], int nsites);

void rec_delaunay(
    delaunay_site_t sites[], /* The sites, sorted left to right. */
    int sl,                  /* Index of first site to consider. */
    int sh,                  /* Index of first site not to consider. */
    quad_arc_t *le,          /* Output: leftmost edge. */
    quad_arc_t *re,          /* Output: rightmost edge. */
    int indent               /* Indentation of debugging messages. */
  );
  /* Recursively create the Delaunay triangulation of {sites[sl..sh-1]}.
    Returns the leftmost edge {le} and the rightmost edge {re} of the 
    trangulation. */

void delaunay_find_common_tangent(quad_arc_t *ldiP, quad_arc_t *rdiP, int indent);
  /* Finds the lower common tangent {T} between two Delaunay diagrams {L,R}.
    Assumes that they are separated by a straight vertical line.
    On input and output, the right face of {*ldiP} and the left face of
    {*rdiP} are the exterior faces of the respective diagrams.
    
    On input, {ORG(*ldiP)} should be the rightmost site of {L};
    On input, {ORG(*rdiP)} should be the leftmost site of {R}.
    
    On output, {ORG(*ldiP)} is the rightmost site of {L} on {T}, and
    {ORG(*rdiP)} is the leftmost site of {R} on {T}. Thus the line
    segment from {ORG(*ldiP)} to {ORG(*rdiP)} is an outer edge of the
    new hull.
    
    This procedure assumes that there are no coincident sites. If
    there are {L} and {R} sites on the separating line (with the same
    {X} coordinate), all those {L} sites must have lower {Y}
    coordinate than tose {R} sites. This is equivalent to assume that
    the sites are perturbed by an infinitesimal horizontal shear {X ->
    X+eps*Y}, or that the separating line is rotated clockwise by an
    infinitesimal angle. Thus, "leftmost" means "highest among the leftmost"
    and "rightmost" means "lowest among the rightmost". */

/* Main procedure: */

quad_arc_t delaunay_build(delaunay_site_t sites[], int nsites)
  {
    quad_arc_t le, re;
    bool_t verbose = DEBUG;
    demand(nsites > 1, "cannot build the delaunay for a single site");
    if (verbose) { fprintf(stderr, "Sorting sites...\n"); }
    sort_sites(sites, nsites);
    if (verbose) { fprintf(stderr, "Recursive build...\n"); }
    rec_delaunay(sites, 0, nsites, &le, &re, 1);
    return (le);
  }

/* Shell-sort the sites into x order, breaking ties by y: */

sign_t delaunay_cmp_sites(delaunay_site_t *a, delaunay_site_t *b)
  {
    int64_t a0 = (int64_t)a->pt.c.c[0];
    int64_t a1 = (int64_t)a->pt.c.c[1];
    int64_t b0 = (int64_t)b->pt.c.c[0];
    int64_t b1 = (int64_t)b->pt.c.c[1];
    int64_t dx = a1*b0 - b1*a0;
    sign_t s = sign_int64(dx);
    if (s == 0)
      { int64_t a2 = (int64_t)a->pt.c.c[2];
        int64_t b2 = (int64_t)b->pt.c.c[2];
        int64_t dy = a2*b0 - b2*a0;
        s = sign_int64(dy);
      }
    return s;
  }

void sort_sites(delaunay_site_t sites[], int nsites)
  {
    int gap;
    for (gap = nsites/2; gap > 0; gap /= 2)
      { int i;
        for (i = gap; i < nsites; i++)
          { int j;
            for (
                j = i-gap; 
                (j >= 0) && (delaunay_cmp_sites(&(sites[j]), &(sites[j+gap])) > 0);
                j -= gap
              ) 
              {
                delaunay_site_t tmp = sites[j]; sites[j] = sites[j+gap]; sites[j+gap] = tmp;
              }
          }
      }
  }

/* Connect two vertices with a new edge: */

quad_arc_t delaunay_make_edge(void)
  { 
    if (DEBUG)
      { return deldebug_make_edge(); }
    else
      { return quad_make_edge(); }
  }

void delaunay_destroy_edge(quad_arc_t e)
  {
    if (DEBUG)
      { deldebug_destroy_edge(e); }
    else
      { quad_destroy_edge(e); }
  }

quad_arc_t delaunay_connect(quad_arc_t a, quad_arc_t b)
  {
    quad_arc_t e0 = delaunay_make_edge();
    quad_arc_t e1 = quad_sym(e0);
    quad_arc_t c = quad_lnext(a);
    SET_ORG(e0, ORG(c));
    SET_ORG(e1, ORG(b));
    quad_splice(e0, c);
    quad_splice(e1, b);
    return e0;
  }

void rec_delaunay(
    delaunay_site_t sites[],
    int sl, int sh,
    quad_arc_t *le, quad_arc_t *re,
    int indent
  )
  {
    bool_t verbose = DEBUG;
    if (verbose) { fprintf(stderr, "%*s+ rec_delaunay [%d..%d]...\n", indent, "", sl, sh-1); }
    assert(sh > sl+1);
    if (sh == sl+2) 
      {
	quad_arc_t a = delaunay_make_edge();
	SET_ORG(a, &sites[sl]); 
        SET_DST(a, &sites[sl+1]);
	*le = a; *re = quad_sym(a);
      }
    else if (sh == sl+3) 
      {
	quad_arc_t a = delaunay_make_edge();
	quad_arc_t b = delaunay_make_edge();
	int ct = delaunay_orient(&sites[sl], &sites[sl+1], &sites[sl+2]);
	quad_splice(quad_sym(a), b);
	SET_ORG(a, &sites[sl]); 
        SET_DST(a, &sites[sl+1]);
	SET_ORG(b, &sites[sl+1]);  
        SET_DST(b, &sites[sl+2]);
	if (ct == 0.0) 
	  { *le = a; *re = quad_sym(b); }
	else 
	  { quad_arc_t c = delaunay_connect(b, a);
	    if (ct > 0.0) 
	      { *le = a; *re = quad_sym(b); }
	    else 
	      { *le = quad_sym(c); *re = c; }
	  }
      }
    else
      {
	quad_arc_t ldo, ldi, rdi, rdo;
	quad_arc_t basel, lcand, rcand;

        int sm = (sl+sh)/2;

        rec_delaunay(sites, sl, sm, &ldo, &ldi, indent+2);
	rec_delaunay(sites, sm, sh, &rdi, &rdo, indent+2);
        
        delaunay_find_common_tangent(&ldi, &rdi, indent);

	basel = delaunay_connect(quad_sym(rdi), ldi);
        if (ORG(ldi) == ORG(ldo)) ldo = quad_sym(basel);
	if (ORG(rdi) == ORG(rdo)) rdo = basel;

	if (verbose) { fprintf(stderr, "%*s  stitching...\n", indent, ""); }
	while (1) 
          {
	    lcand = quad_onext(quad_sym(basel));
	    if (DEBUGGING_QUAD) { deldebug_check_quad_edge(lcand); }
            if (delaunay_rightof(DST(lcand), basel))
	      while (delaunay_incircle(DST(basel), ORG(basel), DST(lcand), DST(quad_onext(lcand)))) 
                { quad_arc_t t = quad_onext(lcand); 
                  delaunay_destroy_edge(lcand); 
                  lcand = t;
                  if (DEBUGGING_QUAD) { deldebug_check_quad_edge(lcand); }
                }

	    rcand = quad_oprev(basel);
	    if (DEBUGGING_QUAD) { deldebug_check_quad_edge(rcand); }
            if (delaunay_rightof(DST(rcand), basel))
	      while (delaunay_incircle(DST(basel), ORG(basel), DST(rcand), DST(quad_oprev(rcand)))) 
                { quad_arc_t t = quad_oprev(rcand); 
                  delaunay_destroy_edge(rcand);
                  rcand = t;
                  if (DEBUGGING_QUAD) { deldebug_check_quad_edge(rcand); }
                }

	    if ((!delaunay_rightof(DST(lcand), basel)) && (!delaunay_rightof(DST(rcand), basel)))
              break;

	    if ( ( ! delaunay_rightof(DST(lcand), basel) ) ||
		 ( delaunay_rightof(DST(rcand), basel) && 
                   delaunay_incircle(DST(lcand), ORG(lcand), ORG(rcand), DST(rcand))
                 )
               )
	      basel = delaunay_connect(rcand, quad_sym(basel));
	    else
	      basel = delaunay_connect(quad_sym(basel), quad_sym(lcand));
            if (DEBUGGING_QUAD) { deldebug_check_quad_edge(basel); }
            if (DEBUGGING_DELAUNAY) 
              { deldebug_set_all_left_faces(basel, NULL);
                deldebug_check_left_triangle(basel, 1);
              }
	  }
	*le = ldo; *re = rdo;
      }
    if (DEBUGGING_QUAD) { deldebug_check_quad_all_edges(); }
    if (verbose) { fprintf(stderr, "%*s- rec_delaunay...\n", indent, ""); }
  }

void delaunay_find_common_tangent(quad_arc_t *ldiP, quad_arc_t *rdiP, int indent)
  {
    quad_arc_t ldi = (*ldiP);
    quad_arc_t rdi = (*rdiP);
    bool_t verbose = DEBUG;
    if (verbose) { fprintf(stderr, "%*s  finding common tangent...\n", indent, ""); }
    while (1) 
      {
        if (verbose) 
          { fprintf(stderr, "%*s", indent+2, ""); deldebug_print_edge("ldi", ldi); fprintf(stderr, "\n");
            fprintf(stderr, "%*s", indent+2, ""); deldebug_print_edge("rdi", rdi); fprintf(stderr, "\n\n");
          }
        if (delaunay_leftof(ORG(rdi), ldi)) ldi = quad_lnext(ldi);
        else if (delaunay_rightof(ORG(ldi), rdi)) rdi = quad_onext(quad_sym(rdi));
        else break;
      }
    (*ldiP) = ldi;
    (*rdiP) = rdi;
  }

r2_t delaunay_r2_from_hi2(hi2_point_t *ph)
  {  
    double h0 = (double)(ph->c.c[0]);
    double h1 = (double)(ph->c.c[1]);
    double h2 = (double)(ph->c.c[2]);
    return (r2_t) {{ h1/h0, h2/h0 }};
  }

hi2_point_t delaunay_hi2_from_r2(r2_t *pc)
  {  
    double px = pc->c[0];
    double py = pc->c[1];
    double pm = fmax(fabs(px),fabs(py));
    int32_t hw;
    if (pm < 1.0)
      { hw = MC; }
    else 
      { hw = (int)floor(((double)MC)/pm);
        demand(hw > 0, "cartesian coordinates are too large");
      }
    int32_t hx = (int)(floor(px*(double)hw + 0.5 + MC)) - MC; assert(abs(hx) <= MC);
    int32_t hy = (int)(floor(py*(double)hw + 0.5 + MC)) - MC; assert(abs(hy) <= MC);
    return (hi2_point_t){{{ hw, hx, hy }}};
  }


bool_t delaunay_rightof(delaunay_site_t *s, quad_arc_t e)
  {
    return delaunay_orient(s, DST(e), ORG(e)) > 0;
  }

bool_t delaunay_leftof(delaunay_site_t *s, quad_arc_t e)
  {
    return delaunay_orient(s, ORG(e), DST(e)) > 0;
  }

sign_t delaunay_orient(delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c)
  {
    return hi2_orient(&(a->pt), &(b->pt), &(c->pt));
  }

bool_t delaunay_incircle(delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c, delaunay_site_t *d)
  {
    if ((a == b) || (a == c) || (a == d) || (b == c) || (b == d) || (c == d)) { return FALSE; }
    return (hi2_in_circle(&(a->pt), &(b->pt), &(c->pt), &(d->pt)) > 0);
  }

/*
** Copyright notice:
**
** Copyright 1996 Institute of Computing, Unicamp.
**
** Permission to use this software for any purpose is hereby granted,
** provided that any substantial copy or mechanically derived version
** of this file that is made available to other parties is accompanied
** by this copyright notice in full, and is distributed under these same
** terms. 
**
** NOTE: this copyright notice does not claim to supersede any copyrights
** that may apply to the original DEC implementation of the quad-edge
** data structure.
**
** DISCLAIMER: This software is provided "as is" with no explicit or
** implicit warranty of any kind.  Neither the authors nor their
** employers can be held responsible for any losses or damages
** that might be attributed to its use.
**
** End of copyright notice.
*/
 
