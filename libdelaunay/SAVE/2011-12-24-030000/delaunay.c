/* Delaunay triangulation by straightline divide-and-conquer. */
/* Last edited on 2011-12-24 02:29:15 by stolfilocal */ 

/* 
** Written by J. Stolfi on april 1993, based on an original
** implementation by Jim Roth (DEC CADM Advanced Group, May 1986).  
** See the copyright notice at the end of this file.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
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

#define MC (delaunay_MAX_COORD)

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

void delaunay_check_edge(quad_arc_t e);
  /* Performs consistency checks on the fields of {quad-edge(e)}. */

void delaunay_check_all_edges(void);
  /* Applies {delaunay_check_edge} on all edges created so far. */

quad_arc_t delaunay_debug_make_edge(void);
  /* Same as {quad_make_edge} but also saves the edge for checking. */

void delaunay_debug_destroy_edge(quad_arc_t e);
  /* Dymmy {quad_destroy_edge} that does not actually free {e}. */

/* Main procedure: */

quad_arc_t delaunay_build(delaunay_site_t sites[], int nsites)
  {
    quad_arc_t le, re;
    bool_t debug = FALSE;
    demand(nsites > 1, "cannot build the delaunay for a single site");
    if (debug) { fprintf(stderr, "Sorting sites...\n"); }
    sort_sites(sites, nsites);
    if (debug) { fprintf(stderr, "Recursive build...\n"); }
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

/* Debugging: */

static int nec = 0; /* Number of edges created so far. */
static quad_arc_vec_t ec; /* Edges created so far are {ex.e[0..nec-1]}. */

quad_arc_t delaunay_debug_make_edge(void)
  { 
    quad_arc_t e = quad_make_edge();
    if (nec == 0) { ec = quad_arc_vec_new(100); }
    quad_arc_vec_expand(&ec, nec);
    ec.e[nec] = e;
    nec++;
    return e;
  }

void delaunay_debug_destroy_edge(quad_arc_t e)
  {
    delaunay_check_all_edges();
    return;
    quad_destroy_edge(e);
    int i;
    int m = 0;
    for (i = 0; i < nec; i++) 
      { if (quad_edge(e) == quad_edge(ec.e[i])) 
          { ec.e[i] = quad_arc_NULL; m++; }
      }
    assert(m == 1);
    delaunay_check_all_edges();
    return;
  }

/* Connect two vertices with a new edge: */

quad_arc_t connect(quad_arc_t a, quad_arc_t b)
  {
    quad_arc_t e;

    e = delaunay_debug_make_edge();
    SET_ORG(e, DST(a));
    SET_DST(e, ORG(b));
    quad_arc_t ta = quad_lnext(a);
    quad_arc_t te = quad_sym(e);
    quad_splice(e, ta);
    quad_splice(te, b);
    return e;
  }

void rec_delaunay(
    delaunay_site_t sites[],
    int sl, int sh,
    quad_arc_t *le, quad_arc_t *re,
    int indent
  )
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "%*s+ rec_delaunay [%d..%d]...\n", indent, "", sl, sh-1); }
    assert(sh > sl+1);
    if (sh == sl+2) 
      {
	quad_arc_t a = delaunay_debug_make_edge();
	SET_ORG(a, &sites[sl]); 
        SET_DST(a, &sites[sl+1]);
	*le = a; *re = quad_sym(a);
      }
    else if (sh == sl+3) 
      {
	quad_arc_t a = delaunay_debug_make_edge();
	quad_arc_t b = delaunay_debug_make_edge();
	int ct = delaunay_orient(&sites[sl], &sites[sl+1], &sites[sl+2]);
	quad_splice(quad_sym(a), b);
	SET_ORG(a, &sites[sl]); 
        SET_DST(a, &sites[sl+1]);
	SET_ORG(b, &sites[sl+1]);  
        SET_DST(b, &sites[sl+2]);
	if (ct == 0.0) 
	  { *le = a; *re = quad_sym(b); }
	else 
	  { quad_arc_t c = connect(b, a);
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
        if (debug) { fprintf(stderr, "%*s  finding common tangent...\n", indent, ""); }
	while (1) 
          {
	    if (debug) { fprintf(stderr, "%*s", indent, ""); delaunay_debug_edge("ldi", ldi); fprintf(stderr, "\n"); }
            if (debug) { fprintf(stderr, "%*s", indent, ""); delaunay_debug_edge("rdi", rdi); fprintf(stderr, "\n\n"); }
            if (delaunay_leftof(ORG(rdi), ldi)) ldi = quad_lnext(ldi);
	    else if (delaunay_rightof(ORG(ldi), rdi)) rdi = quad_onext(quad_sym(rdi));
	    else break;
	  }

	basel = connect(quad_sym(rdi), ldi);
	delaunay_check_edge(basel);
        if (ORG(ldi) == ORG(ldo)) ldo = quad_sym(basel);
	if (ORG(rdi) == ORG(rdo)) rdo = basel;

	if (debug) { fprintf(stderr, "%*s  stitching...\n", indent, ""); }
	while (1) 
          {
	    lcand = quad_onext(quad_sym(basel));
	    delaunay_check_edge(lcand);
            if (delaunay_rightof(DST(lcand), basel))
	      while (delaunay_incircle(DST(basel), ORG(basel), DST(lcand), DST(quad_onext(lcand)))) 
                { quad_arc_t t = quad_onext(lcand); 
                  delaunay_debug_destroy_edge(lcand); 
                  lcand = t;
                  delaunay_check_edge(lcand);
                }

	    rcand = quad_oprev(basel);
	    delaunay_check_edge(rcand);
            if (delaunay_rightof(DST(rcand), basel))
	      while (delaunay_incircle(DST(basel), ORG(basel), DST(rcand), DST(quad_oprev(rcand)))) 
                { quad_arc_t t = quad_oprev(rcand); 
                  delaunay_debug_destroy_edge(rcand);
                  rcand = t;
                  delaunay_check_edge(rcand);
                }

	    if (!delaunay_rightof(DST(lcand), basel) && !delaunay_rightof(DST(rcand), basel)) break;

	    if ( !delaunay_rightof(DST(lcand), basel) ||
		 ( delaunay_rightof(DST(rcand), basel) && 
                   delaunay_incircle(DST(lcand), ORG(lcand), ORG(rcand), DST(rcand))
                 )
               )
	      basel = connect(rcand, quad_sym(basel));
	    else
	      basel = connect(quad_sym(basel), quad_sym(lcand));
            delaunay_check_edge(basel);
	  }
	*le = ldo; *re = rdo;
      }
    delaunay_check_all_edges();
    if (debug) { fprintf(stderr, "%*s- rec_delaunay...\n", indent, ""); }
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
    if ((a == b) || (a == c) || (a == d) || (b == c) || (b == d) || (c == d)) 
      { return FALSE; }
    return hi2_in_circle(&(a->pt), &(b->pt), &(c->pt), &(d->pt));
  }

/* Debugging: */

void delaunay_debug_site (char *msg, delaunay_site_t *a)
  { fprintf(stderr, "%s site[%d] =", msg, a->index);
    i3_gen_print(stderr, &(a->pt.c), ("%+" delaunay_COORD_FMT), "[ ", " ", " ]");
  }

void delaunay_debug_edge (char *msg, quad_arc_t e)
  { fprintf(stderr, "%s ", msg);
    quad_write_arc(stderr, e, 1);
    delaunay_debug_site("  ORG =", ORG(e)); 
    delaunay_debug_site("  DST =", DST(e)); 
  }

void delaunay_check_edge(quad_arc_t e)
  {
    int k;
    for (k = 0; k < 4; k++)
      { assert(! quad_arc_is_null(quad_onext(e))); 
        e = quad_rot(e);
      }
  }

void delaunay_check_all_edges(void)
  {
    int i;
    for (i = 0; i < nec; i++)
      { quad_arc_t e = ec.e[i];
        if (! quad_arc_is_null(e)) { delaunay_check_edge(e); }
      }
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
 
