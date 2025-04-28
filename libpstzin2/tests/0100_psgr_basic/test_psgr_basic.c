#define PROG_NAME "test_gr_basic"
#define PROG_DESC "checks the routines in {psgr.h}, {psgr_plot.h}, {psgr_read.h}, {psgr_write.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-04-27 21:24:04 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_gr_basic_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
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
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2024-12-23 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_gr_basic_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {psgr}."

#define PROG_INFO_OPTS \
  ""

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>
#include <epswr.h>
#include <i2.h>
#include <r2.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_unsafe.h>
#include <psgr_plot.h>
#include <psgr_read.h>
#include <psgr_write.h>
#include <psgr_path.h>
#include <psgr_test.h>
#include <psgr_copy.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED

typedef struct options_t
  { int32_t DUMMY;  /* Placeholder. */
  } options_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

void write_graph(char *fname, psgr_t *gr);
  /* Writes {gr} to file "{fname}". */
  
psgr_t* read_graph(char *fname);
  /* Reads back a graph from file "{fname}". */

psgr_t* test__psgr_test_make_graph__psgr_new__psgr_check_consistency(void);
void test__psgr_num_vertices__psgr_num_edges__psgr_image_size(psgr_t *gr);
void test__psgr_vertex_pixel__psgr_vertex_mark__psgr_vertex_set_mark(psgr_t *gr);
void test__psgr_edge_mark__psgr_edge_set_mark(psgr_t *gr);
void test__psgr_orient_edge__psgr_arc_edge__psgr_arc_dir_bit__psgr_arc_sym(psgr_t *gr);
void test__psgr_arc_org__psgr_arc_dst(psgr_t *gr);
void test__psgr_arc_onext__psgr_arc_oprev__psgr_arc_dnext__psgr_arc_dprev(psgr_t *gr);
void test__psgr_arc_lnext__psgr_arc_lprev(psgr_t *gr);
void test__psgr_arc_split__psgr_arc_unsafe_split(psgr_t *gr);
void test__psgr_arc_delta__psgr_arc_weight__psgr_arc_set_delta__psgr_arc_set_weight(psgr_t *gr);
void test__psgr_vertex_out_arc__psgr_vertex_outdegree(psgr_t *gr);
void test__psgr_arc_path__psgr_arc_starting_dir(psgr_t *gr);
void test__psgr_find_nearest_vertex(psgr_t *gr);
void test__psgr_arc_left_face_curl__psgr_arc_left_face_barycenter(psgr_t *gr);
void test__psgr_copy__psgr_equal__psgr_free(psgr_t *gr);
void test__psgr_add_vertex__psgr_add_edge__psgr_edge_remove__psgr_vertex_remove(psgr_t *gr);
void test__psgr_write_named__psgr_write_file__psgr_read_named__psgr_read_file(psgr_t *gr);
void test__psgr_plot_named__psgr_plot__psgr_plot_path__psgr_plot_coords_from_path(psgr_t *gr);

/* 
psgr_arc_print
psgr_curl_image_fill
psgr_edge_mark_throw
psgr_find_enclosing_face
psgr_get_connecting_arc
psgr_plot_coords_from_pixel
psgr_vertex_from_pixel
psgr_vertex_mark_throw
*/ 

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    /* Create a test graph: */
    psgr_t *gr = test__psgr_test_make_graph__psgr_new__psgr_check_consistency();

    test__psgr_num_vertices__psgr_num_edges__psgr_image_size(gr);
    test__psgr_vertex_pixel__psgr_vertex_mark__psgr_vertex_set_mark(gr);
    test__psgr_edge_mark__psgr_edge_set_mark(gr);
    test__psgr_orient_edge__psgr_arc_edge__psgr_arc_dir_bit__psgr_arc_sym(gr);
    test__psgr_arc_org__psgr_arc_dst(gr);
    test__psgr_arc_onext__psgr_arc_oprev__psgr_arc_dnext__psgr_arc_dprev(gr);
    test__psgr_arc_lnext__psgr_arc_lprev(gr);
    test__psgr_arc_split__psgr_arc_unsafe_split(gr);
    test__psgr_arc_delta__psgr_arc_weight__psgr_arc_set_delta__psgr_arc_set_weight(gr);
    test__psgr_vertex_out_arc__psgr_vertex_outdegree(gr);
    test__psgr_arc_path__psgr_arc_starting_dir(gr);
    test__psgr_find_nearest_vertex(gr);
    test__psgr_arc_left_face_curl__psgr_arc_left_face_barycenter(gr);
    test__psgr_copy__psgr_equal__psgr_free(gr);
    test__psgr_add_vertex__psgr_add_edge__psgr_edge_remove__psgr_vertex_remove(gr);
    test__psgr_write_named__psgr_write_file__psgr_read_named__psgr_read_file(gr);
    test__psgr_plot_named__psgr_plot__psgr_plot_path__psgr_plot_coords_from_path(gr);
    /* Cleanup: */
    psgr_free(gr);
    free(o); o = NULL;
    return 0;
  }
    
psgr_t* test__psgr_test_make_graph__psgr_new__psgr_check_consistency(void)
  { 
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t nf = 2; /* Vertex layers in frame. */
    uint32_t nh = 3; /* Verted rows and cols in hole. */
    psgr_t *gr = psgr_test_make_graph(nf, nh);

    psgr_check_consistency(gr);
    return gr;
  }
    
void test__psgr_plot_named__psgr_plot__psgr_plot_path__psgr_plot_coords_from_path(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    double pixelSize = 5.0;    /* Pixel size in mm. */
    double vertexRadius = 0.5; /* Vertex radius in mm. */
    psgr_plot_named("out/plot-lab0.eps", gr, pixelSize, 0.0, vertexRadius);
    psgr_plot_named("out/plot-lab1.eps", gr, pixelSize, 3.0, vertexRadius);
    /* Calling {psgr_plot_named} will call {psgr_plot} too. */
    /* Calling {psgr_plot} will call {psgr_plot_path} too. */
    /* Calling {psgr_plot_path} will call {psgr_plot_coords_from_path} too. */
  }
  
void test__psgr_num_vertices__psgr_num_edges__psgr_image_size(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    demand(NE == gr->NE, "psgr_num_edges bug");
    
    uint32_t NV = psgr_num_vertices(gr);
    demand(NV == gr->NV, "psgr_num_vertices bug");

    i2_t sz = psgr_image_size(gr);
    int32_t NX = sz.c[0], NY = sz.c[1];
    demand((NX >= 0) && (NX <= psgr_MAX_IMG_SIZE), "invalid image X size");
    demand((NY >= 0) && (NY <= psgr_MAX_IMG_SIZE), "invalid image Y size");
  }
  
void test__psgr_vertex_pixel__psgr_vertex_mark__psgr_vertex_set_mark(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NV = psgr_num_vertices(gr);
    i2_t sz = psgr_image_size(gr);
    int32_t NX = sz.c[0], NY = sz.c[1];
    
    /* Check vertex pixel indices and marks: */
    for (uint32_t kv = 0; kv < NV; kv++)
      { /* Vertex should not be {DELETED} currently: */
        psgr_vertex_t v = (psgr_vertex_t){ kv };
        psgr_mark_t mkv = psgr_vertex_mark(gr, v);
        demand(mkv != DELETED, "vertex is marked as deleted");
        
        /* Test setting and reading the mark: */
        psgr_vertex_set_mark(gr, v, DELETED);
        demand(psgr_vertex_mark(gr, v) == DELETED, "vertex set_mark bug (1)");
        /* Restore the original mark: */
        psgr_vertex_set_mark(gr, v, mkv);
        demand(psgr_vertex_mark(gr, v) == mkv, "vertex set_mark bug (2)");
        
        /* Chech the pixel indices: */
        i2_t ipk = psgr_vertex_pixel(gr, v);
        demand(i2_eq(&ipk, &(gr->vdata[kv].pixel)), "psgr_vertex_pixel bug");
        int32_t x = ipk.c[0], y = ipk.c[1];
        demand((x >= 0) && (x <= NX-1), "bad X pixel index");
        demand((y >= 0) && (y <= NY-1), "bad Y pixel index");
      }
  }
  
void test__psgr_edge_mark__psgr_edge_set_mark(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    
    /* Check edge pixel indices and marks: */
    for (uint32_t ke = 0; ke < NE; ke++)
      { /* Edge should not be {DELETED} currently: */
        psgr_edge_t e = (psgr_edge_t){ ke };
        psgr_mark_t mke = psgr_edge_mark(gr, e);
        demand(mke != DELETED, "edge is marked as deleted");
        
        /* Test setting and reading the mark: */
        psgr_edge_set_mark(gr, e, DELETED);
        demand(psgr_edge_mark(gr, e) == DELETED, "edge set_mark bug (1)");
        /* Restore the original mark: */
        psgr_edge_set_mark(gr, e, mke);
        demand(psgr_edge_mark(gr, e) == mke, "edge set_mark bug (2)");
      }
  }

void test__psgr_orient_edge__psgr_arc_edge__psgr_arc_dir_bit__psgr_arc_sym(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        demand(psgr_edge_mark(gr, e) != DELETED, "edge is marked as deleted");
      
        psgr_arc_t a = psgr_orient_edge(e, 0);
        demand(psgr_arc_edge(a).ix == e.ix, "psgr_arc_edge bug");
        demand(psgr_arc_dir_bit(a) == 0, "psgr_arc_dir_bit bug");

        psgr_arc_t b = psgr_arc_sym(a);
        demand(a.ix != b.ix, "psgr_arc_sym has a fixed arg");
        demand(psgr_arc_edge(b).ix == e.ix, "psgr_arc_sym edge bug");
        demand(psgr_arc_dir_bit(b) == 1, "psgr_arc_sym dir bit bug");
        demand(psgr_arc_sym(b).ix == a.ix, "psgr_arc_sym not involution");
      }
  }
        
void test__psgr_arc_org__psgr_arc_dst(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NV = psgr_num_vertices(gr);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        demand(psgr_edge_mark(gr, e) != DELETED, "edge is marked as deleted");
        psgr_arc_t a = psgr_orient_edge(e, 0);
        psgr_arc_t b = psgr_arc_sym(a);
        
        psgr_vertex_t org = psgr_arc_org(gr, a); demand(org.ix < NV, "invalid org");
        psgr_vertex_t dst = psgr_arc_dst(gr, a); demand(dst.ix < NV, "invalid dst");
        demand(org.ix != dst.ix, "graph has a loop");
        demand(psgr_arc_org(gr, b).ix == dst.ix, "org(sym) bug");
        demand(psgr_arc_dst(gr, b).ix == org.ix, "dst(sym) bug");
      }
  }
        
void test__psgr_arc_onext__psgr_arc_oprev__psgr_arc_dnext__psgr_arc_dprev(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        demand(psgr_edge_mark(gr, e) != DELETED, "edge is marked as deleted");
        psgr_arc_t a = psgr_orient_edge(e, 0);
        psgr_arc_t b = psgr_arc_sym(a);
        
        psgr_arc_t an = psgr_arc_onext(gr, a);
        psgr_arc_t ap = psgr_arc_oprev(gr, a);
        psgr_arc_t bn = psgr_arc_onext(gr, b);
        psgr_arc_t bp = psgr_arc_oprev(gr, b);

        demand(psgr_arc_oprev(gr, an).ix == a.ix, "oprev(onext()) bug");
        demand(psgr_arc_onext(gr, ap).ix == a.ix, "onext(oprev()) bug");
        demand(psgr_arc_oprev(gr, bn).ix == b.ix, "oprev(onext()) bug");
        demand(psgr_arc_onext(gr, bp).ix == b.ix, "onext(oprev()) bug");

        demand(psgr_arc_dprev(gr, a).ix == psgr_arc_sym(bp).ix, "dprev():sym(oprev(sym())) bug");
        demand(psgr_arc_dnext(gr, a).ix == psgr_arc_sym(bn).ix, "dnext():sym(onext(sym())) bug");
        demand(psgr_arc_dprev(gr, b).ix == psgr_arc_sym(ap).ix, "dprev(sym()):sym(oprev()) bug");
        demand(psgr_arc_dnext(gr, b).ix == psgr_arc_sym(an).ix, "dnext(sym()):sym(onext()) bug");
      }
  }
        
void test__psgr_arc_lnext__psgr_arc_lprev(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        demand(psgr_edge_mark(gr, e) != DELETED, "edge is marked as deleted");
        psgr_arc_t a = psgr_orient_edge(e, 0);
        psgr_arc_t b = psgr_arc_sym(a);
        
        psgr_arc_t aon = psgr_arc_onext(gr, a);
        psgr_arc_t aop = psgr_arc_oprev(gr, a);
        psgr_arc_t aln = psgr_arc_lnext(gr, a);
        psgr_arc_t alp = psgr_arc_lprev(gr, a);
        
        psgr_arc_t bon = psgr_arc_onext(gr, b);
        psgr_arc_t bop = psgr_arc_oprev(gr, b);
        psgr_arc_t bln = psgr_arc_lnext(gr, b);
        psgr_arc_t blp = psgr_arc_lprev(gr, b);
        
        demand(aln.ix == bop.ix, "oprev/lnext bug (1)");
        demand(bln.ix == aop.ix, "oprev/lnext bug (2)");
        
        demand(alp.ix == psgr_arc_sym(aon).ix, "onext/lprev bug (1)");
        demand(blp.ix == psgr_arc_sym(bon).ix, "onext/lprev bug (2)");
      }
  }
        
void test__psgr_arc_split__psgr_arc_unsafe_split(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        demand(psgr_edge_mark(gr, e) != DELETED, "edge is marked as deleted");
        psgr_arc_t a = psgr_orient_edge(e, 0);
        psgr_arc_t b = psgr_arc_sym(a);
        
        auto void do_test_split(psgr_arc_t at, psgr_dir_bit_t st);
          /* Tests {} on {at}.  Expects the edge to be {e}, the edge data record to be {ed},
            and the dir bit to be {st}. */

        do_test_split(a, 0);
        do_test_split(b, 1);
        
        void do_test_split(psgr_arc_t at, psgr_dir_bit_t st)
          { 
            psgr_edge_t eat;
            psgr_dir_bit_t sat;
            psgr_arc_split(at, &eat, &sat);  
            demand(eat.ix == e.ix, "split bug: edge");
            demand(sat == st, "split bug: dir bit");

            psgr_edge_t eau;
            psgr_dir_bit_t sau;
            psgr_edge_data_t *edau;
            psgr_arc_unsafe_split(gr, at, &eau, &sau, &edau);  
            demand(eau.ix == e.ix, "unsafe_split bug: edge");
            demand(sau == st, "unsafe_split bug: dir bit");
            demand(edau == &(gr->edata[e.ix]), "unsafe_split bug: data record");
          }
       }
  }
        
void test__psgr_arc_delta__psgr_arc_weight__psgr_arc_set_delta__psgr_arc_set_weight(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        demand(psgr_edge_mark(gr, e) != DELETED, "edge is marked as deleted");
        psgr_arc_t a = psgr_orient_edge(e, 0);
        psgr_arc_t b = psgr_arc_sym(a);
        
        double d1 = 10*dabrandom(-1.0,+1.0);
        psgr_arc_set_delta(gr, a, d1);
        demand(psgr_arc_delta(gr, a) == d1, "delta bug (1)");
        demand(psgr_arc_delta(gr, b) == -d1, "delta/sym bug (1)");
        
        double d2 = 10*dabrandom(-1.0,+1.0);
        psgr_arc_set_delta(gr, b, d2);
        demand(psgr_arc_delta(gr, b) == d2, "delta bug (2)");
        demand(psgr_arc_delta(gr, a) == -d2, "delta/sym bug (2)");
        
        double w1 = dabrandom(0.0001, 9.9999);
        psgr_arc_set_weight(gr, a, w1);
        demand(psgr_arc_weight(gr, a) == w1, "weight bug (1)");
        demand(psgr_arc_weight(gr, b) == w1, "weight/sym bug (1)");
        
        double w2 = dabrandom(0.0001, 9.9999);
        psgr_arc_set_weight(gr, b, w2);
        demand(psgr_arc_weight(gr, b) == w2, "weight bug (2)");
        demand(psgr_arc_weight(gr, a) == w2, "weight/sym bug (2)");
      }
  }
        
void test__psgr_vertex_out_arc__psgr_vertex_outdegree(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NV = psgr_num_vertices(gr);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t kv = 0; kv < NV; kv++)
      { psgr_vertex_t v = (psgr_vertex_t){ kv };
        demand(psgr_vertex_mark(gr, v) != DELETED, "vertex is marked as deleted");
        
        psgr_arc_t a = psgr_vertex_out_arc(gr, v);
        if (a.ix != ANONE.ix) { demand(psgr_arc_org(gr, a).ix == v.ix, "aut_arc/org bug"); }
        
        /* Check outdegree the dumb way: */
        uint32_t deg = psgr_vertex_outdegree(gr, v);
        uint32_t deg_exp = 0;
        if (a.ix != ANONE.ix)
          { demand(psgr_edge_mark(gr, psgr_arc_edge(a)) != DELETED, "edge is marked as deleted");
            for (uint32_t je = 0; je < NE; je++)
              { psgr_edge_t ej = (psgr_edge_t){ je };
                demand(psgr_edge_mark(gr, ej) != DELETED, "edge is marked as deleted");
                psgr_arc_t aj = psgr_orient_edge(ej, 0);
                psgr_arc_t bj = psgr_arc_sym(aj);
                if (psgr_arc_org(gr, aj).ix == v.ix) { deg_exp++; }
                if (psgr_arc_org(gr, bj).ix == v.ix) { deg_exp++; }
              }
          }
        demand(deg == deg_exp, "outdegree bug");
     }
  }
        
void test__psgr_arc_path__psgr_arc_starting_dir(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        psgr_arc_t a = psgr_orient_edge(e, 0);
        psgr_arc_t b = psgr_arc_sym(a);
        
        psgr_path_t aP = psgr_arc_path(gr, a);
        psgr_path_t bP = psgr_arc_path(gr, b);
        demand(aP.n == bP.n, "path/sym bug (n)");
        demand(aP.reverse == ! bP.reverse, "path/sym bug (reverse)");
        demand(aP.v == bP.v, "path/sym bug (v)");

        psgr_vertex_t org = psgr_arc_org(gr, a);
        uint32_t dorg = psgr_vertex_outdegree(gr, org);
        if (dorg >= 3)
          { /* Check ordering of edges out of {org(a)}: */
            psgr_arc_t aop = psgr_arc_oprev(gr, a);
            psgr_arc_t aon = psgr_arc_onext(gr, a);
            assert((aop.ix != a.ix) && (a.ix != aon.ix) && (aon.ix != aop.ix));
            r2_t uj = psgr_arc_starting_dir(gr, aop);
            r2_t ui = psgr_arc_starting_dir(gr, a); 
            r2_t uk = psgr_arc_starting_dir(gr, aon);
            demand(fabs(r2_norm(&ui) - 1.0) < 1.0e-7, "arc_dir_vector not normalized");
            demand(r2_cyclic_order(&uj, &ui, &uk) == +1, "arcs out of vertex are not sorted");
          }
      }
  }
        
void test__psgr_find_nearest_vertex(psgr_t *gr)
  { 
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NV = psgr_num_vertices(gr);
    i2_t sz = psgr_image_size(gr);
    uint32_t NT = 100; /* Number of trials */
    for (uint32_t kt = 0; kt < NT; kt++)
      { r2_t p = (r2_t){{ drandom()*(sz.c[0]-1), drandom()*(sz.c[1]-1) }};
        psgr_vertex_t v_cmp = psgr_find_nearest_vertex(gr, &p);
        demand(v_cmp.ix != VNONE.ix, "nearest_vertex returned null vertex");
        /* Look for the nearest vertex by hand: */
        psgr_vertex_t v_exp = VNONE;
        double dmin = +INF;
        for (uint32_t jv = 0; jv < NV; jv++)
          { psgr_vertex_t v = (psgr_vertex_t){ jv };
            i2_t iq = psgr_vertex_pixel(gr, v);
            r2_t q = (r2_t){{ iq.c[0], iq.c[1] }};
            double dpq = r2_dist(&p, &q);
            if (dpq < dmin) { v_exp = v; dmin = dpq; }
          }
        assert(isfinite(dmin));
        demand(v_cmp.ix == v_exp.ix, "psgr_find_nearest_vertex error");
      }
  }
  
void test__psgr_arc_left_face_curl__psgr_arc_left_face_barycenter(psgr_t *gr)
  { fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NE = psgr_num_edges(gr);
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        psgr_arc_t a = psgr_orient_edge(e, 0);
        double curl_cmp_i = psgr_arc_left_face_curl(gr, a);
        r2_t bar_cmp_i = psgr_arc_left_face_barycenter(gr, a);
        demand(isfinite(curl_cmp_i), "non-finite curl"); 
        psgr_arc_t aj = psgr_arc_lnext(gr, a);
        while (aj.ix != a.ix) 
          { double curl_cmp_j = psgr_arc_left_face_curl(gr, aj);
            if (fabs(curl_cmp_i - curl_cmp_j) > 1.0e-13)
              { fprintf(stderr, "  curl = %24.16e", curl_cmp_i);
                fprintf(stderr, "  a = "); psgr_arc_print(stderr, gr, a); fputc('\n',stderr);
                fprintf(stderr, "  curl = %24.16e", curl_cmp_j);
                fprintf(stderr, "  aj = "); psgr_arc_print(stderr, gr, aj); fputc('\n',stderr);
                fprintf(stderr, "  difference = %24.16e\n", curl_cmp_i - curl_cmp_j);
                demand(FALSE, "curl/lnext error");
              }
            r2_t bar_cmp_j = psgr_arc_left_face_barycenter(gr, aj);
            demand(r2_dist(&bar_cmp_i, &bar_cmp_j) < 1.0e-7, "barycenter/lnex error");
            aj = psgr_arc_lnext(gr, aj);
          }
      }
    fprintf(stderr, "!! psgr_arc_left_face_curl VALUE NOT CHECKED\n");
    fprintf(stderr, "!! psgr_arc_left_face_barycenter VALUE NOT CHECKED\n");
  }

void test__psgr_copy__psgr_equal__psgr_free(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    uint32_t NV_gr = psgr_num_vertices(gr);
    uint32_t NE_gr = psgr_num_edges(gr);
    
    /* Assumes that no vertices or edges are marked {DELETED}: */
    psgr_t *hr = psgr_copy(gr, NULL, NULL);
    psgr_check_consistency(hr);
    demand(psgr_equal(gr, hr), "copy bug");
    uint32_t NV_hr = psgr_num_vertices(hr); assert(NV_hr == NV_gr);
    uint32_t NE_hr = psgr_num_edges(hr); assert(NE_hr == NE_gr);
    
    /* Mark some vertices of {hr} deleted: */
    uint32_t NV_keep = NV_gr;
    uint32_t NE_keep = NE_gr;
    for (uint32_t kdv = 0; kdv < 5; kdv++)
      { psgr_vertex_t v = (psgr_vertex_t){ 7 + 3*kdv };
        assert(v.ix < NV_gr);
        assert(psgr_vertex_mark(hr, v) != DELETED);
        psgr_arc_t a = psgr_vertex_out_arc(hr, v);
        if (a.ix != ANONE.ix)
          { /* Mark incident edges as {DELETED}: */
            psgr_arc_t b = a;
            do 
              { psgr_edge_t eb = psgr_arc_edge(b);
                if (psgr_edge_mark(hr, eb) != DELETED)
                  { psgr_edge_set_mark(hr, eb, DELETED);
                    assert(NE_keep > 0);
                    NE_keep--;
                  }
                b = psgr_arc_onext(hr, b);
              }
            while (b.ix != a.ix);
          }
            
        psgr_vertex_set_mark(hr, v, DELETED);
        assert(NV_keep > 0);
        NV_keep--;
      }
    psgr_check_consistency(hr);
      
    /* Make copy {fr} of {hr}: */
    int32_t *vfr_from_vhr = talloc(NV_hr, int32_t);
    int32_t *efr_from_ehr = talloc(NE_hr, int32_t);
    psgr_t *fr = psgr_copy(hr, vfr_from_vhr, efr_from_ehr);
    psgr_check_consistency(fr);
    uint32_t NV_fr = psgr_num_vertices(fr);
    uint32_t NE_fr = psgr_num_edges(fr);
    demand(NV_fr == NV_keep, "copy NV bug (1)");
    demand(NE_fr == NE_keep, "copy NE bug (1)");
    
    /* Compare vertices of {hr,fr}: */
    for (uint32_t kv_hr = 0; kv_hr < NV_hr; kv_hr++)
      { psgr_vertex_t v_hr = (psgr_vertex_t){ kv_hr };
        int32_t kv_fr = vfr_from_vhr[kv_hr];
        if (kv_fr == VNONE.ix)
          { /* Vertex {v_hr} must have been deleted: */
            demand(psgr_vertex_mark(hr, v_hr) == DELETED, "copy bug: removed undeleted vertex");
          }
        else
          { demand((kv_fr >= 0) && (kv_fr < NV_fr), "copy bug: bad entry in vertex map");
            psgr_vertex_t v_fr = (psgr_vertex_t){ (uint32_t)kv_fr };
            demand(psgr_vertex_mark(fr, v_fr) == psgr_vertex_mark(hr, v_hr), "copy bug: mark changed");
            /* Compare vertices: */
            i2_t ip_hr = psgr_vertex_pixel(hr, v_hr);
            i2_t ip_fr = psgr_vertex_pixel(fr, v_fr);
            demand((ip_hr.c[0] == ip_fr.c[0]) && (ip_hr.c[1] == ip_fr.c[1]), "copy bug: diff pixels");
            /* The {.aout} field should have been checked by {psgr_check_consistency}. */
          }
      }
            
    /* Compare edges of {hr,fr}: */
    for (uint32_t ke_hr = 0; ke_hr < NE_hr; ke_hr++)
      { psgr_edge_t e_hr = (psgr_edge_t){ ke_hr };
        int32_t ke_fr = efr_from_ehr[ke_hr];
        if (ke_fr == ENONE.ix)
          { /* Edge {e_hr} must have been deleted: */
            demand(psgr_edge_mark(hr, e_hr) == DELETED, "copy bug: removed undeleted edge");
          }
        else
          { demand((ke_fr >= 0) && (ke_fr < NE_fr), "copy bug: bad entry in edge map");
            psgr_edge_t e_fr = (psgr_edge_t){ (uint32_t)ke_fr };
            demand(psgr_edge_mark(fr, e_fr) == psgr_edge_mark(hr, e_hr), "copy bug: edge mark changed");
            psgr_arc_t a_hr = psgr_orient_edge(e_hr, 0);
            psgr_vertex_t org_hr = psgr_arc_org(hr, a_hr);
            assert((org_hr.ix >= 0) && (org_hr.ix < NV_hr)); /* Input reqmt. */
            psgr_vertex_t dst_hr = psgr_arc_dst(hr, a_hr);
            assert((dst_hr.ix >= 0) && (dst_hr.ix < NV_hr)); /* Input reqmt. */

            psgr_arc_t a_fr = psgr_orient_edge(e_fr, 0);
            psgr_vertex_t org_fr = psgr_arc_org(fr, a_fr); 
            demand((org_fr.ix >= 0) && (org_fr.ix < NV_fr), "copy bug: invalid edge org");
            demand(org_fr.ix == vfr_from_vhr[org_hr.ix], "copy bug: wrong org");
            psgr_vertex_t dst_fr = psgr_arc_dst(fr, a_fr); 
            demand((dst_fr.ix >= 0) && (dst_fr.ix < NV_fr), "copy bug: invalid edge dst");
            demand(dst_fr.ix == vfr_from_vhr[dst_hr.ix], "copy bug: wrong dst");
            
            /* The {.enext} fields should have been checked by {psgr_check_consistency}. */
            demand(psgr_arc_delta(fr, a_fr) == psgr_arc_delta(hr, a_hr), "copy bug: delta");
            demand(psgr_arc_weight(fr, a_fr) == psgr_arc_weight(hr, a_hr), "copy bug: weight");
            /* Compare paths: */
            psgr_path_t P_hr = psgr_arc_path(hr, a_hr);
            psgr_path_t P_fr = psgr_arc_path(fr, a_fr);
            if (P_hr.v != NULL) { demand(P_fr.v != P_hr.v, "copy bug: shared path nodes"); }
            demand(psgr_path_equal(P_fr, P_hr), "copy bug: paths differ");
          }
      }
            
    psgr_free(hr);
    psgr_free(fr);
    free(vfr_from_vhr);
    free(efr_from_ehr);
  }

void test__psgr_add_vertex__psgr_add_edge__psgr_edge_remove__psgr_vertex_remove(psgr_t *gr)
  { 
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    
    /* We use {psgr_add_vertex} and {psgr_add_edge} to create a 
      graph that is a circle of vertices connected to a central 
      vertex.  Then we use {psgr_edge_remove} and {psgr_vertex_remove}
      to dismantle it. */
    
    /* Assumes that no vertices or edges are marked {DELETED}: */
    uint32_t NX = 640;
    uint32_t NY = 480;
    psgr_t *hr = psgr_new(100, 300, NX,NY);
    psgr_check_consistency(hr);
     
    auto psgr_vertex_t add_and_check_vertex(i2_t *ip);
      /* Adds a vertex at {ip}, checks consistency. */
     
    auto psgr_arc_t add_and_check_edge(psgr_vertex_t vo, psgr_vertex_t vd);
      /* Adds new arc from {vo} to {vd}, checks consistency. */
   
    /* Create a vertex of degree zero at the center: */
    i2_t ictr = (i2_t){{ (int32_t)NX/2, (int32_t)NY/2 }};  /* Center pixel. */
    double rad = 0.40*fmin(NX, NY);      /* Radius of circle (px). */
    psgr_vertex_t vctr = add_and_check_vertex(&ictr);
    
    /* Add rim vertices and radial edges: */
    uint32_t NR = 5; /* Vertices in circle. */
    psgr_vertex_t vrim[NR]; /* Rim vertices, clockwise. */
    psgr_arc_t arad[NR]; /* Arcs from {vctr} to {vrim[0..NR-1]}. */
    for (int32_t kr = 0; kr <= NR; kr++)
      { double ang = (2*M_PI*kr)/NR;
        int32_t x = ictr.c[0] + (int32_t)(floor(rad*sin(ang) + 0.5));
        int32_t y = ictr.c[1] + (int32_t)(floor(rad*cos(ang) + 0.5));
        i2_t ipk = (i2_t){{ x, y }};
        vrim[kr] = add_and_check_vertex(&ipk);
        arad[kr] = add_and_check_edge(vctr, vrim[kr]);
      }
        
    /* Add rim edges: */
    psgr_arc_t arim[NR]; /* Arcs from {vrim[kr]} to {vrim[kr+1]} (mod NR). */
    for (uint32_t kr = 0; kr <= NR; kr++)
      { psgr_vertex_t vo = vrim[kr];
        psgr_vertex_t vd = vrim[(kr + 1) % NR];
        arim[kr] = add_and_check_edge(vo, vd);
      }
      
    /* Check onext links: */
    for (uint32_t kro = 0; kro <= NR; kro++)
      { uint32_t krm = (kro + NR - 1) % NR;
        uint32_t krp = (kro + 1) % NR;
        demand(psgr_arc_onext(hr, arad[kro]).ix == arad[krp].ix, "add_edge bug: onext (arad)");
        demand(psgr_arc_dnext(hr, arad[kro]).ix == psgr_arc_sym(arim[krp]).ix, "add_edge bug: dnext (arad)");
        demand(psgr_arc_dprev(hr, arad[kro]).ix == arim[krm].ix, "add_edge bug: dprev (arad)");
        demand(psgr_arc_onext(hr, arim[kro]).ix == psgr_arc_sym(arim[krm]).ix, "add_edge bug: onext (arim)");
      }    
      
    auto void delete_edge_and_check(psgr_arc_t a);
      /* Deletes the edge of arc{a} and runs consistency checks. */
    
    auto void delete_vertex_and_check(psgr_vertex_t v);
      /* Deletes the vetex {v} (which must have degree zero)
        and runs consistency checks. */
    
    /* Delete the edges: */
    for (int32_t kr = 0; kr <= NR; kr++)
      { delete_edge_and_check(arad[kr]);
        delete_edge_and_check(arim[kr]);
      }
      
    /* Delete the vertices: */
    for (int32_t kr = 0; kr <= NR; kr++)
      { delete_vertex_and_check(vrim[kr]); }
    delete_vertex_and_check(vctr);
    
    demand(psgr_num_vertices(hr) == 0, "remove_vertex bug: not empty");
    demand(psgr_num_edges(hr) == 0, "remove_edge bug: not empty");
    
    psgr_free(hr);
    return;
        
    psgr_vertex_t add_and_check_vertex(i2_t *ip)
      { uint32_t NV_cur = psgr_num_vertices(hr);
        uint32_t NE_cur = psgr_num_edges(hr);
        psgr_vertex_t v = psgr_add_vertex(hr, ip);
        psgr_check_consistency(hr);
        demand(psgr_num_vertices(hr) == NV_cur + 1,  "add_vertex bug: vertex count");
        demand(psgr_num_edges(hr) == NE_cur,  "add_vertex bug: edge count");
        demand(v.ix == NV_cur, "add_vertex bug: vertex index");
        demand(psgr_vertex_out_arc(hr, v).ix == ANONE.ix, "add_vertex bug: out arc not null");
        i2_t ip_check = psgr_vertex_pixel(hr, v);
        demand(i2_eq(ip, &ip_check), "add_vertex bug: pixel");
        psgr_vertex_t v_check = psgr_vertex_from_pixel(hr, ip);
        demand(v_check.ix == v.ix, "add_vertex bug: from_pixel");
        demand(psgr_vertex_mark(hr, v) == CLEAR, "add_vertex bug: mark");
        return v;
      }
    
    psgr_arc_t add_and_check_edge(psgr_vertex_t vo, psgr_vertex_t vd)
      { uint32_t NV_cur = psgr_num_vertices(hr);
        uint32_t NE_cur = psgr_num_edges(hr);
        i2_t ipo = psgr_vertex_pixel(hr, vo);
        i2_t ipd = psgr_vertex_pixel(hr, vd);
        psgr_path_t P = psgr_path_throw(&ipo, 3, &ipd); /* Hope does not cross other paths. */
        double d = dabrandom(-5.0, +5.0);
        double w = dabrandom(0.001, 0.999);
        psgr_arc_t a = psgr_add_edge(hr, vo, vd, d, w, P);
        psgr_check_consistency(hr);
        demand(psgr_num_vertices(hr) == NV_cur,  "add_edge bug: vertex count");
        demand(psgr_num_edges(hr) == NE_cur + 1,  "add_edge bug: edge count");
        psgr_edge_t e = psgr_arc_edge(a);
        demand(psgr_arc_dir_bit(a) == 0, "add_edge bug: dir bit");
        demand(e.ix == NE_cur, "add_edge bug: edge index");
        demand(psgr_arc_org(hr, a).ix == vo.ix, "add_edge bug: org");
        demand(psgr_arc_dst(hr, a).ix == vd.ix, "add_edge bug: dst");
        demand(psgr_arc_delta(hr, a) == d, "add_edge bug: delta");
        demand(psgr_arc_weight(hr, a) == w, "add_edge bug: weight");
        demand(psgr_path_equal(psgr_arc_path(hr, a), P), "add_edge bug: path");
        return a;
      }
      
    auto void delete_edge_and_check(psgr_arc_t a)
      { uint32_t NV_cur = psgr_num_vertices(hr);
        uint32_t NE_cur = psgr_num_edges(hr);
        assert(a.ix != ANONE.ix);
        psgr_edge_t e = psgr_arc_edge(a);
        assert(psgr_edge_mark(hr, e) == CLEAR);
        psgr_arc_t a_onext = psgr_arc_onext(hr, a);
        psgr_arc_t a_oprev = psgr_arc_oprev(hr, a);
        psgr_arc_t a_dnext = psgr_arc_dnext(hr, a);
        psgr_arc_t a_dprev = psgr_arc_dprev(hr, a);
        psgr_vertex_t org = psgr_arc_org(hr, a);
        psgr_vertex_t dst = psgr_arc_dst(hr, a);
        psgr_edge_remove(hr, e, 2);
        psgr_check_consistency(hr);
        demand(psgr_num_vertices(hr) == NV_cur,  "edge_remove bug: vertex count");
        demand(psgr_num_edges(hr) == NE_cur,  "edge_remove bug: edge count");
        demand(psgr_edge_mark(hr, e) == DELETED, "edge_remove bug: mark");
        
        psgr_arc_t org_aout = psgr_vertex_out_arc(hr, org);
        if (a_onext.ix == a.ix)
          { /* Only edge out of {org}: */
            assert(a_oprev.ix == a.ix);
            demand(org_aout.ix == ANONE.ix, "edge_remove bug: last of org");
          }
        else
          { demand(psgr_arc_onext(hr, a_oprev).ix == a_onext.ix, "edge_remove bug: onext ring");
            demand(org_aout.ix != a.ix, "edge_remove bug: org_aout is deleted");
            demand(org_aout.ix != ANONE.ix, "edge_remove bug: org_aout null");
          }
         
        psgr_arc_t dst_aout = psgr_vertex_out_arc(hr, dst);
        if (a_dnext.ix == a.ix)
          { /* Only edge out of {dst}: */
            assert(a_dprev.ix == a.ix);
            demand(dst_aout.ix == ANONE.ix, "edge_remove bug: last of dst");
          }
        else
          { demand(psgr_arc_dnext(hr, a_dprev).ix == a_dnext.ix, "edge_remove bug: dnext ring");
            demand((dst_aout.ix != a.ix) && (dst_aout.ix != psgr_arc_sym(a).ix), "edge_remove bug: dst_aout is deleted");
            demand(dst_aout.ix != ANONE.ix, "edge_remove bug: dst_aout null");
          }
      }
    
    auto void delete_vertex_and_check(psgr_vertex_t v)
      { uint32_t NV_cur = psgr_num_vertices(hr);
        uint32_t NE_cur = psgr_num_edges(hr);
        assert(v.ix != VNONE.ix);
        assert(psgr_vertex_mark(hr, v) == CLEAR);
        psgr_arc_t v_aout = psgr_vertex_out_arc(hr, v);
        assert(v_aout.ix == ANONE.ix); /* Checked earlier. */
        psgr_vertex_remove(hr, v, 2);
        psgr_check_consistency(hr);
        demand(psgr_num_vertices(hr) == NV_cur,  "vertex_remove bug: vertex count");
        demand(psgr_num_edges(hr) == NE_cur,  "vertex_remove bug: edge count");
        demand(psgr_vertex_mark(hr, v) == DELETED, "vertex_remove bug: mark");
      }
  }

void test__psgr_write_named__psgr_write_file__psgr_read_named__psgr_read_file(psgr_t *gr)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    char *fname = "out/graph_g.txt";
    psgr_write_named(fname, gr, TRUE);
    psgr_t *hr = psgr_read_named(fname, TRUE);
    demand(psgr_equal(gr, hr), "read not equal to write");
  }

options_t *tb_parse_options(int32_t argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* PARSE KEYWORD ARGUMENTS: */
    
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
