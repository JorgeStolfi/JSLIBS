#define PROG_NAME "test_gr_basic"
#define PROG_DESC "checks the routines in {pst_gr.h}, {pst_gr_plot.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-03-14 19:32:52 by stolfi */
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
  "  The program checks the functions {pst_gr}."

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

#include <pst_gr.h>
#include <pst_gr_plot.h>
#include <pst_gr_read.h>
#include <pst_gr_write.h>
#include <pst_gr_path.h>
#include <pst_gr_test.h>

#define NONE pst_gr_NONE 

typedef struct options_t
  { int32_t DUMMY;  /* Placeholder. */
  } options_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

void write_graph(char *fname, pst_gr_t *gr);
  /* Writes {gr} to file "{fname}". */
  
pst_gr_t* read_graph(char *fname);
  /* Reads back a graph from file "{fname}". */

void test_pst_gr_arc_moves(pst_gr_t *gr);
void test_pst_gr_path_funcs(pst_gr_t *gr);
void test_pst_gr_find_nearest_vertex(pst_gr_t *gr);
void test_pst_gr_compute_left_face_properties(pst_gr_t *gr);

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    /* Create a test graph: */
    fprintf(stderr, "--- {pst_gr_test_make_graph: pst_gr_new,pst_gr_add_vertex,pst_gr_add_edge} ---\n");
    uint32_t nf = 2; /* Vertex layers in frame. */
    uint32_t nh = 3; /* Verted rows and cols in hole. */
    pst_gr_t *gr = pst_gr_test_make_graph(nf, nh);

    fprintf(stderr, "--- {pst_gr_check_consistency} ---\n");
    pst_gr_check_consistency(gr);

    fprintf(stderr, "--- {pst_gr_plot_named} ---\n");
    pst_gr_plot_named("out/plot-lab0.eps", gr, 0, 1.0);
    pst_gr_plot_named("out/plot-lab1.eps", gr, 8, 1.0);
    
    test_pst_gr_arc_moves(gr);
    test_pst_gr_path_funcs(gr);
    test_pst_gr_find_nearest_vertex(gr);
    test_pst_gr_compute_left_face_properties(gr);

    char *fname = "out/graph_g.txt";
    write_graph(fname, gr);
    pst_gr_t *h = read_graph(fname);
    pst_gr_equal(gr, h);

    /* Cleanup: */
    pst_gr_free(gr);
    free(o); o = NULL;
    return 0;
  }
  
void test_pst_gr_arc_moves(pst_gr_t *gr)
  {
    fprintf(stderr, "--- {pst_gr_arc_{orient_edge,edge,dir_bit,sym,org,dst,oprev,onext,lnext,weight,delta}} ---\n");
    for (pst_gr_edge_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_arc_t ai = pst_gr_orient_edge(ei, 0);
        demand(pst_gr_arc_edge(ai) == ei, "pst_gr_arc_edge bug");
        demand(pst_gr_arc_dir_bit(ai) == 0, "pst_gr_arc_dir_bit bug");
        pst_gr_arc_t bi = pst_gr_arc_sym(ai);
        demand(ai != bi, "pst_gr_arc_sym has a fixed arg");
        demand(pst_gr_arc_edge(bi) == ei, "pst_gr_arc_sym edge bug");
        demand(pst_gr_arc_dir_bit(bi) == 1, "pst_gr_arc_sym dir bit bug");
        demand(pst_gr_arc_sym(bi) == ai, "pst_gr_arc_sym not involution");
        
        pst_gr_vertex_t org = pst_gr_arc_org(gr, ai); demand(org < gr->NV, "invalid org");
        pst_gr_vertex_t dst = pst_gr_arc_dst(gr, ai); demand(dst < gr->NV, "invalid dst");
        demand(org != dst, "graph has a loop");
        demand(pst_gr_arc_org(gr, bi) == dst, "org(sym) bug");
        demand(pst_gr_arc_dst(gr, bi) == org, "dst(sym) bug");
        
        pst_gr_arc_t ani = pst_gr_arc_onext(gr, ai);
        pst_gr_arc_t api = pst_gr_arc_oprev(gr, ai);
        demand(pst_gr_arc_oprev(gr, ani) == ai, "oprev(onext()) bug");
        demand(pst_gr_arc_onext(gr, api) == ai, "onext(oprev()) bug");
        
        pst_gr_arc_t bni = pst_gr_arc_onext(gr, bi);
        pst_gr_arc_t bpi = pst_gr_arc_oprev(gr, bi);
        demand(pst_gr_arc_oprev(gr, bni) == bi, "oprev(onext()) bug");
        demand(pst_gr_arc_onext(gr, bpi) == bi, "onext(oprev()) bug");
        
        demand(pst_gr_arc_lnext(gr, ai) == bpi, "onext/lnext bug (1)");
        demand(pst_gr_arc_lnext(gr, bi) == api, "onext/lnext bug (2)");
        
        demand(pst_gr_arc_weight(gr, bi) == pst_gr_arc_weight(gr, ai), "weight/sym bug");
        demand(pst_gr_arc_delta(gr, bi) == - pst_gr_arc_delta(gr, ai), "delta/sym bug");
        
        uint32_t dorg = pst_gr_outdegree(gr, org);
        pst_gr_arc_t aj = pst_gr_arc_onext(gr, ai);
        for (int32_t k = 1; k < dorg; k++)
          { demand(aj != ai, "outdegree too big");
            aj = pst_gr_arc_onext(gr, aj);
          }
        demand(aj == ai, "outdegree too low");
       }
  }
    
void test_pst_gr_path_funcs(pst_gr_t *gr)
  { fprintf(stderr, "--- {pst_gr_arc_path,pst_gr_arc_start_dir} ---\n");
    for (uint32_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_arc_t ai = pst_gr_orient_edge(ei, 0);
        pst_gr_arc_t bi = pst_gr_arc_sym(ai);
        pst_gr_path_t aP = pst_gr_arc_path(gr, ai);
        pst_gr_path_t bP = pst_gr_arc_path(gr, bi);
        demand(aP.n == bP.n, "path/sym bug (n)");
        demand(aP.reverse == ! bP.reverse, "path/sym bug (reverse)");
        demand(aP.v == bP.v, "path/sym bug (v)");
        pst_gr_vertex_t org = pst_gr_arc_org(gr, ai);
        uint32_t dorg = pst_gr_outdegree(gr, org);
        if (dorg >= 3)
          { /* Check ordering of edges out of {org(ai)}: */
            pst_gr_arc_t aj = pst_gr_arc_oprev(gr, ai);
            pst_gr_arc_t ak = pst_gr_arc_onext(gr, ai);
            assert((aj != ai) && (ai != ak) && (ak != aj));
            r2_t uj = pst_gr_arc_start_dir(gr, aj);
            r2_t ui = pst_gr_arc_start_dir(gr, ai); 
            r2_t uk = pst_gr_arc_start_dir(gr, ak);
            demand(fabs(r2_norm(&ui) - 1.0) < 1.0e-7, "arc_dir_vector not normalized");
            demand(r2_cyclic_order(&uj, &ui, &uk) == +1, "arcs out of vertex are not sorted");
          }
      }
  }
        
void test_pst_gr_find_nearest_vertex(pst_gr_t *gr)
  { 
    fprintf(stderr, "--- {pst_gr_find_nearest_vertex} ---\n");
    for (int32_t vi = 0; vi < gr->NV; vi++)
      { pst_gr_vertex_data_t *vdi = &(gr->vdata[vi]);
        r2_t pi = vdi->coords;
        r2_t di = (r2_t){{ 28*cos(vi), 32*sin(vi) }};
        r2_add(&pi, &di, &pi);
        pst_gr_vertex_t vj_cmp = pst_gr_find_nearest_vertex(gr, &pi);
        pst_gr_vertex_t vj_exp = NONE;
        double dmin = +INF;
        for (pst_gr_vertex_t vk = 0; vk < gr->NV; vk++)
          { pst_gr_vertex_data_t *vdk = &(gr->vdata[vk]);
            r2_t *pk = &(vdk->coords);
            double dik = r2_dist(&pi, pk);
            if (dik < dmin) { vj_exp = vk; dmin = dik; }
          }
        assert(isfinite(dmin));
        demand(vj_cmp == vj_exp, "pst_gr_find_nearest_vertex error");
      }
  }
  
void test_pst_gr_compute_left_face_properties(pst_gr_t *gr)
  { fprintf(stderr, "--- {pst_gr_left_face_curl,pst_gr_left_face_barycenter} ---\n");
    for (uint32_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_arc_t ai = pst_gr_orient_edge(ei, 0);
        double curl_cmp_i = pst_gr_compute_left_face_curl(gr, ai);
        r2_t bar_cmp_i = pst_gr_left_face_barycenter(gr, ai);
        demand(isfinite(curl_cmp_i), "non-finite curl"); 
        pst_gr_arc_t aj = pst_gr_arc_lnext(gr, ai);
        while (aj != ai) 
          { double curl_cmp_j = pst_gr_compute_left_face_curl(gr, aj);
            if (fabs(curl_cmp_i - curl_cmp_j) > 1.0e-13)
              { fprintf(stderr, "  curl = %24.16e", curl_cmp_i);
                fprintf(stderr, "  ai = "); pst_gr_arc_print(stderr, gr, ai); fputc('\n',stderr);
                fprintf(stderr, "  curl = %24.16e", curl_cmp_j);
                fprintf(stderr, "  aj = "); pst_gr_arc_print(stderr, gr, aj); fputc('\n',stderr);
                fprintf(stderr, "  difference = %24.16e\n", curl_cmp_i - curl_cmp_j);
                demand(FALSE, "curl/lnext error");
              }
            r2_t bar_cmp_j = pst_gr_left_face_barycenter(gr, aj);
            demand(r2_dist(&bar_cmp_i, &bar_cmp_j) < 1.0e-7, "barycenter/lnex error");
            aj = pst_gr_arc_lnext(gr, aj);
          }
        fprintf(stderr, "pst_gr_compute_left_face_curl VALUE NOT CHECKED\n");
        fprintf(stderr, "pst_gr_compute_left_face_barycenter VALUE NOT CHECKED\n");
      }
  }

void write_graph(char *fname, pst_gr_t *gr)
  { 
    pst_gr_write_named(fname, gr, TRUE);
  }
    
pst_gr_t* read_graph(char *fname)
  { 
    pst_gr_t *gr = pst_gr_read_named(fname, TRUE);
    return gr;
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
