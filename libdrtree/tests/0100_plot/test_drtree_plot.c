#define PROG_NAME "test_drtree_plot"
#define PROG_VERS "2.0"

/* Draws the evolution of asexual traits. */

/* Last edited on 2023-06-23 11:21:44 by stolfi */

#define test_drtree_plot_C_COPYRIGHT \
  "???"

#define PROG_HELP \
  " " PROG_NAME " \\\n" \
  drtree_test_options_HELP
 
#define PROG_INFO \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Tests the {drtree_compact.h}, {drtree_planar.h}, and {drtree_plot.h} functions. The test uses a randomly" \
  " generated evolution history.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The main output files are \"{outPrefix}-{type}-{tag}-plot.eps\" where {type} is \"comp\" or \"plan\" and {tag} is \"A\", \"B\", etc..\n" \
  "\n" \
  "OPTIONS\n" \
  drtree_test_options_INFO "\n" \
  "" \
  "???"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <cmp.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <vec.h> 
#include <epswr.h> 

#include <drtree.h>
#include <drtree_test.h>
#include <drtree_compact.h>
#include <drtree_planar.h>
#include <drtree_plot.h>

typedef drtree_test_options_t tdp_options_t;

tdp_options_t *tdp_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void tdp_do_test
  ( tdp_options_t *o,
    bool_t planar,
    char *tag,
    int32_t t0,
    int32_t t1,
    int32_t ni,
    drtree_node_t dt[]
  );
  /* Creates an EPS file "{o.outPrefix}-{tag}.eps" with the compact diagram
    of the nodes {dt[0..ni-1}} in the time interval {t0..t1}.
    
    Reference lines are drawn at times {o->tStart,o->tStop,o->tRef} if
    they are in the grids's column range. */

void tdp_plot_named
  ( char *outPrefix, 
    char *tag,
    int32_t tMin,
    int32_t tMax,
    int32_t ncols,
    int32_t nrows,
    int32_t ni, 
    drtree_node_t dt[],
    int32_t rdr[]
  );       
  /* Creates an EPS file "{outPrefix}-{tag}.eps" and draws the diagram
    of the nodes {dt[0..ni-1}} on rows {rdr[0..ni-1]} of the plot grid
    spanning times {tMin..tMax}.
    
    Each time {t} corresponds to column {t-tMin} of the plot grid. 
    Each individual {dt[iq]} is drawn on rows {rdr[iq]}.
    If the life span {dt[iq].tbr..dt[iq].tdt} is disjoint from 
    {0..ncols-1}, the individual is not drawn and {rdr[iq]} is 
    ignored. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */
  
int32_t main(int32_t argc, char **argv)
  { 
    tdp_options_t *o = tdp_parse_options(argc, argv);
    int32_t ni = o->nIndivs; /* number of indivs actually created. */
    
    drtree_node_t dt[ni];
    drtree_test_create_individuals(ni, o->nRoots,o->tStart,o->tStop,o->orphans,o->ageMax,o->nchMax, dt);
    
    /* Determine full time range {tMax..tMin}: */
    int32_t tMin = INT32_MAX;
    int32_t tMax = INT32_MIN;
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_node_t *q = &(dt[iq]);
        assert(! drtree_node_is_null(q));
        if (q->tbr < tMin) { tMin = q->tbr; }
        if (q->tdt > tMax) { tMax = q->tdt; }
      }

    drtree_check_nodes(ni, dt, tMin, tMax);
      
    tdp_do_test(o, FALSE,  "comp-A", tMin, tMax, ni, dt);
    tdp_do_test(o, FALSE,  "comp-B", o->tStart, o->tStop, ni, dt);
    tdp_do_test(o, TRUE,   "plan-A", tMin, tMax, ni, dt);
    tdp_do_test(o, TRUE,   "plan-B", o->tStart, o->tStop, ni, dt);
    
    return 0;
  } 
  
void tdp_do_test
  ( tdp_options_t *o,
    bool_t planar,
    char *tag,
    int32_t t0,
    int32_t t1,
    int32_t ni,
    drtree_node_t dt[]
  )
  {
    int32_t ncols, nrows;
    int32_t rdr[ni];
    drtree_node_t *dc = drtree_clip_time_range(ni, dt, t0, t1);
    drtree_check_nodes(ni, dc, t0, t1);
    if (planar)
      { drtree_planar_arrange(ni, dc, t0, t1, rdr, &ncols, &nrows); }
    else
      { drtree_compact_arrange(ni, dc, t0, t1, rdr, &ncols, &nrows); }
    fprintf(stderr, "  ncols = %d (times {%d .. %d}) nrows = %d\n", ncols, t0, t0+ncols-1, nrows);
    tdp_plot_named(o->outPrefix, tag, t0, t1, ncols, nrows, ni, dc, rdr);
  }

void tdp_plot_named
  ( char *outPrefix, 
    char *tag,
    int32_t tMin,
    int32_t tMax,
    int32_t ncols,
    int32_t nrows,
    int32_t ni, 
    drtree_node_t dt[],
    int32_t rdr[]
  )
  {
    char *name = jsprintf("%s-%s", outPrefix, tag);
    
    double Xstep = 1.0; /* Cell width (mm). */
    double Ystep = 1.0; /* Cell height (mm). */
    
    fprintf(stderr, "  ncols = %d (times {%d .. %d}) nrows = %d\n", ncols, tMin, tMax, nrows);
    assert(tMax == tMin + ncols - 1);
    
    epswr_figure_t *eps = drtree_plot_create_eps_figure(name, ncols, nrows, Xstep, Ystep);
    free(name);
    
    /* Decide which birth dots should be filled: */
    bool_t fill[ni];
    for (int32_t iq = 0; iq < ni; iq++) { fill[iq] = ((iq % 3) != 0); }

    drtree_plot_individuals(eps, tMin, ncols, nrows, Xstep, Ystep, ni, dt, rdr, fill, NULL);

    epswr_set_pen(eps, 0,0,0, 0.25, 0,0);
    epswr_frame(eps);
  
    epswr_end_figure(eps);
  }

tdp_options_t *tdp_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    tdp_options_t *o = drtree_test_parse_options(pp);
    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
