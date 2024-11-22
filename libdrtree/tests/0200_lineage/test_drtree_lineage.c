#define PROG_NAME "test_drtree_lineage"
#define PROG_VERS "2.0"

/* Draws the evolution of asexual traits. */

/* Last edited on 2023-06-25 18:11:00 by stolfi */

#define test_drtree_lineage_C_COPYRIGHT \
  "???"

#define PROG_HELP \
  " " PROG_NAME " \\\n" \
  drtree_test_options_HELP
 
#define PROG_INFO \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Tests the {drtree_lineage.h} and {drtree_plot.h} functions. The test uses a randomly" \
  " generated evolution history.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The main output files are \"{outPrefix}-{tag}-plot.eps\" where {tag} is \"A\", \"B\", etc..\n" \
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
#include <drtree_lineage.h>

typedef drtree_test_options_t tdl_options_t;

tdl_options_t *tdl_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void tdl_do_test
  ( tdl_options_t *o,
    char *tag,
    int32_t ni,
    drtree_node_t dt[],
    int32_t tRef
  );
  /* Creates an EPS file "{o.outPrefix}-{tag}.eps" with the compact diagram
    of the nodes {dt[0..ni-1}} in the time interval {tMin..tMax}
    that comprises all life spans.
    
    Reference lines are drawn at times {o->tStart,o->tStop,o->tRef}.  
    Lineages {L(q|o->tRef|o->tStop)} are highlighted with different colors. */

void tdl_plot_named
  ( char *outPrefix, 
    char *tag,
    int32_t tMin,
    int32_t tMax,
    int32_t ni, 
    drtree_node_t dt[],
    int32_t ncols,
    int32_t nrows,
    int32_t rdr[],
    int32_t tStart,
    int32_t tStop,
    int32_t tRef
  );       
  /* Creates an EPS file "{outPrefix}-{tag}.eps" and draws the diagram
    of the nodes {dt[0..ni-1}} on rows {rdr[0..ni-1]} of the plot grid
    spanning times {tMin..tMax}.
    
    Each time {t} corresponds to column {t-tMin} of the plot grid. 
    Each individual {dt[iq]} is drawn on rows {rdr[iq]}.
    If the life span {dt[iq].tbr..dt[iq].tdt} is disjoint from 
    {0..ncols-1}, the individual is not drawn and {rdr[iq]} is 
    ignored.
    
    Reference lines are drawn at times {tStart,tStop,tRef}
    if they are in the grids's column range.  Also highlights 
    lineages that started at {tRef} and survived until {tStop}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */
  
int32_t main(int32_t argc, char **argv)
  { 
    tdl_options_t *o = tdl_parse_options(argc, argv);
    int32_t ni = o->nIndivs; /* number of indivs actually created. */
    
    drtree_node_t dt[ni];
    drtree_test_create_individuals(ni, o->nRoots,o->tStart,o->tStop,o->orphans,o->ageMax,o->nchMax, dt);
      
    tdl_do_test(o, "A", ni, dt, o->tRef);
    tdl_do_test(o, "B", ni, dt, o->tRef + 3*o->ageMax + 1);
    
    return 0;
  } 
  
void tdl_do_test
  ( tdl_options_t *o,
    char *tag,
    int32_t t0,
    drtree_node_t dt[],
    int32_t tRef
  )
  {
    bool_t planar = TRUE;
    
    int32_t ni = o->nIndivs;
    
    /* Determine full time range: */
    int32_t tMin = INT32_MAX;
    int32_t tMax = INT32_MIN;
    for (int32_t iq = 0; iq < ni; iq++)
      { assert(dt[iq].tbr <= dt[iq].tdt);
        if (dt[iq].tbr < tMin) { tMin = dt[iq].tbr; }
        if (dt[iq].tdt > tMax) { tMax = dt[iq].tdt; }
      }

    int32_t ncols, nrows;
    int32_t rdr[ni];
    if (planar)
      { drtree_planar_arrange(ni, dt, tMin, tMax, rdr, &ncols, &nrows); }
    else
      { drtree_compact_arrange(ni, dt, tMin, tMax, rdr, &ncols, &nrows); }
    fprintf(stderr, "  ncols = %d (times {%d .. %d}) nrows = %d\n", ncols, tMin, tMin+ncols-1, nrows);
    assert(tMax == tMin+ncols-1);
    tdl_plot_named(o->outPrefix, tag, tMin, tMax, ni, dt, ncols, nrows, rdr, o->tStart, o->tStop, tRef);
  }

void tdl_plot_named
  ( char *outPrefix, 
    char *tag,
    int32_t tMin,
    int32_t tMax,
    int32_t ni, 
    drtree_node_t dt[],
    int32_t ncols,
    int32_t nrows,
    int32_t rdr[],
    int32_t tStart,
    int32_t tStop,
    int32_t tRef
  )
  {
    char *name = jsprintf("%s-%s", outPrefix, tag);
    
    double Xstep = 1.0; /* Cell width (mm). */
    double Ystep = 1.0; /* Cell height (mm). */
    
    epswr_figure_t *eps = drtree_plot_create_eps_figure(name, ncols, nrows, Xstep, Ystep);
    free(name);
    
    /* Fing surviving lineages and their founders: */
    int32_t *fnd = drtree_lineage_collect_surviving(ni, dt, tRef, tStop);
    
    /* Founders birth dots are filled: */
    bool_t fill[ni];
    for (int32_t iq = 0; iq < ni; iq++) { fill[iq] = (fnd[iq] == iq); }

    frgb_t rgbStart = (frgb_t){{ 1.000f, 0.200f, 0.000f }};
    drtree_plot_time_line(eps, tMin, ncols, nrows, Xstep, Ystep, &(rgbStart), tStart);
    
    frgb_t rgbStop = (frgb_t){{ 0.200f, 0.200f, 1.000f }};
    drtree_plot_time_line(eps, tMin, ncols, nrows, Xstep, Ystep, &(rgbStop), tStop);
    
    frgb_t rgbRef = (frgb_t){{ 0.200f, 0.700f, 0.200f }};
    drtree_plot_time_line(eps, tMin, ncols, nrows, Xstep, Ystep, &(rgbRef), tRef);

    drtree_plot_individuals(eps, tMin, ncols, nrows, Xstep, Ystep, ni, dt, rdr, fill, fnd);

    epswr_set_pen(eps, 0,0,0, 0.25, 0,0);
    epswr_frame(eps);
  
    epswr_end_figure(eps);
  }

tdl_options_t *tdl_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    tdl_options_t *o = drtree_test_parse_options(pp);
    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
