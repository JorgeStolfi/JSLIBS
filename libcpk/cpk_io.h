#ifndef cpk_io_H
#define cpk_io_H

/* Input/output routines for circle packing. */
/* Last edited on 2025-01-01 02:45:42 by stolfi */ 

#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <interval.h>
#include <r2.h>

#include <cpk_basic.h>
#include <cpk_main.h>
#include <cpk_eps.h>

/* Formats for EUTM X/Y coordinate printout: */
#define XY_FMT "%12.1f"
#define XY_FFMT "%.1f"

/* Formats for Lon/Lat coordinate printout: */
#define LL_FMT "%20.15f"
#define LL_FFMT "%.15f"

/* Formats for candidate vertex weight: */
#define WT_FMT "%12.11f"
#define WT_FFMT "%.11f"

/* Formats for vertex index printout: */
#define VTX_FMT "%7d"
#define VTX_FFMT "%d"

void cpk_plot_solution
  ( char *outDir,
    char *solTag,
    char *demTag,
    interval_t B[],
    cpk_domain_t *C,
    cpk_policy_t *P,
    r2_vec_t *V,
    double_vec_t *W,
    uint32_vec_t *J
  );
  /* Creates a Postscript file "{outDir}/{solTag}-{demTag}.eps"
    showing the problem data and the solution {J}.
    
    The other parameters are: the rectangle to plot {B}, the set of
    valid auction points {C}, the demand points {A}, their radius of
    adequacy {rAdeq}, the EUTM coordinates {V} of the candidate
    auction points (grid points), their weights {W}, the radius of
    auction areas {rAuc}, and the nominal station broadcast radius
    {rSRange}. */

/* 
  READING GEOMETRY FILES 
  
  Each of the procedures {cpk_read_points} and {cpk_read_polys} below
  reads geometric data (points or polylines, respectively) from a file
  named is "{inDir}/{tag}.txt".
  
  Both procedures ignore leading and trailing blanks in each line of
  the file, as well as lines that are entirely blank. Any line whose
  first non-blank character is [0-9.+-] is assumed to be a data line;
  it must contain two numbers, separated by white space. The two
  numbers are assumed to be Longitude and Latitude of a point, in that
  order, expressed as fractional numbers of degrees.
  
  Lines whose first non-blank character is anything other than [0-9.+-]
  is assumed to be a /separator/, whis is treated as described below
  for each routine. */

r2_vec_t cpk_read_points(char *inDir, char *tag);
  /* Reads a list of points from the file "{inDir}/{tag}.txt".
    Ignores separator lines (whose first non-blank character is other
    than [0-9.+-]). */

r2_vec_t cpk_read_polys(char *inDir, char *tag, bool_t closed);
  /* Reads a list of points from the file "{inDir}/{tag}.txt",
    interprets them as the vertices of zero or more polylines. 

    If {closed=TRUE}, the points are interpreted as vertices of a
    polygon. The polygon may have two or more components (islands or
    holes); each of them must be closed, i.e. the last vertex must be
    identical to the first one. If {closed = FALSE}, the points are
    assumed to be vertices of a polyline. A polyline too may have
    several components, but they not have to be closed.

    In either case, the end of each component is defined either by the
    repetition of its first vertex, or by a separator line, or by
    end-of-file. Extra separator lines may be present before or after
    each component. In the returned list, one infinite point {# = (INF,INF)}
    is inserted between every two consecutive components. */

void cpk_pick_ref_coords(r2_vec_t *P, double *refLon, double *refY);
  /* Given a list of Lat/Lon pairs {P}, computes a mean point 
    suitable to use as a reference for Lat/Lon to EUTM conversion,
    and returns its longitude {refLon} and its EUTM Y-coordinate {refY}.
    Any infinite points in {P} are ignored. */

r2_vec_t cpk_LL_to_EUTM(r2_vec_t *P, double refLon, double refY, double magnify);
  /* Converts all finite points of {P} from Lat/Lon to EUTM X/Y, using
    {refLon} as the reference meridian, and {refY} as the reference
    Y-coordinate. All finite XY coordinates are multiplied by
    {magnify} after the conversion. */

r2_vec_t cpk_EUTM_to_LL(r2_vec_t *P, double refLon, double refY, double magnify);
  /* Converts all finite points of {P} from EUTM to Lat/Lon, using
    {refLon} as the reference meridian, and {refY} as the reference
    Y-coordinate. All finite XY coordinates are divided by {magnify}
    before the conversion. */

void cpk_plot_domain(epswr_figure_t *eps, cpk_domain_t *C, cpk_policy_t *P);
  /* Plot the domain {C}, including urban area {Urb}, existing/planned
    stations {Exs}, predefined auctions {Auc}, municipality borders
    {Mun}, and national borders {Nat}. */
  
void cpk_plot_stations
  ( epswr_figure_t *eps, 
    r2_vec_t *P, 
    double rUnc,
    double rRad,
    frgb_t *color
  );
  /* Plots a set of stations showing for each the position {P[k]}, 
    an uncertainty circle (with radius {rUnc}), and 
    the possible broadcast coverage area (with radius {rUnc+rRad}),
    labeled with the index {k}.  Typically one would use {rUnc==0} 
    for existing or planned stations with known coordinates,
    and {rUnc=rAuc} for predetermined auction areas
    whose winner is not yet known. */
   
void cpk_plot_demand(epswr_figure_t *eps, r2_vec_t *P, double rDem);
  /* Plots a set of demand points showing for each the position {P[k]}
    and the area with radius {rDem} where an auction center would be
    interesting to this demand.  Usually {rDem == rAuc}. */
  
void cpk_plot_candidates(epswr_figure_t *eps, r2_vec_t *V, double_vec_t *W);
  /* Plot the candidates {V[i]} and their weights {W[i]} as dots on {eps}. */

void cpk_plot_proposed_auctions
  ( epswr_figure_t *eps, 
    uint32_vec_t *J, 
    r2_vec_t *V, 
    double rAuc, 
    double rRad
  );
  /* Plots a set of proposed auctions as pairs of circles,
    with centers {V[J[j]]} and radii {rAuc} (auctioned area)
    and {rAuc + rRad} (possible broadcast range area). */

void cpk_plot_solution_attributes
  ( epswr_figure_t *eps, 
    interval_t B[], 
    cpk_policy_t *P,
    double_vec_t *W,
    uint32_vec_t *J
  );
  /* Prints the number of auctions {J.nel} and the total weight
    of the auction set {SUM{W[J[i]] : i=0..J.ne-1}}, near
    the upper left corner of the box {B}. */

/* PRINTING */

void cpk_ui2_print(FILE *wr, ui2_t *x, char *fmt);
  /* Prints {x} to {wr} in a standard format, without end-of-line.
    Each component is printed with format {fmt}. */
    
#endif
