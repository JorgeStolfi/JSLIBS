/* pnmift_arc_cost_fn.h - some arc cost functions for PNM segmentation */
/* Last edited on 2010-06-07 13:21:04 by stolfi */

#ifndef pnmift_arc_cost_fn_H
#define pnmift_arc_cost_fn_H

#include <frgb.h>

#include <ift.h>

typedef double pnmift_arc_cost_t;

typedef pnmift_arc_cost_t pnmift_arc_cost_fn_t(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns);
  /* Type of a function that returns an arc cost given the RGB values of two pixels {p,q}
    and the relative arc {ra} connecting them. */

pnmift_arc_cost_fn_t *pnmift_arc_cost_fn_from_name(char *name);
  /* Returns the function {pnmift_arc_cost_fn_{name}} from the list below. */
  
#define pnmift_arc_cost_fn_INFO \
  "        ediff_rgb\n" \
  "        " pnmift_arc_cost_fn_ediff_rgb_INFO "\n" \
  "\n" \
  "        ediff_yuv\n" \
  "        " pnmift_arc_cost_fn_ediff_yuv_INFO "\n" \
  "\n" \
  "        lum\n" \
  "        " pnmift_arc_cost_fn_lum_INFO 

pnmift_arc_cost_t pnmift_arc_cost_fn_ediff_rgb(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns);

#define pnmift_arc_cost_fn_ediff_rgb_INFO \
  "  If {chns=3} the result is the Euclidean difference between the RGB values of {p} and" \
  " {q} divided by {sqrt(3)} so that it ranges in [0_1].  If {chns=1} assumes that {p} and {q} are grays and returns" \
  " the difference between the luminances {Y}.  Undefined for other values of {chns}."

pnmift_arc_cost_t pnmift_arc_cost_fn_ediff_yuv(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns);

#define pnmift_arc_cost_fn_ediff_yuv_INFO \
  "  If {chns=3} the result is the Euclidean difference between the values of {p} and {q} converted" \
  " from RGB to yuv color space, divided by {sqrt(3)} so that it ranges in [0_1].  The {yuv}" \
  " coordiantes are the {YUV} coordintes scaled by {(1+e)/(Y+e)} where {e = 0.01}.  If {chns=1} assumes" \
  " that {p} and {q} are grays and returns the difference between the scaled luminance" \
  " {y = Y*(1+e)/(Y+e)}.  Undefined for other values of {chns}."

pnmift_arc_cost_t pnmift_arc_cost_fn_lum(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns);

#define pnmift_arc_cost_fn_lum_INFO \
  "  The arc cost is the luminance of destination vertex {q}."

#endif
