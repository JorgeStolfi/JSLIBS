/* pnmift_root_cost_fn.h - some arc cost functions for PNM segmentation */
/* Last edited on 2010-06-06 16:54:11 by stolfi */

#ifndef pnmift_root_cost_fn_H
#define pnmift_root_cost_fn_H

#include <frgb.h>

#include <ift.h>

typedef double pnmift_root_cost_t;

typedef pnmift_root_cost_t pnmift_root_cost_fn_t(frgb_t q, int chns);
  /* Type of a function that returns the cost of a trivial path, given the RGB value {q}
    of the path's origin. */

pnmift_root_cost_fn_t *pnmift_root_cost_fn_from_name(char *name);
  /* Returns the function {pnmift_root_cost_fn_{name}} from the list below. */
  
#define pnmift_root_cost_fn_INFO \
  "        zero\n" \
  "        " pnmift_root_cost_fn_zero_INFO "\n" \
  "\n" \
  "        lum\n" \
  "        " pnmift_root_cost_fn_lum_INFO 

pnmift_root_cost_t pnmift_root_cost_fn_zero(frgb_t q, int chns);

#define pnmift_root_cost_fn_zero_INFO \
  "  The root cost is zero."

pnmift_root_cost_t pnmift_root_cost_fn_lum(frgb_t q, int chns);

#define pnmift_root_cost_fn_lum_INFO \
  "  The root cost is the luminance of root pixel {q}."

#endif
