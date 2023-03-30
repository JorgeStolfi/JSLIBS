/* See pst_scaling.h */
/* Last edited on 2023-03-19 15:29:03 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <argparser.h>
#include <affirm.h> 

#include <pst_scaling.h>
#include <pst_basic.h>
#include <pst_argparser.h>

/* INTERNAL PROTOTYPES */

void pst_scaling_debug_vecs
  ( char *label, 
    int32_vec_t *channel, 
    double_vec_t *min,
    double_vec_t *max,
    double_vec_t *ctr,
    double_vec_t *wid
  );

void pst_scaling_debug_params
  ( char *label, 
    int cin, 
    int cot,
    double *min,
    double *max,
    double *ctr,
    double *wid
  );

void pst_scaling_debug_param(char *label, double *par);

/* IMPLEMENTATIONS */

double_vec_t pst_scaling_parse_range_option(argparser_t *pp, char *key, int *NC)
  { if (argparser_keyword_present(pp, key))
      { return pst_double_vec_parse(pp, NC); }
    else
      { return double_vec_new(0); }
  }

#define KNOWN(x) (fabs(x) != INF)
  /* TRUE iff scaling argument {x} has a known value. */

void pst_scaling_complete_params
  ( double *min, 
    double *max, 
    double *ctr, 
    double *wid
  )
  { int nknown = KNOWN(*min) + KNOWN(*max) + KNOWN(*ctr) + KNOWN(*wid);
    demand(nknown == 2, "not enough parameters to define the sample scaling");
    if (KNOWN(*min) && KNOWN(*max))
      { /* Compute {wid} and {ctr} from {min} and {max}: */
        (*ctr) = ((*min) + (*max))/2;
        (*wid) = (*max) - (*min);
      }
    else if (KNOWN(*ctr) && KNOWN(*wid))
      { /* Compute {min} and {max} from {wid} and {ctr}: */
        (*min) = (*ctr) - (*wid)/2;
        (*max) = (*ctr) + (*wid)/2;
      }
    else if (KNOWN(*min) && KNOWN(*wid))
      { /* Compute {max} and {ctr} from {min} and {wid}: */
        (*max) = (*min) + (*wid);
        (*ctr) = (*min) + (*wid)/2;
      }
    else if (KNOWN(*max) && KNOWN(*wid))
      { /* Compute {min} and {ctr} from {max} and {wid}: */
        (*min) = (*max) - (*wid);
        (*ctr) = (*max) - (*wid)/2;
      }
    else if (KNOWN(*min) && KNOWN(*ctr))
      { /* Compute {max} and {wid} from {min} and {ctr}: */
        (*wid) = 2*((*ctr) - (*min));
        (*max) = (*ctr) + (*wid)/2;
      }
    else  if (KNOWN(*max) && KNOWN(*ctr))
      { /* Compute {min} and {wid} from {max} and {ctr}: */
        (*wid) = 2*((*max) - (*ctr));
        (*min) = (*ctr) - (*wid)/2;
      }
    else
      { assert(FALSE); }
    
    /* Adjust {min,max} to avoid divide by zero: */
    if ((*min) == (*max))
      { (*wid) = 1.0e-6 * fabs(*ctr) + 1.0e-20; 
        (*min) = (*ctr) - (*wid)/2;
        (*max) = (*ctr) + (*wid)/2;
        assert((*min) != (*max));
      }
  }

void pst_scaling_use_actual_range
  ( float smin,
    float smax,
    double *min, 
    double *max, 
    double *ctr, 
    double *wid
  ) 
  { 
    demand(KNOWN(smin) && KNOWN(smax), "input sample range is indeterminate");
    int nknown = KNOWN(*min) + KNOWN(*max) + KNOWN(*ctr) + KNOWN(*wid);
    demand(nknown <= 2, "should specify at most two scaling parameters");
    if (nknown < 2)
      { /* At most one argument is known. */
        /* Must determine actual min and max sample values: */
        /* Now use them to complete two known values: */
        if (KNOWN(*min))
          { (*max) = (smax > (*min) ? smax : (*min)); }
        else if (KNOWN(*max))
          { (*min) = (smin < (*max) ? smin : (*max)); }
        else if (KNOWN(*ctr))
          { double dmin = (*ctr) - (double)smin;
            double dmax = (double)smax - (*ctr);
            (*wid) = 2 * (dmin > dmax ? dmin : dmax);
          }
        else if (KNOWN(*wid))
          { (*ctr) = ((double)smax + (double)smin)/2; }
        else 
          { /* No args given, use natural range. */
            (*min) = (double)smin;
            (*max) = (double)smax;
          }
      }
  }

bool_t pst_scaling_parse_uniform(argparser_t *pp, bool_t next)
  { 
    return pst_keyword_present(pp, "-uniform", next);
  }

void pst_scaling_fix_channels(int NC, int32_vec_t *channel)
  { int c;
    if ((channel->ne == 0) && (NC > 0))
      { /* Provide default channel indices: */
        int32_vec_expand(channel, NC-1);
        for (c = 0; c < NC; c++) { channel->e[c] = c; }
        int32_vec_trim(channel, NC);
      }
    else if ((channel->ne == 1) && (NC > 1))
      { /* Replicate the selected channel: */
        int32_vec_expand(channel, NC-1);
        for (c = 1; c < NC; c++) { channel->e[c] = channel->e[0]; }
        int32_vec_trim(channel, NC);
      }
    /* Checking (always, just for paranoia): */
    demand(channel->ne == NC, "inconsistent number of channels");
  }  

void pst_scaling_fix_params
  ( int NC,
    bool_t uniform,
    double_vec_t *min,
    double_vec_t *max,
    double_vec_t *ctr,
    double_vec_t *wid,
    float_image_t *fim, 
    int32_vec_t *channel
  )
  { bool_t debug = FALSE;
    demand((channel == NULL) || (channel->ne == NC), "inconsistent channel map");
    int k;
    if (uniform)
      { /* Uniformize all scaling parameters across all channels, shrink vectors to 1 elem: */
        pst_double_vec_uniformize(min, +INF); 
        pst_double_vec_uniformize(max, +INF); 
        pst_double_vec_uniformize(ctr, +INF); 
        pst_double_vec_uniformize(wid, +INF); 
        if (debug) { pst_scaling_debug_vecs("uniformized ranges", NULL, min, max, ctr, wid); }
        /* Check for over-specification: */
        double *min0 = &(min->e[0]);
        double *max0 = &(max->e[0]);
        double *ctr0 = &(ctr->e[0]);
        double *wid0 = &(wid->e[0]);
        int nknown = KNOWN(*min0) + KNOWN(*max0) + KNOWN(*ctr0) + KNOWN(*wid0);
        demand(nknown <= 2, "should specify at most two scaling parameters");
        /* Obtain scaling data from image if under-specified: */
        if (nknown < 2) 
          { if (fim != NULL)
              { /* Get the actual sample range {smin,smax} over all relevant channels: */
                float smin = +INF; float smax = -INF;
                for (k = 0; k < NC; k++) 
                  { int c = (channel == NULL ? k : channel->e[k]);
                    float_image_update_sample_range(fim, c, &smin, &smax);
                  }
                /* Use the actual sample range {smin,smax} to define the scaling: */
                pst_scaling_use_actual_range(smin, smax, min0, max0, ctr0, wid0);
                if (debug) { pst_scaling_debug_vecs("unif. fixed ranges", NULL, min, max, ctr, wid); }
              }
            else
              { demand(FALSE, "cannot uniformize the sample scaling without the image"); }
          }
        /* Now we must have exactly two known parameters: */
        nknown = KNOWN(*min0) + KNOWN(*max0) + KNOWN(*ctr0) + KNOWN(*wid0);
        assert(nknown == 2);
        /* Complete the missing parameters: */
        pst_scaling_complete_params(min0, max0, ctr0, wid0);
        if (debug) { pst_scaling_debug_vecs("unif. compl ranges", NULL, min, max, ctr, wid); }
        nknown = KNOWN(*min0) + KNOWN(*max0) + KNOWN(*ctr0) + KNOWN(*wid0);
        assert(nknown == 4);
      }
      
    /* Make sure that all scaling arg tuples have {NC} elements: */
    pst_double_vec_regularize(min, NC, +INF);   
    pst_double_vec_regularize(max, NC, +INF);   
    pst_double_vec_regularize(ctr, NC, +INF);
    pst_double_vec_regularize(wid, NC, +INF); 
    if (debug) { pst_scaling_debug_vecs("regularized ranges", channel, min, max, ctr, wid);}

    /* Complete the scalings for each channel: */
    for (k = 0; k < NC; k++)
      { /* Get scaling params for output channel {k}: */
        int c = (channel == NULL ? k : channel->e[k]); 
        double *mink = &(min->e[k]);
        double *maxk = &(max->e[k]);
        double *ctrk = &(ctr->e[k]);
        double *widk = &(wid->e[k]);
        if (debug) { pst_scaling_debug_params("chan orig ranges", c, k, mink, maxk, ctrk, widk); }
        /* Check for over-specification: */
        int nknown = KNOWN(*mink) + KNOWN(*maxk) + KNOWN(*ctrk) + KNOWN(*widk);
        if (uniform) 
          { /* We must have completed and replicated: */
            assert(nknown == 4); 
          }
        else
          { /* Check for over-specified scaling: */
            demand(nknown <= 2, "should specify at most two scaling parameters");
            /* Obtain scaling data from image if under-specified: */
            if (nknown < 2) 
              { if (fim != NULL)
                  { /* Get the actual sample range {smin,smax} over input channel {c}: */
                    float smin = +INF; float smax = -INF;
                    float_image_update_sample_range(fim, c, &smin, &smax);
                    pst_scaling_use_actual_range(smin, smax, mink, maxk, ctrk, widk);
                    if (debug) { pst_scaling_debug_params("chan fixd ranges", c, k, mink, maxk, ctrk, widk); }
                  }
                else
                  { demand(FALSE, "cannot determine the sample scaling without the image"); }
              }
            /* Now we must have exactly two known parameters for this channel: */
            nknown = KNOWN(*mink) + KNOWN(*maxk) + KNOWN(*ctrk) + KNOWN(*widk);
            assert(nknown == 2); 
            /* Complete the missing parameters: */
            pst_scaling_complete_params(mink, maxk, ctrk, widk);
          }
        if (debug) { pst_scaling_debug_params("chan compl ranges", c, k, mink, maxk, ctrk, widk); }
        nknown = KNOWN(*mink) + KNOWN(*maxk) + KNOWN(*ctrk) + KNOWN(*widk);
        assert(nknown == 4);   
      }
  }

void pst_scaling_debug_vecs
  ( char *label, 
    int32_vec_t *channel, 
    double_vec_t *min,
    double_vec_t *max,
    double_vec_t *ctr,
    double_vec_t *wid
  )
  { if (label != NULL) { fprintf(stderr, "%s\n", label); }
    fprintf(stderr, "  lengths:");
    if (channel != NULL) { fprintf(stderr, "  channel = %d", channel->ne); }
    if (min != NULL) { fprintf(stderr, "  min = %d", min->ne); }
    if (max != NULL) { fprintf(stderr, "  max = %d", max->ne); }
    if (ctr != NULL) { fprintf(stderr, "  ctr = %d", ctr->ne); }
    if (wid != NULL) { fprintf(stderr, "  wid = %d", wid->ne); }
    fprintf(stderr, "\n");
    int k = 0;
    while (TRUE)
      { int c = ((channel != NULL) && (k < channel->ne) ? channel->e[k] : -9); 
        double *mink = ((min != NULL) && (k < min->ne) ? &(min->e[k]) : NULL);
        double *maxk = ((max != NULL) && (k < max->ne) ? &(max->e[k]) : NULL);
        double *ctrk = ((ctr != NULL) && (k < ctr->ne) ? &(ctr->e[k]) : NULL);
        double *widk = ((wid != NULL) && (k < wid->ne) ? &(wid->e[k]) : NULL);
        if ((c < 0) && (mink == NULL) && (maxk == NULL) && (ctrk == NULL) && (widk == NULL)) { break; }
        pst_scaling_debug_params("  ", c, k, mink, maxk, ctrk, widk);
        k++;
      }
    fprintf(stderr, "\n");
  }

void pst_scaling_debug_params
  ( char *label, 
    int cin, 
    int cot,
    double *min,
    double *max,
    double *ctr,
    double *wid
  )
  {
    if (label != NULL) { fprintf(stderr, "%s", label); }
    fprintf(stderr, "  cin = %2d", cin);
    fprintf(stderr, "  cot = %2d", cot);
    pst_scaling_debug_param("min", min); 
    pst_scaling_debug_param("max", max); 
    pst_scaling_debug_param("ctr", ctr); 
    pst_scaling_debug_param("wid", wid); 
    fprintf(stderr, "\n");
  }

void pst_scaling_debug_param(char *label, double *par)
  { 
    fprintf(stderr, "  %s = ", label);
    if (par == NULL) 
      { fprintf(stderr, "------------------------"); } 
    else if (KNOWN(*par)) 
      { fprintf(stderr, "%24.16e", *par); } 
    else
      { fprintf(stderr, "????????????????????????"); } 
  }
