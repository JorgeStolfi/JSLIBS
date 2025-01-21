/* See sample_scaling.h */
/* Last edited on 2025-01-21 18:27:19 by stolfi */

#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <float_image.h>
#include <argparser.h>
#include <affirm.h> 

#include <sample_scaling.h>
#include <pst_basic.h>
#include <pst_argparser.h>

/* INTERNAL PROTOTYPES */

double_vec_t sample_scaling_parse_range_option(argparser_t *pp, char *key, int32_t *NC);
  /* If the keyword {key} is present, marks it as parsed, then 
    parses the following arguments as a tuple of one more numbers,
    using {pst_double_vec_parse}.
    
    The parameter {NC} has the same meaning as in {sample_scaling_parse_options_any}. */

bool_t sample_scaling_parse_uniform(argparser_t *pp, bool_t next);
  /* Parses the \"-uniform\" switch; returns TRUE if 
     present, FALSE otherwise.  
     
     If {next} is TRUE, the procedure looks for the keyword only at
     the next argument; otherwise it looks among all arguments that
     are still unparsed.  */
    
void sample_scaling_use_actual_range
  ( float smin,
    float smax,
    double *min_P, 
    double *max_P, 
    double *ctr_P, 
    double *wid_P
  );
  /* Makes sure that no more than two of the scaling range parameters
    {*min_P,*max_P,*ctr_P,*wid_P} are unspecified ({±INF} or {NAN}). The needed
    defaults are obtained from the given range [{smin} _ {smax}].
    
    If exactly two of the arguments are already defined, no
    adjustments are made. Fails if three or more of the arguments are
    already defined. The other cases are described in
    {sample_scaling_set_defaults_from_actual_range_INFO}. */

void sample_scaling_complete_params
  ( double *min_P, 
    double *max_P, 
    double *ctr_P, 
    double *wid_P
  );  
  /* Provides values for any scaling range arguments
    {*min,*max,*ctr,*wid} that are unspecified ({±INF} or {NAN}),
    by using the fundamental equations.  Exactly two of those
    four values must have been specified; the procedure fails 
    otherwise.

    The missing values for the arguments {VMIN}, {VMAX}, {VCTR} and
    {VWID} are chosen so as to satisfy the equations {VWID = VMAX -
    VMIN} and {VCTR = (VMIN + VMAX)/2}. Thus, if any two of those four
    arguments are given, the other two are determined from these
    equations. */

void sample_scaling_debug_vecs
  ( char *label, 
    int32_vec_t *channel, 
    sample_scaling_options_t *sop 
  );

void sample_scaling_debug_params
  ( char *label, 
    int32_t cin, 
    int32_t cot,
    double *min,
    double *max,
    double *ctr,
    double *wid
  );

uint32_t sample_scaling_n_known(double a, double b, double c, double d);
  /* Return the number (in {0..4}) of arguments that are neither {±INF} nor {NAN}. */
  
void sample_scaling_debug_param(char *label, double *par);

/* IMPLEMENTATIONS */

sample_scaling_options_t sample_scaling_parse_options_any(argparser_t *pp, int32_t *NC_P)
  { sample_scaling_options_t sop;
    sop.min = sample_scaling_parse_range_option(pp, "-min",    NC_P);
    sop.max = sample_scaling_parse_range_option(pp, "-max",    NC_P);
    sop.ctr = sample_scaling_parse_range_option(pp, "-center", NC_P);
    sop.wid = sample_scaling_parse_range_option(pp, "-width",  NC_P);

    /* Uniform scaling option: */
    sop.uniform = sample_scaling_parse_uniform(pp, FALSE);
    return sop;
  }


double_vec_t sample_scaling_parse_range_option(argparser_t *pp, char *key, int32_t *NC)
  { if (argparser_keyword_present(pp, key))
      { return pst_double_vec_parse(pp, NC); }
    else
      { return double_vec_new(0); }
  }
  
uint32_t sample_scaling_n_known(double a, double b, double c, double d)
  { uint32_t n = 0;
    if (isfinite(a)) { n++; }
    if (isfinite(b)) { n++; }
    if (isfinite(c)) { n++; }
    if (isfinite(d)) { n++; }
    return n;
  }

void sample_scaling_complete_params
  ( double *min_P, 
    double *max_P, 
    double *ctr_P, 
    double *wid_P
  )
  { uint32_t nknown = sample_scaling_n_known(*min_P, *max_P, *ctr_P, *wid_P);
    demand(nknown == 2, "not enough parameters to define the sample scaling");
    if (isfinite(*min_P) && isfinite(*max_P))
      { /* Compute {wid} and {ctr} from {min} and {max}: */
        (*ctr_P) = ((*min_P) + (*max_P))/2;
        (*wid_P) = (*max_P) - (*min_P);
      }
    else if (isfinite(*ctr_P) && isfinite(*wid_P))
      { /* Compute {min} and {max} from {wid} and {ctr}: */
        (*min_P) = (*ctr_P) - (*wid_P)/2;
        (*max_P) = (*ctr_P) + (*wid_P)/2;
      }
    else if (isfinite(*min_P) && isfinite(*wid_P))
      { /* Compute {max} and {ctr} from {min} and {wid}: */
        (*max_P) = (*min_P) + (*wid_P);
        (*ctr_P) = (*min_P) + (*wid_P)/2;
      }
    else if (isfinite(*max_P) && isfinite(*wid_P))
      { /* Compute {min} and {ctr} from {max} and {wid}: */
        (*min_P) = (*max_P) - (*wid_P);
        (*ctr_P) = (*max_P) - (*wid_P)/2;
      }
    else if (isfinite(*min_P) && isfinite(*ctr_P))
      { /* Compute {max} and {wid} from {min} and {ctr}: */
        (*wid_P) = 2*((*ctr_P) - (*min_P));
        (*max_P) = (*ctr_P) + (*wid_P)/2;
      }
    else  if (isfinite(*max_P) && isfinite(*ctr_P))
      { /* Compute {min} and {wid} from {max} and {ctr}: */
        (*wid_P) = 2*((*max_P) - (*ctr_P));
        (*min_P) = (*ctr_P) - (*wid_P)/2;
      }
    else
      { assert(FALSE); }
    
    /* Adjust {min,max} to avoid divide by zero: */
    if ((*min_P) == (*max_P))
      { (*wid_P) = 1.0e-6 * fabs(*ctr_P) + 1.0e-20; 
        (*min_P) = (*ctr_P) - (*wid_P)/2;
        (*max_P) = (*ctr_P) + (*wid_P)/2;
        assert((*min_P) != (*max_P));
      }
  }

void sample_scaling_use_actual_range
  ( float smin,
    float smax,
    double *min_P, 
    double *max_P, 
    double *ctr_P, 
    double *wid_P
  ) 
  { 
    demand(isfinite(smin) && isfinite(smax), "input sample range is indeterminate");
    uint32_t nknown = sample_scaling_n_known(*min_P, *max_P, *ctr_P, *wid_P);
    demand(nknown <= 2, "should specify at most two scaling parameters");
    if (nknown < 2)
      { /* At most one argument is known. */
        /* Must determine actual min and max sample values: */
        /* Now use them to complete two known values: */
        if (isfinite(*min_P))
          { (*max_P) = (smax > (*min_P) ? smax : (*min_P)); }
        else if (isfinite(*max_P))
          { (*min_P) = (smin < (*max_P) ? smin : (*max_P)); }
        else if (isfinite(*ctr_P))
          { double dmin = (*ctr_P) - (double)smin;
            double dmax = (double)smax - (*ctr_P);
            (*wid_P) = 2 * (dmin > dmax ? dmin : dmax);
          }
        else if (isfinite(*wid_P))
          { (*ctr_P) = ((double)smax + (double)smin)/2; }
        else 
          { /* No args given, use natural range. */
            (*min_P) = (double)smin;
            (*max_P) = (double)smax;
          }
      }
  }

bool_t sample_scaling_parse_uniform(argparser_t *pp, bool_t next)
  { 
    return pst_keyword_present(pp, "-uniform", next);
  }

void sample_scaling_fix_channels(int32_t NC, int32_vec_t *channel)
  { int32_t c;
    if ((channel->ne == 0) && (NC > 0))
      { /* Provide default channel indices: */
        int32_vec_expand(channel, (vec_index_t)NC-1);
        for (c = 0; c < NC; c++) { channel->e[c] = c; }
        int32_vec_trim(channel, (vec_size_t)NC);
      }
    else if ((channel->ne == 1) && (NC > 1))
      { /* Replicate the selected channel: */
        int32_vec_expand(channel, (vec_index_t)NC-1);
        for (c = 1; c < NC; c++) { channel->e[c] = channel->e[0]; }
        int32_vec_trim(channel, (vec_size_t)NC);
      }
    /* Checking (always, just for paranoia): */
    demand(channel->ne == NC, "inconsistent number of channels");
  }  

void sample_scaling_fix_params
  ( sample_scaling_options_t *sop, 
    int32_vec_t *channel, 
    int32_t NC,
    float_image_t *fim
  )
  { bool_t debug = FALSE;
    demand((channel == NULL) || (channel->ne == NC), "inconsistent channel map");
    if (sop->uniform)
      { /* Uniformize all scaling parameters across all channels, shrink vectors to 1 elem: */
        pst_double_vec_uniformize(&(sop->min), +INF); 
        pst_double_vec_uniformize(&(sop->max), +INF); 
        pst_double_vec_uniformize(&(sop->ctr), +INF); 
        pst_double_vec_uniformize(&(sop->wid), +INF); 
        if (debug) { sample_scaling_debug_vecs("uniformized ranges", NULL, sop); }
        /* Check for over-specification: */
        double *min0 = &(sop->min.e[0]);
        double *max0 = &(sop->max.e[0]);
        double *ctr0 = &(sop->ctr.e[0]);
        double *wid0 = &(sop->wid.e[0]);
        uint32_t nknown = sample_scaling_n_known(*min0, *max0, *ctr0, *wid0);
        demand(nknown <= 2, "should specify at most two scaling parameters");
        /* Obtain scaling data from image if under-specified: */
        if (nknown < 2) 
          { if (fim != NULL)
              { /* Get the actual sample range {smin,smax} over all relevant channels: */
                float smin = +INF; float smax = -INF;
                for (int32_t k = 0; k < NC; k++) 
                  { int32_t c = (channel == NULL ? k : channel->e[k]);
                    float_image_update_sample_range(fim, (int32_t)c, &smin, &smax);
                  }
                /* Use the actual sample range {smin,smax} to define the scaling: */
                sample_scaling_use_actual_range(smin, smax, min0, max0, ctr0, wid0);
                if (debug) { sample_scaling_debug_vecs("unif. fixed ranges", NULL, sop); }
              }
            else
              { demand(FALSE, "cannot uniformize the sample scaling without the image"); }
          }
        /* Now we must have exactly two known parameters: */
        nknown = sample_scaling_n_known(*min0, *max0, *ctr0, *wid0);
        assert(nknown == 2);
        /* Complete the missing parameters: */
        sample_scaling_complete_params(min0, max0, ctr0, wid0);
        if (debug) { sample_scaling_debug_vecs("unif. compl ranges", NULL, sop); }
        nknown = sample_scaling_n_known(*min0, *max0, *ctr0, *wid0);
        assert(nknown == 4);
      }
      
    /* Make sure that all scaling arg tuples have {NC} elements: */
    pst_double_vec_regularize(&(sop->min), NC, +INF);   
    pst_double_vec_regularize(&(sop->max), NC, +INF);   
    pst_double_vec_regularize(&(sop->ctr), NC, +INF);
    pst_double_vec_regularize(&(sop->wid), NC, +INF); 
    if (debug) { sample_scaling_debug_vecs("regularized ranges", channel, sop); }

    /* Complete the scalings for each channel: */
    for (int32_t k = 0; k < NC; k++)
      { /* Get scaling params for output channel {k}: */
        int32_t c = (channel == NULL ? k : channel->e[k]); 
        double *mink = &(sop->min.e[k]);
        double *maxk = &(sop->max.e[k]);
        double *ctrk = &(sop->ctr.e[k]);
        double *widk = &(sop->wid.e[k]);
        if (debug) { sample_scaling_debug_params("chan orig ranges", c, k, mink, maxk, ctrk, widk); }
        /* Check for over-specification: */
        uint32_t nknown = sample_scaling_n_known(*mink, *maxk, *ctrk, *widk);
        if (sop->uniform) 
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
                    float_image_update_sample_range(fim, (int32_t)c, &smin, &smax);
                    sample_scaling_use_actual_range(smin, smax, mink, maxk, ctrk, widk);
                    if (debug) { sample_scaling_debug_params("chan fixd ranges", c, k, mink, maxk, ctrk, widk); }
                  }
                else
                  { demand(FALSE, "cannot determine the sample scaling without the image"); }
              }
            /* Now we must have exactly two known parameters for this channel: */
            nknown = sample_scaling_n_known(*mink, *maxk, *ctrk, *widk);
            assert(nknown == 2); 
            /* Complete the missing parameters: */
            sample_scaling_complete_params(mink, maxk, ctrk, widk);
          }
        if (debug) { sample_scaling_debug_params("chan compl ranges", c, k, mink, maxk, ctrk, widk); }
        nknown = sample_scaling_n_known(*mink, *maxk, *ctrk, *widk);
        assert(nknown == 4);   
      }
  }

void sample_scaling_debug_vecs
  ( char *label, 
    int32_vec_t *channel, 
    sample_scaling_options_t *sop
  )
  { if (label != NULL) { fprintf(stderr, "%s\n", label); }
    fprintf(stderr, "  lengths:");
    if (channel != NULL) { fprintf(stderr, "  channel = %d", channel->ne); }
    fprintf(stderr, "  min = %d", sop->min.ne);
    fprintf(stderr, "  max = %d", sop->max.ne);
    fprintf(stderr, "  ctr = %d", sop->ctr.ne);
    fprintf(stderr, "  wid = %d", sop->wid.ne);
    fprintf(stderr, "\n");
    int32_t k = 0;
    while (TRUE)
      { double *mink = (k < sop->min.ne ? &(sop->min.e[k]) : NULL);
        double *maxk = (k < sop->max.ne ? &(sop->max.e[k]) : NULL);
        double *ctrk = (k < sop->ctr.ne ? &(sop->ctr.e[k]) : NULL);
        double *widk = (k < sop->wid.ne ? &(sop->wid.e[k]) : NULL);
        if ((mink == NULL) && (maxk == NULL) && (ctrk == NULL) && (widk == NULL)) { break; }
        if ((channel == NULL) || (k >= channel->ne)) { break; } 
        int32_t c = channel->e[k];
        sample_scaling_debug_params("  ", c, k, mink, maxk, ctrk, widk);
        k++;
      }
    fprintf(stderr, "\n");
  }

void sample_scaling_debug_params
  ( char *label, 
    int32_t cin, 
    int32_t cot,
    double *min,
    double *max,
    double *ctr,
    double *wid
  )
  {
    if (label != NULL) { fprintf(stderr, "%s", label); }
    fprintf(stderr, "  cin = %2d", cin);
    fprintf(stderr, "  cot = %2d", cot);
    sample_scaling_debug_param("min", min); 
    sample_scaling_debug_param("max", max); 
    sample_scaling_debug_param("ctr", ctr); 
    sample_scaling_debug_param("wid", wid); 
    fprintf(stderr, "\n");
  }

void sample_scaling_debug_param(char *label, double *par)
  { 
    fprintf(stderr, "  %s = ", label);
    if (par == NULL) 
      { fprintf(stderr, "------------------------"); } 
    else if (isfinite(*par)) 
      { fprintf(stderr, "%24.16e", *par); } 
    else
      { fprintf(stderr, "????????????????????????"); } 
  }
