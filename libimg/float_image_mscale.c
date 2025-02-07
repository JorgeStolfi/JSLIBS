/* See {float_image_mscale.h}. */
/* Last edited on 2025-01-12 15:33:51 by stolfi */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>
#include <r2.h>
#include <jsfile.h>
#include <wt_table.h>
#include <wt_table_binomial.h>
#include <float_image.h>

#include <float_image_mscale.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

float_image_t *float_image_mscale_shrink(float_image_t *A, float_image_t *M, int32_t NXR, int32_t NYR, int32_t dx, int32_t dy, int32_t nw)
  { 
    /* Generate the 1D weight mask: */
    demand(nw > 0, "invalid filter width");
    double wt[nw];
    wt_table_binomial_fill((uint32_t)nw, wt, NULL);
    wt_table_normalize_sum((uint32_t)nw, wt);
    
    /* Get the image dimensions: */
    int32_t NC = (int32_t)A->sz[0];
    int32_t NXA = (int32_t)A->sz[1];
    int32_t NYA = (int32_t)A->sz[2];

    /* Check the mask dimensions: */
    if (M != NULL)
      { demand((int32_t)M->sz[0] == 1, "weight mask should be monochromatic");
        demand((int32_t)M->sz[1] == NXA, "mask width mismatch");
        demand((int32_t)M->sz[2] == NYA, "mask height mismatch");
      }
    
    /* Create the output image {R}: */
    float_image_t *R = float_image_new(NC, NXR, NYR);
    
    /* Fill the pixels of {R}: */
    int32_t xR, yR, c;
    for(yR = 0; yR < NYR; yR++)
      { for(xR = 0; xR < NXR; xR++)
          { /* Accumulate the weighted pixel sum over the window: */
            double sum_w = 0;    /* Sum of weights. */
            double sum_w_v[NC];  /* Sum of weighted pixels. */
            for (c = 0; c < NC; c++) { sum_w_v[c] = 0; }
            int32_t xD, yD;
            for(yD = 0; yD < nw; yD++)
              { int32_t yA = 2*yR+yD-dy;
                double wty = ((yA < 0) || (yA >= NYA) ? 0.0 : wt[yD]);
                for(xD = 0; xD < nw; xD++)
                  { int32_t xA = 2*xR+xD-dx;
                    double wtx = ((xA < 0) || (xA >= NXA) ? 0.0 : wt[xD]);
                    double w = wtx*wty;
                    if (w > 0)
                      { /* Multiply by the mask weight, if any: */
                        if (M != NULL) { w *= float_image_get_sample(M, 0, xA, yA); }
                        sum_w += w; 
                        for (c = 0; c < NC; c++)
                          { double v = float_image_get_sample(A, c, xA, yA);
                            sum_w_v[c] += w*v;
                          }
                      }
                  }
              }
            /* Store the weighted average (maybe NAN) in the {R} pixel: */
            if (sum_w == 0) { sum_w = 1; }
            for (c = 0; c < NC; c++)
              { double v = sum_w_v[c]/sum_w; 
                float_image_set_sample(R, c, xR, yR, (float)v);
              }
          }
      }
    return R;
  }
  
float_image_t *float_image_mscale_mask_shrink(float_image_t *M, int32_t NXR, int32_t NYR, int32_t dx, int32_t dy, int32_t nw, bool_t harm)
  { /* Generate the 1D weight mask: */
    demand(nw > 0, "invalid filter width");
    double wt[nw];
    wt_table_binomial_fill((uint32_t)nw, wt, NULL);
    wt_table_normalize_sum((uint32_t)nw, wt);
    
    /* Get the image dimensions: */
    int32_t NC = (int32_t)M->sz[0];
    int32_t NXM = (int32_t)M->sz[1];
    int32_t NYM = (int32_t)M->sz[2];
        
    /* Create the output image {R}: */
    float_image_t *R = float_image_new(NC, NXR, NYR);
    
    /* Fill the pixels of {R}: */
    int32_t xR, yR, c;
    for(yR = 0; yR < NYR; yR++)
      { for(xR = 0; xR < NXR; xR++)
          { for (c = 0; c < NC; c++)
              { /* Accumulate the sample sums over the window: */
                double sum_w = 0;    /* Sum of weights. */
                double sum_w_m = 0;  /* Sum of weighted pixels. */
                int32_t xD, yD;
                for(yD = 0; yD < nw; yD++)
                  { int32_t yM = 2*yR+yD-dy;
                    double wty = ((yM < 0) || (yM >= NYM) ? 0.0 : wt[yD]);
                    for(xD = 0; xD < nw; xD++)
                      { int32_t xM = 2*xR+xD-dx;
                        double wtx = ((xM < 0) || (xM >= NXM) ? 0.0 : wt[xD]);
                        double w = wtx*wty;
                        if (w > 0)
                          { sum_w += w;
                            double m = float_image_get_sample(M, c, xM, yM);
                            /* !!! Must find a probabilistic justification for this: */
                            sum_w_m += (harm ? w/m : w*m);
                          }
                      }
                  }
                double m = (harm ? sum_w/sum_w_m : sum_w_m/sum_w);
                assert(! isnan(m));
                float_image_set_sample(R, c, xR, yR, (float)m);
              }
          }
      }
    return R;
  }

r2_t float_image_mscale_point_shrink(r2_t *pA, int32_t dx, int32_t dy, int32_t nw)
  {
    double hw = 0.5*(nw - 1);
    return (r2_t){{ 0.5*(pA->c[0] + dx - hw), 0.5*(pA->c[1] + dy - hw) }};
  }

r2_t float_image_mscale_point_expand(r2_t *pR, int32_t dx, int32_t dy, int32_t nw)
  {
    double hw = 0.5*(nw - 1);
    return (r2_t){{ 2.0*pR->c[0] + hw - dx , 2.0*pR->c[1] + hw - dy }};
  }

int32_t float_image_mscale_rounding_bias(int32_t n)
  {
    int32_t b = 0;
    while (n & 1) { b = 1-b; n /= 2; }
    return b;
  }
  
char *float_image_mscale_file_name(char *filePrefix, int32_t level, int32_t iter, char *tag, char *ext)
  { 
    char *fileName = NULL;
    if (strcmp(filePrefix,"-") == 0)
      { fileName = jsprintf("-"); }
    else
      { /* Typeset {level} and {iter}; leave "" to omit. */
        char *xlevel = (level >= 0 ? jsprintf("-%02d", level) : "");
        char *xiter = (iter >= 0 ? jsprintf("-%09d", iter) : "");
        /* Typeset the filename: */
        if (tag == NULL) { tag = ""; }
        char *tagsep = (tag[0] == 0 ? "" : "-");
        fileName = jsprintf("%s%s%s%s%s.%s", filePrefix, xlevel, xiter, tagsep, tag, ext);
        if (level >= 0) { free(xlevel); }
        if (iter >= 0) { free(xiter); }
      }
    return fileName;
  }

void float_image_mscale_write_file(float_image_t *M, char *filePrefix, int32_t level, int32_t iter, char *tag)
  { char *fileName = float_image_mscale_file_name(filePrefix, level, iter, tag, "fni");
    int32_t indent = (level < 0 ? 0 : 2*level);
    fprintf(stderr, "%*sWriting %s ...", indent, "", fileName);
    FILE* wr = open_write(fileName, FALSE);
    float_image_write(wr, M);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
    free(fileName);
  }

