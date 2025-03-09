/* See {float_image_mscale.h}. */
/* Last edited on 2025-03-02 01:28:40 by stolfi */

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

float_image_t *float_image_mscale_shrink
  ( float_image_t *A,
    int32_t wch,
    bool_t harm,
    int32_t NXR,
    int32_t NYR,
    int32_t dx,
    int32_t dy
  )
  { 
    /* Get the input image dimensions: */
    int32_t NC = (int32_t)A->sz[0];
    int32_t NXA = (int32_t)A->sz[1];
    int32_t NYA = (int32_t)A->sz[2];
    
    /* Create the output image {R}: */
    float_image_t *R = float_image_new(NC, NXR, NYR);
    
    /* Fill the pixels of {R}: */
    for (int32_t yR = 0; yR < NYR; yR++)
      { for (int32_t xR = 0; xR < NXR; xR++)
          { /* Accumulate the weighted pixel sum over the window: */
            /* For {c!=wch}, the value of {sum_w_vA[c]} is the sum of the sample {A[c,X',Y']}
              times the weight of pixel {A[X',Y']} times {A[c.X',Y']},
              for all {X',Y'} in the {2×2} block that maps to {R[X,Y]}.
              If {wch} is a valid channel index, {sum_w_vA[wch]}
              is the simple (unweighted) sum of the weights
              of those four pixels, if {harm} is false, or of 
              their inverses, if {harm} is true. */
            
            double sum_wA = 0;    /* Sum of pixel weights in block. */
            double sum_wA_vA[NC]; /* Sum of weighted pixel values in block. */
            for (int32_t c = 0; c < NC; c++) { sum_wA_vA[c] = 0; }
            for (int32_t s = 0; s < 2; s++)
              { int32_t yA = 2*yR - dy + s;
                for (int32_t r = 0; r < 2; r++)
                  { int32_t xA = 2*xR - dx + r;
                    if ((xA >= 0) && (xA < NXA) && (yA >= 0) && (yA < NYA))
                      { /* Get the reliability weight, if any: */
                        double wA = ((wch >= 0) && (wch < NC) ? float_image_get_sample(A, wch, xA, yA) : 1.0); 
                        demand(isfinite(wA) && (wA >= 0), "invalid weight value");
                        if (wA > 0)
                          { sum_wA += wA; 
                            float vA[NC];
                            float_image_get_pixel(A, xA, yA, vA);
                            for (int32_t c = 0; c < NC; c++)
                              { if (c == wch)
                                  { sum_wA_vA[c] += (harm ? 1/vA[c] : vA[c]); }
                                else
                                  { sum_wA_vA[c] += wA*vA[c]; }
                              }
                          }
                        else
                          { if (harm) { sum_wA_vA[wch] = +INF; } }
                      }
                  }
              }
            /* Store the weighted average (maybe NAN) in the {R} pixel: */
            if (sum_wA == 0) { sum_wA = 1; }
            float vR[NC];
            bool_t bad = FALSE;
            for (int32_t c = 0; c < NC; c++)
              { double vRc;
                if (c == wch)
                  { vRc = sum_wA_vA[c]/4.0;
                    if (harm) { vRc = 1.0/vRc; }
                  }
                else
                  { vRc = sum_wA_vA[c]/sum_wA; }
                vR[c] = (float)vRc;
                if (! isfinite(vR[c])) { bad = TRUE; }
              }
            if (bad) { for (int32_t c = 0; c < NC; c++) { vR[c] = (c == wch ? 0.0 : NAN); }}
            float_image_set_pixel(R, xR, yR, vR);
          }
      }
    return R;
  }

float_image_t *float_image_mscale_expand
  ( float_image_t *R,
    int32_t NXA,
    int32_t NYA,
    int32_t dx,
    int32_t dy
  )
  {
    /* Get the input image dimensions: */
    int32_t NC = (int32_t)R->sz[0];
    int32_t NXR = (int32_t)R->sz[1];
    int32_t NYR = (int32_t)R->sz[2];
    
    /* Create the output image {R}: */
    float_image_t *A = float_image_new(NC, NXA, NYA);
    float_image_fill(A, 0.0); /* In case there are undefined margins. */
   
    /* Expand pixels of {R} into {2×2} blocks of {A}. */
    float vR[NC];
    for(int32_t yR = 0; yR < NYR; yR++)
      { for(int32_t xR = 0; xR < NXR; xR++)
          { float_image_get_pixel(R, xR, yR, vR);
            for(int32_t s = 0; s < 2; s++)
              { int32_t yA = 2*yR - dy + s;
                for(int32_t r = 0; r < 2; r++)
                  { int32_t xA = 2*xR - dx + r;
                    if ((xA >= 0) && (xA < NXA) && (yA >= 0) && (yA < NYA))
                      { float_image_set_pixel(A, xA, yA, vR); }
                  }
              }
          }
      }
    return A;
  }

r2_t float_image_mscale_point_shrink(r2_t *pA, int32_t dx, int32_t dy)
  {
    return (r2_t){{ 0.5*(pA->c[0] + dx - 0.5), 0.5*(pA->c[1] + dy - 0.5) }};
  }

r2_t float_image_mscale_point_expand(r2_t *pR, int32_t dx, int32_t dy)
  {
    return (r2_t){{ 2.0*pR->c[0] + 0.5 - dx , 2.0*pR->c[1] + 0.5 - dy }};
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
        char *xiter = (iter == 0 ? "-beg" : (iter < 0 ? "-end" : jsprintf("-%09d", iter)));
        /* Typeset the filename: */
        if (tag == NULL) { tag = ""; }
        char *tagsep = (tag[0] == 0 ? "" : "-");
        fileName = jsprintf("%s%s%s%s%s.%s", filePrefix, xlevel, xiter, tagsep, tag, ext);
        if (level >= 0) { free(xlevel); }
        if (iter > 0) { free(xiter); }
      }
    return fileName;
  }

void float_image_mscale_write_file(float_image_t *M, char *filePrefix, int32_t level, int32_t iter, char *tag)
  { char *fileName = float_image_mscale_file_name(filePrefix, level, iter, tag, "fni");
    int32_t indent = (level < 0 ? 0 : 2*level+2);
    fprintf(stderr, "%*swriting %s ...", indent, "", fileName);
    FILE* wr = open_write(fileName, FALSE);
    float_image_write(wr, M);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
    free(fileName);
  }

