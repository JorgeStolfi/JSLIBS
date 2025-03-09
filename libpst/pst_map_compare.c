/* See pst_map_compare.h */
/* Last edited on 2025-03-04 18:32:40 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_mscale.h>
#include <pst_map_shift_values.h>

#include <pst_map_compare.h>

void pst_compare_check_sizes
  ( float_image_t *A,
    int32_t *NC_A,
    float_image_t *B,
    int32_t *NC_B,
    int32_t NC,
    int32_t *NX,
    int32_t *NY
  )
  { demand(A != NULL, "image {A} must be not null");
    demand(B != NULL, "image {B} must be not null");
    
    float_image_get_size(A, NC_A, NX, NY); 
    
    int32_t NX_B, NY_B;
    float_image_get_size(B, NC_B, &NX_B, &NY_B);
    demand(((*NC_A) == NC) || ((*NC_A) == NC-1), "map {A} has wrong channel count");
    
    demand((NX_B == (*NX)) && (NY_B == (*NY)), "map size mismatch");
    demand(((*NC_B) == NC) || ((*NC_B) == NC-1), "map {B} has wrong channel count");
  }

float_image_t *pst_map_compare
  ( float_image_t *A,
    float_image_t *B,
    int32_t NC,
    double avgA[],
    double devA[],
    double avgB[],
    double devB[],
    uint32_t undef[],
    double avgE[],
    double devE[],
    double devRelE[]
  )
  { 
    int32_t wch = NC-1; /* Index of weight channel. */
    
    int32_t NC_A, NC_B, NX, NY;
    pst_compare_check_sizes(A, &NC_A, B, &NC_B, NC, &NX, &NY);
    assert((NC_A == NC) || (NC_A == NC-1));
    assert((NC_B == NC) || (NC_B == NC-1));

    float_image_t *E = float_image_new(NC, NX, NY);
    
    /* Compute the mean values of {A,B,E}: */
    double sum_AW[NC];
    double sum_BW[NC];
    double sum_EW[NC];
    double sum_W[NC];
    for (int32_t c = 0; c < NC; c++)
      { sum_AW[c] = sum_BW[c] = sum_EW[c] = sum_W[c] = 0.0; 
        undef[c] = (uint32_t)(NX*NY); /* To be decremented. */
      }
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { /* Fetch the pixels {vA,vB} from {A,B}, with default weight 1.0: */
            float vA[NC], vB[NC];
            float_image_get_pixel(A, x, y, vA);
            if (NC_A < NC) { assert(wch == NC_A); vA[wch] = 1.0; }
            float_image_get_pixel(B, x, y, vB);
            if (NC_B < NC) { assert(wch == NC_B); vB[wch] = 1.0; }
            
            /* Get weight {w} of those pixels: */
            double wA = vA[wch];
            demand(isfinite(wA) && wA >= 0, "invalid weight value in {A}");
            double wB = vB[wch];
            demand(isfinite(wB) && wB >= 0, "invalid weight value in {B}");
            
            /* Get weight {wE} of difference: */
            double wE = 1.0/(1.0/wA + 1.0/wB);
            if (! isfinite(wE)) { wE = 0.0; }

            /* Compute pixel {vE} of {E}: */
            float vE[NC];
            for (int32_t c = 0; c < NC; c++) 
              { if (c == wch)
                  { vE[c] = (float)wE; }
                else 
                  { assert ((c < NC_A) && (c < NC_B));
                    if (isfinite(vA[c]) && isfinite(vB[c]) && (wE > 0))
                      { vE[c] = vA[c] - vB[c];
                        if (! isfinite(vE[c])) { vE[c] = NAN; }
                      }
                    else 
                      { vE[c] = NAN; }
                  }
              }
            if (wE > 0)
              { /* Accumulate weighted sums of samples: */
                for (int32_t c = 0; c < NC; c++) 
                  { if (isfinite(vE[c]))
                      { /* Accumulate values: */
                        if (c < NC_A) { sum_AW[c] += vA[c]*wE; }
                        if (c < NC_B) { sum_BW[c] += vB[c]*wE; }
                        sum_EW[c] += vE[c]*wE;
                        sum_W[c] += wE;
                        assert(undef[c] > 0);
                        undef[c]--;
                      }
                  }
              }
            /* Store difference in error image: */
            float_image_set_pixel(E, x, y, vE);
          }
      }

    /* Compute sample means: */
    for (int32_t c = 0; c < NC; c++)
      { avgA[c] = (sum_W[c] == 0 ? 0.5 : sum_AW[c]/sum_W[c]);
        avgB[c] = (sum_W[c] == 0 ? 0.5 : sum_BW[c]/sum_W[c]);
        avgE[c] = (sum_W[c] == 0 ? 0.0 : sum_EW[c]/sum_W[c]);
      }
      
    /* Compute sample deviations: */
    double sum_dA2W[NC];
    double sum_dB2W[NC];
    double sum_dE2W[NC];
    for (int32_t c = 0; c < NC; c++)
      { sum_dA2W[c] = sum_dB2W[c] = sum_dE2W[c] = 0.0; }
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float vE[NC];
            float_image_get_pixel(E, x, y, vE);
            /* Get error sample and weight: */
            double wE = ((wch >= 0) && (wch < NC) ? vE[wch] : 1.0);
            assert(isfinite(wE) && (wE >= 0));
            if (wE > 0)
              { /* Fetch the pixels {vA,vB} from {A,B}: */
                float vA[NC], vB[NC];
                float_image_get_pixel(A, x, y, vA);
                float_image_get_pixel(B, x, y, vB);
                for (int32_t c = 0; c < NC; c++) 
                  { if (isfinite(vA[c]) && isfinite(vB[c]))
                      { /* Accumulate squares: */
                        double dA = vA[c] - avgA[c]; sum_dA2W[c] += dA*dA*wE;
                        double dB = vB[c] - avgB[c]; sum_dB2W[c] += dB*dB*wE;
                        double dE = vE[c] - avgE[c]; sum_dE2W[c] += dE*dE*wE;
                      }
                  }
              }
          }
      }
    /* Compute the deviations: */
    for (int32_t c = 0; c < NC; c++)
      { devA[c] = (sum_W[c] == 0 ? 0.0 : sqrt(sum_dA2W[c]/sum_W[c]));
        devB[c] = (sum_W[c] == 0 ? 0.0 : sqrt(sum_dB2W[c]/sum_W[c]));
        devE[c] = (sum_W[c] == 0 ? 0.0 : sqrt(sum_dE2W[c]/sum_W[c]));
        double devAB = hypot(devA[c],devB[c])/M_SQRT2;
        devRelE[c] = ((! isfinite(devAB)) || (devAB == 0) ? 0.0 : devE[c]/devAB);
      }

    return E;
  }

void pst_map_compare_analyze_and_write
  ( float_image_t *A,
    float_image_t *B, 
    int32_t NC,
    double change[],
    char *filePrefix,
    int32_t level,
    int32_t iter,
    char *tagA,
    char *tagB,
    char *tagE,
    bool_t verbose
  )
  { bool_t debug = FALSE;
    int32_t indent = (level <= -1 ? 0 : 2*level+2);
    if (debug) { fprintf(stderr, "%*sentering {%s}\n", indent, "", __FUNCTION__); }
    
    demand(A != NULL, "image {A} must be not null");

    int32_t NX, NY;
    float_image_get_size(A, NULL, &NX, &NY);
    int32_t wch = NC-1; /* Index of weight channel. */
      
    if (tagA != NULL)
      { if (verbose) { fprintf(stderr, "%*swriting the map {%s} ...\n", indent, "", tagA); }
        float_image_mscale_write_file(A, filePrefix, level, iter, tagA);
      }
      
    if (B != NULL) 
      { if (tagB != NULL)
          { if (verbose) { fprintf(stderr, "%*swriting the map {%s} ...\n", indent, "", tagB); }
            float_image_mscale_write_file(B, filePrefix, level, iter, tagB);
          }
        double avgA[NC], devA[NC], avgB[NC], devB[NC], avgE[NC], devE[NC], devRelE[NC];
        uint32_t undef[NC];
        float_image_t *E = pst_map_compare(A, B, NC, avgA, devA, avgB, devB, undef, avgE, devE, devRelE);
        
        if (verbose) 
          { fprintf(stderr, "%*siteration %5d: change per channel", indent, "", iter);
            for (int32_t c = 0; c < NC; c++) { fprintf(stderr, " [%d] = %14.11f", c, change[c]); }
            fprintf(stderr, "\n");
            if (debug)
              { for (int32_t c = 0; c < NC; c++)
                  { fprintf(stderr, "%*s  A[%d] = %14.11f ± %14.11f\n", indent, "", c, avgA[c], devA[c]);
                    fprintf(stderr, "%*s  B[%d] = %14.11f ± %14.11f\n", indent, "", c, avgB[c], devB[c]);
                    fprintf(stderr, "%*s  E[%d] = %14.11f ± %14.11f", indent, "", c, avgE[c], devE[c]);
                    fprintf(stderr, " (rel %14.11f)  undef = %6d\n", devRelE[c], undef[c]); 
                  }
              }
          }

        if (tagE != NULL)
          { /* Write the error map {E}: */
            if (debug) { fprintf(stderr, "%*sshifting the error map {E} to zero mean ...\n", indent, ""); }
            pst_map_shift_values(E, wch, -1, avgE);
            if (verbose) { fprintf(stderr, "%*swriting the height error map {E} ...\n", indent, ""); }
            float_image_mscale_write_file(E, filePrefix, level, iter, tagE);
            
            /* Open the stats file: */
            char *fileName = float_image_mscale_file_name(filePrefix, level, iter, tagE, "txt");
            if (verbose) { fprintf(stderr, "%*swriting %s ...\n", indent, "", fileName); }
            FILE* wr = fopen(fileName, "wt");
            
            free(fileName);
            fprintf(wr, "%2d %6d %6d %9d", level, NX, NY, iter);
            fprintf(wr, " %3d", NC);
            for (int32_t c = 0; c < NC; c++)
              { fprintf(wr, "  %d", c);
                fprintf(wr, "  %14.11f", change[c]);
                fprintf(wr, "  %14.11f %14.11f", avgA[c], devA[c]);
                fprintf(wr, "  %14.11f %14.11f", avgB[c], devB[c]);
                fprintf(wr, "  %14.11f %14.11f %14.11f", avgE[c], devE[c], devRelE[c]);
              }
            if (wr == stdout) { fflush(wr); } else { fclose(wr); }
          }
      }
    if (debug) { fprintf(stderr, "%*sexiting {pst_map_compare_level_analyze_and_write}\n", indent, ""); }
  }
