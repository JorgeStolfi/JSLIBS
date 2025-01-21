/* See pst_height_map.h */
/* Last edited on 2025-01-19 22:28:49 by stolfi */

#include <math.h>
#include <affirm.h>
#include <assert.h>
#include <stdint.h>

#include <pst_height_map.h>
#include <pst_interpolate.h>
#include <float_image.h>
#include <float_image_mscale.h>

float_image_t *pst_height_map_shrink(float_image_t *IZ)
  {
    int32_t NCI, NXI, NYI;
    float_image_get_size(IZ, &NCI, &NXI, &NYI);
    demand((NCI == 1) || (NCI == 2), "height map must have 1 or 2 channels");

    int32_t NXJ = NXI/2+1;
    int32_t NYJ = NYI/2+1;
    float_image_t *JZ = float_image_new(NCI, (int32_t)NXJ, (int32_t)NYJ);
    
    for (int32_t yJ = 0; yJ < NYJ; yJ++)
      { for (int32_t xJ = 0; xJ < NXJ; xJ++)
          { double sum_wz = 0;
            double sum_w = 0; 
            float min_w = +INF;
            for (int32_t dx = 0; dx <= 1; dx++)
              { for (int32_t dy = 0; dy <= 1; dy++)
                  { int32_t xI = 2*xJ + dx;
                    int32_t yI = 2*yJ + dy;
                    if ((xI < NXI) && (yI < NYI))
                      { float v[NCI];
                        float_image_get_pixel(IZ, xI, yI, v);
                        double w = (NCI == 2 ? v[1] : 1.0);
                        sum_wz += w*v[0];
                        sum_w += w;
                        min_w = fminf(min_w, v[1]);
                      }
                  }
              }
            float vJ[2];
            vJ[0] = (float)(sum_w == 0 ? 0.0 : 0.5*sum_wz/sum_w);
            if (NCI == 2) { vJ[1] = min_w; }
            float_image_set_pixel(JZ, xJ, yJ, vJ);
          }
      }
    return JZ;
  }

float_image_t *pst_height_map_expand(float_image_t *JZ, int32_t NX, int32_t NY)
  { int32_t NCJ, NXJ, NYJ;
    float_image_get_size(JZ, &NCJ, &NXJ, &NYJ);
    demand((NCJ == 1) || (NCJ == 2), "input height map must have 1 or 2 channels");
    
    int32_t NXI = NX;
    int32_t NYI = NY;
    demand(NXJ == NXI/2 + 1, "inconsistent {X} sizes");
    demand(NYJ == NYI/2 + 1, "inconsistent {Y} sizes");
    float_image_t *IZ = float_image_new(NCJ, (int32_t)NXI, (int32_t)NYI);
    
    for (int32_t yJ = 0; yJ < NYJ; yJ++)
      { for (int32_t xJ = 0; xJ < NXJ; xJ++)
          { float vJ[2];
            float_image_get_pixel(JZ, xJ, yJ, vJ);
            vJ[0] = 2*vJ[0];
            for (int32_t dx = 0; dx <= 1; dx++)
              { for (int32_t dy = 0; dy <= 1; dy++)
                  { int32_t xI = 2*xJ - dx;
                    int32_t yI = 2*yJ - dy;
                    if ((xI >= 0) && (xI < NXI) && (yI >= 0) && (yI < NYI))
                      { float_image_set_pixel(IZ, xI, yI, vJ); }
                  }
              }
          }
      }
    return IZ;
  }

float_image_t *pst_height_map_compare
  ( float_image_t *AZ,
    float_image_t *BZ,
    bool_t zero_mean,
    double *sAZP,
    double *sBZP,
    double *sEZP,
    double *sreP
  )
  { int32_t NCA, NX, NY;
    float_image_get_size(AZ, &NCA, &NX, &NY);
    demand((NCA == 1) || (NCA == 2), "input height map {AZ} must have 1 or 2 channels");
    float_image_check_size(BZ, -1, NX, NY, "height map size mismatch");
    int32_t NCB = (int32_t)BZ->sz[0];
    demand((NCB == 1) || (NCB == 2), "input height map {BZ} must have 1 or 2 channels");
    
       
    /* Compute the mean values of {AZ,BZ,EZ}: */
    double sum_AZW = 0;
    double sum_BZW = 0;
    double sum_EZW = 0;
    double sum_W = 0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { /* Get relevant samples from original image: */
            double vA = float_image_get_sample(AZ, 0, x, y);
            double wA = (NCA == 2 ? float_image_get_sample(AZ, 1, x, y) : 1.0);
            double vB = float_image_get_sample(BZ, 0, x, y);
            double wB = (NCB == 2 ? float_image_get_sample(BZ, 1, x, y) : 1.0);
            double w = wA*wB;
            double vE = vA - vB;
            /* Accumulate values: */
            sum_AZW += vA*w;
            sum_BZW += vB*w;
            sum_EZW += vE*w;
            sum_W += w;
          }
      }
    double avgA = (sum_W == 0 ? 0.5 : sum_AZW/sum_W);
    double avgB = (sum_W == 0 ? 0.5 : sum_BZW/sum_W);
    double avgE = (sum_W == 0 ? 0.0 : sum_EZW/sum_W);
      
    /* Fill {EZ} and RMS values: */
    int32_t NCE = 2;
    float_image_t *EZ = float_image_new(NCE, NX, NY);
    double sum_AdZ2W = 0.0;
    double sum_BdZ2W = 0.0;
    double sum_EZ2W = 0.0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { /* Get relevant samples from original image: */
            double vA = float_image_get_sample(AZ, 0, x, y);
            double wA = (NCA == 2 ? float_image_get_sample(AZ, 1, x, y) : 1.0);
            double vB = float_image_get_sample(BZ, 0, x, y);
            double wB = (NCB == 2 ? float_image_get_sample(BZ, 1, x, y) : 1.0);
            double w = wA*wB;
            double vE = vA - vB;
            if (zero_mean) { vE = vE - avgE; }
            /* Store difference in error image: */
            float_image_set_sample(EZ, 0, x, y, (float)vE);
            float_image_set_sample(EZ, 1, x, y, (float)w);
            /* Accumuate squares: */
            vA = vA - avgA;
            vB = vB - avgB;
            sum_AdZ2W += vA*vA*w;
            sum_BdZ2W += vB*vB*w;
            sum_EZ2W += vE*vE*w;
          }
      }
    /* Compute the RMS values and errors: */
    double sAZ = (sum_W == 0 ? 0.0 : sqrt(sum_AdZ2W/sum_W));
    double sBZ = (sum_W == 0 ? 0.0 : sqrt(sum_BdZ2W/sum_W));
    double sEZ = (sum_W == 0 ? 0.0 : sqrt(sum_EZ2W/sum_W));
    double sMZ = hypot(sAZ,sBZ)/M_SQRT2;
    double sre = sEZ/sMZ;
    
    /* Return results: */
    if (sAZP != NULL) { (*sAZP) = sAZ; }
    if (sBZP != NULL) { (*sBZP) = sBZ; }
    if (sEZP != NULL) { (*sEZP) = sEZ; }
    if (sreP != NULL) { (*sreP) = sre; }

    return EZ;
  }

void pst_height_map_level_analyze_and_write
  ( char *filePrefix,
    int32_t level,
    int32_t iter,
    double change,
    float_image_t *CZ,
    float_image_t *RZ, 
    bool_t writeImages,
    bool_t writeError
  )
  {
    int32_t NCC, NX, NY;
    float_image_get_size(CZ, &NCC, &NX, &NY);
    demand((NCC == 1) || (NCC == 2), "input height map must have 1 or 2 channels");
    if (RZ != NULL) { float_image_check_size(RZ, -1, NX, NY, "height map size mismatch"); }
    int32_t NCR = (int32_t)RZ->sz[0];
    demand((NCR == 1) || (NCR == 2), "reference height map must have 1 or 2 channels");
    
    int32_t indent = (level <= -1 ? 0 : 2*level+2);
      
    if (writeImages)
      { float_image_mscale_write_file(CZ, filePrefix, level, iter, "Z"); }
      
    if (RZ != NULL)
      { double sAZ, sBZ, sEZ, sre;
        float_image_t *EZ = pst_height_map_compare(CZ, RZ, TRUE, &sAZ, &sBZ, &sEZ, &sre);
        if (writeImages)
          { float_image_mscale_write_file(EZ, filePrefix, level, iter, "eZ"); }
        if (writeError)
          { char *fileName = float_image_mscale_file_name(filePrefix, level, iter, "eZ", "txt");
            fprintf(stderr, "%*swriting %s ...", indent, "", fileName);
            FILE* wr = fopen(fileName, "wt");
            assert(wr != NULL); 
            fprintf
              ( wr, "%2d %6d %6d %9d %14.11f  %14.11f %14.11f  %14.11f %14.11f\n",
                level, NX, NY, iter, change, sAZ, sBZ, sEZ, sre
              );
            if (wr == stdout) { fflush(wr); } else { fclose(wr); }
            fprintf(stderr, "\n");
            free(fileName);
          }
      }
  }

void pst_height_map_shift_to_zero_mean(float_image_t *Z)
  { 
    int32_t NC, NX, NY;
    float_image_get_size(Z, &NC, &NX, &NY);
    demand((NC == 1) || (NC == 2), "height map must have 1 or 2 channels");
    
    /* Compute the mean {Z} value: */
    double sum_wz = 0;
    double sum_w = 0;
    for (int32_t y = 0; y <= NY; y++)
      { for (int32_t x = 0; x <= NX; x++)
          { double z = float_image_get_sample(Z, 0, x, y);
            double w = (NC == 2 ? float_image_get_sample(Z, 1, x, y) : 1.0);
            if (! isfinite(z)) 
              { w = 0.0;
                if (NC == 2) { float_image_set_sample(Z, 1, x, y, 0.0); }
              }
            sum_wz += w*z;
            sum_w += w;
          }
      }
    double mean_z = (sum_w == 0 ? 0.0 : sum_wz/sum_w);
    for (int32_t y = 0; y <= NY; y++)
      { for (int32_t x = 0; x <= NX; x++)
          { float *smp = float_image_get_sample_address(Z, 0, x, y);
            (*smp) = (float)((*smp) - mean_z); 
          }
      }
  }    
 
