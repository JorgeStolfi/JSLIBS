/* See pst_height_map.h */
/* Last edited on 2024-12-23 06:00:41 by stolfi */

#include <math.h>
#include <affirm.h>
#include <assert.h>
#include <stdint.h>

#include <pst_height_map.h>
#include <pst_interpolate.h>
#include <float_image.h>
#include <float_image_mscale.h>

void pst_height_map_expand
  ( float_image_t *JZ,
    float_image_t *JW,
    float_image_t *IZ,
    float_image_t *IW
  )
  { 
    assert(JZ->sz[0] == 1);
    uint32_t NXJ = (uint32_t)JZ->sz[1];
    uint32_t NYJ = (uint32_t)JZ->sz[2];
    
    if (JW != NULL) { assert(JW->sz[0] == 1); assert(JW->sz[1] == NXJ); assert(JW->sz[2] == NYJ); }
    
    assert(IZ->sz[0] == 1);
    uint32_t NXI = (uint32_t)IZ->sz[1]; assert(NXJ == NXI/2 + 1);
    uint32_t NYI = (uint32_t)IZ->sz[2]; assert(NYJ == NYI/2 + 1);
    
    if (IW != NULL) { assert(IW->sz[0] == 1); assert(IW->sz[1] == NXI); assert(IW->sz[2] == NYI); }
    
    for (int32_t yI = 0; yI < NYI; yI++)
      { for (int32_t xI = 0; xI < NXI; xI++)
          { 
	    int32_t xJ = xI/2;
	    int32_t yJ = yI/2;
	    double v,w;
	    
            if ((xI%2) == 0)
              { if ( (yI%2) == 0 )
                  { v = float_image_get_sample (JZ, 0, xJ,yJ);
                    w = (JW == NULL ? 1.0 : float_image_get_sample(JZ,0,xJ,yJ));
                  }
                else
                  { pst_interpolate_four_samples(JZ,JW, 0, xJ,yJ, xJ,yJ+1, &v,&w); }
              }
	    else
              { if ( (yI%2) == 0 )
                  { pst_interpolate_four_samples(JZ,JW, 0, xJ,yJ, xJ+1,yJ, &v,&w); }
                else
                  { double wa,wb;
                    double va,vb;
                    pst_interpolate_four_samples(JZ,JW, 0, xJ,yJ,   xJ+1,yJ+1, &va,&wa);
                    pst_interpolate_four_samples(JZ,JW, 0, xJ,yJ+1, xJ+1,yJ,   &vb,&wb);
                    if ((wa == 0) && (wb == 0))
                      { w = 0; v = (va+vb)/2.0; }
                    else
                      { w = wa+wb; v = (wa*va + wb*vb)/w; }
                  }
              }
            /* Dont forget to multiply the interpolation by two ! */
            float_image_set_sample(IZ, 0, xI, yI, (float)(2*v));
            if (IW != NULL) { float_image_set_sample(IW, 0, xI, yI, (float)w); }
         }
      }
  }

float_image_t *pst_height_map_shrink(float_image_t *IZ, uint32_t avgWidth)
  {
    int32_t NX_JZ = (int32_t)IZ->sz[1]/2+1;
    int32_t NY_JZ = (int32_t)IZ->sz[2]/2+1;
    int32_t dxy = (int32_t)(avgWidth-1)/2;
    float_image_t *JZ = float_image_mscale_shrink(IZ, NULL, NX_JZ, NY_JZ, dxy, dxy, (int32_t)avgWidth);
    float_image_rescale_samples(JZ, 0, 0.0, 1.0, 0.0, 0.5);
    return JZ;
  }

float_image_t *pst_height_map_compare
  ( float_image_t *AZ,
    float_image_t *BZ,
    float_image_t *W,
    bool_t zero_mean,
    double *sAZP,
    double *sBZP,
    double *sEZP,
    double *sreP
  )
  { 
    assert(AZ->sz[0] == 1);
    assert(BZ->sz[0] == 1);
    int32_t NX = (int32_t)AZ->sz[1]; assert(BZ->sz[1] == NX);
    int32_t NY = (int32_t)AZ->sz[2]; assert(BZ->sz[2] == NY);
    
    if (W != NULL)
      { assert(W->sz[0] == 1);
        assert(W->sz[1] == NX);
        assert(W->sz[2] == NY);
      }
      
    /* Compute the mean values of {AZ,BZ,EZ}: */
    double sum_AZW = 0;
    double sum_BZW = 0;
    double sum_EZW = 0;
    double sum_W = 0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { /* Get relevant samples from original image: */
            double vA = float_image_get_sample(AZ, 0, x, y);
            double vB = float_image_get_sample(BZ, 0, x, y);
            double vW = float_image_get_sample(W, 0, x, y);
            double vE = vA - vB;
            /* Accumulate values: */
            sum_AZW += vA*vW;
            sum_BZW += vB*vW;
            sum_EZW += vE*vW;
            sum_W += vW;
          }
      }
    double avgA = (sum_W == 0 ? 0.5 : sum_AZW/sum_W);
    double avgB = (sum_W == 0 ? 0.5 : sum_BZW/sum_W);
    double avgE = (sum_W == 0 ? 0.0 : sum_EZW/sum_W);
      
    /* Fill {EZ} and RMS values: */
    float_image_t *EZ = float_image_new(1, NX, NY);
    double sum_AdZ2W = 0.0;
    double sum_BdZ2W = 0.0;
    double sum_EZ2W = 0.0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { /* Get relevant samples from original image: */
            double vA = float_image_get_sample(AZ, 0, x, y);
            double vB = float_image_get_sample(BZ, 0, x, y);
            double vW = float_image_get_sample(W, 0, x, y);
            double vE = vA - vB;
            if (zero_mean) { vE = vE - avgE; }
            /* Store difference in error image: */
            float_image_set_sample(EZ, 0, x, y, (float)vE);
            /* Accumuate squares: */
            vA = vA - avgA;
            vB = vB - avgB;
            sum_AdZ2W += vW*vA*vA;
            sum_BdZ2W += vW*vB*vB;
            sum_EZ2W += vW*vE*vE;
          }
      }
    /* Compute the RMS values and errors: */
    double sAZ = sqrt(sum_AdZ2W/sum_W);
    double sBZ = sqrt(sum_BdZ2W/sum_W);
    double sEZ = sqrt(sum_EZ2W/sum_W);
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
    uint32_t level,
    bool_t levelTag,
    uint32_t iter,
    bool_t iterTag,
    double change,
    float_image_t *CZ,
    float_image_t *RZ,
    float_image_t *U, 
    bool_t writeImages,
    bool_t writeError,
    uint32_t indent
  )
  {
    demand(CZ->sz[0] == 1, "bad CZ channels");
    int32_t NX = (int32_t)CZ->sz[1]; 
    int32_t NY = (int32_t)CZ->sz[2]; 
    
    if (RZ != 0)
      { demand(RZ->sz[0] == 1, "bad RZ channels");
        demand(NX == RZ->sz[1], "bad RZ cols"); 
        demand(NY == RZ->sz[2], "bad RZ rows"); 
      }
    
    if (U != 0)
      { demand(U->sz[0] == 1, "bad U channels");
        demand(NX == U->sz[1], "bad U cols"); 
        demand(NY == U->sz[2], "bad U rows"); 
      }
      
    int32_t levelQ = (levelTag ? (int32_t)level : -1);
    int32_t iterQ = (iterTag ? (int32_t)iter : -1);
    
    if (writeImages)
      { float_image_mscale_write_file(CZ, filePrefix, levelQ, iterQ, "Z", (int32_t)indent); }
      
    if (RZ != NULL)
      {
        double sAZ, sBZ, sEZ, sre;
        float_image_t *EZ = pst_height_map_compare(CZ, RZ, U, TRUE, &sAZ, &sBZ, &sEZ, &sre);
        if (writeImages)
          { float_image_mscale_write_file(EZ, filePrefix, levelQ, iterQ, "eZ", (int32_t)indent); }
        if (writeError)
          { char *fileName = float_image_mscale_file_name(filePrefix, levelQ, iterQ, "eZ", "txt");
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
