/* See pst_weight_map.h */
/* Last edited on 2024-12-24 19:13:48 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <float_image.h>
#include <float_image_mscale.h>
#include <wt_table.h>
#include <wt_table_binomial.h>

#include <pst_basic.h>
#include <pst_weight_map.h>
 
/* IMPLEMENTATIONS */

float_image_t *pst_weight_map_shrink(float_image_t *IW, bool_t harmonic, uint32_t avgWidth)
  { 
    demand(IW->sz[0] == 1, "weight map is not mono");
    int32_t NX_JW = (int32_t)((IW->sz[1]+1)/2);
    int32_t NY_JW = (int32_t)((IW->sz[2]+1)/2);
    uint32_t dxy = (avgWidth-1)/2;
    return float_image_mscale_mask_shrink(IW, (int32_t)NX_JW, (int32_t)NY_JW, (int32_t)dxy, (int32_t)dxy, (int32_t) avgWidth, harmonic);
  }
  
  
float_image_t *pst_weight_map_expand_height_weights(float_image_t *IW)
  {
    int32_t NX,NY;
    NX = (int32_t)IW->sz[1];
    NY = (int32_t)IW->sz[2];
    float_image_t* HW = float_image_new(1,(int32_t)NX+1,(int32_t)NY+1);
    for (int32_t x = 0; x < NX+1; x++){
      for (int32_t y = 0; y < NY+1; y++){
        double c00,c01,c10,c11;

        c00 = ( (x  > 0) && (y > 0) ? 1: 0 );
        c10 = ( (x  < NX) && (y > 0) ? 1: 0 );
        c01 = ( (x  > 0) && (y < NY) ? 1: 0 );
        c11 = ( (x  < NX) && (y < NY) ? 1 : 0 );

        double Sc = c00+c01+c10+c11;
        assert(Sc > 0);

        double w00,w01,w10,w11;

        w00 = ( (x  > 0) && (y > 0) ? float_image_get_sample(IW,0,x-1,y-1): 0 );
        w10 = ( (x  < NX) && (y > 0) ? float_image_get_sample(IW,0,x,y-1): 0 );
        w01 = ( (x  > 0) && (y < NY) ? float_image_get_sample(IW,0,x-1,y): 0 );
        w11 = ( (x  < NX) && (y < NY) ? float_image_get_sample(IW,0,x,y): 0 );

    /*    double f00,f01,f10,f11;

        f00 = (  (w00 > 0) ? 1/w00: 0);
        f10 = (  (w10 > 0) ? 1/w10: 0);
        f01 = (  (w01 > 0) ? 1/w01: 0);
        f11 = (  (w11 > 0) ? 1/w11: 0);

        double Sf = f00+f10+f01+f11;
        double Sw = w00+w01+w10+w11;
        double w = 0;
        if (Sf !=  0) {
          w = Sw/Sf;
        }
    */    
        double w = (w00+w01+w10+w11)/(Sc);

        float_image_set_sample(HW,0,x,y,(float)w);
      }
    }

    return HW;
  }

float_image_t *pst_weight_map_slope_to_height(float_image_t *W, bool_t harmonic, int32_t NXV, int32_t NYV)
  { 
    demand(W->sz[0] == 1, "weight map is not mono");
    int32_t NXW = (int32_t)W->sz[1];
    int32_t NYW = (int32_t)W->sz[2];
    
    /* Decide the window width {nw}, either 2 or 3: */
    uint32_t nw;
    if ((NXV == NXW) && (NYV == NYW))
      { /* Same size, average {3x3} blocks: */
        nw = 3;
      }
    else if ((NXV == NXW + 1) && (NYV == NYW + 1))
      { /* New maks is one larger, average {2x2} blocks: */
        nw = 2;
      }
    else
      { demand(FALSE, "wrong weight map size"); }

    /* Get the weight table: */
    double wt[nw];
    uint32_t stride;
    wt_table_binomial_fill(nw, wt, &stride);
    wt_table_normalize_sum(nw, wt);
    assert(stride == (nw - 1)/2);

    /* Compute displacement of window: */
    int32_t dxy = ((int32_t)nw - 1)/2;
    
    /* Alocate the output image: */
    float_image_t *V = float_image_new(1, (int32_t)NXV, (int32_t)NYV);
    
    for (int32_t y = 0; y < NYV; y++)
      { for (int32_t x = 0; x < NXV; x++)
          { double sum_cw = 0.0; /* Sum of window weights times {W} samples. */
            double sum_c = 0.0;  /* Sum of window weights. */
            for (uint32_t r = 0; r < nw; r++)
              { for (uint32_t s = 0; s < nw; s++) 
                  { /* Get the sample from {W}: */
                    int32_t xx = x + (int32_t)r - dxy;
                    int32_t yy = y + (int32_t)s - dxy;
                    float wrs = 
                      ( (xx >= 0) && (xx < NXW) && (yy >= 0) && (yy < NYW) ?
                        float_image_get_sample(W, 0, xx, yy) :
                        0.0f
                      );
                    demand(wrs >= 0.0, "negative weight");
                    demand(! isnan(wrs), "weight is NAN");
                    if (harmonic) { wrs = (float)(1.0/wrs); }
                    /* Get the window weight {c}: */
                    double c = wt[r]*wt[s];  /* Window weight. */
                    /* Accumulate: */
                    sum_cw += c*wrs;
                    sum_c += c;
                  }
              }
            double avg = (harmonic ? sum_c/sum_cw : sum_cw/sum_c);
            float_image_set_sample(V, 0, x, y, (float)avg);
          }
      }
    return V;
  }
  
float_image_t *pst_weight_map_heights_from_slopes(int32_t NX_Z, int32_t NY_Z,float_image_t *GW)
  {
    float_image_t* EW = float_image_new(1, (int32_t)NX_Z, (int32_t)NY_Z);
    if ( GW != NULL ){ assert( (GW->sz[1] == (NX_Z-1)) && (GW->sz[2] == (NY_Z -1)) ); }
    for (int32_t y = 0 ; y < NY_Z ; y++)
      { for (int32_t x = 0 ; x < NX_Z ; x++ )
          {
            double wmm = ( (x  > 0 ) && (y > 0) ? (GW != NULL ? float_image_get_sample(GW,0,x-1,y-1) : 1): 0);
            double wmp = ( (x  > 0 ) && (y < (NY_Z -1)) ? (GW != NULL ? float_image_get_sample(GW,0,x-1,y) : 1): 0);
            double wpm = ( (x  < (NX_Z -1) ) && (y > 0) ? (GW != NULL ? float_image_get_sample(GW,0,x,y-1) : 1): 0);
            double wpp = ( (x  < (NX_Z -1) ) && (y < (NY_Z -1) ) ? (GW != NULL ? float_image_get_sample(GW,0,x,y) : 1): 0);

            double w = 16.0/(1/wmm + 1/wpm + 1/wmp + 1/wpp);
            float_image_set_sample(EW,0,x,y,(float)w);

          }
      }
    return EW;
  }
