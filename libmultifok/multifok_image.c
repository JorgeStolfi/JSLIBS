/* See {multifok_image.h}. */
/* Last edited on 2025-04-11 08:55:11 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <i2.h>
#include <affirm.h>
#include <jsprintf.h>

#include <float_image.h>
#include <float_image_paint.h>

#include <multifok_image.h>

float_image_t *multifok_image_set_weight_channel(float_image_t *img, float_image_t *wht, int32_t chns, bool_t clear)
  { int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    demand((NC == chns) || (NC == chns+1), "arg {chsn} inconsistent with channel count");
  
    float_image_t *timg = float_image_new(chns+1, NX, NY);
    
    /* Copy the data channels of {img} to {timg}: */
    for (int32_t c = 0; c < chns; c++)
      { float_image_assign_channel_rectangle(timg, c, 0,NX-1, 0,NY-1, img, c, 0,0); }
      
    if (wht != NULL)
      { /* Set the weight channel from {wht}: */
        float_image_check_size(wht, 1, NX, NY, "inconsistent sizes {img,wht}");
        float_image_assign_channel_rectangle(timg, chns, 0,NX-1, 0,NY-1, wht, 0, 0,0);
      }
    else if (NC == chns+1)
      { /* Copy weight channel from {img} too: */
        float_image_assign_channel_rectangle(timg, chns, 0,NX-1, 0,NY-1, img, chns, 0,0);
      }
    else
      { /* Fill weight channel with ones: */
        float_image_fill_channel(timg, chns, 1.0);
      }
      
    if (clear)
      { for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { float wxy = float_image_get_sample(timg, chns, x, y);
                if (wxy == 0)
                  { for (int32_t c = 0; c < chns; c++) 
                      { float_image_set_sample(timg, c, x, y, NAN); }
                  }
              }
          }
      }
                      
    return timg;
  }

void multifok_image_draw_crosses(float_image_t *img, int32_t ch, uint32_t NQ, i2_t pix[], float val)
  { 
    int32_t NC = (int32_t)img->sz[0];
    
    /* Draw crosses over the image: */
    double rad = 6.0;
    double hwd = 0.5;
    bool_t empty = TRUE;
    bool_t diagonal = TRUE;
    for (uint32_t kq = 0;  kq < NQ; kq++)
      { double xctr = pix[kq].c[0] + 0.5;
        double yctr = pix[kq].c[1] + 0.5;
        for (int32_t ic = 0;  ic < NC; ic++)
          { float_image_paint_cross(img, ic, xctr, yctr, rad, empty, 2.0*hwd, diagonal, 0.0f, 3);
            float valk = (ic == ch ? val : 0.0f);
            float_image_paint_cross(img, ic, xctr, yctr, rad, empty, hwd, diagonal, valk, 3);
          }
      }
  }

#define multifok_image_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

