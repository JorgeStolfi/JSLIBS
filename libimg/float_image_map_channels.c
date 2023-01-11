/* See float_image_map_channels.h */
/* Last edited on 2023-01-07 14:12:35 by stolfi */ 

#define _GNU_SOURCE
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <float_image.h>

#include <float_image_map_channels.h>

/* INTERNAL PROTOTYPES */

void float_image_map_channels
  ( float_image_t *imgA, 
    float_image_t *imgB, 
    float_image_map_channels_proc_t map 
  )
  {
    int32_t NCA, NCB, NX, NY;
    float_image_get_size(imgA, &NCA, &NX, &NY);
    float_image_get_size(imgB, &NCB, NULL, NULL);
    float_image_check_size(imgB, -1, NX, NY);
    float vA[NCA];
    float vB[NCB];
    for (int32_t iy = 0; iy < NY; iy++) 
      { for (int32_t ix = 0; ix < NX; ix++) 
          { float_image_get_pixel(imgA, ix, iy, vA);
            map(NCA, vA, NCB, vB);
            float_image_set_pixel(imgB, ix, iy, vB);
          }
      }
  }

void float_image_map_channels_RGB_to_YUV
  ( float_image_t *imgA, 
    float_image_t *imgB
  )
  {
    int32_t NCA, NCB;
    float_image_get_size(imgA, &NCA, NULL, NULL);
    demand(NCA >= 3, "input image must have at least 3 channels");
    float_image_get_size(imgB, &NCB, NULL, NULL);
        
    auto void map_RGB_to_YUV(int32_t NCA, float vA[], int32_t NCB, float vB[]);
      /* Assumes the first 3 samples {vA[0..2]} are an RGB triple,
        masp it to European TV YUV coordinates, and stores in {vB[0..2]}.
        If {NCA} and {NCB} are 4, copies {vA[3]} to {vB[3]}. */
    
    float_image_map_channels(imgA, imgB, map_RGB_to_YUV);
    
    return;
    
    void map_RGB_to_YUV(int32_t NCA, float vA[], int32_t NCB, float vB[])
      { frgb_t p = (frgb_t){{ vA[0], vA[1], vA[2] }};
        frgb_to_YUV(&p);
        for (int32_t ic = 0; ic < NCB; ic++)
          { vB[ic] = (ic < 3 ? p.c[ic] : (ic < NCA ? vA[ic] : 0.0f)); }
      }
  }
