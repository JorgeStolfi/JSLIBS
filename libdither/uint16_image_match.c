/* See {uint16_image_match.h}.  */
/* Last edited on 2024-12-26 19:34:28 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <jspnm.h>
#include <jsmath.h>
#include <affirm.h>
#include <bool.h>
#include <yuvhacks.h>
#include <uint16_color_table.h>

#include <uint16_image_match.h>

int32_t uint16_image_match_rgb
  ( int32_t r, int32_t g, int32_t b, 
    uint16_t maxval,
    uint16_color_table_t *chv, 
    uint16_t mapmaxval
  )
  { uint32_t NH = chv->ne;
    demand(NH > 0, "empty histogram");
    float s = (float)(1.0/maxval);
    float y1 = s*(float)r;
    float u1 = s*(float)g;
    float v1 = s*(float)b;
    float dist = +INF;;
    s = (float)(1.0/mapmaxval);
    int32_t ind = -1;
    for (int32_t i = 0; ((dist > 0) && (i < NH)); i++)
      { ppm_pixel_t *qq = &(chv->e[i].color);
        float y2 = s * qq->c[0];
        float u2 = s * qq->c[1];
        float v2 = s * qq->c[2];
        float newdist = (y1-y2)*(y1-y2) + (u1-u2)*(u1-u2) + (v1-v2)*(v1-v2);
        if (newdist < dist) { ind = i; dist = newdist; }
      }
    assert((ind >= 0) && (ind < NH));
    return ind;
  }

int32_t uint16_image_match_yuv 
  ( int32_t r, int32_t g, int32_t b, 
    uint16_t maxval,
    uint16_color_table_t *chv, 
    uint16_t mapmaxval
  )
  { uint32_t NH = chv->ne;
    demand(NH > 0, "empty histogram");
    float y1, u1, v1;
    rgb_to_yuv(r, g, b, maxval, &y1, &u1, &v1);
    float dist = +INF;
    int32_t ind = -1;
    for (int32_t i = 0; ((dist > 0) && (i < NH)); i++)
      { ppm_pixel_t *qq = &(chv->e[i].color);
        int32_t ri = qq->c[0];
        int32_t gi = qq->c[1];
        int32_t bi = qq->c[2];
        float y2, u2, v2;
        rgb_to_yuv(ri, gi, bi, mapmaxval, &y2, &u2, &v2);
        float newdist = (y1-y2)*(y1-y2) + (u1-u2)*(u1-u2) + (v1-v2)*(v1-v2);
        if (newdist < dist) { ind = i; dist = newdist; }
      }
    assert((ind >= 0) && (ind < NH));
    return ind;
  }
