/* See {test_voxb_mark}.h  */
/* Last edited on 2022-10-20 05:49:17 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <r3.h>
#include <r3x3.h>
#include <r3_motion.h>
#include <ppv_array.h>

#include <voxb_obj.h>
#include <voxb_splat.h>

#include <test_voxb_mark.h>


void test_voxb_mark_corners(ppv_array_t *A, r3_t *ctr, r3_t *rad)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);
    r3_gen_print (stderr, ctr, "%.2f", "  ctr = ( ", " ", " )\n");
    r3_gen_print (stderr, rad, "%.2f", "  rad = ( ", " ", " )\n");

    double ballR = 5.0; /* Ball radius in voxels. */
    
    auto bool_t ball(r3_t *p);
      /* Indicator function for a ball of radius {ballR}. */
    
    /* Array dimensions: */
    int32_t ix, iy, iz;
    for (ix = -1; ix <= +1; ix += 2)
      { for (iy = -1; iy <= +1; iy += 2)
          { for (iz = -1; iz <= +1; iz += 2)
              { /* Pick a corner of the domain: */
                r3_motion_state_t S; /* Ball position and pose. */
                S.p = (r3_t){{ ctr->c[0] + ix*rad->c[0], ctr->c[1] + iy*rad->c[1], ctr->c[2] + iz*rad->c[2] }};
                r3x3_ident(&(S.M)); 
                voxb_splat_object(A, ball, &S, ballR, voxb_op_OR, FALSE);
              }
          }
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "exit %s\n", __FUNCTION__);

    return;
    
    /* INTERNAL IMPLEMENTATIONS */

    bool_t ball(r3_t *p)
      { return voxb_obj_ball(p, ballR); }
  }
    
void test_voxb_mark_edges(ppv_array_t *A, r3_t *ctr, r3_t *rad, int32_t ax)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);
    r3_gen_print (stderr, ctr, "%.2f", "  ctr = ( ", " ", " )\n");
    r3_gen_print (stderr, rad, "%.2f", "  rad = ( ", " ", " )\n");

    int32_t xax = (ax + 1) % 3; /* The actual "{X}" axis. */
    int32_t yax = (ax + 2) % 3; /* The actual "{Y}" axis. */
    int32_t zax = (ax + 3) % 3; /* The actual "{Z}" axis. */

    double rodR = 3.0; /* Rod radius in voxels. */
    double rodH = rad->c[zax]; /* Rod half-height in voxels. */
    double rodF = 0.0; /* Radius of rod fillet. */
    fprintf(stderr, "  rodH = %.2f\n", rodH);
        
    auto bool_t rod(r3_t *p);
      /* Indicator function for a canonical rod of radius {rodR}. */
    
    /* Splat the rods: */
    int32_t ix, iy;
    for (ix = -1; ix <= +1; ix += 2)
      { for (iy = -1; iy <= +1; iy += 2)
          { /* Pick a corner of the domain: */
            r3_motion_state_t S; /* Rod position and pose. */
            S.p.c[xax] = ctr->c[xax] + ix*rad->c[xax];
            S.p.c[yax] = ctr->c[yax] + iy*rad->c[yax];
            S.p.c[zax] = ctr->c[zax];
            r3x3_zero(&(S.M)); 
            S.M.c[0][xax] = 1.0;
            S.M.c[1][yax] = 1.0;
            S.M.c[2][zax] = 1.0;
            voxb_splat_object(A, rod, &S, hypot(rodH,rodR), voxb_op_OR, FALSE);
          }
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "exit %s\n", __FUNCTION__);

    return;
    
    /* INTERNAL IMPLEMENTATIONS */

    bool_t rod(r3_t *p)
      { return voxb_obj_rod(p, rodH, rodR, rodF); }
  }

