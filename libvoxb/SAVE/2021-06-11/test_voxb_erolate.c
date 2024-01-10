/* See {test_voxb_erolate.h} */
/* Last edited on 2021-06-11 13:22:38 by jstolfi */

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

#include <test_voxb_erolate.h>

#define NAX (ppv_array_NAXES)

void test_voxb_erolate_obj(ppv_array_t *a, int32_t isec, double smr);
  /* Splats a test object at the center of {a}, then applies erosion,
    dilation, or open-close op selected by the section number {isec},
    with smoothing radius {smr}. */

void test_voxb_erolate(ppv_array_t *a, r3_t *ctr, r3_t *rad, double smr)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);
    
    ppv_size_t NZ = a->size[0];
    ppv_size_t NY = a->size[1];
    ppv_size_t NX = a->size[2];
    
    
    /* Assume that {NX} is the largest, and split the tomogram in {X}. */
    int32_t nsec = 4;            /* Number of sections. */
    ppv_size_t szsec = NX/nsec;  /* Size of each section. */
    for (int32_t isec = 0; isec < nsec; isec++)
      { ppv_size_t xskip = isec*szsec;
        ppv_size_t yskip = (NY - szsec)/2;
        ppv_size_t zskip = (NZ - szsec)/2;
        ppv_array_t asec = (*a);
        ppv_crop(&asec, 0, zskip, szsec);
        ppv_crop(&asec, 1, yskip, szsec);
        ppv_crop(&asec, 2, xskip, szsec);
        
        test_voxb_erolate_obj(&asec, isec, smr);
      }
    return;     
  }

void test_voxb_erolate_obj(ppv_array_t *a, int32_t isec, double smr)
  { 
    ppv_size_t NZ = a->size[0];
    ppv_size_t NY = a->size[1];
    ppv_size_t NX = a->size[2];
    
    /* Center and radius of array: */
    r3_t ctr = (r3_t){{ 0.5*(double)NX, 0.5*(double)NY, 0.5*(double)NZ }};
    r3_t rad = (r3_t){{ 0.5*(double)NX, 0.5*(double)NY, 0.5*(double)NZ }};
    double rad_min = fmin(rad.c[0],fmin(rad.c[1],rad.c[2]));
    
    double cubeA = M_PI/12;       /* Tilting angle. */
    double cubeH = 0.25*rad_min;  /* Half-side of cube. */

    fprintf(stderr, "{test_voxb_erolate_obj}: [ NZ NY NX ] = [ %lu %lu %lu ]\n", NZ, NY, NX);
    r3_gen_print(stderr, &ctr, "%7.3f", "ctr = ( ", " ", " )\n");
    r3_gen_print(stderr, &rad, "%7.3f", "rad = ( ", " ", " )\n");
    fprintf(stderr, "cubeH = %7.3f\n", cubeH);

    r3x3_t R; /* Initial rotation matrix.*/
    r3_t u = (r3_t){{ 1.0, 0.0, 0.0 }};
    r3_t v = (r3_t){{ cos(cubeA), 0.0, sin(cubeA) }};
    r3x3_u_v_rotation(&u, &v, &R); 
    
    r3_motion_state_t state;
    state.p = ctr;
    state.M = R; 
    
    auto bool_t cube(r3_t *p);
      /* Indicator function for a canonical cube with side {2*cubeH} and sharp edges. */
       
    voxb_splat_object(a, cube, &state, sqrt(3)*cubeH, voxb_op_OR);
    /* ppv_index_t ixctr[NAX] = { 30,30,30,0,0,0 }; */
    /* ppv_set_sample(a, ixctr, 1); */
    
    /* Squeeze and stretch cube to make the pin: */
    r3x3_t L;
    r3x3_ident(&L);
    L.c[0][0] = 0.25;
    L.c[1][1] = 2.00;
    L.c[2][2] = 0.25;
    r3x3_mul(&L, &R, &(state.M));
    voxb_splat_object(a, cube, &state, sqrt(3)*cubeH, voxb_op_OR);
    
    /* Stretch and squeeze cube to make the hole: */
    r3x3_t K;
    r3x3_ident(&K);
    K.c[0][0] = 0.5;
    K.c[1][1] = 0.5;
    K.c[2][2] = 2.0;
    r3x3_mul(&K, &R, &(state.M));
    voxb_splat_object(a, cube, &state, sqrt(3)*cubeH, voxb_op_SUB);

    /* Now apply erosion/dilation: */
     if (isec == 0)
      { /* No smoothing. */ }
    else if (isec == 1)
      { voxb_erolate_with_ball(a, +smr); }
    else if (isec == 2)
      { voxb_erolate_with_ball(a, +smr);
        voxb_erolate_with_ball(a, -2*smr);
        voxb_erolate_with_ball(a, +smr);
      }
    else if (isec == 3)
      { voxb_erolate_with_ball(a, -smr); }
    else
      { assert(FALSE); }
   
    return;
      
    bool_t cube(r3_t *p)
      { return voxb_obj_cube(p, cubeH, 0.0); }
  }    

