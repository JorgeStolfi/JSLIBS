/* See {test_voxb_obj.h} */
/* Last edited on 2025-01-04 23:56:42 by stolfi */

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

#include <test_voxb_objs.h>

void test_voxb_objs(ppv_array_t *A, r3_t *ctr, r3_t *rad)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);

    /* Get the half-size {orad} of largest cube that fits in the box: */
    double orad = fmin(rad->c[0], fmin(rad->c[1], rad->c[2])); 

    test_voxb_objs_ball         (A, ctr, orad);
    test_voxb_objs_donut        (A, ctr, orad);
    test_voxb_objs_rod          (A, ctr, orad);
    test_voxb_objs_tube         (A, ctr, orad);
    test_voxb_objs_cube_hole    (A, ctr, orad);
    test_voxb_objs_box          (A, ctr, orad);
    test_voxb_objs_rounded_box  (A, ctr, orad);
    test_voxb_objs_cup          (A, ctr, orad);
    
    fprintf(stderr, "\n");
    fprintf(stderr, "exit %s\n", __FUNCTION__);
    return;
  }

void test_voxb_objs_ball(ppv_array_t *A, r3_t *ctr, double rad)
  { 
    double ballR = 0.05*rad;  /* Radius of ball. */
    double ballX = 0.80*rad;  /* {X}-displacement. */
    
    r3_motion_state_t state;
    state.p = (*ctr); state.p.c[0] += ballX;
    r3x3_ident(&(state.M)); 
    
    auto bool_t ball(r3_t *p);
      /* Indicator function for a canonical ball of radius {ballR}. */
       
    voxb_splat_object(A, ball, &state, ballR, voxb_op_OR, FALSE);
    
    return;

    bool_t ball(r3_t *p)
      { return voxb_obj_ball(p, ballR); }
  }      

void test_voxb_objs_donut(ppv_array_t *A, r3_t *ctr, double rad)
  { 
    double majR = 0.25*rad;  /* Radius of donut midline. */
    double minR = majR/3; /* Radius of dough. */
    
    r3_motion_state_t state;
    state.p = (*ctr);
    r3x3_ident(&(state.M)); 
    
    auto bool_t donut(r3_t *p);
      /* Indicator function for a canonical donut of major radius {majR},
        minor radius {minR}. */
       
    voxb_splat_object(A, donut, &state, majR + minR, voxb_op_OR, FALSE);
    
    return;

    bool_t donut(r3_t *p)
      { return voxb_obj_donut(p, minR, majR, 2); }
  }      

void test_voxb_objs_rod(ppv_array_t *A, r3_t *ctr, double rad)
  { 
    double rodH = 0.20*rad;   /* Half-height of rod. */
    double rodR = 0.10*rad;   /* Radius of rod. */
    double rodF = 0.25*rodR;  /* Fillet radius. */
    double rodA = M_PI/12;    /* Tilting angle. */
    
    r3_motion_state_t state;
    state.p = (*ctr);
    r3_t u = (r3_t){{ 0.0, 0.0, 1.0 }};
    r3_t v = (r3_t){{ sin(rodA), 0.0, cos(rodA) }};
    r3x3_u_to_v_rotation(&u, &v, &(state.M)); 
    
    auto bool_t rod(r3_t *p);
      /* Indicator function for a canonical rod with height {2*rodH},
        radius {rodR}, fillet radius {rodF}. */
       
    voxb_splat_object(A, rod, &state, hypot(rodH, rodR), voxb_op_OR, FALSE);
    
    return;
      
    bool_t rod(r3_t *p)
      { return voxb_obj_rod(p, rodH, rodR, rodF); }
  }    
  
void test_voxb_objs_tube(ppv_array_t *A, r3_t *ctr, double rad)
  { 
    double tubeH = 0.20*rad;    /* Half-height of tube. */
    double tubeRi = 0.10*rad;   /* Radius of tube. */
    double tubeRo = 0.30*rad;   /* Radius of tube. */
    double tubeF = 0.25*tubeRi; /* Fillet radius. */
    double tubeA = M_PI/12;     /* Tilting angle. */
    
    double tubeX = -0.70*rad; /* X displacement. */
    double tubeY = +0.20*rad; /* Y displacement. */

    r3_motion_state_t state;
    state.p = (*ctr); state.p.c[0] += tubeX; state.p.c[1] += tubeY;
    r3_t u = (r3_t){{ 0.0, 0.0, 1.0 }};
    r3_t v = (r3_t){{ sin(tubeA), 0.0, cos(tubeA) }};
    r3x3_u_to_v_rotation(&u, &v, &(state.M)); 
    
    auto bool_t tube(r3_t *p);
      /* Indicator function for a canonical tube with height {2*tubeH},
        radii {tubeRi,tubeRo}, fillet radius {tubeF}. */
       
    voxb_splat_object(A, tube, &state, hypot(tubeH, tubeRo), voxb_op_OR, FALSE);
    
    return;
      
    bool_t tube(r3_t *p)
      { return voxb_obj_tube(p, tubeH, tubeRi,tubeRo, tubeF); }
  }    
  
void test_voxb_objs_cube_hole(ppv_array_t *A, r3_t *ctr, double rad)
  { 
    double cubeZ = 0.50*rad;   /* Distance above array center of cube center and hole center. */
    double cubeA = M_PI/12;    /* Tilting angle. */
    
    r3_motion_state_t state;
    state.p = (*ctr); state.p.c[2] += cubeZ;
    r3_t u = (r3_t){{ 1.0, 0.0, 0.0 }};
    r3_t v = (r3_t){{ cos(cubeA), 0.0, sin(cubeA) }};
    r3x3_u_to_v_rotation(&u, &v, &(state.M)); 

    double cubeH = 0.15*rad;   /* Half-side of cube. */
    double cubeF = 0.25*cubeH; /* Fillet radius. */
    
    auto bool_t cube(r3_t *p);
      /* Indicator function for a canonical cube with side {2*cubeH},
        fillet radius {cubeF}. */
       
    voxb_splat_object(A, cube, &state, sqrt(3)*cubeH, voxb_op_OR, FALSE);
    
    /* Stretch and squeeze cube to make the hole: */
    r3x3_t K;
    r3x3_ident(&K);
    K.c[0][0] = 0.5;
    K.c[1][1] = 0.5;
    K.c[2][2] = 2.0;
    r3x3_mul(&K, &(state.M), &(state.M));
    voxb_splat_object(A, cube, &state, sqrt(3)*cubeH, voxb_op_SUB, FALSE);
    
    return;
      
    bool_t cube(r3_t *p)
      { return voxb_obj_cube(p, cubeH, cubeF); }
  }    

void test_voxb_objs_box(ppv_array_t *A, r3_t *ctr, double rad)
  {
    double boxZ = -0.40*rad;  /* {Z}-distance from array center to box center. */
    double boxA = -M_PI/12;   /* Tilting angle. */
    
    r3_motion_state_t state;
    state.p = (*ctr); state.p.c[2] += boxZ;
    r3_t u = (r3_t){{ 1.0, 0.0, 0.0 }};
    r3_t v = (r3_t){{ cos(boxA), 0.0, sin(boxA) }};
    r3x3_u_to_v_rotation(&u, &v, &(state.M)); 

    double boxRX = 0.18*rad;   /* Half {X}-side of box. */
    double boxRY = 0.12*rad;   /* Half {Y}-side of box. */
    double boxRZ = 0.06*rad;   /* Half {Z}-side of box. */
    double boxF =  0.03*rad;   /* Fillet radius. */
    
    double boxRXYZ = hypot(hypot(boxRX, boxRY), boxRZ);
    
    auto bool_t box(r3_t *p);
      /* Indicator function for a canonical box with sides {2*boxRX}, 
        {2*boxRY}, {2*boxRZ}, fillet radius {boxF}. */
       
    voxb_splat_object(A, box, &state, boxRXYZ, voxb_op_OR, FALSE);
    
    return;
    
    bool_t box(r3_t *p)
      { return voxb_obj_box(p, boxRX, boxRY, boxRZ, boxF); }
  }

void test_voxb_objs_rounded_box(ppv_array_t *A, r3_t *ctr, double rad)
  {
    double boxZ = -0.65*rad;  /* {Z}-distance from array center to box center. */
    double boxA = -M_PI/12;   /* Tilting angle. */
    
    r3_motion_state_t state;
    state.p = (*ctr); state.p.c[2] += boxZ;
    r3_t u = (r3_t){{ 1.0, 0.0, 0.0 }};
    r3_t v = (r3_t){{ cos(boxA), 0.0, sin(boxA) }};
    r3x3_u_to_v_rotation(&u, &v, &(state.M)); 

    double boxRX = 0.24*rad;    /* Half {X}-side of box. */
    double boxRY = 0.18*rad;    /* Half {Y}-side of box. */
    double boxRZ = 0.12*rad;    /* Half {Z}-side of box. */
    double boxRoundR = boxRY/2; /* The base corner rounding radius. */ 
    double boxF =  0.03*rad;    /* Fillet radius. */
    
    double boxRXYZ = hypot(hypot(boxRX, boxRY), boxRZ);
    
    auto bool_t box(r3_t *p);
      /* Indicator function for a canonical box with sides {2*boxRX}, 
        {2*boxRY}, {2*boxRZ}, fillet radius {boxF}. */
       
    voxb_splat_object(A, box, &state, boxRXYZ, voxb_op_OR, FALSE);
    
    return;
    
    bool_t box(r3_t *p)
      { return voxb_obj_rounded_box(p, boxRX, boxRY, boxRZ, boxRoundR, boxF); }
  }

void test_voxb_objs_cup(ppv_array_t *A, r3_t *ctr, double rad)
  {
    /* Basic cup parameters: */
    double cupX = 0.65*rad;
    double cupY = 0.65*rad;
    double cupZ = 0.65*rad;

    double cupH = 0.30*rad;  /* Half-height of cup. */
    
    /* Derived cup parameters: */
    double cupR = 1.5*cupH;        /* Radius of cup. */
    double cupT = 0.20*cupR;       /* Thickness of cup wall. */
    double cupF = 2*cupT;          /* Radius of base filet. */
    
    /* Position and orientation of cup: */
    r3_motion_state_t state;
    state.p = (r3_t){{ cupX, cupY, cupZ }};
    r3x3_ident(&(state.M));
    
    auto bool_t cup(r3_t *p);
      /* Indicator function for the cup. */
       
    fprintf(stderr, "splatting round cup: halfH = %.2f  R = %.2f  thk = %.2f  fillR = %.2f\n", cupH, cupR, cupT, cupF);
    voxb_splat_object(A, cup, &state, hypot(cupH, cupR), voxb_op_OR, FALSE);
    
    return;
      
    bool_t cup(r3_t *p)
      { return voxb_obj_round_cup(p, cupH, cupR, cupT, cupF); }
  }
  
