/* See {test_voxm_tubes.h} */
/* Last edited on 2021-06-22 13:47:52 by jstolfi */

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

#include <r3_path.h>
#include <voxm_splat.h>
#include <voxm_splat_tube.h>

#include <test_voxm_tubes.h>
#include <test_voxm_mark.h>

void test_voxm_tubes(ppv_array_t *A, r3_t *ctr, r3_t *rad, double fuzzR)
  { 
    int NT = 3; /* Number of independent tests. */
    
    /* Which tubes to splat: */
    bool_t doit[NT];
    doit[0] =  TRUE; /* {test_voxm_tubes_helix}. */
    doit[1] =  TRUE; /* {test_voxm_tubes_segment}. */
    doit[2] =  TRUE; /* {test_voxm_tubes_bezier}. */
    
    /* Divide box into equal sub_boxes: */
    int k; 
    int NB = 0; /* Number of active tests: */
    for (k = 0; k < NT; k++) { if (doit[k]) { NB++; } } 
    r3_t rsu = (*rad); rsu.c[0] /= NB;  /* Half-size of each sub-box. */
    r3_t csu; r3_sub(ctr, rad, &csu); r3_add(&csu, &rsu, &csu);  /* Center of next sub-box. */
    
    /* Splat the tests: */
    for (k = 0; k < NT; k++)
      { if (doit[k])
          { /* Mark the {Y} and {Z} edges of the sub-box: */
            test_voxm_mark_edges(A, &csu, &rsu, 1, fuzzR);
            test_voxm_mark_edges(A, &csu, &rsu, 2, fuzzR);

            /* Splat the tubes: */
            switch (k)
              { case 0: test_voxm_tubes_helix  (A, &csu, &rsu, fuzzR); break;
                case 1: test_voxm_tubes_segment(A, &csu, &rsu, fuzzR); break;
                case 2: test_voxm_tubes_bezier (A, &csu, &rsu, fuzzR); break;
                default: assert(FALSE);
              }

            /* Advance to the next box: */
            csu.c[0] += 2*rsu.c[0];
          }
      }
  }
  
void test_voxm_tubes_helix(ppv_array_t *A, r3_t *ctr, r3_t *rad, double fuzzR)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);

    int N = 3; /* Number of tubes. */
    r3_motion_state_t S[N]; /* Placement states of tubes. */
    double t0[N]; /* Start time. */
    double t1[N]; /* End time. */
    double len[N]; /* Representative arc lengths. */
    double ang[N]; /* Representative turn angles. */
    double hht[N]; /* Representative advances. */
    
    /* Define the end states assuming that the domain has radius 1 and center at {(0,0,0)}: */
    int k;
    for (k = 0; k < N; k++)
      { double pY = 2*(k + 0.5)/N - 1.0;
        double pZ = -0.50;
        S[k].p = (r3_t){{ 00.00, pY, pZ }};
        r3x3_ident(&(S[k].M));
      }
    /* Tilt the middle spiral: */
    int j;
    double ca = cos(M_PI/6);
    double sa = sin(M_PI/6);
    for (j = 0; j < 3; j++)
      { double Xj = S[1].M.c[0][j];
        double Zj = S[1].M.c[2][j];
        S[1].M.c[0][j] = + ca*Xj - sa*Zj;
        S[1].M.c[2][j] = + sa*Xj + ca*Zj;
      }
    
    /* Tube parameters: */
    double inR = 4.0;
    double otR = 7.5;
    
    /* Scale and translate the end states for center {ctr} and radius {rad}: */
    double trad = fmin(rad->c[0], fmin(rad->c[1], rad->c[2]));  /* Half-side of inscribed cube. */
    double scale = trad - otR;
    r3_t shift = (*ctr);
    for (k = 0; k < N; k++)
      { test_voxm_rescale_r3_motion_state(&(S[k]), scale, &shift); }

    /* Define the times, lengths, angles, and steps (already scaled): */
    double hrad = 0.25*trad; /* Typical helix radius. */
    t0[0] = 00.000; t1[0] = +1.000; len[0] = 1.00*hrad*M_PI; ang[0] = 2*M_PI; hht[0] = 0.50*hrad+ 2*otR;
    t0[1] = -1.000; t1[1] = +2.000; len[1] = 1.00*hrad*M_PI; ang[1] = 2*M_PI; hht[1] = 0.35*hrad+ 2*otR;
    t0[2] = 00.000; t1[2] = +1.000; len[2] = 1.50*hrad*M_PI; ang[2] = 2*M_PI; hht[2] = 0.25*hrad+ 2*otR;

    /* Splat the tubes, then clear the bores: */
    int isub;
    for (isub = 0; isub < 2; isub++)
      { bool_t sub = (isub == 1); /* FALSE to splat the tubes, TRUE to clear the hole. */
        for (k = 0; k < N; k++)
          { 
            fprintf(stderr, "  helix %d  t = [ %.3f _ %.3f ]  len = %.2f", k, t0[k], t1[k], len[k]);
            fprintf(stderr, "  ang = %.2f deg  hht = %.2f\n", ang[k]*180/M_PI, hht[k]);
            r3_motion_state_debug(stderr, &(S[k]), "  ", "placement");
            r3_motion_state_t S0, S1;
            voxm_splat_tube_round_helix(A, t0[k], t1[k], &(S[k]), len[k], ang[k], hht[k], inR, otR, fuzzR, sub, &S0, &S1);
            r3_motion_state_debug(stderr, &(S0), "  ", "initial");
            r3_motion_state_debug(stderr, &(S1), "  ", "final");
          }
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "exit %s\n", __FUNCTION__);
    
    return;
  }

void test_voxm_tubes_segment(ppv_array_t *A, r3_t *ctr, r3_t *rad, double fuzzR)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);

    int N = 4; /* Number of tubes. */
    r3_path_state_t S[N], T[N]; /* Initial and final states of tubes. */
    
    /* Define the end states assuming that the domain has radius 1 and center at {(0,0,0)}: */
    /* A straight tube parallel to the {X}-axis: */
    S[0].p = (r3_t){{ -1.00, -0.50, +0.50 }}; S[0].v = (r3_t){{ +1.0, 00.0, 00.0 }};
    T[0].p = (r3_t){{ 00.00, -0.50, +0.50 }}; T[0].v = (r3_t){{ +1.0, 00.0, 00.0 }};
    
    /* A nearly half-circular U-turn extension of the former: */
    double v1 = sqrt(0.50)*M_PI; 
    S[1].p = T[0].p;                          S[1].v = (r3_t){{ +v1, 00.00, 00.00 }};
    T[1].p = (r3_t){{ 00.00, +0.50, -0.50 }}; T[1].v = (r3_t){{ -v1, 00.00, 00.00 }};
    
    /* An S-shaped path extending the former: */
    double v2 = 1.00*M_PI; 
    S[2].p = T[1].p;                          S[2].v = (r3_t){{ -v2, 00.00, 00.00 }};
    T[2].p = (r3_t){{ 00.00, -0.50, -0.50 }}; T[2].v = (r3_t){{ -v2, 00.00, 00.00 }};
    
    /* A fork from the former: */
    double v3 = 0.25*M_PI; 
    S[3].p = T[1].p;                          S[3].v = (r3_t){{ -v3,   00.00, 00.00 }};
    T[3].p = (r3_t){{ -0.50, +1.00, -0.50 }}; T[3].v = (r3_t){{ 00.00, +v3,   00.00 }};

    /* Tube parameters: */
    double inR = 4.0;
    double otR = 7.5;
    
    /* Scale and translate the end states for center {ctr} and radius {rad}: */
    double trad = fmin(rad->c[0], fmin(rad->c[1], rad->c[2]));  /* Half-side of inscribed cube. */
    double scale = trad - otR;
    r3_t shift = (*ctr); /* shift.c[0] += 00.00*trad; */
    int k;
    for (k = 0; k < N; k++)
      { S[k].t = 0.00;
        test_voxm_rescale_r3_path_state(&(S[k]), scale, &shift);
        T[k].t = 1.00;
        test_voxm_rescale_r3_path_state(&(T[k]), scale, &shift);
      }
    
    /* Splat the tubes, then clear the bores: */
    int isub;
    for (isub = 0; isub < 2; isub++)
      { bool_t sub = (isub == 1); /* FALSE to splat the tubes, TRUE to clear the hole. */
        for (k = 0; k < N; k++)
          { voxm_splat_tube_round_segment(A, &(S[k]), &(T[k]), inR, otR, fuzzR, sub); }
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "exit %s\n", __FUNCTION__);
    
    return;
  }
  
void test_voxm_tubes_bezier(ppv_array_t *A, r3_t *ctr, r3_t *rad, double fuzzR)
  { 
    fprintf(stderr, "enter %s\n", __FUNCTION__);
    
    int N = 3; /* Number of tubes. */
    r3_t p0[N], p1[N], p2[N], p3[N]; /* Bezier controlpoints of segments. */
    
    /* Define the end states assuming that the domain has radius 1 and center at {(0,0,0)}: */
    p0[0] = (r3_t){{ -1.00, 00.00, -1.00 }}; 
    p1[0] = (r3_t){{ -1.00, 00.00, +1.00 }};
    p2[0] = (r3_t){{ +1.00, -1.00, -1.00 }};
    p3[0] = (r3_t){{ +1.00, 00.00, 00.00 }};
    
    p0[1] = p3[0]; 
    r3_mix(-1.0, &(p2[0]), +2.0, &(p3[0]), &(p1[1]));
    p2[1] = (r3_t){{ -1.00, +1.00, +1.00 }};
    p3[1] = (r3_t){{ -1.00, 00.00, 00.00 }};
    
    p0[2] = (r3_t){{ 00.00, -0.60, +0.50 }}; 
    p1[2] = (r3_t){{ 00.00, -0.20, +0.50 }};
    p2[2] = (r3_t){{ 00.00, +0.20, +0.50 }};
    p3[2] = (r3_t){{ 00.00, +0.60, +0.50 }};
    
    /* Tube parameters: */
    double inR = 4.0;
    double otR = 7.5;
    
    /* Scale and translate the control points center {ctr} and radius {rad}: */
    double trad = fmin(rad->c[0], fmin(rad->c[1], rad->c[2]));  /* Half-side of inscribed cube. */
    double scale = trad - otR;
    r3_t shift = (*ctr);
    int k;
    for (k = 0; k < N; k++)
      { test_voxm_rescale_r3(&(p0[k]), scale, &shift);
        test_voxm_rescale_r3(&(p1[k]), scale, &shift);
        test_voxm_rescale_r3(&(p2[k]), scale, &shift);
        test_voxm_rescale_r3(&(p3[k]), scale, &shift);
      }

    /* Splat the tubes, then clear the bores: */
    int isub;
    for (isub = 0; isub < 2; isub++)
      { bool_t sub = (isub == 1); /* FALSE to splat the tubes, TRUE to clear the hole. */
        for (k = 0; k < N; k++)
          { voxm_splat_tube_round_bezier(A, &(p0[k]), &(p1[k]), &(p2[k]), &(p3[k]), inR, otR, fuzzR, sub); }
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "exit %s\n", __FUNCTION__);
  }
    
void test_voxm_rescale_r3(r3_t *p, double scale, r3_t *shift)
  { r3_mix(1.0, shift, scale, p, p); }
    
void test_voxm_rescale_r3_path_state(r3_path_state_t *P, double scale, r3_t *shift)
  { r3_mix(scale, &(P->p), 1.0, shift, &(P->p));
    r3_scale(scale, &(P->v), &(P->v));
  }

void test_voxm_rescale_r3_motion_state(r3_motion_state_t *P, double scale, r3_t *shift)
  { r3_mix(scale, &(P->p), 1.0, shift, &(P->p)); }
