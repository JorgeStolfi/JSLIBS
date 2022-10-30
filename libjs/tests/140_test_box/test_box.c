#define PROG_NAME "test_box"
#define PROG_DESC "tests the ordered table search procedure"
#define PROG_VERS "1.1"

/* Last edited on 2022-10-30 12:11:52 by stolfi */
/* Created on 2021-09-25 or earler by J. Stolfi, UNICAMP */

#define test_box_COPYRIGHT \
  "Copyright © 2021  by the State University of Campinas (UNICAMP)"

#define _ISOC99_SOURCE 1
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include <jsrandom.h>
#include <affirm.h>
#include <interval.h>
#include <interval_io.h>
#include <bool.h>

#include <box.h>

#define GLOB_LO (-10.0)
#define GLOB_HI (+10.0)

/* PROTOTYPES */

int32_t main (int32_t argc, char **argv);
void do_tests(box_dim_t d);
void test_attribs(box_dim_t d, interval_t B[]);
void test_equal(box_dim_t d, interval_t B[]);
void test_include_point(box_dim_t d, interval_t B[]);
void test_join_meet(box_dim_t d, interval_t B[]);
void test_split(box_dim_t d, interval_t B[]);
void test_widen_round(box_dim_t d, interval_t B[]);
void test_box_map(box_dim_t d, interval_t B[]);
void test_face(box_dim_t d, interval_t B[]);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    srand(4615);
    
    int32_t nt = 10000;
    
    for (int32_t it = 0; it < nt; it++)
      { /* Choose the dimension {d}: */
        box_dim_t d = (box_dim_t)int32_abrandom(0, box_MAX_DIM);
        do_tests(d);
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
  
void do_tests(box_dim_t d)
  { 
    /* TESTING: void box_throw(box_dim_t d, double elo, double ehi, double p_empty, double p_single, interval_t B[]); */
    /* TESTING: void box_is_empty(box_dim_t d, interval_t B[]); */
    /* TESTING: void box_is_empty(box_dim_t d, interval_t B[]); */

    interval_t B[d]; box_throw(d, GLOB_LO, GLOB_HI, 0.100, 0.100, B);
  
    test_attribs(d, B);
    
    test_equal(d, B);
    
    test_include_point(d, B);
    test_join_meet(d, B);
    test_widen_round(d, B);
    
    test_split(d, B);
            
    test_box_map(d, B);   
    test_face(d, B);   
  }
  
void test_attribs(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_{lo_corner,hi_corner,corner,center,half_widths,widths,max_width,radius}} ---\n");
    
    double ctr[d], h[d], w[d], plo[d], phi[d], pdir[d];
    double er, mw;

    interval_side_t dir[d];
    for (int32_t i = 0; i < d; i++) { dir[i] = (interval_side_t)int32_abrandom(0, 1); }
    
    /* TESTING: void box_lo_corner(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_hi_corner(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_corner(box_dim_t d, interval_t B[], interval_side_t dir[], double p[]); */
    /* TESTING: void box_center(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_half_widths(box_dim_t d, interval_t B[], double h[]); */
    /* TESTING: void box_widths(box_dim_t d, interval_t B[], double w[]); */
    /* TESTING: double box_max_width(box_dim_t d, interval_t B[]); */
    /* TESTING: double box_radius(box_dim_t d, interval_t B[]); */
    
    if (! box_is_empty(d, B))
      { box_lo_corner(d, B, plo);
        box_hi_corner(d, B, phi);
        box_corner(d, B, dir, pdir);
        box_center(d, B, ctr);
      }
    box_half_widths(d, B, h);
    box_widths(d, B, w);
    mw = box_max_width(d, B);
    er = box_radius(d, B);
    
    double mw_cmp = 0;
    double sum2 = 0;
    for (int32_t i = 0; i < d; i++)
      { interval_t Bi = B[i];
        if (! box_is_empty(d, B))
          { assert(LO(Bi) == plo[i]);
            assert(HI(Bi) == phi[i]);
            assert(Bi.end[dir[i]] == pdir[i]);
            assert(interval_mid(&Bi) == ctr[i]);
          }
        assert(interval_rad(&Bi) == h[i]);
        assert(fabs(w[i] - 2*h[i]) < 1.0e-13);
        mw_cmp = fmax(mw_cmp, w[i]);
        sum2 += h[i]*h[i];
      }
    assert(mw_cmp == mw);
    assert(fabs(er - sqrt(sum2)) < 1.0e-12);
  }
     
void test_equal(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_equal} ---\n");

    interval_t A[d]; box_throw(d, GLOB_LO, GLOB_HI, 0.100, 0.100, A);

    assert(box_equal(d, B, B));
    if (box_is_empty(d, A) && box_is_empty(d, B))
      { assert(box_equal(d, A, B)); }
    else if (box_is_empty(d, A) || box_is_empty(d, B))
      { assert(! box_equal(d, A, B)); }
    else
      { bool_t eq = TRUE;
        for (int32_t i = 0; i < d; i++)
          { interval_t Ai = A[i];
            interval_t Bi = B[i];
            if (LO(Ai) != LO(Bi)) { eq = FALSE; }
            if (HI(Ai) != HI(Bi)) { eq = FALSE; }
          }
        assert(eq == box_equal(d, A, B));
      }
  }
  
void test_include_point(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_include_point} ---\n");

    /* TESTING: void box_include_point(box_dim_t d, interval_t B[], double p[], interval_t C[]); */

    interval_t C[d];
    double p[d];
    for (int32_t i = 0; i < d; i++) { p[i] = dabrandom(GLOB_LO, GLOB_HI); }

    box_include_point(d, B, p, C);
    
    for (int32_t i = 0; i < d; i++)
      { interval_t Bi = B[i];
        interval_t Ci = C[i];
        if (box_is_empty(d, B))
          { assert(LO(Ci) == p[i]);
            assert(HI(Ci) == p[i]);
          }
        else
          { assert(LO(Ci) == fmin(LO(Bi), p[i]));
            assert(HI(Ci) == fmax(HI(Bi), p[i]));
          }
      }
  }

void test_join_meet(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_join,box_meet} ---\n");

    /* TESTING: void box_join(box_dim_t d, interval_t A[], interval_t B[], interval_t C[]); */
    /* TESTING: void box_meet(box_dim_t d, interval_t A[], interval_t B[], interval_t C[]); */

    interval_t A[d]; box_throw(d, GLOB_LO, GLOB_HI, 0.100, 0.100, A);

    interval_t CJ[d], CM[d];
    box_join(d, A, B, CJ);
    box_meet(d, A, B, CM);
    
    if (box_is_empty(d, A))
      { assert(box_equal(d, B, CJ));
        assert(box_is_empty(d, CM));
      }
    else if (box_is_empty(d, B))
      { assert(box_equal(d, A, CJ));
        assert(box_is_empty(d, CM));
      }
    else
      { /* Check join: */
        assert(! box_is_empty(d, CJ));
        for (int32_t i = 0; i < d; i++)
          { interval_t Ai = A[i];
            interval_t Bi = B[i];
            interval_t CJi = CJ[i];
            assert(LO(CJi) == fmin(LO(Ai), LO(Bi)));
            assert(HI(CJi) == fmax(HI(Ai), HI(Bi)));
          }
        /* Check whether meet should be empty: */
        double mety = FALSE;
        for (int32_t i = 0; i < d; i++)
          { interval_t Ai = A[i];
            interval_t Bi = B[i];
            double cmlo = fmax(LO(Ai), LO(Bi));
            double cmhi = fmin(HI(Ai), HI(Bi));
            if (cmlo > cmhi) { mety = TRUE; }
          }
        if (mety)
          { assert(box_is_empty(d, CM)); }
        else
          { for (int32_t i = 0; i < d; i++)
              { interval_t Ai = A[i];
                interval_t Bi = B[i];
                interval_t CMi = CM[i];
                assert(LO(CMi) == fmax(LO(Ai), LO(Bi)));
                assert(HI(CMi) == fmin(HI(Ai), HI(Bi)));
              }
          }
      }
  }

void test_split(box_dim_t d, interval_t B[])
  {
    bool_t debug = FALSE;
    
    fprintf(stderr, "--- testing {box_split} ---\n");
    if (d == 0) { return; }

    interval_t BLO[d], BMD[d], BHI[d];
    box_axis_t a = (box_axis_t)int32_abrandom(0,d-1);
    double x;
    bool_t ety = FALSE;
    if (box_is_empty(d, B))
      { x = 3.1415926; ety = TRUE; }
    else 
      { double toss = drandom();
        double alo = LO(B[a]), ahi = HI(B[a]);
        if (toss < 0.100)
          { /* Pick a value below the lower bound: */
            x = dabrandom(GLOB_LO - 1, alo);
          }
        else if (toss > 0.900)
          { /* Pick a value bove the upper bound: */
            x = dabrandom(ahi, GLOB_HI + 1);
          }
        else if (toss < 0.300)
          { /* Pick the lower bound: */
            x = alo;
          }
        else if (toss > 0.700)
          { /* Pick the upper bound: */
            x = ahi;
          }
        else if (alo== ahi)
          { /* Pick the singleton value: */
            x = alo;
          }
        else
          { /* Pick a value inside the range: */
            do 
              { x = dabrandom(alo, ahi); } 
            while ((x != alo) && (x != ahi));
          }
      }
      
    if (debug)
      { for (int32_t i = 0; i < d; i++)
          { interval_t Bi = B[i];
            fprintf(stderr, "  B[%02d] = ", i);
            interval_gen_print(stderr, &Bi,   "%12.6f", "( ", " ", " )");
            if (i == a) { fprintf(stderr, "  x = %12.6f", x); }
            fprintf(stderr, "\n");
          }
      }

    box_split(d, B, a, x, BLO, BMD, BHI);
    
    if (ety)
      {
        assert(box_is_empty(d, BLO));
        assert(box_is_empty(d, BMD));
        assert(box_is_empty(d, BHI));
      }
    else
      { interval_t Ba = B[a];
        assert(LO(Ba) <= HI(Ba));
        bool_t ety_BLO = (x <= LO(Ba));
        bool_t ety_BHI = (x >= HI(Ba));
        bool_t ety_BMD_sing = (LO(Ba) == HI(Ba)) && (x != LO(Ba));
        bool_t ety_BMD_open = (LO(Ba) != HI(Ba)) && ((x <= LO(Ba)) || (x >= HI(Ba)));
        bool_t ety_BMD = (ety_BMD_sing || ety_BMD_open);
        if (debug) 
          { fprintf(stderr, "  ety BLO = %d BMD = %d BHI = %d\n", ety_BLO, ety_BMD, ety_BHI); }
      
        for (int32_t i = 0; i < d; i++)
          { interval_t Bi = B[i];
            interval_t BLOi = BLO[i];
            interval_t BMDi = BMD[i];
            interval_t BHIi = BHI[i];
            if (debug)
              { fprintf(stderr, "  splitting axis a = %d by x = %12.6f\n", a, x);
                interval_gen_print(stderr, &Bi,   "%12.6f", "    Bi =   ( ", " ", " )\n");
                interval_gen_print(stderr, &BLOi, "%12.6f", "    BLOi = ( ", " ", " )\n");
                interval_gen_print(stderr, &BMDi, "%12.6f", "    BMDi = ( ", " ", " )\n");
                interval_gen_print(stderr, &BHIi, "%12.6f", "    BHIi = ( ", " ", " )\n");
              }
            assert(LO(Bi) <= HI(Bi));
            if (ety_BLO) 
              { assert(LO(BLOi) > HI(BLOi)); }
            else 
              { assert(LO(BLOi) == LO(Bi)); 
                if (i == a)
                  { assert(HI(BLOi) == fmin(x, HI(Bi))); }
                else
                  { assert(HI(BLOi) == HI(Bi)); }
              }
            
            if (ety_BHI)
              { assert(LO(BHIi) > HI(BHIi)); }
            else 
              { if (i == a)
                  { assert(LO(BHIi) == fmax(x, LO(Bi))); }
                else
                  { assert(LO(BHIi) == LO(Bi)); }
                assert(HI(BHIi) == HI(Bi)); 
              }
            
            if (ety_BMD)
              { assert(LO(BMDi) > HI(BMDi)); }
            else 
              { if (i == a)
                  { assert(LO(BMDi) == x); assert(HI(BMDi) == x); }
                else
                  { assert(LO(BMDi) == LO(Bi)); assert(HI(BMDi) == HI(Bi)); }
              }
          }
      }
  }

void test_widen_round(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_{widen,widen_var,round}} ---\n");

    /* TESTING: void box_widen(box_dim_t d, interval_t B[], double margin, interval_t C[]); */
    /* TESTING: void box_widen_var(box_dim_t d, interval_t B[], double mrglo[], double mrghi[], interval_t C[]); */
    /* TESTING: void box_round(box_dim_t d, interval_t B[], double unit, bool_t out, interval_t C[]); */

    double margin = 2.718281828;

    interval_t CWP[d], CWM[d];
    box_widen(d, B, +margin, CWP);
    box_widen(d, B, -margin, CWM);

    interval_t CWV[d];
    double mrglo[d], mrghi[d];
    box_widen_var(d, B, mrglo, mrghi, CWV);

    interval_t CRO[d], CRI[d];
    double unit = 1.42;
    box_round(d, B, unit, TRUE,  CRO);
    box_round(d, B, unit, FALSE, CRI);
                
    fprintf(stderr, "!!! NOT TESTED !!!\n");
  }
                  
void test_box_map(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_{point,box}_{mpa,unmap)}} ---\n");
      
    /* TESTING: void box_point_map(box_dim_t d, double z[], interval_t B[], double x[]); */
    /* TESTING: void box_point_unmap(box_dim_t d, double x[], interval_t B[], double z[]); */
    /* TESTING: void box_box_map(box_dim_t d, interval_t Z[], interval_t B[], interval_t X[]); */
    /* TESTING: void box_box_unmap(box_dim_t d, interval_t X[], interval_t B[], interval_t Z[]); */
                
    fprintf(stderr, "!!! NOT TESTED !!!\n");

  }
                  
void test_face(box_dim_t d, interval_t B[])
  {
    fprintf(stderr, "--- testing {box_face_{dimension,position,signature,rel_box,print}} ---\n");
      
    /* TESTING: box_dim_t box_face_dimension(box_dim_t m, box_face_index_t fi); */
    /* TESTING: box_signed_dir_t box_face_position(box_face_index_t fi, box_axis_index_t j); */
    /* TESTING: void box_face_signature(box_dim_t m, box_face_index_t fi, box_signed_dir_t dir[]); */
    /* TESTING: void box_face_rel_box(box_dim_t m, box_face_index_t fi, interval_t F[]); */
    /* TESTING: void box_face_print(FILE *wr, box_dim_t m, box_face_index_t fi); */
                
    fprintf(stderr, "!!! NOT TESTED !!!\n");

  }
