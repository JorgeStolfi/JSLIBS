#define PROG_NAME "test_box"
#define PROG_DESC "tests the ordered table search procedure"
#define PROG_VERS "1.1"

/* Last edited on 2025-02-07 18:18:43 by stolfi */
/* Created on 2021-09-25 or earler by J. Stolfi, UNICAMP */

#define test_box_COPYRIGHT \
  "Copyright � 2021  by the State University of Campinas (UNICAMP)"

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
#include <jsmath.h>

#include <box.h>

#define GLOB_LO (-10.0)
#define GLOB_HI (+10.0)

/* PROTOTYPES */

int32_t main (int32_t argc, char **argv);
void tbx_do_tests(box_dim_t d, bool_t verbose);
void tbx_print_point(FILE *wr, box_dim_t d, double p[], char *pref, char *fmt, char *sep, char *suff);
  /* Prints {p[0..d-1]} on {wr} preceded by {pref} and followed by {suff}.
    Each element is printed with {fmt} and elements are separated by {sep}. */

void test_print(box_dim_t d, interval_t B[], bool_t verbose);
void test_attribs(box_dim_t d, interval_t B[], bool_t verbose);
void test_dimension__measure(box_dim_t d, interval_t B[], bool_t verbose);
void test_equal(box_dim_t d, interval_t B[], bool_t verbose);
void test_disjoint__contained(box_dim_t d, interval_t B[], bool_t verbose);
void test_has_point(box_dim_t d, interval_t B[], bool_t verbose);
void test_from_center_and_radii(box_dim_t d, interval_t B[], bool_t verbose);
void test_include_point(box_dim_t d, interval_t B[], bool_t verbose);
void test_join__meet(box_dim_t d, interval_t B[], bool_t verbose);
void test_split(box_dim_t d, interval_t B[], bool_t verbose);
void test_widen__round(box_dim_t d, interval_t B[], bool_t verbose);
void test_shift__unshift(box_dim_t d, interval_t B[], bool_t verbose);
void test_box_map(box_dim_t d, interval_t B[], bool_t verbose);
void test_face(box_dim_t d, interval_t B[], bool_t verbose);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    srand(4615*417+2);
    
    int32_t nt = 10000;
    
    box_dim_t d = 0;
    for (uint32_t it = 0;  it < nt; it++)
      { bool_t verbose = (it < 10);
        tbx_do_tests(d, verbose);
        if (it < 5)
          { d++; }
        else 
          { d = (box_dim_t)int32_abrandom(5, box_MAX_DIM); }
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
  
void tbx_do_tests(box_dim_t d, bool_t verbose)
  { 
    /* TESTING: void box_throw(box_dim_t d, double elo, double ehi, double p_empty, double p_single, interval_t B[]); */
    /* TESTING: void box_is_empty(box_dim_t d, interval_t B[]); */
    /* TESTING: void box_is_empty(box_dim_t d, interval_t B[]); */

    interval_t B[d]; box_throw(d, GLOB_LO, GLOB_HI, 0.100, 0.100, B);
  
    test_print(d, B, verbose);
  
    test_attribs(d, B, verbose);
    test_dimension__measure(d, B, verbose);
    
    test_equal(d, B, verbose);
    test_disjoint__contained(d, B, verbose);
    test_has_point(d, B, verbose);

    test_from_center_and_radii(d, B, verbose);
    test_include_point(d, B, verbose);
    test_join__meet(d, B, verbose);
    test_widen__round(d, B, verbose);
    test_shift__unshift(d, B, verbose);
    test_box_map(d, B, verbose); 
    
    test_split(d, B, verbose);
            
    test_box_map(d, B, verbose);   
    test_face(d, B, verbose);   
  }
  
void test_print(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if ((! verbose) || (d > 8)) { return; }
    if (verbose) { fprintf(stderr, "--- testing {box_{print,gen_print,face_print}} d = %d ---\n", d); }
    
    /* TESTING: void box_lo_corner(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_hi_corner(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_corner(box_dim_t d, interval_t B[], interval_side_t dir[], double p[]); */
    /* TESTING: void box_center(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_radii(box_dim_t d, interval_t B[], double h[]); */
    /* TESTING: void box_widths(box_dim_t d, interval_t B[], double w[]); */
    /* TESTING: double box_max_width(box_dim_t d, interval_t B[], bool_t verbose); */
    /* TESTING: double box_radius(box_dim_t d, interval_t B[], bool_t verbose); */
    
    fprintf(stderr, "testing {box_print}: ");
    box_print(stderr, d, B);
    fprintf(stderr, "\n");
    
    fprintf(stderr, "testing {box_gen_print}:\n");
    box_gen_print(stderr, d, B, "%+.4f", "B = ", " � ", "\n");
    
    if (d <= 4)
      { fprintf(stderr, "testing {box_face_print}:\n");
        uint32_t nf = (uint32_t)ipow(3,d);
        for (uint32_t fi = 0; fi < nf; fi++)
          { fprintf(stderr, "face %3d = ", fi); 
            box_face_print(stderr, d, fi);
            fprintf(stderr, "\n"); 
          }
      }
  }
  
void test_attribs(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_{lo_corner,hi_corner,corner,center,radii,widths,max_width,radius}} d = %d ---\n", d); }
    
    double ctr[d], rad[d], wid[d], plo[d], phi[d], pdir[d];
    double er, mw;

    interval_side_t dir[d];
    for (uint32_t i = 0;  i < d; i++) { dir[i] = (interval_side_t)int32_abrandom(0, 1); }
    
    /* TESTING: void box_lo_corner(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_hi_corner(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_corner(box_dim_t d, interval_t B[], interval_side_t dir[], double p[]); */
    /* TESTING: void box_center(box_dim_t d, interval_t B[], double p[]); */
    /* TESTING: void box_radii(box_dim_t d, interval_t B[], double rad[]); */
    /* TESTING: void box_widths(box_dim_t d, interval_t B[], double wid[]); */
    /* TESTING: double box_max_width(box_dim_t d, interval_t B[], bool_t verbose); */
    /* TESTING: double box_radius(box_dim_t d, interval_t B[], bool_t verbose); */
    
    if (! box_is_empty(d, B))
      { box_lo_corner(d, B, plo);
        box_hi_corner(d, B, phi);
        box_corner(d, B, dir, pdir);
        box_center(d, B, ctr);
      }
    box_radii(d, B, rad);
    box_widths(d, B, wid);
    mw = box_max_width(d, B);
    er = box_radius(d, B);
    
    if ((d == 0) || box_is_empty(d, B))
      { for (int32_t i = 0;  i < d; i++)
          { assert(wid[i] == -INF);
            assert(rad[i] == -INF);
            assert(mw == -INF);
          }
      }
    else
      { double mw_cmp = -INF;
        double sum2 = 0;
        for (uint32_t i = 0;  i < d; i++)
          { interval_t Bi = B[i];
            assert(LO(Bi) == plo[i]);
            assert(HI(Bi) == phi[i]);
            assert(Bi.end[dir[i]] == pdir[i]);
            assert(interval_mid(&Bi) == ctr[i]);
            assert(interval_rad(&Bi) == rad[i]);
            assert(fabs(wid[i] - 2*rad[i]) < 1.0e-13);
            mw_cmp = fmax(mw_cmp, wid[i]);
            sum2 += rad[i]*rad[i];
          }
        assert(fabs(mw_cmp - mw) < 1.0e-12);
        assert(fabs(er - sqrt(sum2)) < 1.0e-12);
      }
  }
     
void test_dimension__measure(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_dimension,box_measure} d = %d ---\n", d); }

    uint32_t dim = 0;
    double mes_hi = 1.0;
    double mes_lo = 1.0;
    double eps = 1.0e-14;
    double tiny = 1.0e-320;
    for (uint32_t i = 0; i < d; i++)
      { double lo = LO(B[i]), hi = HI(B[i]);
        if (lo > hi) { dim = 0; mes_lo = 0.0; mes_hi = 0.0; break; }
        if (lo < hi)
          { dim++; 
            mes_lo = mes_lo * (1-eps)*(hi-lo); 
            mes_hi = mes_hi * (1+eps)*(hi-lo); 
            if (mes_hi < tiny) { mes_hi = tiny; }
          }
      }
    demand(dim == box_dimension(d, B), "{box_dimension} failed");
    assert(mes_lo <= mes_hi);
    if (mes_hi == 0)
      { /* Box must be empty: */
        assert(box_is_empty(d, B));
        demand(box_measure(d, B) == 0.0, "empty box has nonzero {box_measure}");
      }
    else
      { /* Box must be non-empty: */
        assert(! box_is_empty(d, B));
        double mes = box_measure(d, B);
        demand(mes != 0, "non-empty box has zero {box_measure}");
        demand((mes_lo <= mes) && (mes <= mes_hi), "{box_measure} failed");
      }
  }

void test_equal(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_equal} d = %d ---\n", d); }

    interval_t A[d]; box_throw(d, GLOB_LO, GLOB_HI, 0.100, 0.100, A);

    assert(box_equal(d, B, B));
    if (box_is_empty(d, A) && box_is_empty(d, B))
      { assert(box_equal(d, A, B)); }
    else if (box_is_empty(d, A) || box_is_empty(d, B))
      { assert(! box_equal(d, A, B)); }
    else
      { bool_t eq = TRUE;
        for (uint32_t i = 0;  i < d; i++)
          { interval_t Ai = A[i];
            interval_t Bi = B[i];
            if (LO(Ai) != LO(Bi)) { eq = FALSE; }
            if (HI(Ai) != HI(Bi)) { eq = FALSE; }
          }
        assert(eq == box_equal(d, A, B));
      }
  }
  
void test_from_center_and_radii(box_dim_t d, interval_t B[], bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- testing {box_from_center_and_radii(} d = %d ---\n", d); }

    /* TESTING: void box_from_center_and_radii(box_dim_t d, double ctr[], double rad[], interval_t B[]);  */
    double ctr[d], rad[d];
    bool_t ety = FALSE;
    for (uint32_t i = 0;  i < d; i++)
      { ctr[i] = drandom();
        if (drandom() < 0.1)
          { rad[i] = -1; ety = TRUE; }
        else if (drandom() < 0.1)
          { rad[i] = 0; }
        else if (drandom() < 0.1)
          { rad[i] = +INF; }
        else
          { rad[i] = fabs(drandom()); }
      }
    box_from_center_and_radii(d, ctr, rad, B);
    
    for (uint32_t i = 0;  i < d; i++)
      { interval_t Bi = B[i];
        if (ety)
          { demand(interval_is_empty(&Bi), "box should be empty"); }
        else
          { double Bctri, Bradi;
            interval_mid_rad(&Bi, &Bctri, &Bradi);
            if ((! isnan(Bctri)) && (! isnan(Bradi)))
              { demand ((Bradi == +INF) == (rad[i] == +INF), "inconsistent infinities");
                if (isfinite(Bradi))
                  { double mag = 1.0e-200;
                    mag += fabs(Bradi) + fabs(rad[i]);
                    mag += fabs(Bctri) + fabs(ctr[i]);
                    demand(fabs(Bradi - rad[i]) <= 1.0e-8*mag, "radius mismatch");
                    demand(fabs(Bctri - ctr[i]) <= 1.0e-8*mag, "center mismatch");
                  }
              }
          }
      }
  }
  
void test_include_point(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_include_point} d = %d ---\n", d); }

    /* TESTING: void box_include_point(box_dim_t d, interval_t B[], double p[], interval_t C[]); */

    interval_t C[d];
    double p[d];
    for (uint32_t i = 0;  i < d; i++) { p[i] = dabrandom(GLOB_LO, GLOB_HI); }

    box_include_point(d, B, p, C);
    
    for (uint32_t i = 0;  i < d; i++)
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

void test_join__meet(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_join,box_meet} d = %d ---\n", d); }

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
        for (uint32_t i = 0;  i < d; i++)
          { interval_t Ai = A[i];
            interval_t Bi = B[i];
            interval_t CJi = CJ[i];
            assert(LO(CJi) == fmin(LO(Ai), LO(Bi)));
            assert(HI(CJi) == fmax(HI(Ai), HI(Bi)));
          }
        /* Check whether meet should be empty: */
        double mety = FALSE;
        for (uint32_t i = 0;  i < d; i++)
          { interval_t Ai = A[i];
            interval_t Bi = B[i];
            double cmlo = fmax(LO(Ai), LO(Bi));
            double cmhi = fmin(HI(Ai), HI(Bi));
            if (cmlo > cmhi) { mety = TRUE; }
          }
        if (mety)
          { assert(box_is_empty(d, CM)); }
        else
          { for (uint32_t i = 0;  i < d; i++)
              { interval_t Ai = A[i];
                interval_t Bi = B[i];
                interval_t CMi = CM[i];
                assert(LO(CMi) == fmax(LO(Ai), LO(Bi)));
                assert(HI(CMi) == fmin(HI(Ai), HI(Bi)));
              }
          }
      }
  }

void test_disjoint__contained(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_disjoint,box_contained} d = %d ---\n", d); }

    /* TESTING: void box_disjoint(box_dim_t d, interval_t A[], interval_t B[], interval_t C[]); */
    /* TESTING: void box_contained(box_dim_t d, interval_t A[], interval_t B[], interval_t C[]); */
    
    auto void throw_cute_interval(double loR, double hiR, double *loS_P, double *hiS_P, bool_t *disj_P, bool_t *cont_P);
      /* Given a non-empty interval {R=[loR _ hiR]}, throws a non-empty interval {S=[loS _ hiS]}
        that is likely to have interesting relation to {R}, such as overlap, containment,
        and coincident borders.  Also returns booleans {disj} that tells if {R} and {S} are disjoint and
        {cont} that tells if {S} is contained in {R}. (These take into account the convention
        that a non-empty interval is closed iff it is a singleton.) */
    
    uint32_t nt = 100; /* Number of cases to try. */
    for (uint32_t it = 0; it < nt; it++)
      { /* Create box {A} that is in interesting postion relative to {B}: */
        interval_t A[d];
        bool_t disjoint;  /* True if {A,B} are disjoint along some axis. */
        bool_t contained;  /* True if {A} is contained in {B} along every axis. */
        if (d == 0)
          { /* The boxes are always identical and non-empty: */
            disjoint = FALSE;
            contained = TRUE;
          }
        else if (drandom() < 0.10)
          { /* Set {A} to an empty box: */
            box_empty(d, A);
            /* No matter what {B} is: */
            disjoint = TRUE;
            contained = TRUE;
          } 
        else if (box_is_empty(d, B))
          { /* Set {A} to random non-empty box: */
            for (int32_t i = 0; i < d; i++) 
              { LO(A[i]) = 3.14; HI(A[i]) = LO(A[i]) + (drandom() < 0.5 ? 0.0 : 1.0); }
            /* Since {B} is empty and {A} is non-empty: */
            disjoint = TRUE;
            contained = FALSE;
          }
        else 
          { /* Set {A} to random non-empty box likely in special relation to {B}: */
            disjoint = FALSE;
            contained = TRUE;
            for(int32_t i = 0; i < d; i++)
              { bool_t disj_i, cont_i;
                throw_cute_interval(LO(B[i]), HI(B[i]), &(LO(A[i])), &(HI(A[i])), &disj_i, &cont_i);
                if (disj_i) { disjoint = TRUE; }
                if (! cont_i) { contained = FALSE; }
              }
          }
        if (verbose && (it < 3))
          { fprintf(stderr, "testing with");
            box_gen_print(stderr, d, A, "%+24.16e", " A = ", " x ", "");
            box_gen_print(stderr, d, B, "%+24.16e", " B = ", " x ", "");
            fprintf(stderr, " expecting disjoint = %c  contained = %c\n", "FT"[disjoint], "FT"[contained]);
          }
        demand(disjoint == box_disjoint(d, A, B), "box_disjoint failed");
        demand(contained == box_contained(d, A, B), "box_contained failed");
      }
    
    void throw_cute_interval(double loR, double hiR, double *loS_P, double *hiS_P, bool_t *disj_P, bool_t *cont_P)
      { assert(loR <= hiR);
        /* Get some interesting values for {loS,hiS}: */
        uint32_t nx = 0;
        double x[8]; /* Interesting values for {loS,hiS} are {x[0..nc-1]} */
        x[nx] = loR - 2; nx++;
        x[nx] = loR - 1; nx++;
        x[nx] = loR; nx++;
        if (loR < hiR) 
          { x[nx] = 0.75*loR + 0.25*hiR; nx++;
            x[nx] = 0.25*loR + 0.75*hiR; nx++;
            x[nx] = hiR; nx++;
          }
        x[nx] = hiR + 1; nx++;
        x[nx] = hiR + 2; nx++;
        assert(nx <= 8);
        
        /* Now choose randomly a pair of interesting values */
        uint32_t i = uint32_abrandom(0, nx-1);
        uint32_t j = uint32_abrandom(0, nx-1);
        if (i > j) { uint32_t t = i; i = j; j = t; }
        double loS = x[i];
        double hiS = x[j];
        assert(loS <= hiS);
        bool_t disj;
        bool_t cont;
        if (loR < hiR)
          { /* {R} is open */
            if (loS < hiS)
              { /* {R} and {S} are open: */
                disj = (hiS <= loR) || (loS >= hiR);
                cont = (loS >= loR) && (hiS <= hiR);
              }
            else
              { /* {R} is open, {S} is singleton: */
                disj = (hiS <= loR) || (loS >= hiR);
                cont = (loS > loR) && (hiS < hiR);
              }
          }
        else
          { /* {R} is singleton: */
            if (loS < hiS)
              { /* {R} is singleton, {S} is open: */
                disj = (hiS <= loR) || (loS >= hiR);
                cont = FALSE;
              }
            else
              { /* {R} and {S} are both singletons: */
                disj = (loS != loR);
                cont = (loS == loR);
              }
          }
        (*loS_P) = loS;
        (*hiS_P) = hiS;
        (*disj_P) = disj;
        (*cont_P) = cont;
      }
  }

void test_has_point(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_has_point} d = %d ---\n", d); }
    
    auto void throw_cute_value(double loR, double hiR, double *v_P, bool_t *inside_P);
      /* Given a non-empty interval {R=[loR _ hiR]}, throws a value {v}
        that is likely to have interesting relation to {R}, such as
        overlap, containment, and coincident borders. Also returns
        boolean {inside} that tells if {v} is inside {R}, taking into
        account the convention that a non-empty interval is closed iff
        it is a singleton. */
    
    uint32_t nt = 100; /* Number of cases to try. */
    for (uint32_t it = 0; it < nt; it++)
      { /* Create point {p} that is in interesting postion relative to {B}: */
        double p[d];
        bool_t inside;  /* True if {p} is inside {B}. */
        if (d == 0)
          { /* The box {B} is not empty and contains {p}: */
            inside = TRUE;
          }
        else if (box_is_empty(d, B))
          { /* The box {B} is empty: */
            inside = FALSE;
          }
        else 
          { /* Set {p} to random point likely in special relation to {B}: */
            inside = TRUE;
            for(int32_t i = 0; i < d; i++)
              { bool_t inside_i;
                throw_cute_value(LO(B[i]), HI(B[i]), &(p[i]), &inside_i);
                if (! inside_i) { inside = FALSE; }
              }
          }
        if (verbose && (it < 3))
          { fprintf(stderr, "testing with");
            tbx_print_point(stderr, d, p, " p = ( ", "%+24.16e", " ", " )");
            box_gen_print(stderr, d, B, "%+24.16e", " B = ", " x ", "");
            fprintf(stderr, " expecting inside = %c\n", "FT"[inside]);
          }
        demand(inside == box_has_point(d, B, p), "box_has_point failed");
      }
    
    void throw_cute_value(double loR, double hiR, double *v_P, bool_t *inside_P)
      { assert(loR <= hiR);
        /* Get some interesting values for {v}: */
        uint32_t nx = 0;
        double x[5]; /* Interesting values for {loS,hiS} are {x[0..nx-1]} */
        x[nx] = loR - 1; nx++;
        x[nx] = loR; nx++;
        if (loR < hiR) 
          { x[nx] = 0.75*loR + 0.25*hiR; nx++;
            x[nx] = hiR; nx++;
          }
        x[nx] = hiR + 1; nx++;
        assert(nx <= 5);
        
        /* Now choose randomly an interesting value */
        uint32_t i = uint32_abrandom(0, nx-1);
        double v = x[i];
        bool_t inside;
        if (loR < hiR)
          { /* {R} is open */
            inside = ((loR < v) && (v < hiR));
          }
        else if (loR == hiR)
          { /* {R} is singleton: */
            inside = (v == loR);
          }
        else
          { inside = FALSE; }
        (*v_P) = v;
        (*inside_P) = inside;
      }
  }

void test_split(box_dim_t d, interval_t B[], bool_t verbose)
  {
    bool_t debug = FALSE;
    
    if (verbose || debug) { fprintf(stderr, "--- testing {box_split} d = %d ---\n", d); }
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
        else if (alo == ahi)
          { /* Pick the singleton value: */
            x = alo;
          }
        else
          { /* Pick a value inside the range: */
            x = dabrandom(alo, ahi);
          }
      }
      
    if (debug)
      { box_gen_print(stderr, d, B, "%+.6f", "B = \n  ", "\n  ", "\n");
        fprintf(stderr, "  splitting axis a = %d by x = %12.6f\n", a, x);
      }

    box_split(d, B, a, x, BLO, BMD, BHI);
    
    if (ety)
      { assert(box_is_empty(d, BLO));
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
      
        for (uint32_t i = 0;  i < d; i++)
          { interval_t Bi = B[i];
            interval_t BLOi = BLO[i];
            interval_t BMDi = BMD[i];
            interval_t BHIi = BHI[i];
            if (debug)
              { interval_gen_print(stderr, &Bi,   "%12.6f", "    Bi =   ( ", " ", " )\n");
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

void test_widen__round(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_{widen,widen_var,round}} d = %d ---\n", d); }

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
                
    if (verbose) { fprintf(stderr, "!!! NOT TESTED !!!\n"); }
  }
                  
void test_shift__unshift(box_dim_t d, interval_t B[], bool_t verbose)
  { if (verbose) { fprintf(stderr, "--- testing {box_{shift,unshift}} d = %d ---\n", d); }
    interval_t C[d];
    interval_t D[d];
    double v[d];
    for (int32_t i = 0; i < d; i++) { v[i] = dabrandom(-5.0, +5.0); }
    box_shift(d, B, v, C);
    box_unshift(d, C, v, D);
    if (box_is_empty(d, B))
      { demand(box_is_empty(d, C), "box_shift of empty box is not empty");
        demand(box_is_empty(d, D), "box_unshift of empty box is not empty");
      }
    else
      { demand(! box_is_empty(d, C), "box_shift of non-empty box is empty");
        demand(! box_is_empty(d, D), "box_unshift of non-empty box is empty");
        for (int32_t i = 0; i < d; i++)
          { for (int32_t ke = 0; ke <= 1; ke++)
              { demand(C[i].end[ke] == B[i].end[ke] + v[i], "box_shift bug");
                demand(D[i].end[ke] == C[i].end[ke] - v[i], "box_unshift bug");
              }
          }
      } 
  }

void test_box_map(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_{point,box}_{mpa,unmap)}} d = %d ---\n", d); }
      
    /* TESTING: void box_point_map(box_dim_t d, double z[], interval_t B[], double x[]); */
    /* TESTING: void box_point_unmap(box_dim_t d, double x[], interval_t B[], double z[]); */
    /* TESTING: void box_box_map(box_dim_t d, interval_t Z[], interval_t B[], interval_t X[]); */
    /* TESTING: void box_box_unmap(box_dim_t d, interval_t X[], interval_t B[], interval_t Z[]); */
                
    if (verbose) { fprintf(stderr, "!!! NOT TESTED !!!\n"); }

  }
                  
void test_face(box_dim_t d, interval_t B[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- testing {box_face_{dimension,position,signature,rel_box,print}} d = %d ---\n", d); }
      
    /* TESTING: box_dim_t box_face_dimension(box_dim_t m, box_face_index_t fi); */
    /* TESTING: box_signed_dir_t box_face_position(box_face_index_t fi, box_axis_index_t j); */
    /* TESTING: void box_face_signature(box_dim_t m, box_face_index_t fi, box_signed_dir_t dir[]); */
    /* TESTING: void box_face_rel_box(box_dim_t m, box_face_index_t fi, interval_t F[]); */
    /* TESTING: void box_face_print(FILE *wr, box_dim_t m, box_face_index_t fi); */
                
    if (verbose) { fprintf(stderr, "!!! NOT TESTED !!!\n"); }

  }

void tbx_print_point(FILE *wr, box_dim_t d, double p[], char *pref, char *fmt, char *sep, char *suff)
  { 
    if (sep == NULL) { sep = " "; }
    if (pref != NULL) { fputs(pref, wr); }
    for (uint32_t i = 0; i < d; i++)
      { if (i > 0) { fputs(sep, wr); }
        fprintf(wr, fmt, p[i]);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }
        
