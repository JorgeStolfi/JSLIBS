/* See box.h */
/* Last edited on 2024-11-15 19:49:02 by stolfi */ 

/* We need to set these in order to get {asinh}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <stdint.h>
#include <stdlib.h>
#include <fenv.h>
#include <fpu_control.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <interval_io.h>
#include <set32.h>
#include <jsrandom.h>
#include <box.h>

/* BOXES */

bool_t box_is_empty(box_dim_t d, interval_t B[])
  {
    for (uint32_t i = 0;  i < d; i++) 
      { if (interval_is_empty(&(B[i]))) { return TRUE; } }
    return FALSE;
  }

void box_lo_corner(box_dim_t d, interval_t B[], double p[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { interval_t *Bi = &(B[i]);
        demand(! interval_is_empty(Bi), "empty box");
        p[i] = LO(*Bi);
      }
  }

void box_hi_corner(box_dim_t d, interval_t B[], double p[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { interval_t *Bi = &(B[i]);
        demand(! interval_is_empty(Bi), "empty box");
        p[i] = HI(*Bi); 
      }
  }

void box_corner(box_dim_t d, interval_t B[], interval_side_t dir[], double p[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { interval_t *Bi = &(B[i]);
        demand(! interval_is_empty(Bi), "empty box");
        p[i] = B[i].end[dir[i]];
      }
  }

void box_center(box_dim_t d, interval_t B[], double p[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { interval_t *Bi = &(B[i]);
        demand(! interval_is_empty(Bi), "empty box");
        p[i] = interval_mid(Bi);
      }
  }

void box_half_widths(box_dim_t d, interval_t B[], double h[])
  { 
    if (box_is_empty(d, B))
      { for (uint32_t i = 0;  i < d; i++) { h[i] = 0.0; } }
    else
      { for (uint32_t i = 0;  i < d; i++) 
          { interval_t *Bi = &(B[i]);
            h[i] = interval_rad(Bi);
          }
      }
  }

void box_widths(box_dim_t d, interval_t B[], double w[])
  { if (box_is_empty(d, B))
      { for (uint32_t i = 0;  i < d; i++) { w[i] = 0.0; } }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        for (uint32_t i = 0;  i < d; i++) { w[i] = HI(B[i]) - LO(B[i]); }
        fesetround(oround);
      }
  }

double box_max_width(box_dim_t d, interval_t B[])
  { if ((d == 0) || (box_is_empty(d, B)))
      { return 0.0; }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        double sz = 0;
        for (uint32_t i = 0;  i < d; i++) 
          { double szi = HI(B[i]) - LO(B[i]);
            if (szi > sz) { sz = szi; }
          }
        fesetround(oround);
        return sz;
      }
  }

double box_radius(box_dim_t d, interval_t B[])
  { if ((d == 0) || (box_is_empty(d, B)))
      { return 0.0; }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        double sum2 = 0;
        for (uint32_t i = 0;  i < d; i++) 
          { double rdi = (HI(B[i]) * 0.5) - (LO(B[i]) * 0.5);
            sum2 += rdi*rdi;
          }
        double rad = sqrt(sum2);
        fesetround(oround);
        return rad;
      }
  }

bool_t box_equal(box_dim_t d, interval_t A[], interval_t B[])
  { if (box_is_empty(d, A))
      { return box_is_empty(d, B); }
    else if (box_is_empty(d, B))
      { return box_is_empty(d, A); }
    else 
      { for (uint32_t i = 0;  i < d; i++) 
          { interval_t *Ai = &(A[i]);
            interval_t *Bi = &(B[i]);
            if ((LO(*Ai) != LO(*Bi)) || (HI(*Ai) != HI(*Bi)))
              { return FALSE; }
          }
        return TRUE;
      }
  }

/* BOX CREATION AND MODIFICATION */

void box_empty(box_dim_t d,  interval_t C[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { C[i] = (interval_t){{ +INF, -INF }}; }
  }

void box_include_point(box_dim_t d, interval_t B[], double p[], interval_t C[])
  { 
     if (box_is_empty(d, B))
      { for (uint32_t i = 0;  i < d; i++) { C[i] = (interval_t){{ p[i], p[i] }}; } }
    else
      { for (uint32_t i = 0;  i < d; i++) 
          { interval_t *Bi = &(B[i]);
            interval_t *Ci = &(C[i]);
            double z = p[i];
            LO(*Ci) = fmin(LO(*Bi), z);
            HI(*Ci) = fmax(HI(*Bi), z);
          }
      }
    return;
  }
  
void box_join(box_dim_t d, interval_t A[], interval_t B[], interval_t C[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { interval_t *Ai = &(A[i]);
        if (LO(*Ai) > HI(*Ai)) { for (uint32_t i = 0;  i < d; i++) { C[i] = B[i]; } return ;}
        interval_t *Bi = &(B[i]);
        if (LO(*Bi) > HI(*Bi)) { for (uint32_t i = 0;  i < d; i++) { C[i] = A[i]; } return ;}
        interval_t *Ci = &(C[i]);
        LO(*Ci) = fmin(LO(*Ai), LO(*Bi));
        HI(*Ci) = fmax(HI(*Ai), HI(*Bi));
      }
    return;
  }

void box_meet(box_dim_t d, interval_t A[], interval_t B[], interval_t C[])
  { 
    for (uint32_t i = 0;  i < d; i++) 
      { interval_t *Ai = &(A[i]);
        if (LO(*Ai) > HI(*Ai)) { box_empty(d, C); return ;}
        interval_t *Bi = &(B[i]);
        if (LO(*Bi) > HI(*Bi)) { box_empty(d, C); return ;}
        interval_t *Ci = &(C[i]);
        LO(*Ci) = fmax(LO(*Ai), LO(*Bi));
        HI(*Ci) = fmin(HI(*Ai), HI(*Bi));
        if (LO(*Ci) > HI(*Ci)) { box_empty(d, C); return; }
      }
    return;
  }

void box_widen(box_dim_t d, interval_t B[], double margin, interval_t C[])
  {
    if (box_is_empty(d, B))
      { box_empty(d, C); }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        for (uint32_t i = 0;  i < d; i++) 
          { interval_t *Bi = &(B[i]);
            interval_t *Ci = &(C[i]);
            double t = -LO(*Bi) + margin; /* To round down. */
            LO(*Ci) = -t;
            HI(*Ci) = HI(*Bi) + margin;
          }
        fesetround(oround);
      }
    return;
  }
  
void box_widen_var(box_dim_t d, interval_t B[], double mrglo[], double mrghi[], interval_t C[])
  {
    if (box_is_empty(d, B))
      { box_empty(d, C); }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        for (uint32_t i = 0;  i < d; i++) 
          { interval_t *Bi = &(B[i]);
            interval_t *Ci = &(C[i]);
            double t = -LO(*Bi) + mrglo[i]; /* To round down. */
            LO(*Ci) = -t;
            HI(*Ci) = HI(*Bi) + mrghi[i];
          }
        fesetround(oround);
      }
    return;
  }

void box_round(box_dim_t d, interval_t B[], double unit, bool_t out, interval_t C[])
  {
    demand(unit > 0, "invalid {unit}");
    if (box_is_empty(d, B))
      { box_empty(d, C); }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        for (uint32_t i = 0;  i < d; i++) 
          { interval_t *Bi = &(B[i]);
            interval_t *Ci = &(C[i]);
            if (out)
              { LO(*Ci) = -(unit*ceil((-LO(*Bi))/unit)); /* To get proper rounding. */
                HI(*Ci) = unit*ceil(LO(*Bi)/unit);
                assert(LO(*Ci) <= HI(*Ci));
              }
            else
              { LO(*Ci) = -(unit*floor((-LO(*Bi))/unit)); /* To get proper rounding. */
                HI(*Ci) = unit*floor(LO(*Bi)/unit);
                if (LO(*Ci) > HI(*Ci)) 
                  { box_empty(d, C);  break; }
              }
          }
        fesetround(oround);
      }
    return;
  }

void box_split
  ( box_dim_t d, 
    interval_t B[],  
    box_axis_t a,  
    double x,  
    interval_t BLO[],  
    interval_t BMD[],  
    interval_t BHI[] 
  )
  { bool_t debug = FALSE;
    
    demand((a >= 0) && (a < d), "invalid axis");
    if (box_is_empty(d, B))
      { if (BLO != NULL) { box_empty(d, BLO); }
        if (BMD != NULL) { box_empty(d, BMD); }
        if (BHI != NULL) { box_empty(d, BHI); }
      }
    else
      { /* Decide which parts are not empty: */
        interval_t *Ba = &(B[a]);
        assert(LO(*Ba) <= HI(*Ba));
        if ((LO(*Ba) < x) && (x < HI(*Ba)))
          { /* All three parts are non-empty: */
            for (uint32_t i = 0;  i < d; i++) 
              { if (i == a)
                  { if (BLO != NULL) { LO(BLO[i]) = LO(B[i]); HI(BLO[i]) = x; assert(LO(BLO[i]) < HI(BLO[i])); }
                    if (BMD != NULL) { LO(BMD[i]) = x; HI(BMD[i]) = x; }
                    if (BHI != NULL) { LO(BHI[i]) = x; HI(BHI[i]) = HI(B[i]); assert(LO(BHI[i]) < HI(BHI[i])); }
                  }
                else
                  { if (BLO != NULL) { BLO[i] = B[i]; }
                    if (BMD != NULL) { BMD[i] = B[i]; } 
                    if (BHI != NULL) { BHI[i] = B[i]; }
                  }
              }
          }
        else if ((LO(*Ba) == x) && (HI(*Ba) == x))
          { /* {BLO} and {BHI} are empty, {BMD} is {B}: */
            if (BLO != NULL) { box_empty(d, BLO); }
            if (BMD != NULL) { for (uint32_t i = 0;  i < d; i++) { BMD[i] = B[i]; } }
            if (BHI != NULL) { box_empty(d, BHI); }
          }
        else if (x <= LO(*Ba))
          { /* {BLO} and {BMD} are empty, {BHI} is {B}: */
            if (BLO != NULL) { box_empty(d, BLO); }
            if (BMD != NULL) { box_empty(d, BMD); }
            if (BHI != NULL) { for (uint32_t i = 0;  i < d; i++) { BHI[i] = B[i]; } }
          }
        else if (x >= HI(*Ba))
          { /* {BMD} {BHI} are empty, {BLO} is {B}: */
            if (BLO != NULL) { for (uint32_t i = 0;  i < d; i++) { BLO[i] = B[i]; } }
            if (BMD != NULL) { box_empty(d, BMD); }
            if (BHI != NULL) { box_empty(d, BHI); }
          }
        else
          { assert(FALSE); }
      }
    if (debug)
      { box_gen_print(stderr, d, B, "%12.6f", "B =   ", " × ", "\n");
        if (BLO != NULL) { box_gen_print(stderr, d, BLO, "%12.6f", "BLO = ", " × ", "\n"); }
        if (BMD != NULL) { box_gen_print(stderr, d, BMD, "%12.6f", "BMD = ", " × ", "\n"); }
        if (BHI != NULL) { box_gen_print(stderr, d, BHI, "%12.6f", "BHI = ", " × ", "\n"); }
      }
  }

void box_throw(box_dim_t d, double elo, double ehi, double p_empty, double p_single, interval_t B[])
  {
    if (d == 0) { return; }
    if (drandom() < p_empty)
      { /* Generate an empty box: */
        box_empty(d, B); 
        assert(box_is_empty(d, B));
      }
    else
      { for (uint32_t i = 0;  i < d; i++)
          { if (drandom() < p_single)
              { /* Generate a singleton interval: */
                double x = NAN;
                do
                  { x = dabrandom(elo, ehi); }
                while ((x == elo) || (x == ehi));
                B[i] = (interval_t){{ x, x }};
              }
            else
              { /* Generate an open interval: */
                double x = NAN, y = NAN;
                do
                  { x = dabrandom(elo, ehi);
                    y = dabrandom(elo, ehi);
                    if (x > y) { double t = x; x = y; y = t; }
                  }
                while (x == y);
                B[i] = (interval_t){{ x, y }};
              }
          }
        assert(! box_is_empty(d, B));
      }
  }


void box_point_map(box_dim_t d, double z[], interval_t B[], double x[])
  { 
    int32_t j = 0;
    for (uint32_t i = 0;  i < d; i++)
      { double lo = LO(B[i]);
        double hi = HI(B[i]);
        if (lo < hi)
          { double wd = hi - lo;
            x[i] = lo + wd*z[j];
            j++;
          }
        else
          { demand(lo == hi, "empty box"); 
            x[i] = lo;
          }
      }
  }

void box_point_unmap(box_dim_t d, double x[], interval_t B[], double z[])
  { int32_t j = 0;
    for (uint32_t i = 0;  i < d; i++)
      { double lo = LO(B[i]);
        double hi = HI(B[i]);
        if (lo < hi)
          { double wd = hi - lo;
            z[j] = (x[i] - lo)/wd;
            j++;
          }
        else
          { demand(lo == hi, "empty box"); }
      }
  }

void box_box_map(box_dim_t d, interval_t Z[], interval_t B[], interval_t X[])
  { int32_t i, j;
    for (i = 0, j = 0; i < d; i++)
      { double lo = LO(B[i]);
        double hi = HI(B[i]);
        if (lo < hi)
          { double wd = hi - lo;
            LO(X[i]) = lo + wd*LO(Z[j]);
            HI(X[i]) = lo + wd*HI(Z[j]);
            j++;
          }
        else
          { demand(lo == hi, "empty box"); 
            LO(X[i]) = HI(X[i]) = lo;
          }
      }
  }

void box_box_unmap(box_dim_t d, interval_t X[], interval_t B[], interval_t Z[])
  { int32_t i, j;
    for (i = 0, j = 0; i < d; i++)
      { double lo = LO(B[i]);
        double hi = HI(B[i]);
        if (lo < hi)
          { double wd = hi - lo;
            LO(Z[j]) = (LO(X[i]) - lo)/wd;
            HI(Z[j]) = (HI(X[i]) - lo)/wd;
            j++;
          }
        else
          { affirm(lo == hi, "empty box"); }
      }
  }

/* RELATIVE FACE INDICES */

box_dim_t box_face_dimension(box_dim_t m, box_face_index_t fi)
  { while (fi != 0)
      { box_signed_dir_t dirj = ((fi + 1) % 3) - 1; 
        if (dirj != 0) { m--; }
        fi /= 3;
      }
    return m;
  }

box_signed_dir_t box_face_position(box_face_index_t fi, box_axis_index_t j)
  { while (j > 0) { j--; fi /= 3; }
    return ((fi + 1) % 3) - 1; 
  }

void box_face_signature(box_dim_t m, box_face_index_t fi, box_signed_dir_t dir[])
  { box_axis_index_t j;
    for (j = 0; j < m; j++)
      { dir[j] = ((fi + 1) % 3) - 1; 
        fi /= 3;
      }
  }
 
void box_face_rel_box(box_dim_t m, box_face_index_t fi, interval_t F[])
  { box_axis_index_t j;
    for (j = 0; j < m; j++)
      { box_signed_dir_t dirj = ((fi + 1) % 3) - 1; 
        double lo = 0.0, hi = 1.0;
        if (dirj == -1)
          { hi = 0.0; }
        else if (dirj == +1)
          { lo = 1.0; }
        F[j] = (interval_t){{lo,hi}};
        fi /= 3;
      }
  }

void box_print(FILE *wr, box_dim_t d, interval_t B[])
  { 
    box_gen_print(wr, d, B, NULL, NULL, NULL, NULL); 
  }

void box_gen_print
  ( FILE *wr, 
    box_dim_t d, 
    interval_t B[], 
    char *fmt, 
    char *pref, 
    char *sep, 
    char *suff
  )
  { if (fmt == NULL) { fmt = "%16.8e"; }
    if (sep == NULL) { sep = " "; }

    if (pref != NULL) { fputs(pref, wr); }
    bool_t ety = box_is_empty(d, B);
    for (uint32_t i = 0;  i < d; i++) 
      { if (i > 0) { fputs(sep, wr); }
        if (ety) 
          { fputs("[]", wr); }
        else
          { interval_t *Bi = &(B[i]);
            if (LO(*Bi) == HI(*Bi))
              { /* Singleton: */
                fputs("{", wr);
                fprintf(wr, fmt, LO(*Bi));
                fputs("}", wr);
              }
            else
              { /* Open interval: */
                interval_gen_print(wr, Bi, fmt, "(", "_", ")");
              }
          }
      }
    if (suff != NULL) { fputs(suff, wr); }
  }

void box_face_print(FILE *wr, box_dim_t m, box_face_index_t fi)
  { box_axis_index_t j;
    for (j = 0; j < m; j++)
      { uint32_t dj = (fi % 3); 
        fprintf(stderr, "%c", ("0+-")[dj]);  
        fi /= 3;
      }
  }
