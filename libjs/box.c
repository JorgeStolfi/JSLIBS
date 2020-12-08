/* See box.h */
/* Last edited on 2017-08-02 10:46:38 by jstolfi */ 

#include <affirm.h>
#include <interval.h>
#include <set32.h>
#include <box.h>

#include <stdlib.h>
#include <fenv.h>
#include <fpu_control.h>
#include <math.h>

/* BOXES */

void box_lo_corner(box_dim_t d, interval_t B[], double p[])
  { int i;
    for (i = 0; i < d; i++) { p[i] = LO(B[i]); }
  }

void box_hi_corner(box_dim_t d, interval_t B[], double p[])
  { int i;
    for (i = 0; i < d; i++) { p[i] = HI(B[i]); }
  }

void box_corner(box_dim_t d, interval_t B[], interval_side_t dir[], double p[])
  { int i;
    for (i = 0; i < d; i++) { p[i] = B[i].end[dir[i]]; }
  }

void box_center(box_dim_t d, interval_t B[], double p[])
  { int i;
    for (i = 0; i < d; i++) { p[i] = interval_mid(&(B[i])); }
  }

void box_half_widths(box_dim_t d, interval_t B[], double h[])
  { int i;
    for (i = 0; i < d; i++) { h[i] = interval_rad(&(B[i])); }
  }

void box_widths(box_dim_t d, interval_t B[], double w[])
  { int oround = fegetround();
    fesetround(FE_UPWARD);
    int i;
    for (i = 0; i < d; i++) { w[i] = HI(B[i]) - LO(B[i]); }
    fesetround(oround);
  }

double box_max_width(box_dim_t d, interval_t B[])
  { int oround = fegetround();
    fesetround(FE_UPWARD);
    int i;
    double sz = 0;
    for (i = 0; i < d; i++) 
      { double szi = HI(B[i]) - LO(B[i]);
        if (szi > sz) { sz = szi; }
      }
    fesetround(oround);
    return sz;
  }

double box_radius(box_dim_t d, interval_t B[])
  { int oround = fegetround();
    fesetround(FE_UPWARD);
    int i;
    double sum2 = 0;
    for (i = 0; i < d; i++) 
      { double rdi = (HI(B[i]) * 0.5) - (LO(B[i]) * 0.5);
        sum2 += rdi*rdi;
      }
    double rad = sqrt(sum2);
    fesetround(oround);
    return rad;
  }

void box_join(box_dim_t d, interval_t A[], interval_t B[], interval_t C[])
  { int i;
    for (i = 0; i < d; i++) { C[i] = interval_join(&(A[i]), &(B[i])); }
  }

void box_meet(box_dim_t d, interval_t A[], interval_t B[], interval_t C[])
  { int i;
    for (i = 0; i < d; i++) { C[i] = interval_meet(&(A[i]), &(B[i])); } 
  }

void box_split
  ( box_dim_t d, 
    interval_t B[],  
    box_axis_t a,  
    double x,  
    interval_t BLO[],  
    interval_t BHI[] 
  )
  { int i;
    for (i = 0; i < d; i++) 
      { if (i == a)
          { LO(BLO[i]) = LO(B[i]); 
            HI(BLO[i]) = x; 
            LO(BHI[i]) = x;
            HI(BHI[i]) = HI(B[i]);
          } 
        else 
          { BLO[i] = B[i]; BHI[i] = B[i]; }
      }
  }

void box_point_map(box_dim_t d, double z[], interval_t B[], double x[])
  { int i, j;
    for (i = 0, j = 0; i < d; i++)
      { double lo = LO(B[i]);
        double hi = HI(B[i]);
        if (lo < hi)
          { double wd = hi - lo;
            x[i] = lo + wd*z[j];
            j++;
          }
        else
          { affirm(lo == hi, "empty box"); 
            x[i] = lo;
          }
      }
  }

void box_point_unmap(box_dim_t d, double x[], interval_t B[], double z[])
  { int i, j;
    for (i = 0, j = 0; i < d; i++)
      { double lo = LO(B[i]);
        double hi = HI(B[i]);
        if (lo < hi)
          { double wd = hi - lo;
            z[j] = (x[i] - lo)/wd;
            j++;
          }
        else
          { affirm(lo == hi, "empty box"); }
      }
  }

void box_box_map(box_dim_t d, interval_t Z[], interval_t B[], interval_t X[])
  { int i, j;
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
          { affirm(lo == hi, "empty box"); 
            LO(X[i]) = HI(X[i]) = lo;
          }
      }
  }

void box_box_unmap(box_dim_t d, interval_t X[], interval_t B[], interval_t Z[])
  { int i, j;
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
 
void box_face_rel_box(box_dim_t m, box_face_index_t fi, interval_t box[])
  { box_axis_index_t j;
    for (j = 0; j < m; j++)
      { box_signed_dir_t dirj = ((fi + 1) % 3) - 1; 
        double lo = 0.0, hi = 1.0;
        if (dirj == -1)
          { hi = 0.0; }
        else if (dirj == +1)
          { lo = 1.0; }
        box[j] = (interval_t){{lo,hi}};
        fi /= 3;
      }
  }

void box_face_print(FILE *wr, box_dim_t m, box_face_index_t fi)
  { box_axis_index_t j;
    for (j = 0; j < m; j++)
      { int dj = (fi % 3); fprintf(stderr, "%c", ("0+-")[dj]);  fi /= 3; }
  }
