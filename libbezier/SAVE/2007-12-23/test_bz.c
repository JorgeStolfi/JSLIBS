/* Last edited on 2007-12-23 22:44:45 by stolfi */
/* Test of bz_basic.h and bz_patch.h routines. */

#include <bz_basic.h>
#include <bz_patch.h>
#include <affirm.h>
#include <interval.h>
#include <box.h>
#include <pswr.h>
#include <stdio.h>
#include <math.h>

/* Dimension of grid's domain: */
#define D (2)  

/* Dimension of grid's range space {(X,Y,F(X,Y))}: */
#define BD (3)

/* Number of faces in a cell of dimension {D}, namely {3^D}: */
#define NF (9)

/* Flatness tolerance (mm): */
#define TOL_MM (0.5)

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

void test_bezier(bz_degree_t g, double ratio);
  /* Tests of plotting routines, for patches of degree {g} and
    splitting at the given {ratio}. */
    
void twiddle_bezier(bz_patch_t *b);
  /* Modifies all control points of {b} by different amounts. */
   
void print_point(FILE *wr, char *pref, double *x, bz_patch_ddim_t n, char *fmt, char *suff);
  /* Prints {x[0..d-1]} with format {fmt}, all bracketed by {pref} and {suff}. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { fprintf(stderr, "=== EASY TEST ===\n");
    test_bezier(1, 0.5);
    fprintf(stderr, "=== HARD TEST ===\n");
    test_bezier(2, 0.5);
    fprintf(stderr, "=== DONE ===\n");
    return(0);
  }

void test_bezier(bz_degree_t g, double ratio)
  {
    bz_patch_t b = bz_patch_make_test(D, BD, g);
    interval_t bbox[BD];
    box_axis_t ax;
    double x[D];
    double bx[BD];
    int i, j, f;
    int NS = 3; /* Sampling points per axis. */

    affirm(b.m == 2, "wrong dimension");
    bz_patch_print(stderr, &(b), "%6.2f");
    fprintf(stderr, "\n");

    bz_patch_compute_bbox(&b, bbox);
    
    fprintf(stderr, "\n--- testing bz_patch_eval ---\n\n");
    for (i = 0; i < NS; i++)
      { for (j = 0; j < NS; j++)
          { x[0] = ((double)i)/((double)NS-1);
            x[1] = ((double)j)/((double)NS-1);
            bz_patch_eval(&b, x, bx);
            print_point(stderr, "  x = ", x, D, "%6.2f", "");
            print_point(stderr, "  bx = ", bx, BD, "%6.2f", "\n");
          }
      }
    
    fprintf(stderr, "\n--- testing bz_patch_get_face ---\n\n");
    for (f = 0; f < NF; f++)
      { int fd = box_face_dimension(D, f);
        bz_patch_t t = bz_patch_uniform_new(fd, BD, g);
        box_signed_dir_t dir[D];
        fprintf(stderr, "extracting face %d (", f);
        box_face_print(stderr, D, f);
        fprintf(stderr, ") dimension = %d\n", fd);
        box_face_signature(D, f, dir);
        bz_patch_get_face(&b, dir, &t);
        bz_patch_print(stderr, &t, "%6.2f");
        fprintf(stderr, "\n");
      }
    
    fprintf(stderr, "\n--- testing bz_set_face ---\n\n");
    for (f = 0; f < NF; f++)
      { int fd = box_face_dimension(D, f);
        bz_patch_t t = bz_patch_uniform_new(fd, BD, g);
        box_signed_dir_t dir[D];
        fprintf(stderr, "extracting face %d (", f);
        box_face_print(stderr, D, f);
        fprintf(stderr, ") dimension = %d\n", fd);
        box_face_signature(D, f, dir);
        bz_patch_get_face(&b, dir, &t);
        bz_patch_print(stderr, &t, "%6.2f");
        fprintf(stderr, "setting face to its current shape (no-op), result =\n");
        bz_patch_set_face(&b, dir, &t);
        bz_patch_print(stderr, &b, "%6.2f");
        twiddle_bezier(&t);
        fprintf(stderr, "setting face to new shape =\n");
        bz_patch_print(stderr, &t, "%6.2f");
        fprintf(stderr, " result =\n");
        bz_patch_set_face(&b, dir, &t);
        bz_patch_print(stderr, &b, "%6.2f");
        fprintf(stderr, "\n");
      }
    
    fprintf(stderr, "\n--- testing bz_multiaffine_approx ---\n\n");
    { bz_patch_t t = bz_patch_uniform_new(D, BD, 1);
    double err;
      fprintf(stderr, "original patch = \n");
      bz_patch_print(stderr, &b, "%6.2f");
      fprintf(stderr, "\n");
      bz_patch_multiaffine_approx(&b, &t);
      fprintf(stderr, "approximation = \n");
      bz_patch_print(stderr, &t, "%6.2f");
      fprintf(stderr, "\n");
      err = bz_patch_multiaffine_error(&b, &t);
      fprintf(stderr, "error = %8.6f\n", err);
    }
    
    fprintf(stderr, "\n--- testing bz_split (ratio = %6.4f) ---\n\n", ratio);
    { fprintf(stderr, "original =\n");
      bz_patch_print(stderr, &b, "%6.2f");

      for (ax = 0; ax < D; ax++)
        { bz_patch_t bLO = bz_patch_uniform_new(D, BD, g);
          bz_patch_t bHI = bz_patch_uniform_new(D, BD, g);
          fprintf(stderr, "splitting across axis %d\n", ax);
          bz_patch_split(&b, ax, ratio, &bLO, &bHI);
          fprintf(stderr, "low half =\n");
          bz_patch_print(stderr, &bLO, "%6.2f");
          fprintf(stderr, "high half =\n");
          bz_patch_print(stderr, &bHI, "%6.2f");
          fprintf(stderr, "\n");
        }
    }
  }
  
void twiddle_bezier(bz_patch_t *b)
  {
    int bstep[bz_patch_MAX_DDIM], bsz;
    int j, k, ix;
    double *c = &(b->c[0]);
    bz_patch_compute_steps(b, bstep, &bsz);
    double sc = ((double)bsz*b->n);
    k = 0;
    for (ix = 0; ix < bsz; ix++)
      { for (j = 0; j < b->n; j++)
        { c[k] += (k+0.5)/sc; k++; }
      }
  }

void print_point(FILE *wr, char *pref, double *x, bz_patch_ddim_t n, char *fmt, char *suff)
  { int j;
    fprintf(wr, pref);
    for (j = 0; j < n; j++,x++)
      { fprintf(wr, fmt, (*x)); if (j != 0) { fprintf(wr, " "); } }
    fprintf(wr, suff);
  }
