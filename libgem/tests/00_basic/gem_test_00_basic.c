#define PROG_NAME "gem_test_00_basic"
#define PROG_DESC "basic tests of the {oct.h} procedures"
#define PROG_VERS "1.0"

/* Last edited on 2015-01-01 23:54:43 by stolfilocal */ 

#define PROG_COPYRIGHT \
  "Copyright © 2014  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2014-07-11"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <affirm.h>

#include <gem.h>

void gem_test_new_step_splice_free(int d);
void gem_test_data(int d);
/* void gem_test_traverse(int d); */
void gem_test_write_read(char *name, int d, gem_ref_t m);

gem_ref_t gem_test_make_square(int d, bool_t rev1, bool_t rev2);
  /* Builds a {d}-dimensional cone of a 2-dimensional gem consisting of a barycentric division of a square,
    with opposite sides identified with parallel or antiparallel orientations dependinn on {rev1,rev2}. */

gem_ref_t gem_test_make_cross_polytope(int d);
  /* Builds a {d+1}-dimensional gem that is a {d+1}-dimensional cross
    polytope divided into simplices by coning each facet with the
    origin. It has {2^(d+1)} cells, {2*(d+1)} vertices. The vertex at
    the origin has vertex-color {d+1}. The cell walls with wall-color
    {d+1} are all unattached. 
    
    In particular, if {d} is {-1}, the result has a single unattached
    0-dimensional cell. */

gem_ref_t gem_test_make_star(int d, int n);
  /* Builds a star of {n} {d}-dimensional simplices sharing a {d-2}-face ({n} must be even). */

int main(int argc, char **argv);
void write_gem(char *name, int d, gem_ref_t a);
gem_ref_t read_gem(char *name, int d);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    fprintf(stderr, "sizeof(void *) = %d\n", (int)(sizeof(void*)));
    fprintf(stderr, "sizeof(gem_ref_t) = %d\n", (int)(sizeof(gem_ref_t)));
    
    gem_test_new_step_splice_free(0);
    gem_test_new_step_splice_free(1);
    gem_test_new_step_splice_free(2);
    gem_test_new_step_splice_free(3);
    gem_test_new_step_splice_free(5);
    
    gem_test_data(0);
    gem_test_data(1);
    gem_test_data(2);
    gem_test_data(5);
    
    /* gem_test_traverse(-1); */
    /* gem_test_traverse(0); */
    /* gem_test_traverse(1); */
    /* gem_test_traverse(2); */
    /* gem_test_traverse(3); */
    /* gem_test_traverse(5); */
    
    gem_test_write_read("spher", 2, gem_test_make_square(2, FALSE, FALSE));
    gem_test_write_read("projp", 2, gem_test_make_square(2, TRUE,  TRUE));
    gem_test_write_read("klein", 2, gem_test_make_square(2, FALSE, TRUE));
    /* gem_test_write_read("octah", 2, make_gem_octahedron()); */
    /* gem_test_write_read("star4", 2, make_gem_star(2, 4)); */
    /* gem_test_write_read("star6", 3, make_gem_star(3, 6)); */
    return 0;
  }
    
void gem_test_new_step_splice_free(int d)
  {
    fprintf(stderr, "--- testing {gem_node_new}, {gem_step}, {gem_splice}, {gem_node_free} (d = %d) ---\n", d);
    bool_t debug = FALSE;
    gem_ref_t r = gem_node_new(d); demand(gem_node_dim(r) == d, "bad {r} dimension");
    gem_ref_t s = gem_node_new(d); demand(gem_node_dim(s) == d, "bad {s} dimension");
    int i, j;
    for (i = 0; i <= d; i++) 
      { demand(gem_step(r, i) == r, "{gem_node_new}/{gem_step} error {r}.1");
        demand(gem_step(s, i) == s, "{gem_node_new}/{gem_step} error {s}.1");
      }
    /* Check simple {gem_splice} of two detached or attached faces: */
    for (i = 0; i <= d; i++) 
      { if (debug) { fprintf(stderr, "  splicing i = %d ...", i); }
        gem_splice(r, s, i);
        if (debug) { fprintf(stderr, " checking ..."); }
        for (j = 0; j <= d; j++) 
          { gem_ref_t r0 = gem_step(r, j);
            gem_ref_t s0 = gem_step(s, j);
            gem_ref_t r1 = (i == j ? s : r);
            gem_ref_t s1 = (i == j ? r : s);
            if (debug) { fprintf(stderr, " r:%d=%s", j, (r0 == r ? "r" : (r0 == s ? "s" : "!"))); }
            if (debug) { fprintf(stderr, " s:%d=%s", j, (s0 == r ? "r" : (s0 == s ? "s" : "!"))); }
            demand(r0 == r1, "{gem_splice} error {r}.1");
            demand(s0 == s1, "{gem_splice} error {s}.1");
          }
        if (debug) { fprintf(stderr, " unsplicing ..."); }
        gem_splice(r, s, i);
        for (j = 0; j <= d; j++) 
          { if (debug) { fprintf(stderr, " %d", j); }
            demand(gem_step(r, i) == r, "{gem_splice} error {r}.2");
            demand(gem_step(s, i) == s, "{gem_splice} error {s}.2");
          }
        if (debug) { fprintf(stderr, "\n"); }
      }
    gem_node_free(r);
    gem_node_free(s);
  }

void gem_test_data(int d)
  {
    fprintf(stderr, "--- testing {gem_set_data}, {gem_get_data} (d = %d) ---\n", d);
    gem_ref_t r = gem_node_new(d); demand(gem_node_dim(r) == d, "bad {r} dimension");
    int d1 = 418;
    gem_set_data(r, d1);
    int d2 = gem_get_data(r);
    demand(d1 == d2, "{gem_set_data}/{gem_get_data} error 1");
    gem_node_free(r);
  }
    
// void gem_test_traverse(int d)
//   {
//     fprintf(stderr, "--- testing {gem_residues_enum}, {gem_traverse}, {gem_get_label}, {gem_component_free} (d = %d) ---\n", d);
//     bool_t debug = TRUE;
//     if (debug) { fprintf(stderr, "  creating cross polytope\n"); }
//     gem_ref_t p = gem_test_make_cross_polytope(d);
//     
//     int nc = (1 << (d+1)); /* Number of cells. */
//     gem_ref_vec_t node = gem_ref_vec_new(100);
//     int nn = 0;
//     if (debug) { fprintf(stderr, "  traversing it\n"); }
//     gem_traverse(p, d, &node, &nn);
//     demand(nn == nc, "wrong node count");
//     demand(node.e[0] == p, "root is not the first visited node");
//     int k;
//     for (k = 0; k < nn; k++)
//       { int lab = gem_get_label(node.e[k]);
//         demand(lab == k, "wrong label");
//       }
//     
//     /* void gem_residues_enum(gem_ref_t root, int d, int rcol[], int r, int scol[], int s, int na, gem_ref_t nodes[], int *nnP); */
//     if (debug) { fprintf(stderr, "  recycling it\n"); }
//     gem_component_free(p, d);
//     free(node.e);
//     fprintf(stderr, "!! NOT FULLY TESTED\n");
//   }
    
gem_ref_t gem_test_make_cross_polytope(int d)
  { /* Note that {d} may be -1, in which case we create  */
    bool_t debug = TRUE;
    /* Create the {2^(d+1)} cells: */
    int nc = (1 << (d+1)); /* Number of cells. */
    if (debug) { fprintf(stderr, "    creating %d nodes\n", nc); }
    gem_ref_t *p = notnull(malloc(nc*sizeof(gem_ref_t)), "no mem");
    int k;
    for (k = 0; k < nc; k++)
      { p[k] = gem_node_new(d+1); gem_set_data(p[k], k); }
    /* Glue them. Each node {p[k]} is glued to {p[k1]} if {k1>k} and {k,k1} differ by 1 bit. */
    if (debug) { fprintf(stderr, "    gluing the nodes"); }
    for (k = 0; k < nc; k++)
      { int i;
        for (i = 0; i <= d; i++)
          { int k1 = k ^ (1 << i); 
            if (k1 > k)
              { if (debug) { fprintf(stderr, " %d:%d-%d:%d", k,i,k1,i); }
                gem_splice(p[k], p[k1], i);
              }
          }
      }
    if (debug) { fprintf(stderr, "\n"); }
    gem_ref_t res = p[0];
    free(p);
    if (debug) { fprintf(stderr, "    done\n"); }
    return res;
  }

gem_ref_t gem_test_make_square(int d, bool_t rev1, bool_t rev2)
  { gem_ref_t a[4], b[4];
    int i;
    for (i = 0; i < 4; i++) 
      { a[i] = gem_node_new(d); gem_set_data(a[i], 2*i);
        b[i] = gem_node_new(d); gem_set_data(b[i], 2*i + 1);
      }
    for (i = 0; i < 4; i++) 
      { int i1 = (i + 1) % 4;
        gem_splice(a[i], b[i], 0);
        gem_splice(b[i], a[i1], 1);
      }
    if (rev1)
      { gem_splice(a[0], a[2], 2);
        gem_splice(b[0], b[2], 2); 
      }
    else
      { gem_splice(a[0], b[2], 2);
        gem_splice(b[0], a[2], 2);
      }
    if (rev2)
      { gem_splice(a[1], a[3], 2);
        gem_splice(b[1], b[3], 2);
      }
    else
      { gem_splice(a[1], b[3], 2);
        gem_splice(b[1], a[3], 2);
      }
    return a[0];
  } 
     
void gem_test_write_read(char *name, int d, gem_ref_t m)
  { fprintf(stderr, "--- testing write, read %s (d = %d) ---\n", name, d);
    /* write_gem(name, d, m); */
    /* read_gem(name, d, m); */
    fprintf(stderr, "!! NOT TESTED\n");
  }

void write_gem(char *name, int d, gem_ref_t a)
  { assert(FALSE); }
  
gem_ref_t read_gem(char *name, int d)
  { assert(FALSE); }
