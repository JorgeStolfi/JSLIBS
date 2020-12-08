#define PROG_NAME "g3map_test_00_basic"
#define PROG_DESC "basic tests of the {g3map.h} procedures"
#define PROG_VERS "1.0"

/* Last edited on 2016-05-17 14:49:28 by stolfilocal */ 

#define PROG_COPYRIGHT \
  "Copyright © 2015  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2015-12-21"
  
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
#include <g3map.h>

/* INTENAL PROTOTYPES */

void g3map_test_make_step_splice_vert(void);
void g3map_test_make_step_splice_edge(void);
void g3map_test_make_step_splice_face(int M, int N);

void g3map_test_make_step_splice_cell(bool_t cube);
  /* Tests splicing of cells.  If {cube} is true,
    glues two cubes, else glues two tetrahedra. */

void g3map_test_face_make(int N, g3map_place_t e[]);
  /* Creates a polygon with {N} edges, which are returned in {e[0..N-1]}. Also
   does some consistency checks on {g3map_face_make}. */ 
   
g3map_place_t g3map_test_cube_make(void);
  /* Builds a detached cell with topology of a cube. 
     Returns a place on one of the faces. */
  
g3map_place_t g3map_test_tetra_make(void);
  /* Builds a detached cell with topology of a tertahedron. 
     Returns a place on one of the faces. */

void g3map_test_check_cycle_degree(int N, g3map_place_t a, int c0, int c1);
  /* Walks on the map along the cycle 
    that starts at {a} and alternates with steps {c0}
    and {c1}. Checks whether that cycle has exactly {N} distinct
    steps ({2*N} distinct places). */

/* void g3map_test_data(int d); */
/* void g3map_test_traverse(int d); */
/* void g3map_test_write_read(char *name, int d, g3map_place_t m); */

/* g3map_place_t g3map_test_make_square(int d, bool_t rev1, bool_t rev2); */
  /* Builds a {d}-dimensional cone of a 2-dimensional gem consisting of a barycentric division of a square,
    with opposite sides identified with parallel or antiparallel orientations dependinn on {rev1,rev2}. */

/* g3map_place_t g3map_test_make_cross_polytope(int d); */
  /* Builds a {d+1}-dimensional gem that is a {d+1}-dimensional cross
    polytope divided into simplices by coning each facet with the
    origin. It has {2^(d+1)} cells, {2*(d+1)} vertices. The vertex at
    the origin has vertex-color {d+1}. The cell walls with wall-color
    {d+1} are all unattached. 
    
    In particular, if {d} is {-1}, the result has a single unattached
    0-dimensional cell. */

/* g3map_place_t g3map_test_make_star(int d, int n); */
  /* Builds a star of {n} {d}-dimensional simplices sharing a {d-2}-face ({n} must be even). */

int main(int argc, char **argv);
/* void write_gem(char *name, int d, g3map_place_t a); */
/* g3map_place_t read_gem(char *name, int d); */

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    fprintf(stderr, "sizeof(void *) = %d\n", (int)(sizeof(void*)));
    fprintf(stderr, "sizeof(g3map_place_t) = %d\n", (int)(sizeof(g3map_place_t)));
    
    g3map_test_make_step_splice_vert();
    
    g3map_test_make_step_splice_edge();
    
    g3map_test_make_step_splice_face(5, 7);
    g3map_test_make_step_splice_face(2,2);
    g3map_test_make_step_splice_face(1,1);
    
    g3map_test_make_step_splice_cell(FALSE);
    g3map_test_make_step_splice_cell(TRUE);
    
    
    /* g3map_test_data(0); */
    /* g3map_test_data(1); */
    /* g3map_test_data(2); */
    
    /* g3map_test_write_read("cubbe", 2, g3map_test_cube_make()); */
    /* g3map_test_write_read("tetra", 2, g3map_test_tetra_make()); */

    /* g3map_test_write_read("spher", 2, g3map_test_make_square(2, FALSE, FALSE)); */
    /* g3map_test_write_read("projp", 2, g3map_test_make_square(2, TRUE,  TRUE)); */
    /* g3map_test_write_read("klein", 2, g3map_test_make_square(2, FALSE, TRUE)); */
    /* g3map_test_write_read("octah", 2, make_g3map_octahedron()); */
    /* g3map_test_write_read("star4", 2, make_g3map_star(2, 4)); */
    /* g3map_test_write_read("star6", 3, make_g3map_star(3, 6)); */
    return 0;
  }
    
void g3map_test_make_step_splice_vert(void)
  {
    fprintf(stderr, "--- testing {g3map_vert_make} ---\n");
    /* bool_t debug = FALSE; */
    g3map_place_t r = g3map_vert_make(); demand(gem_node_dim(r) == 3, "bad {r} dimension");
    int i;
    for (i = 0; i <= 3; i++) 
      { demand(gem_step(r, i) == r, "{g3map_vert_make} error {r}.1"); }
  }
    
void g3map_test_make_step_splice_edge(void)
  {
    fprintf(stderr, "--- testing {g3map_edge_make,g3map_edge_splice} ---\n");
    /* bool_t debug = FALSE; */
    int N = 5;
    /* Create {N} unattached edges: */
    int k;
    g3map_place_t u[N]; /* Origins of edges. */
    g3map_place_t v[N]; /* Destinations of edges. */
    g3map_place_t e[N]; /* Edges. */
    for (k = 0; k < N; k++)
      { 
        u[k] = g3map_vert_make(); demand(gem_node_dim(u[k]) == 3, "bad {u} dimension");
        v[k] = g3map_vert_make(); demand(gem_node_dim(v[k]) == 3, "bad {v} dimension");
        e[k] = g3map_edge_make(u[k], v[k]); 
        demand(e[k] == u[k], "{} returned wrong end");
        demand(g3map_step(u[k], 0) == v[k], "{g3map_edge_make}/{gem_step} error {v}.1");
        demand(g3map_step(v[k], 0) == u[k], "{g3map_edge_make}/{gem_step} error {v}.2");
        int i;
        for (i = 1; i <= 3; i++) 
          { demand(gem_step(u[k], i) == u[k], "{g3map_edge_make} error {e}.1");
            demand(gem_step(v[k], i) == v[k], "{g3map_edge_make} error {e}.2");
          }
      }
    /* Splice them into a cycle: */
    for (k = 0; k < N; k++)
      { /* Attach edge {e[k]} to the next one: */
        int k1 = (k+1)%N;
        g3map_edge_splice(g3map_step(e[k],0), e[k1]);
        
        /* Check them: */
        demand(g3map_step(u[k], 0) == v[k], "{g3map_edge_splice}/{g3map_step} error {e}.3");
        demand(g3map_step(v[k], 0) == u[k], "{g3map_edge_splice}/{g3map_step} error {e}.4");
        
        demand(g3map_step(u[k1], 0) == v[k1], "{g3map_edge_splice}/{g3map_step} error {e}.5");
        demand(g3map_step(v[k1], 0) == u[k1], "{g3map_edge_splice}/{g3map_step} error {e}.6");
        
        demand(g3map_step(v[k], 1) == u[k1], "{g3map_edge_splice}/{g3map_step} error {e}.7");
        demand(g3map_step(u[k1], 1) == v[k], "{g3map_edge_splice}/{g3map_step} error {e}.8");
        
        if (k < N-1)
          { demand(g3map_step(v[k1], 1) == v[k1], "{g3map_edge_splice}/{g3map_step} error {e}.9"); 
            demand(g3map_step(u[0], 1) == u[0], "{g3map_edge_splice}/{g3map_step} error {e}.10"); 
          }
        else
          { demand(g3map_step(v[0], 1) == u[1], "{g3map_edge_splice}/{g3map_step} error {e}.11");
            demand(g3map_step(u[1], 1) == v[0], "{g3map_edge_splice}/{g3map_step} error {e}.12");
          }
        int i;
        for (i = 2; i <= 3; i++) 
          { demand(gem_step(v[k], i) == v[k], "{g3map_edge_make} error {e}.13");
            demand(gem_step(u[k1], i) == u[k1], "{g3map_edge_make} error {e}.14");
          }
      }
    /* Unsplice them: */
    for (k = 0; k < N; k++)
      { /* Attach edge {e[k]} to the next one: */
        int k1 = (k+1)%N;
        g3map_edge_splice(g3map_step(e[k],0), e[k1]);
        
        /* Check them: */
        demand(g3map_step(u[k], 0) == v[k], "{g3map_edge_splice}/{g3map_step} error {e}.20");
        demand(g3map_step(v[k], 0) == u[k], "{g3map_edge_splice}/{g3map_step} error {e}.21");
        
        demand(g3map_step(u[k1], 0) == v[k1], "{g3map_edge_splice}/{g3map_step} error {e}.22");
        demand(g3map_step(v[k1], 0) == u[k1], "{g3map_edge_splice}/{g3map_step} error {e}.23");
        
        demand(g3map_step(v[k], 1) == v[k], "{g3map_edge_make}/{g3map_step} error {e}.24");
        demand(g3map_step(u[k1], 1) == u[k1], "{g3map_edge_make}/{g3map_step} error {e}.25");
        
        int i;
        for (i = 2; i <= 3; i++) 
          { demand(gem_step(v[k], i) == v[k], "{g3map_edge_make} error {e}.26");
            demand(gem_step(u[k1], i) == u[k1], "{g3map_edge_make} error {e}.27");
          }
      }
  }
     
void g3map_test_make_step_splice_face(int M, int N)
  {
    fprintf(stderr, "--- testing {g3map_face_make}, {g3map_face_splice} ---\n");
    /* bool_t debug = FALSE; */
    
    /* Create two faces {a,b}: */
    g3map_place_t a[M]; /* Edges of first polygon. */
    g3map_place_t b[N]; /* Edges of second polygon. */
    
    g3map_test_face_make(M, a);
    g3map_test_face_make(N, b);

    /* Glue them by 1 edge: */
    g3map_face_splice(a[0], b[0]);
     
    /* Check: */
    int k;
    for (k = 0; k < M+N; k++)
      { 
        if (k == 0)
          { demand(g3map_step(a[k], 2) == b[k], "{g3map_face_splice}/{g3map_step} error {a}.22");
            demand(g3map_step(b[k], 2) == a[k], "{g3map_face_splice}/{g3map_step} error {b}.23");

            demand(gem_step(a[k], 3) == a[k], "{g3map_face_splice}/{g3map_step} error {a}.24");
            demand(gem_step(b[k], 3) == b[k], "{g3map_face_splice}/{g3map_step} error {b}.25");
          }
        else
          { if (k < M)
              { demand(g3map_step(a[k], 2) == a[k], "{g3map_face_splice}/{g3map_step} error {a}.31");
                demand(gem_step(a[k], 3) == a[k], "{g3map_face_splice}/{g3map_step} error {a}.33");
              }
            if (k < N)
              { demand(g3map_step(b[k], 2) == b[k], "{g3map_face_splice}/{g3map_step} error {b}.41");
                demand(gem_step(b[k], 3) == b[k], "{g3map_face_splice}/{g3map_step} error {b}.43");
              }
          }
      }
      
  }
    
void g3map_test_face_make(int N, g3map_place_t e[])
  { 
    /* Create {N} unattached edges {e[0..M-1]: */
    int k;
    for (k = 0; k < N; k++)
      { g3map_place_t uk = g3map_vert_make();
        g3map_place_t vk = g3map_vert_make();
        e[k] = g3map_edge_make(uk, vk); 
      }
    /* Connect them into a circuit: */
    g3map_place_t r = g3map_face_make(N, e);
    demand(r == e[0], "{g3map_face_make} invalid return value"); 

    /* Ensure that corners are properly attached: */
    for (k = 0; k < N; k++)
      { int k1 = (k+1)%N;
        g3map_place_t v0 = g3map_step(e[k], 0); /* Destintation of {e[k]}. */
        g3map_place_t u1 = e[k1];               /* Origin of {e[k+1]}. */
        demand(g3map_step(v0, 1) == u1, "{g3map_face_make}/{g3map_step} error {e}.7");
        demand(g3map_step(u1, 1) == v0, "{g3map_face_make}/{g3map_step} error {e}.8");
      }
      
    /* Check whether face has {N} sides: */
    g3map_test_check_cycle_degree(N, e[0], 0, 1);
    
    /* Ensure that face is unattached: */
    for (k = 0; k < N; k++)
      { g3map_place_t u0 = e[k];  /* Origin of {e[k]}. */
        g3map_place_t v0 = g3map_step(e[k], 0); /* Destintation of {e[k]}. */
        int i;
        for (i = 2; i <= 3; i++)
          { demand(gem_step(u0, i) == u0, "{g3map_face_make} not free {u0}");
            demand(gem_step(v0, i) == v0, "{g3map_face_make} not free {v0}");
          }
      }
      
    /* !!! Should test unsplice, splice several sides, etc. !!! */
  }
  
void g3map_test_make_step_splice_cell(bool_t cube)
  {
    char *cname = ((char*[2]){ "tetra", "cube" })[cube];
    fprintf(stderr, "--- testing {g3map_cell_splice} of %s ---\n", cname);
    /* bool_t debug = FALSE; */
    
    /* Create two cubes: */
    fprintf(stderr, "  %s a\n", cname);
    g3map_place_t a = (cube ? g3map_test_cube_make() : g3map_test_tetra_make());
    fprintf(stderr, "  %s b\n", cname);
    g3map_place_t b = (cube ? g3map_test_cube_make() : g3map_test_tetra_make());
    
    /* Glue by one face: */
    fprintf(stderr, "  splicing the %ss\n", cname);
    g3map_cell_splice(a, b);
    fprintf(stderr, "  done\n");
    
    /* Check the splicing: */
    g3map_place_t ak = a;
    g3map_place_t bk = b;
    
    fprintf(stderr, "  (");
    int ne = (cube ? 4 : 3);
    int k;
    for (k = 0; k < 2*ne; k++)
      { fprintf(stderr, " %d", k);
        demand(g3map_step(ak, 3) == bk, "{g3map_cell_splice}/{g3map_step} error {a}.1");
        demand(g3map_step(bk, 3) == ak, "{g3map_cell_splice}/{g3map_step} error {b}.1");
        /* Walk around the face: */
        if ((k % 2) == 0)
          { /* Move to the other end of the edge: */
            ak = g3map_step(ak, 0);
            bk = g3map_step(bk, 0);
          }
        else
          { /* Move to the adjacent edge:  */
            ak = g3map_step(ak, 1);
            bk = g3map_step(bk, 1);
          }
      }
    assert(ak == a);
    assert(bk == b);
    fprintf(stderr, " )\n");
      
    /* !!! Should test unsplice, splice several faces, etc. !!! */
  }
  
g3map_place_t g3map_test_cube_make(void)
  {
    /* Create 6 squares {f[0..5]} with 24 unglued edges {e[0..23]}: */
    g3map_place_t f[6]; /* One place in each face. */
    g3map_place_t e[24]; /* Edges of face {k} are {e[4*k..4*k+3]}. */
    int k, s;
    for (k = 0; k < 6; k++)
      { g3map_test_face_make(4, &(e[4*k]));
        f[k] = e[4*k];
      }
      
    /* Splice them, so that {f[k]} is opposite to {f[k+3]}: */
    for (s = 0; s < 2; s++)
      { /* Glue faces {f[3*s..3*s+2]} to tach other: */
        for (k = 0; k < 3; k++)
          { int k1 = (k+1)%3;
            /* Splice {f[k]} to {f[k+1]}: */
            g3map_face_splice(g3map_step(e[4*k+14*s], 1), e[4*k1+14*s]);
          }
      }
    g3map_face_splice(e[ 1], e[16]);
    g3map_face_splice(e[ 2], e[23]);
    g3map_face_splice(e[ 5], e[20]);
    g3map_face_splice(e[ 6], e[15]);
    g3map_face_splice(e[ 9], e[12]);
    g3map_face_splice(e[10], e[19]);
    
    /* Ensure that the facets are all unattached: */
    for (k = 0; k < 6; k++)
      { demand(gem_step(f[k], 3) == f[k], "{g3map_test_cube_make}/{g3map_step} error {f}.1"); }
    
    /* Check whether every face has four sides: */
    for (k = 0; k < 6; k++)
      { g3map_test_check_cycle_degree(4, f[k], 0, 1); }
      
    /* Collect the vertices: */
    g3map_place_t v[8]; /* One place near each vertex. */
    for (s = 0; s < 2; s++)
      { int kf = 3*s; /* Face index (0 or 3). */
        int kv = 4*s; /* Index of first vertex of face {kf} (0 or 4). */
        v[kv] = f[kf]; /* Vertices of face {kf} will vbe {v[0..3]}: */
        for (k = 1; k < 4; k++) { v[kv+k] = g3map_step(g3map_step(v[kv+k-1], 0), 1); }
      }
    
    /* Check whether every vertex has three corners: */
    for (k = 0; k < 8; k++)
      { g3map_test_check_cycle_degree(3, v[k], 2, 1); }
    
    /* !!! Should check more !!! */
    
    /* Return a place on {f[0]} such that the face is CCW seen from outside: */
    return e[0];
  }
  
g3map_place_t g3map_test_tetra_make(void)
  {
    /* Create 4 triangles {f[0..3]} with 12 unglued edges {e[0..11]}: */
    g3map_place_t f[4]; /* One place in each face. */
    g3map_place_t e[11]; /* Edges of face {k} are {e[3*k..3*k+2]}. */
    int k, s;
    for (k = 0; k < 4; k++)
      { g3map_test_face_make(3, &(e[3*k]));
        f[k] = e[3*k];
      }
      
    /* Splice them in pairs: {f[0]} with {f[1]}, {f[2]} with {f[3]): */
    for (s = 0; s < 2; s++)
      { /* Glue faces {f[2*s]} with {f[2*s+1]}: */
        g3map_face_splice(e[6*s], g3map_step(e[6*s + 3],0));
      }
      
    g3map_face_splice(e[ 1], g3map_step(e[10],0));
    g3map_face_splice(e[ 2], g3map_step(e[ 8],0));
    g3map_face_splice(e[ 4], g3map_step(e[ 7],0));
    g3map_face_splice(e[ 5], g3map_step(e[11],0));
    
    /* Ensure that the facets are all unattached: */
    for (k = 0; k < 4; k++)
      { demand(gem_step(f[k], 3) == f[k], "{g3map_test_cube_make}/{g3map_step} error {f}.1"); }
    
    /* Check whether every face has three sides: */
    for (k = 0; k < 4; k++)
      { g3map_test_check_cycle_degree(3, f[k], 0, 1); }
      
    /* Collect the vertices: */
    g3map_place_t v[4]; /* One place near each vertex. */
    for (s = 0; s < 4; s++) { v[s] = e[3*s]; }
    
    /* Check whether every vertex has three corners: */
    for (k = 0; k < 4; k++)
      { g3map_test_check_cycle_degree(3, v[k], 2, 1); }
    
    /* !!! Should check more !!! */
    
    /* Return a place on {f[0]}: */
    return e[0];
  }
  
void g3map_test_check_cycle_degree(int N, g3map_place_t a, int c0, int c1)
  {
    int k, r;
    g3map_place_t u0 = a;  /* Place around the cycle. */
    for (k = 0; k < N; k++)
      { 
        /* At this point, {u0} is place index {k} in the cycle, with same hand as {a}. */
        g3map_place_t v0 = g3map_step(u0, c0); /* The {c0}-opposite of place {k}. */
        g3map_place_t u1 = g3map_step(v0, c1); /* Place {k+1}. */
        demand(g3map_step(v0, c1) == u1, "error {e}.7");
        demand(g3map_step(u1, c1) == v0, "error {e}.8");

        /* Check for repetitions: */
        g3map_place_t x0 = a;
        for (r = 0; r < k; r++)
          { demand(x0 != u0, "repeated place 1");
            demand(x0 != v0, "repeated place 2");
            g3map_place_t y0 = g3map_step(x0, c0); 
            demand(y0 != u0, "repeated place 3");
            demand(y0 != v0, "repeated place 4");
            g3map_place_t x1 = g3map_step(y0, c1); 
            x0 = x1;
          }
        u0 = u1;
      }
    demand(u0 == a, "wrong cycle degree");
  }

// 
// void g3map_test_data(int d)
//   {
//     fprintf(stderr, "--- testing {g3map_set_data}, {g3map_get_data} (d = %d) ---\n", d);
//     g3map_place_t r = g3map_node_new(d); demand(g3map_node_dim(r) == d, "bad {r} dimension");
//     int d1 = 418;
//     g3map_set_data(r, d1);
//     int d2 = g3map_get_data(r);
//     demand(d1 == d2, "{g3map_set_data}/{g3map_get_data} error 1");
//     g3map_node(r);
//   }
//     
// void g3map_test_traverse(int d)
//   {
//     fprintf(stderr, "--- testing {g3map_residues_enum}, {g3map_traverse}, {g3map_get_label}, {g3map_component} (d = %d) ---\n", d);
//     bool_t debug = TRUE;
//     if (debug) { fprintf(stderr, "  creating cross polytope\n"); }
//     g3map_place_t p = g3map_test_make_cross_polytope(d);
//     
//     int nc = (1 << (d+1)); /* Number of cells. */
//     g3map_place_vec_t node = g3map_place_vec_new(100);
//     int nn = 0;
//     if (debug) { fprintf(stderr, "  traversing it\n"); }
//     g3map_traverse(p, d, &node, &nn);
//     demand(nn == nc, "wrong node count");
//     demand(node.e[0] == p, "root is not the first visited node");
//     int k;
//     for (k = 0; k < nn; k++)
//       { int lab = g3map_get_label(node.e[k]);
//         demand(lab == k, "wrong label");
//       }
//     
//     /* void g3map_residues_enum(g3map_place_t root, int d, int rcol[], int r, int scol[], int s, int na, g3map_place_t nodes[], int *nnP); */
//     if (debug) { fprintf(stderr, "  recycling it\n"); }
//     g3map_component(p, d);
//     free(node.e);
//     fprintf(stderr, "!! NOT FULLY TESTED\n");
//   }
//     
// g3map_place_t g3map_test_make_cross_polytope(int d)
//   { /* Note that {d} may be -1, in which case we create  */
//     bool_t debug = TRUE;
//     /* Create the {2^(d+1)} cells: */
//     int nc = (1 << (d+1)); /* Number of cells. */
//     if (debug) { fprintf(stderr, "    creating %d nodes\n", nc); }
//     g3map_place_t *p = notnull(malloc(nc*sizeof(g3map_place_t)), "no mem");
//     int k;
//     for (k = 0; k < nc; k++)
//       { p[k] = g3map_node_new(d+1); g3map_set_data(p[k], k); }
//     /* Glue them. Each node {p[k]} is glued to {p[k1]} if {k1>k} and {k,k1} differ by 1 bit. */
//     if (debug) { fprintf(stderr, "    gluing the nodes"); }
//     for (k = 0; k < nc; k++)
//       { int i;
//         for (i = 0; i <= d; i++)
//           { int k1 = k ^ (1 << i); 
//             if (k1 > k)
//               { if (debug) { fprintf(stderr, " %d:%d-%d:%d", k,i,k1,i); }
//                 g3map_splice(p[k], p[k1], i);
//               }
//           }
//       }
//     if (debug) { fprintf(stderr, "\n"); }
//     g3map_place_t res = p[0];
//     free(p);
//     if (debug) { fprintf(stderr, "    done\n"); }
//     return res;
//   }
// 
// g3map_place_t g3map_test_make_square(int d, bool_t rev1, bool_t rev2)
//   { g3map_place_t a[4], b[4];
//     int i;
//     for (i = 0; i < 4; i++) 
//       { a[i] = g3map_node_new(d); g3map_set_data(a[i], 2*i);
//         b[i] = g3map_node_new(d); g3map_set_data(b[i], 2*i + 1);
//       }
//     for (i = 0; i < 4; i++) 
//       { int i1 = (i + 1) % 4;
//         g3map_splice(a[i], b[i], 0);
//         g3map_splice(b[i], a[i1], 1);
//       }
//     if (rev1)
//       { g3map_splice(a[0], a[2], 2);
//         g3map_splice(b[0], b[2], 2); 
//       }
//     else
//       { g3map_splice(a[0], b[2], 2);
//         g3map_splice(b[0], a[2], 2);
//       }
//     if (rev2)
//       { g3map_splice(a[1], a[3], 2);
//         g3map_splice(b[1], b[3], 2);
//       }
//     else
//       { g3map_splice(a[1], b[3], 2);
//         g3map_splice(b[1], a[3], 2);
//       }
//     return a[0];
//   } 
//      
// void g3map_test_write_read(char *name, int d, g3map_place_t m)
//   { fprintf(stderr, "--- testing write, read %s (d = %d) ---\n", name, d);
//     /* write_gem(name, d, m); */
//     /* read_gem(name, d, m); */
//     fprintf(stderr, "!! NOT TESTED\n");
//   }
// 
// void write_gem(char *name, int d, g3map_place_t a)
//   { assert(FALSE); }
//   
// g3map_place_t read_gem(char *name, int d)
//   { assert(FALSE); }
