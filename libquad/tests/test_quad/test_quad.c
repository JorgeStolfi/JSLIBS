#define PROG_NAME "test_quad"
#define PROG_DESC "basic tests of the {quad.h} procedures"
#define PROG_VERS "1.0"

/* Last edited on 2011-12-22 19:47:06 by stolfilocal */ 

#define PROG_COPYRIGHT \
  "Copyright � 2007  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2009-03-06"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <frgb.h>
#include <argparser.h>
#include <jsfile.h>

#include <quad.h>

quad_arc_t make_map_torus(void);
  /* Builds a Torus bottle with two edges. */

quad_arc_t make_map_star(int n);
  /* Builds a star with {n} edges. */

quad_arc_t make_map_pyramid(int n);
  /* Builds a pyramid with {n} sides. */

int main(int argc, char **argv);
void putwr (quad_arc_t a);
void do_tests(char *name, quad_arc_t m);
void write_map(char *name, quad_arc_t a);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { do_tests("torus", make_map_torus());
    do_tests("star4", make_map_star(4));
    do_tests("pyra4", make_map_pyramid(4));
    return 0;
  }
  
void do_tests(char *name, quad_arc_t m)
  {
    /* Check edge/tumblecode decomposition: */
    quad_edge_t E = quad_edge(m);
    assert(E != NULL);
    fprintf(stderr, "sizeof(void *) = %d\n", (int)(sizeof(void*)));
    fprintf(stderr, "sizeof(quad_edge_t) = %d\n", (int)(sizeof(quad_edge_t)));
    fprintf(stderr, "sizeof(quad_arc_t) = %d\n", (int)(sizeof(quad_arc_t)));
    write_map(name, m);
    int it;
    for (it = 0; it < 4; it++)
      { quad_bits_t t = it;
        quad_arc_t e = quad_orient(E, t);
        assert(t == quad_tumble_code(e));
        assert(E == quad_edge(e));
      
        fprintf(stderr, "  %-10s ", "e =");          quad_write_arc(stderr, e, 1);              fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "rot(e) =");     quad_write_arc(stderr, quad_rot(e), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "sym(e) =");     quad_write_arc(stderr, quad_sym(e), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "tor(e) =");     quad_write_arc(stderr, quad_tor(e), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "onext(e) =");   quad_write_arc(stderr, quad_onext(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "dnext(e) =");   quad_write_arc(stderr, quad_dnext(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "lnext(e) =");   quad_write_arc(stderr, quad_lnext(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "rnext(e) =");   quad_write_arc(stderr, quad_rnext(e), 1);  fprintf(stderr, "\n");

        fprintf(stderr, "Checking effect of tumbling ops on tubmble bits ...\n");
        assert(quad_lon_bit(quad_sym(e))    != quad_lon_bit(e));
        
        assert(quad_dop_bit(quad_sym(e))    == quad_dop_bit(e));
        assert(quad_dop_bit(quad_rot(e))    != quad_dop_bit(e));
        assert(quad_dop_bit(quad_tor(e))    != quad_dop_bit(e));
        
        fprintf(stderr, "Checking tumbling identities ...\n");
        assert(quad_sym(e) != e);
        assert(quad_rot(e) != e);
        assert(quad_tor(e) != e);
        assert(quad_rot(e) != quad_tor(e));
        
        assert(quad_sym(quad_rot(e)) != e);
        
        assert(quad_sym(quad_sym(e)) == e);
        assert(quad_rot(quad_rot(e)) == quad_sym(e));
        assert(quad_rot(quad_rot(e)) == quad_sym(e));
        
        fprintf(stderr, "Checking some walking identities ...\n");
        assert(quad_onext(quad_oprev(e)) == e);
        assert(quad_lnext(quad_lprev(e)) == e);
        assert(quad_dnext(quad_dprev(e)) == e);
        assert(quad_rnext(quad_rprev(e)) == e);
        
        assert(quad_rnext(quad_onext(e)) == quad_sym(e));
        assert(quad_lnext(quad_dnext(e)) == quad_sym(e));
        assert(quad_onext(quad_lnext(e)) == quad_sym(e));
        assert(quad_dnext(quad_rnext(e)) == quad_sym(e));

        assert(quad_rot(quad_onext(e)) == quad_lnext(quad_rot(e)));
        assert(quad_rot(quad_dprev(e)) == quad_rprev(quad_rot(e)));
        /* ... many more ... */
      }
  }
    
quad_arc_t make_map_torus(void)
  { quad_arc_t a, b;
    a = quad_make_edge(); 
    b = quad_make_edge();
    quad_splice(a, b);
    quad_splice(quad_sym(a), a);
    quad_splice(quad_sym(b), a);
    return a;
  } 
     
quad_arc_t make_map_star(int n)
  { quad_arc_t a, b;
    a = quad_make_edge();
    int k;
    for(k = 1; k < n; k++)
      { b = quad_make_edge();
        quad_splice(a, b);
        a = b;
      }
    return a;
  } 
     
quad_arc_t make_map_pyramid(int n)
  { /* Build an {n}-armed star: */
    quad_arc_t a = make_map_star(n);
    /* Build an {n}-sided ring: */
    quad_arc_t b = make_map_star(n);
    b = quad_rot(b);
    /* Stitch them together: */
    quad_arc_t c = quad_sym(a);
    int k;
    for(k = 0; k < n; k++)
      { quad_splice(b, c);
        c = quad_dnext(c);
        b = quad_lnext(b);
      }
    return c;
  } 
 
void write_map(char *name, quad_arc_t a)
  { char *filename = NULL;
    asprintf(&filename, "out/%s.quad", name);
    FILE *wr = open_write(filename, TRUE);
    quad_arc_vec_t root = quad_arc_vec_new(1); /* Root list. */
    root.e[0] = a;
    quad_arc_vec_t et = quad_arc_vec_new(0); /* Edge table. */
    quad_write_map(wr, &root, &et);
    fclose(wr);
    free(filename);
  }

