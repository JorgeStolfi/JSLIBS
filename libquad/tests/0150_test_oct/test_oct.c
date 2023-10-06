#define PROG_NAME "test_oct"
#define PROG_DESC "basic tests of the {oct.h} procedures"
#define PROG_VERS "1.0"

/* Last edited on 2023-10-05 20:21:31 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2007  State University of Campinas (UNICAMP)\n\n" jslibs_copyright
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2009-03-06"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <frgb.h>
#include <argparser.h>
#include <jsfile.h>

#include <oct.h>

oct_arc_t make_map_klein(void);
  /* Builds a Klein bottle with two edges. */

oct_arc_t make_map_star(int32_t n);
  /* Builds a star with {n} edges. */

oct_arc_t make_map_pyramid(int32_t n);
  /* Builds a pyramid with {n} sides. */

int32_t main(int32_t argc, char **argv);
void putwr (oct_arc_t a);
void do_tests(char *name, oct_arc_t m);
void write_map(char *name, oct_arc_t a);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { do_tests("klein", make_map_klein());
    do_tests("star4", make_map_star(4));
    do_tests("pyra4", make_map_pyramid(4));
    return 0;
  }
  
void do_tests(char *name, oct_arc_t m)
  {
    /* Check edge/tumblecode decomposition: */
    oct_edge_t ed = oct_edge(m);
    assert(ed != NULL);
    fprintf(stderr, "sizeof(void *) = %d\n", (int32_t)(sizeof(void*)));
    fprintf(stderr, "sizeof(oct_edge_t) = %d\n", (int32_t)(sizeof(oct_edge_t)));
    fprintf(stderr, "sizeof(oct_arc_t) = %d\n", (int32_t)(sizeof(oct_arc_t)));
    write_map(name, m);
    int32_t it;
    for (it = 0; it < 8; it++)
      { oct_bits_t t = it;
        oct_arc_t e = oct_orient(ed, t);
        assert(t == oct_tumble_code(e));
        assert(ed == oct_edge(e));
      
        fprintf(stderr, "  %-10s ", "e =");          oct_write_arc(stderr, e, 1);             fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "rot(e) =");     oct_write_arc(stderr, oct_rot(e), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "sym(e) =");     oct_write_arc(stderr, oct_sym(e), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "tor(e) =");     oct_write_arc(stderr, oct_tor(e), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "vflip(e) =");   oct_write_arc(stderr, oct_vflip(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "fflip(e) =");   oct_write_arc(stderr, oct_fflip(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "duar(e) =");    oct_write_arc(stderr, oct_duar(e), 1);   fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "dual(e) =");    oct_write_arc(stderr, oct_dual(e), 1);   fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "onext(e) =");   oct_write_arc(stderr, oct_onext(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "dnext(e) =");   oct_write_arc(stderr, oct_dnext(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "lnext(e) =");   oct_write_arc(stderr, oct_lnext(e), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "rnext(e) =");   oct_write_arc(stderr, oct_rnext(e), 1);  fprintf(stderr, "\n");

        fprintf(stderr, "Checking effect of tumbling ops on tubmble bits ...\n");
        assert(oct_lon_bit(oct_sym(e))    != oct_lon_bit(e));
        assert(oct_lon_bit(oct_vflip(e))  != oct_lon_bit(e));
        assert(oct_lon_bit(oct_fflip(e))  == oct_lon_bit(e));
        
        assert(oct_trn_bit(oct_sym(e))    != oct_trn_bit(e));
        assert(oct_trn_bit(oct_vflip(e))  == oct_trn_bit(e));
        assert(oct_trn_bit(oct_fflip(e))  != oct_trn_bit(e));
        
        assert(oct_dop_bit(oct_sym(e))    == oct_dop_bit(e));
        assert(oct_dop_bit(oct_vflip(e))  == oct_dop_bit(e));
        assert(oct_dop_bit(oct_fflip(e))  == oct_dop_bit(e));
        assert(oct_dop_bit(oct_rot(e))    != oct_dop_bit(e));
        assert(oct_dop_bit(oct_tor(e))    != oct_dop_bit(e));
        assert(oct_dop_bit(oct_dual(e))   != oct_dop_bit(e));
        assert(oct_dop_bit(oct_duar(e))   != oct_dop_bit(e));
        
        assert(oct_cir_bit(oct_sym(e))    == oct_cir_bit(e));
        assert(oct_cir_bit(oct_vflip(e))  != oct_cir_bit(e));
        assert(oct_cir_bit(oct_fflip(e))  != oct_cir_bit(e));
        assert(oct_cir_bit(oct_rot(e))    == oct_cir_bit(e));
        assert(oct_cir_bit(oct_tor(e))    == oct_cir_bit(e));
        assert(oct_cir_bit(oct_dual(e))   != oct_cir_bit(e));
        assert(oct_cir_bit(oct_duar(e))   != oct_cir_bit(e));
        
        fprintf(stderr, "Checking tumbling identities ...\n");
        assert(oct_fflip(e) != e);
        assert(oct_vflip(e) != e);
        assert(oct_fflip(e) != oct_vflip(e));
        assert(oct_sym(e) != e);
        assert(oct_rot(e) != e);
        assert(oct_tor(e) != e);
        assert(oct_rot(e) != oct_tor(e));
        
        assert(oct_fflip(oct_rot(e)) != e);
        assert(oct_vflip(oct_rot(e)) != e);
        assert(oct_sym(oct_rot(e)) != e);
        
        assert(oct_fflip(oct_fflip(e)) == e);
        assert(oct_vflip(oct_vflip(e)) == e);
        assert(oct_sym(oct_sym(e)) == e);
        assert(oct_rot(oct_rot(e)) == oct_sym(e));
        assert(oct_vflip(oct_fflip(e)) == oct_sym(e));
        assert(oct_rot(oct_rot(e)) == oct_sym(e));
        
        fprintf(stderr, "Checking some walking identities ...\n");
        assert(oct_onext(oct_oprev(e)) == e);
        assert(oct_lnext(oct_lprev(e)) == e);
        assert(oct_dnext(oct_dprev(e)) == e);
        assert(oct_rnext(oct_rprev(e)) == e);
        
        assert(oct_rnext(oct_onext(e)) == oct_sym(e));
        assert(oct_lnext(oct_dnext(e)) == oct_sym(e));
        assert(oct_onext(oct_lnext(e)) == oct_sym(e));
        assert(oct_dnext(oct_rnext(e)) == oct_sym(e));

        assert(oct_onext(oct_fflip(e)) == oct_fflip(oct_oprev(e)));
        assert(oct_dnext(oct_fflip(e)) == oct_fflip(oct_dprev(e)));
        assert(oct_lnext(oct_vflip(e)) == oct_vflip(oct_lprev(e)));
        assert(oct_rnext(oct_vflip(e)) == oct_vflip(oct_rprev(e)));

        assert(oct_rot(oct_onext(e)) == oct_lnext(oct_rot(e)));
        assert(oct_rot(oct_dprev(e)) == oct_rprev(oct_rot(e)));
        /* ... many more ... */
      }
  }
    
oct_arc_t make_map_klein(void)
  { oct_arc_t a, b;
    a = oct_make_edge(); 
    b = oct_make_edge();
    oct_splice(a, b);
    oct_splice(oct_sym(a), a);
    oct_splice(oct_fflip(oct_sym(b)), a);
    return a;
  } 
     
oct_arc_t make_map_star(int32_t n)
  { oct_arc_t a, b;
    a = oct_make_edge();
    int32_t k;
    for(k = 1; k < n; k++)
      { b = oct_make_edge();
        oct_splice(a, b);
        a = b;
      }
    return a;
  } 
     
oct_arc_t make_map_pyramid(int32_t n)
  { /* Build an {n}-armed star: */
    oct_arc_t a = make_map_star(n);
    /* Build an {n}-sided ring: */
    oct_arc_t b = make_map_star(n);
    b = oct_fflip(oct_rot(b));
    /* Stitch them together: */
    oct_arc_t c = oct_sym(a);
    int32_t k;
    for(k = 0; k < n; k++)
      { oct_splice(b, c);
        c = oct_dnext(c);
        b = oct_lnext(b);
      }
    return c;
  } 
 
void write_map(char *name, oct_arc_t a)
  { char *filename = NULL;
    asprintf(&filename, "out/%s.oct", name);
    FILE *wr = open_write(filename, TRUE);
    oct_arc_vec_t root = oct_arc_vec_new(1); /* Root list. */
    root.e[0] = a;
    oct_arc_vec_t A = oct_arc_vec_new(0); /* Edge table. */
    oct_write_map(wr, &root, &A);
    fclose(wr);
    free(filename);
  }

