#define PROG_NAME "test_haf"
#define PROG_DESC "basic tests of the {haf.h} procedures"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-10 00:05:38 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright � 2023  State University of Campinas (UNICAMP)\n\n" jslibs_copyright
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2023-10-05"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

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
#include <jsprintf.h>

#include <haf.h>
#include <haf_write.h>
#include <haf_read.h>
#include <haf_enum.h>
#include <haf_shapes.h>

int32_t main(int32_t argc, char **argv);
void putwr (haf_arc_t a);

void test_edge_enum(void);
void do_test_edge_enum(uint64_t NR, haf_arc_t root[], uint64_t NE_exp, uint64_t NC_exp);
  /* Tests {haf_enum_edges} on {root[0..NR-1]}.
    Checks whether the edge and component counts match {NE_exp} and {NC_exp},
    respectively. */

void test_basic(char *name, haf_arc_t m);
  /* Tests the basic operation on the component reached from {m}.
    Also creates an additional component {m1}, and calls 
    {test_enum_write_read(name,m,m1,NULL)}. */
    
void test_enum_write_read(char *name, haf_arc_t a, haf_arc_t b, haf_arc_t c);
  /* Writes to file "out/{name}.haf" the components reachable from {a,b,c}
    and reads them back.  Any of {a,b,c} which are {NULL} are ignored. */

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { test_edge_enum();
    test_basic("torus", haf_shapes_torus());
    test_basic("star4", haf_shapes_star(4));
    test_basic("pyra4", haf_shapes_pyramid(4));
    return 0;
  }
  
void test_basic(char *name, haf_arc_t m)
  {
    /* Check edge/tumblecode decomposition: */
    haf_edge_t ed = haf_edge(m);
    assert(ed != NULL);
    fprintf(stderr, "sizeof(void *) = %d\n", (int32_t)(sizeof(void*)));
    fprintf(stderr, "sizeof(haf_edge_t) = %d\n", (int32_t)(sizeof(haf_edge_t)));
    fprintf(stderr, "sizeof(haf_arc_t) = %d\n", (int32_t)(sizeof(haf_arc_t)));
    for (uint8_t it = 0; it < 2; it++)
      { haf_dir_bit_t t = it;
        haf_arc_t a = haf_orient(ed, t);
        assert(t == haf_dir_bit(a));
        assert(ed == haf_edge(a));
        assert(haf_base_arc(ed) == haf_orient(ed, 0));
      
        fprintf(stderr, "  %-10s ", "a =");          haf_write_arc(stderr, a, 1);             fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "sym(a) =");     haf_write_arc(stderr, haf_sym(a), 1);    fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "onext(a) =");   haf_write_arc(stderr, haf_onext(a), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "dnext(a) =");   haf_write_arc(stderr, haf_dnext(a), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "lnext(a) =");   haf_write_arc(stderr, haf_lnext(a), 1);  fprintf(stderr, "\n");
        fprintf(stderr, "  %-10s ", "rnext(a) =");   haf_write_arc(stderr, haf_rnext(a), 1);  fprintf(stderr, "\n");

        fprintf(stderr, "Checking effect of {.sym} on direction bit ...\n");
        assert(haf_dir_bit(haf_sym(a)) != haf_dir_bit(a));
        
        fprintf(stderr, "Checking {.sym} involution identities ...\n");
        assert(haf_sym(a) != a);
        assert(haf_sym(haf_sym(a)) == a);
       
        fprintf(stderr, "Checking some walking identities ...\n");
        assert(haf_onext(haf_oprev(a)) == a);
        assert(haf_lnext(haf_lprev(a)) == a);
        assert(haf_dnext(haf_dprev(a)) == a);
        assert(haf_rnext(haf_rprev(a)) == a);
        
        assert(haf_rnext(haf_onext(a)) == haf_sym(a));
        assert(haf_lnext(haf_dnext(a)) == haf_sym(a));
        assert(haf_onext(haf_lnext(a)) == haf_sym(a));
        assert(haf_dnext(haf_rnext(a)) == haf_sym(a));
        /* ... many more ... */
      } 

    haf_arc_t m1 = haf_shapes_orange(5);
    test_enum_write_read(name, m, m1, NULL);
  }
  
void test_edge_enum(void)
  { fprintf(stderr, "> test_edge_enum\n");
    
    uint32_t m_z = 1; 
    haf_arc_t z = haf_shapes_star(m_z);
    uint64_t NE_exp_z = m_z;
  
    uint32_t m_a = 10; 
    haf_arc_t a = haf_shapes_star(m_a);
    uint64_t NE_exp_a = m_a;

    uint32_t m_b =  7;
    haf_arc_t b = haf_shapes_pyramid(m_b);
    uint64_t NE_exp_b = 2*m_b;

    /* Root pointers and expected edge count: */
    haf_arc_t root[3];
    
    fprintf(stderr, "  -- single edge:\n");
    root[0] = z;
    do_test_edge_enum(1, root, NE_exp_z, 1);

    fprintf(stderr, "  -- star(%u):\n", m_a);
    root[0] = a;
    do_test_edge_enum(1, root, NE_exp_a, 1);

    fprintf(stderr, "  -- pyramid(%u):\n", m_b);
    root[0] = b;
    do_test_edge_enum(1, root, NE_exp_b, 1);

    fprintf(stderr, " star(%u) and pyramid(%u)�2:\n", m_a, m_b);
    root[0] = a;
    root[1] = b;
    root[2] = haf_sym(haf_lnext(haf_lnext(a)));
    do_test_edge_enum(3, root, NE_exp_a + NE_exp_b, 2);
    
    fprintf(stderr, "< test_edge_enum\n");
    
  }
    
void do_test_edge_enum(uint64_t NR, haf_arc_t root[], uint64_t NE_exp, uint64_t NC_exp)
  {
    bool_t debug = FALSE;
    
    if (debug) { fprintf(stderr, "  > {do_test_edge_enum}\n"); }
    
    /* Enumerate edges: */
    haf_edge_id_t eid0 = 1003;
    haf_edge_vec_t A = haf_edge_vec_new(1000);
    haf_edge_id_vec_t oid = haf_edge_id_vec_new(10);
    haf_arc_vec_t C = haf_arc_vec_new(10);
    haf_enum_edges(NR, root, eid0, &A, &oid, &C);
    demand(A.ne == NE_exp, "edge count does not check");
    demand(oid.ne == NE_exp, "edge count does not check");
    demand(C.ne == NC_exp, "edge count does not check");
    
    /* Check that edges have been numbered correctly: */
    for(haf_edge_count_t ke = 0; ke < A.ne; ke++) 
      { haf_edge_t edk = A.e[ke];
        haf_arc_t ak = haf_base_arc(edk);
        demand(ak != NULL, "null arc in input list");
        haf_arc_id_t akid = haf_arc_id(ak);
        demand(akid == 2*(ke + eid0), "bad arc number (even)");
        demand(haf_sym(ak) != NULL, "null {.sym} pointer");
        demand(haf_arc_id(haf_sym(ak)) == akid + 1, "bad arc number (odd)");
      }

    /* Do another enum to check {oid}: */
    haf_edge_id_t eid1 = 1950;
    haf_enum_edges(NR, root, eid1, &A, &oid, &C);
    demand(A.ne == NE_exp, "edge count does not check (A)");
    demand(oid.ne == NE_exp, "edge count does not check (oid)");
    demand(C.ne == NC_exp, "component count does not check");
    for(haf_edge_count_t ke = 0; ke < A.ne; ke++) 
      { haf_edge_t edk = A.e[ke];
        demand(edk != NULL, "null edge in input list");
        demand(haf_edge_id(edk) == ke + eid1, "bad new edge number in table {A}");
        demand(oid.e[ke] == ke + eid0, "bad old edge number in table {oid}");
      }

    if (debug) { fprintf(stderr, "  < {do_test_edge_enum}\n"); }
  }    
 
void test_enum_write_read(char *name, haf_arc_t a, haf_arc_t b, haf_arc_t c)
  { fprintf(stderr, "> test_enum_write_read}\n");
    fprintf(stderr, "  testing enum/write/read name = %s\n", name);
    
    haf_arc_t root_wr[3]; /* Root list. */
    haf_arc_count_t NR_wr = 0;
    if (a != NULL) { root_wr[NR_wr] = a; NR_wr++; }
    if (b != NULL) { root_wr[NR_wr] = b; NR_wr++; }
    if (c != NULL) { root_wr[NR_wr] = c; NR_wr++; }

    fprintf(stderr, "running {haf_enum_edges}...\n");
    haf_edge_id_t eid0_wr = 4615;
    haf_edge_vec_t A_wr = haf_edge_vec_new(0); /* The edges. */
    haf_edge_id_vec_t oid = haf_edge_id_vec_new(10);
    haf_arc_vec_t C_wr = haf_arc_vec_new(0); /* One arc on each component. */
    haf_enum_edges(NR_wr, root_wr, eid0_wr, &A_wr, &oid, &C_wr);
    haf_edge_count_t NE_wr = A_wr.ne;
    fprintf(stderr, "found %u edges\n", A_wr.ne);
    fprintf(stderr, "found %u connected components\n", C_wr.ne);
    assert(C_wr.ne <= NR_wr);
    for (haf_edge_count_t ke = 0; ke < NE_wr; ke++)
      { demand(haf_edge_id(A_wr.e[ke]) == eid0_wr + ke, "edge id not consistent with index"); }

    fprintf(stderr, "running {haf_write_map}...\n");
    char *filename = jsprintf("out/%s.haf", name);
    FILE *wr = open_write(filename, TRUE);
    haf_write_map(wr, NE_wr, A_wr.e, eid0_wr, NR_wr, root_wr);
    fclose(wr);
    
    fprintf(stderr, "running {haf_read_map}...\n");
    FILE *rd = open_read(filename, TRUE);
    haf_edge_vec_t A_rd = haf_edge_vec_new(0); /* The edges read. */
    haf_arc_vec_t R_rd = haf_arc_vec_new(0); /* One arc on each component. */
    haf_edge_id_t eid0_rd = 4634;
    haf_read_map(rd, &A_rd, &eid0_rd, &R_rd);
    fclose(rd);
    
    fprintf(stderr, "comparing map written and map read\n");
    haf_edge_count_t NE_rd = A_rd.ne;
    demand(NE_rd == NE_wr, "edge counts do not match");
    demand(eid0_rd == eid0_wr, "min edge ids do not match");
    for (haf_edge_count_t ke = 0; ke < NE_rd; ke++)
      { haf_edge_id_t eidk_wr = haf_edge_id(A_wr.e[ke]);
        haf_edge_id_t eidk_rd = haf_edge_id(A_rd.e[ke]);
        demand(eidk_rd == eidk_wr, "read and write edge ids do not match");
        demand(eidk_rd == eid0_wr + ke, "read edge id not consistent with index");
      }

    for (haf_edge_count_t ke = 0; ke < NE_rd; ke++)
      { haf_edge_t edk_wr = A_wr.e[ke];
        haf_edge_t edk_rd = A_rd.e[ke];
        for (uint8_t t = 0; t < 2; t++)
          { haf_arc_t akt_wr = haf_orient(edk_wr, t);
            haf_arc_t bkt_wr = haf_lnext(akt_wr);
            
            haf_arc_t akt_rd = haf_orient(edk_rd, t);
            haf_arc_t bkt_rd = haf_lnext(akt_rd);
            
            demand(haf_arc_id(bkt_rd) == haf_arc_id(bkt_wr), "{.lnext} arc ids do not match");
          }
      }

    haf_arc_t *root_rd = R_rd.e;
    haf_arc_count_t NR_rd = R_rd.ne;
    demand(NR_rd == NR_wr, "root counts do not match");
    for (haf_arc_count_t kr = 0; kr < NR_rd; kr++)
      { haf_arc_t rk_rd = root_rd[kr];
        haf_arc_t rk_wr = root_wr[kr];
        demand(haf_arc_id(rk_rd) == haf_arc_id(rk_wr), "roots arc ids do not match");
      }

    free(filename);
    fprintf(stderr, "< test_enum_write_read}\n");
  }
