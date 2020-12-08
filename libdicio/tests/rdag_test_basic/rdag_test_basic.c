/* Basic tests of {rdag.h} and {rdag_io.h} */
/* Last edited on 2009-10-28 22:57:51 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>

#include <rdag.h>
#include <rdag_def.h>
#include <rdag_io.h>

int main (int argc, char **argv);
void rdag_test_show(FILE *wr, char *title, rdag_t *D);

typedef rdag_node_data_t ND; 

int main (int argc, char **argv)
  {
    rdag_t *D = rdag_new(10, 8, 1, 1000);
    rdag_test_show(stderr, "--- newly allocated ---", D);
    
    rdag_node_t n0 = rdag_node_NULL;
    rdag_node_t n1 = rdag_node_from_fields(D, n0, 1, 1, n0);
    rdag_node_t n2 = rdag_node_from_fields(D, n1, 2, 1, n0);
    rdag_node_t n3 = rdag_node_from_fields(D, n2, 3, 1, n1);
    rdag_node_t n4 = rdag_node_from_fields(D, n1, 2, 1, n0);
    assert(n4 == n2);
    assert(n4 != n3);
    rdag_test_show(stderr, "--- with three nodes ---", D);
    
    auto rdag_disp_t enter(uint32_t len, rdag_symbol_t o, rdag_node_t s);
    auto rdag_disp_t exit(uint32_t len, rdag_symbol_t o, rdag_node_t s);
    auto rdag_disp_t push(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t);
    auto rdag_disp_t pop(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t);
    
    rdag_disp_t de = rdag_enum_paths(D, 1, n3, enter, push, pop, exit);
    fprintf(stderr, "result = %u", de);
    
    rdag_disp_t enter(uint32_t len, rdag_symbol_t o, rdag_node_t s)
      { fprintf(stderr, "%*senter(%u,%u,%u)\n", len, "", len, o, s);
        return rdag_disp_FINE;
      }
   
    rdag_disp_t exit(uint32_t len, rdag_symbol_t o, rdag_node_t s)
      { fprintf(stderr, "%*sexit(%u,%u,%u)\n", len, "", len, o, s); 
        return rdag_disp_FINE; 
      }
      
    rdag_disp_t push(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t)
      { fprintf(stderr, "%*spush(%u,%u,%u,%u,%u)\n", len, "", len, s, i, o, t); 
        return rdag_disp_FINE; 
      }
      
    rdag_disp_t pop(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t)
      { fprintf(stderr, "%*spop(%u,%u,%u,%u,%u)\n", len, "", len, s, i, o, t); 
        return rdag_disp_FINE; 
      }

    return 0;
  }
  
void rdag_test_show(FILE *wr, char *title, rdag_t *D)
  {
    fprintf(wr, "%s\n", title);
    fprintf(wr, "  bit sizes: nn = %u  ni = %u  no = %u\n", D->nn, D->ni, D->no);
    fprintf(wr, "  max_node = %u\n", D->max_node);
    fprintf(wr, "  max_alloc_node = %u\n", D->max_alloc_node);
    fprintf(wr, "  hash_valid = %c\n", "FT"[D->hash_valid]);
    fprintf(wr, "  hash_size = %u\n", D->hash_size);
    fprintf(wr, "  masks:");
    fprintf(wr, "  node = %u", D->mask_node);
    fprintf(wr, "  i_symbol = %u", D->mask_i_symbol);
    fprintf(wr, "  o_symbol = %u", D->mask_o_symbol);
    fprintf(wr, "  \n");
    fprintf(wr, "  shifts:");
    fprintf(wr, "  f-link = %u", D->shift_f_link);
    fprintf(wr, "  i-mark = %u", D->shift_i_mark);
    fprintf(wr, "  o-mark = %u", D->shift_o_mark);
    fprintf(wr, "  p-link = %u", D->shift_p_link);
    fprintf(wr, "  \n");
    fprintf(wr, "  node data:\n");
    uint32_t k;
    for (k = 0; k < D->max_node; k++)
      { rdag_node_t s = k + 1;
        rdag_node_data_t dt; rdag_node_data_get(D, s, &dt);
        fprintf(wr, "  ");
        fprintf(wr, "  node = %10u", s);
        fprintf(wr, "  f = %10u", dt.f_link);
        fprintf(wr, "  i = %10u", dt.i_mark);
        fprintf(wr, "  o = %10u", dt.o_mark);
        fprintf(wr, "  p = %10u", dt.p_link);
        if (D->hash_valid) { fprintf(wr, "  hnext = %u", D->hnext[s-1]); }
        fprintf(wr, "  \n");
      }
    if (D->hash_valid)
      {
        fprintf(wr, "  hash bucket heads:\n");
        uint32_t h;
        for (h = 0; h < D->hash_size; h++)
          { if (D->hinit[h] != rdag_node_NULL)
              { fprintf(wr, "  ");
                fprintf(wr, "  hashval = %10u", h);
                fprintf(wr, "  node = %10u", D->hinit[h]);
                fprintf(wr, "  \n");
              }
          }
      }
  }
