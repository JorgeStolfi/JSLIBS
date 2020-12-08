/* Procedures of {raut.h} that require access to the internal rep. */
/* Last edited on 2009-10-31 00:05:37 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <bool.h>

#include <rdag.h>
#include <raut.h>
#include <raut_def.h>

raut_t *raut_new(uint32_t nn, uint32_t ni, rdag_node_t max_alloc_node)
  {
    raut_t *A = (raut_t *)notnull(malloc(sizeof(raut_t)), "no mem");
    A->D = rdag_new(nn, ni, 1, max_alloc_node);
    A->root = raut_state_VOID;
    A->doc = NULL;
    return A;
  }

raut_state_t raut_root_get(raut_t *A)
  { 
    return A->root;
  }

rdag_node_t raut_node_max(raut_t *A)
  { 
    return rdag_node_max(A->D);
  }

rdag_t *raut_dag_get(raut_t *A)
  { 
    return A->D;
  }

char *raut_doc_get(raut_t *A)
  { 
    return A->doc;
  }

void raut_doc_set(raut_t *A, char *doc)
  { 
    A->doc = doc;
  }
  
bool_t raut_state_is_accepting(raut_t *A, raut_state_t u)
  { 
    demand(u.ac <= 1, "invalid accept bit");
    return(u.ac == 1);
  }

bool_t raut_state_has_arcs(raut_t *A, raut_state_t u)
  { 
    demand(u.ac <= 1, "invalid accept bit");
    return (u.nd != rdag_node_NULL);
  }

void raut_root_set(raut_t *A, raut_state_t u)
  { 
    demand(u.ac <= 1, "invalid accept bit");
    demand(u.nd <= rdag_node_max(A->D), "invalid state node");
    A->root = u;
  }

rdag_disp_t raut_enum_suffs(raut_t *A, raut_state_t u, rdag_string_action_t *proc)
  {
    
    rdag_symbol_vec_t str = rdag_symbol_vec_new(50);
    int nstr = 0;
    
    bool_t debug = FALSE;
    
    auto rdag_disp_t enter(uint32_t len, rdag_node_t t, rdag_symbol_t o);
    auto rdag_disp_t push(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t);
    auto rdag_disp_t pop(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t);

    rdag_disp_t enter(uint32_t len, rdag_symbol_t o, rdag_node_t t)
      { if (o == 1)
          { /* Process suffix {str[0..nstr-1]}: */
            rdag_disp_t de = proc(nstr, str.e);
            return de;
          }
        else
          { /* Keep going: */
            return rdag_disp_FINE;
          }
      }
    
    rdag_disp_t push(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t)
      { /* Make sure that there is space in {str}: */
        if (debug) { fprintf(stderr, "  %*spush(%u)\n", 2*len, "", i); }
        rdag_symbol_vec_expand(&str, nstr);
        /* Append {i} to {str}: */
        str.e[nstr] = i;
        nstr++;
        /* Keep going: */
        return rdag_disp_FINE;
      }

    rdag_disp_t pop(uint32_t len, rdag_node_t s, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t t)
      { /* Remove the last letter of {str}: */
        if (debug) { fprintf(stderr, "  %*spop(%u)\n", 2*len, "", i); }
        nstr--;
        /* Keep going: */
        return rdag_disp_FINE;
      }

    
    rdag_t *D = raut_dag_get(A);
    rdag_disp_t dd = rdag_enum_paths(D, u.ac, u.nd, &enter, &push, &pop, NULL);
    return dd;
  }

raut_state_t raut_arc_set(raut_t *A, raut_state_t u, rdag_symbol_t i, raut_state_t v)
  {
    demand(u.ac <= 1, "invalid accept bit");
    demand(u.nd <= raut_node_max(A), "invalid state node");
    
    rdag_t *D = A->D;
    
    auto rdag_node_t do_set(rdag_node_t s);
      /* Returns a node {s'} which is like {s} except that 
        it has a subnode with i-mark {i}, o-mark {v.ac},
        and p-link {v.nd}. */
         
    auto rdag_node_t do_del(rdag_node_t s);
      /* Returns a node {s'} which is like {s} except that 
        it has no subnode with i-mark {i}. */
         
    rdag_node_t do_set(rdag_node_t s)
      {
        if (s == rdag_node_NULL)
          { return rdag_node_from_fields(D, s, i, v.ac, v.nd); }
        else 
          { rdag_node_data_t sdt; /* Data fields of {s}. */
            rdag_node_data_get(A->D, s, &sdt);
            if (sdt.i_mark < i)
              { return rdag_node_from_fields(D, s, i, v.ac, v.nd); }
            else if (sdt.i_mark == i)
              { return rdag_node_from_fields(D, sdt.f_link, i, v.ac, v.nd); }
            else
              { sdt.f_link = do_set(sdt.f_link);
                return rdag_node_from_data(D, &sdt);
              }
          }
      }
      
    rdag_node_t do_del(rdag_node_t s)
      {
        if (s == rdag_node_NULL)
          { return s; }
        else 
          { rdag_node_data_t sdt; /* Data fields of {s}. */
            rdag_node_data_get(A->D, s, &sdt);
            if (sdt.i_mark < i)
              { return s; }
            else if (sdt.i_mark == i)
              { return sdt.f_link; }
            else
              { sdt.f_link = do_del(sdt.f_link);
                return rdag_node_from_data(D, &sdt);
              }
          }
      }
      
    /* Choose between add or delete: */
    rdag_node_t r;
    if ((v.ac == 0) && (v.nd == rdag_node_NULL))
      { r = do_del(u.nd); }
    else
      { r = do_set(u.nd); }
    return (raut_state_t){ .ac = u.ac, .nd = r };
  }

raut_state_t raut_arc_follow(raut_t *A, raut_state_t u, rdag_symbol_t i)
  {
    demand(u.ac <= 1, "invalid accept bit");
    demand(u.nd <= raut_node_max(A), "invalid state node");
    
    rdag_node_t s = rdag_subnode_find(A->D, u.nd, i);
    raut_state_t v;
    if (s == rdag_node_NULL)
      { v = raut_state_VOID; }
    else
      { v.ac = rdag_o_mark(A->D, s); 
        v.nd = rdag_p_link(A->D, s);
        /* Paranoia check of {Suff(u) != {}} invariant: */
        assert((v.ac == 1) || (v.nd != rdag_node_NULL));
      }
    return v;
  }

rdag_node_t raut_node_max_alloc(raut_t *A)
  {
    return rdag_node_max_alloc(A->D);
  }

uint32_t raut_reachable_node_count(raut_t *A, raut_state_t u)
  {
    rdag_node_t root[1] = { u.nd };
    return rdag_reachable_node_count(A->D, 1, root);
  }

void raut_crunch(raut_t *A)
  {
    /* Get the root state of {A}: */
    raut_state_t root_state = raut_root_get(A);
    /* Crunch, retaining the root node: */
    rdag_node_t root_node[1] = { root_state.nd };
    rdag_crunch(A->D, 1, root_node);
    /* Update the node f the root state of {A}: */
    root_state.nd = root_node[0];
    raut_root_set(A, root_state);
  }

void raut_expand(raut_t *A, rdag_node_t max_alloc_node)
  {
    rdag_expand(A->D, max_alloc_node);
  }

void raut_free(raut_t *A)
  { 
    rdag_free(A->D);
    free(A);
  }
