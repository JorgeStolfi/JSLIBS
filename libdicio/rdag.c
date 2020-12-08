/* Procedures of {rdag.h} that do not require access to the internal rep. */
/* Last edited on 2009-10-29 22:04:23 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <bool.h>
#include <vec.h>

#include <rdag.h>

rdag_node_t rdag_subnode_find(rdag_t *D, rdag_node_t s, rdag_symbol_t i)
  {
    /* !!! Assumes that transitions are ordered. Should we remove that invariant? !!! */
    while (s != rdag_node_NULL)
      { rdag_symbol_t j = rdag_i_mark(D, s);
        if (j < i) 
          { return rdag_node_NULL; }
        else if (j > i) 
          { s = rdag_f_link(D, s); }
        else
          { return s; }
      }
    return rdag_node_NULL;
  }

uint32_t rdag_subnode_count(rdag_t *D, rdag_node_t s)
  {
    uint32_t n = 0;
    while (s != rdag_node_NULL) { n++; s = rdag_f_link(D, s); }
    return n;
  }

rdag_node_t rdag_subnode_first(rdag_t *D, rdag_node_t s)
  {
    demand(s != rdag_node_NULL, "node is null");
    rdag_node_t p = rdag_node_NULL;
    while (s != rdag_node_NULL) { p = s; s = rdag_f_link(D, s); }
    assert(p != rdag_node_NULL);
    return p;
  }

rdag_disp_t rdag_enum_paths
  ( rdag_t *D, 
    rdag_symbol_t o_root,
    rdag_node_t t_root,
    rdag_node_action_t *enter,
    rdag_p_step_action_t *push, 
    rdag_p_step_action_t *pop,
    rdag_node_action_t *exit
  )
  {
    demand(t_root <= rdag_node_max(D), "invalid initial symbol");
    demand(o_root <= rdag_o_mark_max(D), "invalid initial node");
   
    uint32_t len = 0;
    
    /* INTERNAL PROTOTYPES */

    auto rdag_disp_t do_enum_paths(rdag_symbol_t o, rdag_node_t s);
      /* Does {rdag_enum_paths} starting at {s} (a generic sucessor of
        {t_root}, reached through the output mark {o}), with the path
        length starting at {len} intead of 0. Returns {rdag_disp_STOP}
        if the enumeration should be aborted, and {rdag_disp_FINE}
        otherwise. */

    auto rdag_disp_t do_enum_rest(rdag_node_t r);
      /* Enumerates subpaths that start at {r}, a subnode of some node
        {s}, and at all subnodes of {r}. Returns {rdag_disp_SKIP}
        if the remaining arcs out of {s} should be skipped;
        {rdag_disp_STOP} if the enumeration should be aborted; and
        {rdag_disp_FINE} otherwise. */

    rdag_disp_t d_paths_root = do_enum_paths(o_root, t_root);
    switch(d_paths_root)
      { case rdag_disp_STOP: return d_paths_root;
        case rdag_disp_SKIP: assert(FALSE);
        case rdag_disp_FINE: break;
        default: assert(FALSE);
      }
    return rdag_disp_FINE;
    
    /* INTERNAL IMPLEMENTATIONS */

    rdag_disp_t do_enum_paths(rdag_symbol_t o, rdag_node_t s)
      { 
        rdag_disp_t d_enter = (enter == NULL ? rdag_disp_FINE : enter(len, o, s));
        switch(d_enter)
          { case rdag_disp_STOP: return d_enter;
            case rdag_disp_SKIP: break;
            case rdag_disp_FINE: 
              { rdag_disp_t d_rest = do_enum_rest(s);
                switch(d_rest)
                  { case rdag_disp_STOP: return d_rest;
                    case rdag_disp_SKIP: break;
                    case rdag_disp_FINE: break;
                    default: assert(FALSE);
                  }
              }
              break;
            default: assert(FALSE);
          }
        rdag_disp_t d_exit = (exit == NULL ? rdag_disp_FINE : exit(len, o, s));
        switch(d_exit)
          { case rdag_disp_STOP: return d_exit;
            case rdag_disp_SKIP: break;
            case rdag_disp_FINE: break;
            default: assert(FALSE);
          }
        
        return rdag_disp_FINE;
      }

    rdag_disp_t do_enum_rest(rdag_node_t r)
      {
        if (r == rdag_node_NULL){ return rdag_disp_FINE; }
        /* Enumerate subpaths that begin with all the other arcs out of {s}: */
        rdag_disp_t d_rest = do_enum_rest(rdag_f_link(D,r));
        switch(d_rest)
          { case rdag_disp_STOP: return d_rest;
            case rdag_disp_SKIP: return d_rest;
            case rdag_disp_FINE: break;
            default: assert(FALSE);
          }

        /* Now enumerate those subpaths that start with this p-step out of {s}: */
        rdag_symbol_t ir = rdag_i_mark(D,r);
        rdag_symbol_t ot = rdag_o_mark(D,r);
        rdag_node_t t = rdag_p_link(D,r);
        rdag_disp_t d_push = (push == NULL ? rdag_disp_FINE : push(len, r, ir, ot, t));
        switch(d_push)
          { case rdag_disp_STOP: return d_push;
            case rdag_disp_SKIP: break;
            case rdag_disp_FINE: 
              { /* Enumerate the subpaths that start beyond the p-step: */
                len++;
                rdag_disp_t d_paths = do_enum_paths(ot, t);
                switch(d_paths)
                  { case rdag_disp_STOP: len--; return d_paths;
                    case rdag_disp_SKIP: assert(FALSE);
                    case rdag_disp_FINE: break;
                    default: assert(FALSE);
                  }
              }
              break;
            default: assert(FALSE);
          }
        len--;
        rdag_disp_t d_pop = (pop == NULL ? rdag_disp_FINE : pop(len, r, ir, ot, t));
        return d_pop;
      }
  }

rdag_node_t *rdag_copy(rdag_t *OLD, rdag_t *NEW, rdag_node_t old_s, rdag_node_t map[])
  {
    /* We could do a linear scan instead of a recursive enumeration, */
    /*   but then time would be proportional to the total DAG size, */
    /*   instead of size of reachable part. */

    /* Make sure that we have a node map: */
    if (map == NULL)
      { /* Map array not given, allocate one and fill it with {rdag_node_NULL}: */
        uint32_t map_size = rdag_node_max(OLD);
        map = (rdag_node_t *)notnull(malloc(map_size * sizeof(rdag_node_t)), "no mem");
        int k;
        for (k = 0; k < map_size; k++) { map[k] = rdag_node_NULL; }
      }
    assert(map != NULL);

    auto rdag_node_t do_copy(rdag_node_t old_t);
      /* Copies the node {old_t} of {OLD} and all its descendants, as needed.
        Returns the corresponding node in {NEW}. */
      
    rdag_node_t do_copy(rdag_node_t old_t)
     {
       if (old_t == rdag_node_NULL) { return rdag_node_NULL; }
       if (map[old_t-1] != rdag_node_NULL) { return map[old_t-1]; }
       
       /* Get the old fields of {old_t}: */
       rdag_node_data_t old_dt;
       rdag_node_data_get(OLD, old_t, &old_dt);
       
       /* Get the data of the corresponding new node: */
       rdag_node_data_t new_dt;
       new_dt.f_link = do_copy(old_dt.f_link);
       new_dt.i_mark = old_dt.i_mark;
       new_dt.o_mark = old_dt.o_mark;
       new_dt.p_link = do_copy(old_dt.p_link);
       
       /* Add the new node to {NEW}: */
       rdag_node_t new_t = rdag_node_from_data(NEW, &new_dt);
       map[old_t-1] = new_t;
       return new_t;
     }

    (void)do_copy(old_s);
    return map;
  }

vec_typeimpl(rdag_symbol_vec_t, rdag_symbol_vec, rdag_symbol_t);
