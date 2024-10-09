/* See {tosl_arc_list.h} */
/* Last edited on 2024-10-08 22:55:46 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <haf.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_arc_list.h>

tosl_arc_id_t tosl_arc_list_merge
  ( tosl_arc_id_t *L0_P,
    tosl_arc_id_t *L1_P,
    tosl_coord_t Zp,
    tosl_mesh_t *mesh
  )
  {
    tosl_arc_id_t A = *L1_P; (*L1_P) = -1;

    while ((*L0_P) != -1)
      { tosl_arc_id_t ia = tosl_arc_list_pop(L0_P, mesh);
        tosl_coord_t Zorg = mesh->Vpos[mesh->Arc[ia].ivorg].c[2]; /* {Z} of origin of {ia}.  */
        tosl_arc_id_t ja = tosl_sym(ia);
        tosl_coord_t Zdst = mesh->Vpos[mesh->Arc[ja].ivorg].c[2]; /* {Z} of destination of {ia}.  */
        assert(Zorg < Zdst); /* Arcs must be up-going. */
        assert(Zorg < Zp); /* Arcs must start below {Zp}. */
        if (Zdst > Zp)
          { tosl_arc_list_add(&A, ia, mesh); }
        else
          { assert(Zdst < Zp); /* No plane through vertex. */ }
      }
    (*L0_P) = -1;
    return A;
  }

int32_t tosl_arc_list_len(tosl_arc_id_t L, tosl_mesh_t *mesh)
  { int32_t N = 0;
    if (L != -1)
      { tosl_arc_id_t ia = L;
        do { N++; ia = mesh->Arc[ia].succ; } while (ia != L);
      }
    return N;
  }

tosl_arc_id_t tosl_arc_list_pop(tosl_arc_id_t *L_P, tosl_mesh_t *mesh)
  { 
    tosl_arc_id_t ia = (*L_P);
    assert(ia != -1);
    tosl_arc_id_t ia_pred = mesh->Arc[ia].pred;
    tosl_arc_id_t ia_succ = mesh->Arc[ia].succ;
    if (ia_pred == ia)
      { /* Last element in list: */
        assert(ia_succ == ia);
        (*L_P) = -1;
      }
    else
      { /* At least two elements: */
        mesh->Arc[ia_pred].succ = ia_succ;
        mesh->Arc[ia_succ].pred = ia_pred;
        mesh->Arc[ia].pred = ia;
        mesh->Arc[ia].succ = ia;
        (*L_P) = ia_succ;
      }
    return ia;
  }
  
void tosl_arc_list_remove(tosl_arc_id_t *L_P, tosl_arc_id_t ia, tosl_mesh_t *mesh)
  { 
    assert(ia != -1);
    tosl_arc_id_t ia_pred = mesh->Arc[ia].pred;
    tosl_arc_id_t ia_succ = mesh->Arc[ia].succ;
    if (ia_pred == ia)
      { /* Last element in list: */
        assert(ia_succ == ia);
        assert((*L_P) == ia);
        (*L_P) = -1;
      }
    else
      { /* At least two elements: */
        mesh->Arc[ia_pred].succ = ia_succ;
        mesh->Arc[ia_succ].pred = ia_pred;
        mesh->Arc[ia].pred = ia;
        mesh->Arc[ia].succ = ia;
        if ((*L_P) == ia) { (*L_P) = ia_succ; }
      }
  }
  
void tosl_arc_list_add(tosl_arc_id_t *L_P, tosl_arc_id_t ia, tosl_mesh_t *mesh)
  { /* Require {ia} to be a singleton list: */
    assert(ia != -1);
    assert(mesh->Arc[ia].pred == ia);
    assert(mesh->Arc[ia].succ == ia);
    
    tosl_arc_id_t ja = (*L_P);
    if (ja == -1)
      { /* List {L} is empty, becomes {ia}: */
        (*L_P) = ia;
      }
    else
      { /* List {L} is not empty, insert before head: */
        tosl_arc_id_t ka = mesh->Arc[ja].pred;
        assert(mesh->Arc[ka].succ == ja);
        mesh->Arc[ka].succ = ia;
        mesh->Arc[ja].pred = ia;
        mesh->Arc[ia].pred = ka;
        mesh->Arc[ia].succ = ja;
      }
  }
