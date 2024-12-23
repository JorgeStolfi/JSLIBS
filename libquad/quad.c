/* See quad.h. */
/* Last edited on 2024-12-22 11:03:27 by stolfi */
  
/* Written by J. Stolfi in Apr/1993, based on an original
  implementation by Jim Roth (DEC CADM Advanced Group, May/1986).
  Changed by J. Stolfi in Dec/2011 to be procedure-based instead of
  macro-based. */ 

#define quad_C_copyright \
  "Copyright © 2011 State University of Campinas (UNICAMP).\n\n" jslibs_copyright "\n\n" \
  "NOTE: this copyright notice does not claim to supersede any copyrights" \
  " that may apply to the original DEC implementation of the quad-edge" \
  " data structure."

#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <vec.h>
#include <bool.h>
#include <jswsize.h>
#include <affirm.h>
#include <jsmath.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>

#include <quad.h>

struct quad_rep_t {};

#define UADDR(x) ((uaddress_t)(x)) 
  /* An {oct_arc_t} as an unsigned integer. */

#define TMASK ((uaddress_t)3u) 
  /* Mask for the tumble code of an {oct_arc_t}. */

#define EMASK (~ ((uaddress_t)3u)) 
  /* Mask for the edge address an {oct_arc_t}. */

#define EDGE(e) ((quad_edge_t)(UADDR(e) & EMASK))
  /* Edge record address from arc {e}. */

quad_edge_t quad_edge(quad_arc_t e)
  { return EDGE(e); }

#define TUMBLECODE(e) (UADDR(e) & TMASK)
  /* The tumble code of an arc {e} (as an address-size word). */

#define BASEARC(ed) ((quad_arc_t)(UADDR(ed) & EMASK))
  /* The base arc of the edge record {ed}. */

#define ORIENT(ed,tc) ((quad_arc_t)(UADDR(ed) | ((tc) & TMASK)))
  /* The base arc of the edge record {ed}. */

quad_bits_t quad_tumble_code(quad_arc_t e)
  { return TUMBLECODE(e); }

quad_arc_t quad_orient(quad_edge_t ed, quad_bits_t tc)
  { return ORIENT(ed,tc); }

bool_t quad_arc_is_null(quad_arc_t e)
  { return (EDGE(e) == NULL); }

#define LONBIT(e) (quad_bits_t)(((UADDR(e) & 2u) >> 1))
#define DOPBIT(e) (quad_bits_t)((UADDR(e) & 1u))
  /*  
    In this implementation, the tumble code {tc} of an arc {e} with base
    arc {b} consists of two bits {s=LONBIT(e)} and {d=DOPBIT(e)}. 
    They are such that {e = b.sym^s.rot^d}, or  {e = b.rot^{2*s+d}}. */

quad_bits_t quad_lon_bit(quad_arc_t e)
  { return LONBIT(e); }

quad_bits_t quad_dop_bit(quad_arc_t e)
  { return DOPBIT(e); }

/* Edge orientation operators: */

#define ROT(e) ((quad_arc_t)((UADDR(e) & EMASK) + ((UADDR(e)+1u) & TMASK)))
#define SYM(e) ((quad_arc_t)((UADDR(e) & EMASK) + ((UADDR(e)+2u) & TMASK)))
#define TOR(e) ((quad_arc_t)((UADDR(e) & EMASK) + ((UADDR(e)+3u) & TMASK)))

quad_arc_t quad_nop(quad_arc_t e)
  { return e; }

quad_arc_t quad_sym(quad_arc_t e)
  { return SYM(e); }

quad_arc_t quad_rot(quad_arc_t e)
  { return ROT(e); }

quad_arc_t quad_tor(quad_arc_t e)
  { return TOR(e); }

typedef struct quad_edge_rec_t { /* Edge record (four arcs). */
    quad_arc_t next[4];  /* Topological links. */
    void *data[4];       /* Client data fields. */
    uint64_t mark;       /* For internal use; set to 0 when record is created. */
    uint64_t eid;        /* Edge record number. */
  } quad_edge_rec_t;

/* Vertex/face walking operators (internal): */
#define ONEXT(e)    (EDGE(e)->next[(UADDR(e)+0u) & TMASK])
#define ROTRNEXT(e) (EDGE(e)->next[(UADDR(e)+1u) & TMASK])
#define SYMDNEXT(e) (EDGE(e)->next[(UADDR(e)+2u) & TMASK])
#define TORLNEXT(e) (EDGE(e)->next[(UADDR(e)+3u) & TMASK])

#define RNEXT(e) (TOR(ROTRNEXT(e)))
#define DNEXT(e) (SYM(SYMDNEXT(e)))
#define LNEXT(e) (ROT(TORLNEXT(e)))

#define OPREV(e) (ROT(ROTRNEXT(e)))
#define DPREV(e) (TOR(TORLNEXT(e)))
#define RPREV(e) (SYMDNEXT(e))
#define LPREV(e) (SYM(ONEXT(e)))

/* Function versions of the walking macros: */

quad_arc_t quad_onext(quad_arc_t e)
  { return ONEXT(e); }

quad_arc_t quad_oprev(quad_arc_t e)
  { return OPREV(e); }

quad_arc_t quad_dnext(quad_arc_t e)
  { return DNEXT(e); }

quad_arc_t quad_dprev(quad_arc_t e)
  { return DPREV(e); }

quad_arc_t quad_lnext(quad_arc_t e)
  { return LNEXT(e); }

quad_arc_t quad_lprev(quad_arc_t e)
  { return LPREV(e); }

quad_arc_t quad_rnext(quad_arc_t e)
  { return RNEXT(e); }

quad_arc_t quad_rprev(quad_arc_t e)
  { return RPREV(e); }

quad_arc_t quad_walk(quad_arc_t e, int32_t r, int32_t n)
  { /* Reduce {r} mod 4: */
    r = (r & 3);
    /* Apply {rot^r}: */
    switch(r)
      { case 0: /* nop */ break;
        case 1: e = ROT(e); break;
        case 2: e = SYM(e); break;
        case 3: e = TOR(e); break;
        default: assert(FALSE);
      }
    /* Apply {onext^n} or {oprev^{-n}}: */
    while(n > 0) { e = ONEXT(e); n--; }
    while(n < 0) { e = OPREV(e); n++; }
    /* Apply {tor^r}: */
    switch(r)
      { case 0: /* nop */ break;
        case 1: e = TOR(e); break;
        case 2: e = SYM(e); break;
        case 3: e = ROT(e); break;
        default: assert(FALSE);
      }
    /* We are done; */
    return e;
  }

/* Data pointers: */

#define ODATA(e) (EDGE(e)->data[(UADDR(e)+0u) & TMASK])
#define RDATA(e) (EDGE(e)->data[(UADDR(e)+1u) & TMASK])
#define DDATA(e) (EDGE(e)->data[(UADDR(e)+2u) & TMASK])
#define LDATA(e) (EDGE(e)->data[(UADDR(e)+3u) & TMASK])

void * quad_odata(quad_arc_t e)
  { return ODATA(e); }

void * quad_rdata(quad_arc_t e)
  { return RDATA(e); }

void * quad_ddata(quad_arc_t e)
  { return DDATA(e); }

void * quad_ldata(quad_arc_t e)
  { return LDATA(e); }

void quad_set_odata(quad_arc_t e, void *p)
  { ODATA(e) = p; }

void quad_set_rdata(quad_arc_t e, void *p)
  { RDATA(e) = p; }

void quad_set_ddata(quad_arc_t e, void *p)
  { DDATA(e) = p; }

void quad_set_ldata(quad_arc_t e, void *p)
  { LDATA(e) = p; }

/* Orientation bits: */

#define SYMBIT(e)   ((UADDR(e) & 2u)>>1)  
#define ROTBIT(e)   (UADDR(e) & 1u)
#define QUADBITS(e) (UADDR(e) & TMASK)
  /* Bits that distinguish the various flavors of the same edge:
      {SYMBIT(e)} distinguishes {e} from {SYM(e)}.
      {ROTBIT(e)} distinguishes {{e,SYM(e)}} from {{ROT(e),TOR(e)}}.
      {QUADBITS(e) = 2*SYMBIT(e) + ROTBIT(e)} distinguishes the four arcs
        {ROT^k(e)} from each other.
    Nothing can be assumed between {QUADBITS(e)} and {QUADBITS(ONEXT(e))}. */

/* Make a new edge: */

#define MARK(e)  (EDGE(e)->mark)
  /* The {mark} field of the edge record of arc [e}. */
  
#define EDGE_ID(e)  (EDGE(e)->eid)
  /* The {eid} field of the edge record of arc [e}. */
  
quad_arc_t quad_make_edge(void)
  {
    quad_arc_t e;

    e = (quad_arc_t) malloc(sizeof(quad_edge_rec_t));
    ONEXT(e) = e;
    SYMDNEXT(e) = SYM(e);
    ROTRNEXT(e) = TOR(e);
    TORLNEXT(e) = ROT(e);
    ODATA(e) = NULL;
    DDATA(e) = NULL;
    LDATA(e) = NULL;
    RDATA(e) = NULL;
    MARK(e) = 0;
    EDGE_ID(e) = 0;
    return e;
  }

/* Delete an edge: */

void quad_destroy_edge(quad_arc_t e)
  {
    quad_arc_t f = SYM(e);
    if (ONEXT(e) != e) quad_splice(e, OPREV(e));
    if (ONEXT(f) != f) quad_splice(f, OPREV(f));  
    free(EDGE(e));
  }

/* Edge numbers: */

quad_edge_id_t quad_edge_id(quad_edge_t ed)
  { 
    return ed->eid;
  }

void quad_set_edge_id(quad_edge_t ed, quad_edge_id_t eid)
  { 
    ed->eid = eid;
  }

/* Splice primitive: */

void quad_splice(quad_arc_t a, quad_arc_t b)
  {
    quad_arc_t ta = ONEXT(a);
    quad_arc_t tb = ONEXT(b);
    
    quad_arc_t u = ROT(ta);
    quad_arc_t v = ROT(tb);

    quad_arc_t tu = ONEXT(u);
    quad_arc_t tv = ONEXT(v);
    
    ONEXT(a) = tb;
    ONEXT(b) = ta;
    ONEXT(u) = tv;
    ONEXT(v) = tu;    
    /* Checking: */
    assert(quad_onext(a) == tb);
    assert(quad_onext(b) == ta);
    assert(quad_onext(u) == tv);
    assert(quad_onext(v) == tu);
  }

/* Enumerate edge quads */

void quad_do_enum(quad_arc_t a, void visit_proc(quad_arc_t e), uint64_t mark);
  /* Enumerates all primal edges reachable from {a}, setting their {mark} fields to {mark}.
    Assumes that any edge that has the {mark} field equal to {mark} has already been 
    visited. */

uint64_t next_mark = 1;  /* An integer greater than any edge's {mark} field. */

void quad_enum(quad_arc_vec_t *root, void visit_proc(quad_arc_t e))
  {
    uint64_t mark = next_mark;
    assert(mark != 0);
    next_mark++;
    for (uint32_t i = 0;  i < root->ne; i++) 
      { quad_do_enum(root->e[i], visit_proc, mark); }
  }

void quad_do_enum (quad_arc_t e, void visit_proc(quad_arc_t e), uint64_t mark)
  {
    while (MARK(e) != mark)
      { visit_proc(e);
        MARK(e) = mark;
        quad_do_enum (ONEXT(SYM(e)), visit_proc, mark);
        e = ONEXT(e);
      }
  }

quad_edge_id_t quad_renumber_edges(quad_arc_vec_t *root, quad_arc_vec_t *A)
  {
    quad_edge_id_t ne = 0;
    
    auto void renumber_edge(quad_arc_t p);
      /* Visit-proc that renumbers the base edge of {p} with {ne}
        and increments {ne}. */

    quad_enum(root, &renumber_edge);
    if (A != NULL) { quad_arc_vec_trim(A, (uint32_t)ne); }
    return ne;

    /* IMPLEMENTATIONS OF LOCAL PROCS */
    
    void renumber_edge(quad_arc_t p)
      { quad_edge_t ed = quad_edge(p);
        ed->eid = ne; 
        if (A != NULL) { quad_arc_vec_expand(A, (int32_t)ne); A->e[ne] = p; }
        ne++;
      }
  }

/* ARC INPUT/OUTPUT */

void quad_write_arc(FILE *wr, quad_arc_t e, uint32_t width)
  { fprintf(wr, 
      ("%*" uint64_u_fmt ":%u%u"), width, 
      EDGE(e)->mark, 
      (uint32_t)SYMBIT(e), (uint32_t)DOPBIT(e));
  }

quad_arc_t quad_read_arc(FILE *rd, quad_arc_vec_t *A)
  { /* Parse the edge number {eid} and get the edge {ed} from the edge table {ET}: */
    uint64_t eid = fget_uint64(rd, 10);
    demand(eid < A->ne, "invalid edge number");
    quad_edge_t ed = quad_edge(A->e[eid]);
    /* Parse ":" and the tumble code {tc} in binary: */
    fget_match(rd, ":");
    uint64_t tc = fget_uint64(rd, 2);
    demand(tc <= 3, "invalid tumble code");
    /* Put it all together: */
    return quad_orient(ed, (quad_bits_t)tc);
  }

/* MAP INPUT/OUTPUT */

#define FILE_TYPE "quad-edge"
#define FILE_VERSION "2011-12-22"

void quad_write_map(FILE *wr, quad_arc_vec_t *root, quad_arc_vec_t *A)
  { 
    /* Write the header line: */
    filefmt_write_header(wr, FILE_TYPE, FILE_VERSION);
    /* Grab the number {nr} of roots and the root index width {wa}: */
    uint32_t nr = root->ne;
    /* Compute the width in digits {wr} of the root index: */
    uint32_t dr = (nr < 10 ? 1 : digits(nr-1));
    /* Make sure that {A} is non-null and points to an empty table: */
    quad_arc_vec_t A_local = quad_arc_vec_new(0); /* A local edge table. */
    if (A == NULL) { A = &A_local; }
    
    /* Renumbers all edge records reachable from {root[0..nr-1]} and saves them in {A}: */
    quad_edge_id_t ne = quad_renumber_edges(root, A);
    assert(ne == A->ne);
    
    /* Determine the width {eE} in digits of the edge number: */
    uint32_t dE = (ne < 10 ? 1 : digits(ne-1));
    /* We should have zero edges if and only if we have zero roots: */
    assert((nr == 0) == (ne == 0)); 
    /* Write the root and edge counts: */
    fprintf(wr, "roots = %d\n", nr);
    fprintf(wr, "edges = %lu\n", ne);
    /* Write the roots, one per line: */
    for (uint32_t i = 0;  i < nr; i++)
      { /* Write {i} and the root arc number {i}: */
        fprintf(wr, "%*u ", dr, i);
        quad_write_arc(wr, root->e[i], dE);
        fputc('\n', wr);
      }
    /* Write the edge records, one per line: */
    for (int64_t eid = 0; eid < ne; eid++)
      { /* Get the reachable edge number {eid}: */
        quad_edge_t ed = quad_edge(A->e[eid]);
        /* Check whether the renumbering worked: */
        assert(ed->eid == eid);
        /* Write the edge's number and its {onext} links: */
        fprintf(wr, ("%*" uint64_u_fmt), dE, ed->eid);
        for (uint32_t r = 0;  r < 4; r++)
          { fputc(' ', wr); quad_write_arc(wr, ed->next[r], dE); }
        fputc('\n', wr);
      }
    /* Write the footer line: */
    filefmt_write_footer(wr, FILE_TYPE);
    /* Reclaim the edge table if local: */
    if (A == &A_local) { quad_arc_vec_trim(&A_local, 0); }
    return;
    
  }

void quad_read_map(FILE *rd, quad_arc_vec_t *root, quad_arc_vec_t *A)
  { /* Check and consume the header line: */
    filefmt_read_header(rd, FILE_TYPE, FILE_VERSION);
    /* Parse the root count {nr} and the edge count {ne}: */
    uint32_t nr = nget_uint32(rd, "roots", 10); fget_eol(rd);
    uint32_t ne = nget_uint32(rd, "edges", 10); fget_eol(rd);
    /* Make sure that {A} is non-null and points to a table with {ne} slots: */
    quad_arc_vec_t A_local = quad_arc_vec_new(0); /* A local edge table. */
    if (A == NULL) { A = &A_local; }
    quad_arc_vec_trim(A, ne);
    /* Create the edge records, save their base arcs in {A->e[0..ne-1]}: */
    for (uint64_t eid = 0; eid < ne; eid++) 
      { quad_arc_t a = quad_make_edge(); 
        EDGE(a)->eid = eid; 
        A->e[eid] = a;
      }
    /* (Re)alocate the root record and read the root arcs, one per line: */
    quad_arc_vec_trim(root, nr);
    for (int64_t i = 0; i < nr; i++)
      { /* Parse the root index {i} and the root arc {root->e[i]}: */
        uint64_t iread = fget_uint64(rd, 10);
        demand(iread == i, "root index mismatch");
        /* Parse the root arc number {i}, save it in {root}: */
        root->e[i] = quad_read_arc(rd, A); 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    /* Read the contents of the edge records {0..ne-1}: */
    for (int64_t eid = 0; eid < ne; eid++) 
      { /* Parse the edge number {eid}: */
        uint64_t eid_read = fget_uint64(rd, 10);
        demand(eid_read == eid, "edge number mismatch");
        /* Get the edge {ed} from the edge table {A}: */
        quad_edge_t ed = quad_edge(A->e[eid]);
        /* Read its links {ed->next[0..3]}: */
        for (uint32_t r = 0;  r < 4; r++) { ed->next[r] = quad_read_arc(rd, A); } 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    
    /* Check and consume the footer line: */
    filefmt_read_footer(rd, FILE_TYPE);
    /* Reclaim the edge table if local: */
    if (A == &A_local) { quad_arc_vec_trim(&A_local, 0); }
  }

vec_typeimpl(quad_arc_vec_t,quad_arc_vec,quad_arc_t);
vec_typeimpl(quad_edge_vec_t,quad_edge_vec,quad_edge_t);
