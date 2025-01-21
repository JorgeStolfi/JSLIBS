/* See oct.h. */
/* Last edited on 2025-01-09 23:01:04 by stolfi */

/* This implementation was originally created by J. Stolfi in Apr/1993.
  It was based on the orientable-map version {quad.c} implemented by
  Jim Roth (DEC CADM Advanced Group) in May/1986, which was the basis for
  the full Modula-3 file {Oct.m3} by Rober M. Rosi in 1993.
  The latter was converted to C by J. Stolfi in 1996.  It was substantially
  revised by J. Stolfi between Jan/2007 and Mar/2009, to add support
  general orbit enumeration and other enhancements. */

#define oct_C_copyright \
  "Copyright © 1996, 2006 State University of Campinas (UNICAMP).\n\n" jslibs_copyright "\n\n" \
  "NOTE: this copyright notice does not claim to supersede any copyrights" \
  " that may apply to the original DEC implementation of the quad-edge" \
  " data structure."

#include <assert.h>
#include <malloc.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <vec.h>
#include <bool.h>
#include <jswsize.h>
#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>
#include <jsmath.h>

#include <oct.h>
#include <oct_enum.h>

#define UADDR(x) ((uaddress_t)(x)) 

#define TMASK ((uaddress_t)7u) 
  /* Mask for the tumble code of an {oct_arc_t}. */

#define EMASK (~ TMASK) 
  /* Mask for the edge address an {oct_arc_t}. */

#define EDGE(e) ((oct_edge_t)(UADDR(e) & EMASK))
  /* Extracts the edge record from an arc reference. */

oct_edge_t oct_edge(oct_arc_t e)
  { return EDGE(e); }

#define TUMBLECODE(e) (UADDR(e) & TMASK)
  /* The tumble code of an arc {e} (as an address-size word). */

#define BASEARC(E) ((oct_arc_t)(UADDR(E) & EMASK))
  /* The base arc of the edge record {ed}. */

#define ORIENT(ed,tc) ((oct_arc_t)(UADDR(ed) | ((tc) & TMASK)))
  /* The arc whose edge is {ed} tumble code {tc}. */

oct_bits_t oct_tumble_code(oct_arc_t e)
  { return TUMBLECODE(e); }

oct_arc_t oct_orient(oct_edge_t ed, oct_bits_t tc)
  { return ORIENT(ed,tc); }

bool_t oct_is_null(oct_arc_t e)
  { return (EDGE(e) == NULL); }

#define SYMBIT(e) ((oct_bits_t)((UADDR(e) & 4u) >> 2))  
#define DOPBIT(e) ((oct_bits_t)((UADDR(e) & 2u) >> 1))
#define FLPBIT(e) ((oct_bits_t)(UADDR(e) & 1u))
  /*  
    In this implementation, the tumble code {tc} of an arc {e} with base
    arc {b} consists of three bits {s=SYMBIT(e)}, {d=DOPBIT(e)}, and
    {f=FLPBIT(e)}. They are such that

      {e = b.sym^s.rot^d.fflip^f}, or 

      {e = b.rot^{2*s+d}.fflip^f}. */

#define ROTBITS(e) ((UADDR(e) >> 1) & 3u)
  /* The bits {s} and {d} together, i.e. the integer {r=2*s+d}.
    If {b} is the base edge of {e}, then {e = b.rot^r.fflip^f}. */

oct_bits_t oct_lon_bit(oct_arc_t e)
  { return SYMBIT(e); }

oct_bits_t oct_trn_bit(oct_arc_t e)
{ return (oct_bits_t)(SYMBIT(e) ^ FLPBIT(e)); }

oct_bits_t oct_dop_bit(oct_arc_t e)
  { return DOPBIT(e); }

oct_bits_t oct_cir_bit(oct_arc_t e)
  { return FLPBIT(e); }

#define SYMFIX(e) ((((UADDR(e) >> 1) ^ UADDR(e)) << 2) & TMASK)
  /* The exclusive-or of {ROTBIT} and {FLPBIT}, positioned at {SYMBIT}.
    A common sub-expression for the formulas below. */

#define NOP(e)   (e)
#define ROT(e)   ((oct_arc_t)(UADDR(e) ^ SYMFIX(e) ^ 2u))
#define SYM(e)   ((oct_arc_t)(UADDR(e) ^ 4u))
#define TOR(e)   ((oct_arc_t)(UADDR(e) ^ SYMFIX(e) ^ 6u))
#define FFLIP(e) ((oct_arc_t)(UADDR(e) ^ 1u))
#define VFLIP(e) ((oct_arc_t)(UADDR(e) ^ 5u))
#define DUAL(e)  ((oct_arc_t)(UADDR(e) ^ SYMFIX(e) ^ 3u))
#define DUAR(e)  ((oct_arc_t)(UADDR(e) ^ SYMFIX(e) ^ 7u))
  /* The basic tumbling operators, assuming that the tumble code is
    defined as above. They act on the tumble code as follows:

      arc              sdf tc | rot    sym    tor    fflip  vflip  dual   duar
      ---------------  --- - | -----  -----  -----  -----  -----  -----  -----
      e = b.nop        000 0 | 010 2  100 4  110 6  001 1  101 5  011 3  111 7     
      e = b.nop.fflip  001 1 | 111 7  101 5  011 3  000 0  100 4  110 6  010 2     
      e = b.rot        010 2 | 100 4  110 6  000 0  011 3  111 7  101 5  001 1     
      e = b.rot.fflip  011 3 | 001 1  111 7  101 5  010 2  110 6  000 0  100 4     
      e = b.sym        100 4 | 110 6  000 0  010 2  101 5  001 1  111 7  011 3     
      e = b.sym.fflip  101 5 | 011 3  001 1  111 7  100 4  000 0  010 2  110 6     
      e = b.tor        110 6 | 000 0  010 2  100 4  111 7  011 3  001 1  101 5     
      e = b.tor.fflip  111 7 | 101 5  011 3  001 1  110 6  101 2  100 4  000 0 

   Apart from column ordering, this is merely the multiplication table
   for the group of tumble operators:
   
            | nop    fflip  rot    dual   sym    vflip  tor    duar
      ----- | -----  -----  -----  -----  -----  -----  -----  -----
      nop   | 000 0  001 1  010 2  011 3  100 4  101 5  110 6  111 7     
      fflip | 001 1  000 0  111 7  110 6  101 5  100 4  011 3  010 2     
      rot   | 010 2  011 3  100 4  101 5  110 6  111 7  000 0  001 1     
      dual  | 011 3  010 2  001 1  000 0  111 7  110 6  101 5  100 4     
      sym   | 100 4  101 5  110 6  111 7  000 0  001 1  010 2  011 3     
      vflip | 101 5  100 4  011 3  010 2  001 1  000 0  111 7  110 6     
      tor   | 110 6  111 7  000 0  001 1  010 2  011 3  100 4  101 5     
      duar  | 111 7  110 6  101 5  100 4  011 3  101 2  001 1  000 0 

   Perhaps a different encoding can make the product operation into a
   single arithmetic/boolean operation... */

oct_arc_t oct_nop(oct_arc_t e)
  { return e; }

oct_arc_t oct_sym(oct_arc_t e)
  { return SYM(e); }

oct_arc_t oct_rot(oct_arc_t e)
  { return ROT(e); }

oct_arc_t oct_tor(oct_arc_t e)
  { return TOR(e); }

oct_arc_t oct_fflip(oct_arc_t e)
  { return FFLIP(e); }

oct_arc_t oct_vflip(oct_arc_t e)
  { return VFLIP(e); }

oct_arc_t oct_dual(oct_arc_t e)
  { return DUAL(e); }

oct_arc_t oct_duar(oct_arc_t e)
  { return DUAR(e); }

oct_arc_t oct_rot_fflip(oct_arc_t e, int32_t r, int32_t f)
  { /* Reduce {r} mod 8 (math mod): */
    r = r & 3;
    /* Apply {ROT} {r} times: */
    if (r == 1)
      { e = ROT(e); }
    else if (r == 3)
      { e = TOR(e); }
    else
      { e = SYM(e); }
    /* Reduce {f} mod 2 (math mod): */
    f = f & 1;
    /* Apply {FFLIP} {f} times: */
    if (f != 0) 
      { e = FFLIP(e); }
    return e;
  }

/* EDGE RECORD FIELDS */

typedef struct oct_edge_rec_t 
  { oct_arc_t next[4];   /* Topological links. */
    uint64_t eid;        /* Edge number (multipurpose). */
  } oct_edge_rec_t;
  /* An {oct_edge_rec_t} describes the connections of an unoriented,
    undirected edge {ed} of the map, as well as of its dual edge. 
    If {b} is the base edge of {ed}, then {next[i] = b.rot^i.onext}.
    That is,
    
      {next[0] = onext(b)}
      {next[1] = onext(rot(b))}
      {next[2] = onext(sym(b))}
      {next[3] = onext(tor(b))}
      
    */

#define NMASK ((((uint64_t)1) << 61) - 1)
  /* Mask for valid edge numbers, {0..2^61-1}. */

uint64_t oct_edge_id(oct_edge_t ed)
  { return ed->eid; }

void oct_set_edge_id(oct_edge_t ed, uint64_t eid)
  { demand((eid & (~ NMASK)) == 0, "edge number is too big"); 
    ed->eid = eid;
  }

uint64_t oct_arc_num(oct_arc_t *e)
  { return ((EDGE(e))->eid << 3) | TUMBLECODE(e); }

/* WALKING  OPERATORS 
  
  To understand the macros below, observe that
  {ROTBITS(FFLIP(e)) == ROTBITS(e)} 
  {ROTBITS(TOR(FFLIP(e))) == ROTBITS(ROT(e))} 
  {ROTBITS(ROT(FFLIP(e))) == ROTBITS(TOR(e))} 
  for any {e}. */

/* Origin ring: */

#define ONEXT(e) \
  ( FLPBIT(e) ? \
    DUAL((EDGE(e))->next[ROTBITS(DUAR(e))]) : \
    (EDGE(e))->next[ROTBITS(e)] \
  )

#define OPREV(e) \
  ( FLPBIT(e) ? \
    FFLIP((EDGE(e))->next[ROTBITS(e)]) : \
    ROT((EDGE(e))->next[ROTBITS(ROT(e))]) \
  )

/* Destination ring: */

#define DNEXT(e) \
  ( FLPBIT(e) ? \
    DUAR((EDGE(e))->next[ROTBITS(ROT(e))]) : \
    SYM((EDGE(e))->next[ROTBITS(SYM(e))]) \
  )

#define DPREV(e) \
  ( FLPBIT(e) ? \
    VFLIP((EDGE(e))->next[ROTBITS(SYM(e))]) : \
    TOR((EDGE(e))->next[ROTBITS(TOR(e))]) \
  )

/* Left face ring: */

#define LNEXT(e) \
  ( FLPBIT(e) ? \
    FFLIP((EDGE(e))->next[ROTBITS(SYM(e))]) : \
    ROT((EDGE(e))->next[ROTBITS(TOR(e))]) \
  )

#define LPREV(e) \
  ( FLPBIT(e) ? \
    ROT(FFLIP((EDGE(e))->next[ROTBITS(TOR(e))])) : \
    SYM((EDGE(e))->next[ROTBITS(e)]) \
  )

/* Right face ring: */

#define RNEXT(e) \
  ( FLPBIT(e) ? \
    VFLIP((EDGE(e))->next[ROTBITS(e)]) : \
    TOR((EDGE(e))->next[ROTBITS(ROT(e))]) \
  )

#define RPREV(e) \
  ( FLPBIT(e) ? \
    DUAL((EDGE(e))->next[ROTBITS(ROT(e))]) : \
    (EDGE(e))->next[ROTBITS(SYM(e))] \
  )

/* Function versions of the walking macros: */

oct_arc_t oct_onext(oct_arc_t e)
  { return ONEXT(e); }

oct_arc_t oct_oprev(oct_arc_t e)
  { return OPREV(e); }

oct_arc_t oct_dnext(oct_arc_t e)
  { return DNEXT(e); }

oct_arc_t oct_dprev(oct_arc_t e)
  { return DPREV(e); }

oct_arc_t oct_lnext(oct_arc_t e)
  { return LNEXT(e); }

oct_arc_t oct_lprev(oct_arc_t e)
  { return LPREV(e); }

oct_arc_t oct_rnext(oct_arc_t e)
  { return RNEXT(e); }

oct_arc_t oct_rprev(oct_arc_t e)
  { return RPREV(e); }

oct_arc_t oct_walk(oct_arc_t e, int32_t r, int32_t n)
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

/* COUNTING INCIDENCES */

uint oct_odegree(oct_arc_t e)
  { uint deg = 0;
    oct_arc_t p = e;
    do { deg++; p = oct_onext(p); } while (p != e);
    return deg;
  }

uint oct_ddegree(oct_arc_t e)
  { uint deg = 0;
    oct_arc_t p = e;
    do { deg++; p = oct_dnext(p); } while (p != e);
    return deg;
  }

uint oct_ldegree(oct_arc_t e)
  { uint deg = 0;
    oct_arc_t p = e;
    do { deg++; p = oct_lnext(p); } while (p != e);
    return deg;
  }

uint oct_rdegree(oct_arc_t e)
  { uint deg = 0;
    oct_arc_t p = e;
    do { deg++; p = oct_rnext(p); } while (p != e);
    return deg;
  }

/* EDGE CREATION AND DESTRUCTION */

#define ESIZE (sizeof(oct_edge_rec_t))
  /* Size of an edge record in bytes. */

oct_arc_t oct_make_edge(void)
  { /* Allocate a new edge record: */
    oct_edge_t ed = (oct_edge_t)notnull(malloc(ESIZE), "no mem");
    /* Get the edge's base arc: */
    oct_arc_t e = BASEARC(ed);
    /* Set the fields for a map on {S^2} with {e} as the single (non-loop) edge {e}: */
    ed->next[0] = NOP(e);
    ed->next[1] = TOR(e);
    ed->next[2] = SYM(e);
    ed->next[3] = ROT(e);
    ed->eid = 0;
    return e;
  }

void oct_destroy_edge(oct_arc_t e)
  { /* Grab the edge record {ed}: */
    oct_edge_t ed = EDGE(e);
    /* Make sure that the edge is isolated: */
    oct_arc_t f = SYM(e);
    if (EDGE(ed->next[0]) != ed) { oct_splice(e, OPREV(e)); }
    if (EDGE(ed->next[2]) != ed) { oct_splice(f, OPREV(f)); }
    /* Paranoia check -- assumes that the oct-edge is well-formed: */
    for (uint32_t i = 0; i < 4; i++) { assert(EDGE(ed->next[i]) == ed); }
    /* OK, reclaim the edge record: */
    free(ed);
  }

/* MAP SPLICING */

void oct_splice(oct_arc_t a, oct_arc_t b)
  {
    demand(b != FFLIP(ONEXT(a)), "invalid splice");
    if (a == b) { return; }
    
    oct_arc_t ta = ONEXT(a), tb = ONEXT(b);
    oct_arc_t c = ROT(ta), d = ROT(tb); 
    oct_arc_t tc = ONEXT(c), td = ONEXT(d);
    
    if (!FLPBIT(a))
      { EDGE(a)->next[ROTBITS(a)] = tb; }
    else
      { EDGE(a)->next[ROTBITS(DUAR(a))] = DUAL(tb); }
    
    if (!FLPBIT(b))
      { EDGE(b)->next[ROTBITS(b)] = ta; }
    else
      { EDGE(b)->next[ROTBITS(DUAR(b))] = DUAL(ta);  }
    
    if (!FLPBIT(c)) 
      { EDGE(c)->next[ROTBITS(c)] = td; }
    else
      { EDGE(c)->next[ROTBITS(DUAR(c))] = DUAL(td); }          
    
    if (!FLPBIT(d))
      { EDGE(d)->next[ROTBITS(d)] = tc; }
    else
      { EDGE(d)->next[ROTBITS(DUAR(d))] = DUAL(tc); }
      
    /* Checking: */
    assert(oct_onext(a) == tb);
    assert(oct_onext(b) == ta);
    assert(oct_onext(c) == td);
    assert(oct_onext(d) == tc);
  }

/* ARC INPUT/OUTPUT */

void oct_write_arc(FILE *wr, oct_arc_t e, uint32_t width)
  { oct_edge_t ed = EDGE(e);
    oct_bits_t tc = TUMBLECODE(e);
    fprintf
      ( wr, ("%*" uint64_u_fmt ":%u%u%u"), width, 
        ed->eid, (uint32_t)((tc & 4u) >> 2), (uint32_t)((tc & 2u) >> 1), (uint32_t)(tc & 1u)
      );
  }

oct_arc_t oct_read_arc(FILE *rd, oct_arc_vec_t *A)
  { /* Parse the edge number {eid} and get the edge {ed} from the edge table {ET}: */
    uint64_t eid = fget_uint64(rd, 10);
    demand(eid < A->ne, "invalid edge number");
    oct_edge_t ed = oct_edge(A->e[eid]);
    /* Parse ":" and the tumble code {tc} in binary: */
    fget_match(rd, ":");
    uint64_t tc = fget_uint64(rd, 2);
    demand(tc <= 7, "invalid tumble code");
    /* Put it all together: */
    return oct_orient(ed, (oct_bits_t)tc);
  }

/* MAP INPUT/OUTPUT */

#define FILE_TYPE "oct-edge"
#define FILE_VERSION "2007-01-16"

void oct_write_map(FILE *wr, oct_arc_vec_t *root, oct_arc_vec_t *A)
  { 
  
    uint32_t NE = 0;  /* Count of reachable octets (edge records). */
    
    auto bool_t renumber_edge(oct_arc_t p);
      /* Visit-proc that renumbers the base edge of {p} with {NE}
        and increments {NE}. */
    
    /* Write the header line: */
    filefmt_write_header(wr, FILE_TYPE, FILE_VERSION);
    /* Grab the number {NR} of roots and the root index width {wa}: */
    uint32_t NR = root->ne;
    /* Compute the width in digits {wr} of the root index: */
    uint32_t dr = (NR < 10 ? 1 : digits(NR-1));
    /* Make sure that {A} is non-null and points to an empty table: */
    oct_arc_vec_t A_local = oct_arc_vec_new(0); /* A local edge table. */
    if (A == NULL) 
      { A = &A_local; }
    else
      { oct_arc_vec_trim(A, 0); }
    /* Renumber all edge records reachable from {root[0..NR-1]}: */
    oct_enum_octets(*root, &renumber_edge, A);
    assert(A->ne == NE);
    /* Determine the width {eE} in digits of the edge number: */
    uint32_t dE = (NE < 10 ? 1 : digits(NE-1));
    /* We should have zero edges if and only if we have zero roots: */
    assert((NR == 0) == (NE == 0)); 
    /* Write the root and edge counts: */
    fprintf(wr, "roots = %d\n", NR);
    fprintf(wr, "edges = %d\n", NE);
    /* Write the roots, one per line: */
    uint32_t i;
    for (i = 0; i < NR; i++)
      { /* Write {i} and the root arc number {i}: */
        fprintf(wr, "%*u ", dr, i);
        oct_write_arc(wr, root->e[i], dE);
        fputc('\n', wr);
      }
    /* Write the edge records, one per line: */
    uint64_t eid;
    for (eid = 0; eid < NE; eid++)
      { /* Get the reachable edge number {eid}: */
        oct_edge_t ed = oct_edge(A->e[eid]);
        /* Check whether the renumbering worked: */
        assert(ed->eid == eid);
        /* Write the edge's number and its {onext} links: */
        fprintf(wr, ("%*" uint64_u_fmt), dE, eid);
        for (uint32_t r = 0; r < 4; r++)
          { fputc(' ', wr); oct_write_arc(wr, ed->next[r], dE); }
        fputc('\n', wr);
      }
    /* Write the footer line: */
    filefmt_write_footer(wr, FILE_TYPE);
    /* Reclaim the edge table if local: */
    if (A == &A_local) { oct_arc_vec_trim(&A_local, 0); }
    return;
    
    /* IMPLEMENTATIONS OF LOCAL PROCS */
    
    bool_t renumber_edge(oct_arc_t p)
      {
        oct_edge_t ed = oct_edge(p);
        ed->eid = NE; NE++;
        return FALSE;
      }
  }

void oct_read_map(FILE *rd, oct_arc_vec_t *root, oct_arc_vec_t *A)
  { /* Check and consume the header line: */
    filefmt_read_header(rd, FILE_TYPE, FILE_VERSION);
    /* Parse the root count {NR} and the edge count {NE}: */
    uint32_t NR = nget_uint32(rd, "roots", 10); fget_eol(rd);
    uint32_t NE = nget_uint32(rd, "edges", 10); fget_eol(rd);
    /* Make sure that {A} is non-null and points to a table with {NE} slots: */
    oct_arc_vec_t A_local = oct_arc_vec_new(0); /* A local edge table. */
    if (A == NULL) { A = &A_local; }
    oct_arc_vec_trim(A, NE);
    /* Create the edge records, save their base arcs in {A->e[0..NE-1]}: */
    uint64_t eid;
    for (eid = 0; eid < NE; eid++) 
      { A->e[eid] = oct_make_edge(); 
        oct_set_edge_id(oct_edge(A->e[eid]), eid);
      }
    /* (Re)alocate the root record and read the root arcs, one per line: */
    oct_arc_vec_trim(root, NR);
    uint64_t i;
    for (i = 0; i < NR; i++)
      { /* Parse the root index {i} and the root arc {root->e[i]}: */
        uint64_t iread = fget_uint64(rd, 10);
        demand(iread == i, "root index mismatch");
        /* Parse the root arc number {i}, save it in {root}: */
        root->e[i] = oct_read_arc(rd, A); 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    /* Read the contents of the edge records {0..NE-1}: */
    for (eid = 0; eid < NE; eid++) 
      { /* Parse the edge number {eid}: */
        uint64_t eid_read = fget_uint64(rd, 10);
        demand(eid_read == eid, "edge number mismatch");
        /* Get the edge {ed} from the edge table {A}: */
        oct_edge_t ed = oct_edge(A->e[eid]);
        /* Read its links {ed->next[0..3]}: */
        for (uint32_t r = 0; r < 4; r++) { ed->next[r] = oct_read_arc(rd, A); } 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    
    /* Check and consume the footer line: */
    filefmt_read_footer(rd, FILE_TYPE);
    /* Reclaim the edge table if local: */
    if (A == &A_local) { oct_arc_vec_trim(&A_local, 0); }
  }

vec_typeimpl(oct_arc_vec_t,oct_arc_vec,oct_arc_t);
vec_typeimpl(oct_edge_vec_t,oct_edge_vec,oct_edge_t);
