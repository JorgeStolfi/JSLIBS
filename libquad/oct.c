/* See oct.h. */
/* Last edited on 2016-04-01 02:58:18 by stolfilocal */

#define oct_C_copyright \
  "Copyright � 1996, 2006 Institute of Computing, Unicamp."

/* AUTHORS

  This implementation was originally created by J. Stolfi in April 1993.
  It was based on the orientable-map version {quad.c} implemented by
  Jim Roth (DEC CADM Advanced Group) in May 1986, which was the basis for
  the full Modula-3 file {Oct.m3} by Rober M. Rosi in 1993.
  The latter was converted to C by J. Stolfi in 1996.  It was substantially
  revised by J. Stolfi between January 2007 and March 2009, to add support
  general orbit enumeration and other enhancements. */

#define _GNU_SOURCE
#include <assert.h>
#include <malloc.h>
#include <stdint.h>

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
  /* The base arc of the edge record {E}. */

#define ORIENT(E,t) ((oct_arc_t)(UADDR(E) | ((t) & TMASK)))
  /* The base arc of the edge record {E}. */

oct_bits_t oct_tumble_code(oct_arc_t e)
  { return TUMBLECODE(e); }

oct_arc_t oct_orient(oct_edge_t E, oct_bits_t t)
  { return ORIENT(E,t); }

bool_t oct_is_null(oct_arc_t e)
  { return (EDGE(e) == NULL); }

#define SYMBIT(e) ((oct_bits_t)((UADDR(e) & 4u) >> 2))  
#define DOPBIT(e) ((oct_bits_t)((UADDR(e) & 2u) >> 1))
#define FLPBIT(e) ((oct_bits_t)(UADDR(e) & 1u))
  /*  
    In this implementation, the tumble code {t} of an arc {e} with base
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

      arc              sdf t | rot    sym    tor    fflip  vflip  dual   duar
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

oct_arc_t oct_rot_fflip(oct_arc_t e, int r, int f)
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
    uint64_t num;        /* Edge number (multipurpose). */
  } oct_edge_rec_t;
  /* An {oct_edge_rec_t} describes the connections of an unoriented,
    undirected edge {E} of the map, as well as of its dual edge. 
    If {b} is the base edge of {E}, then {next[i] = b.rot^i.onext}.
    That is,
    
      {next[0] = onext(b)}
      {next[1] = onext(rot(b))}
      {next[2] = onext(sym(b))}
      {next[3] = onext(tor(b))}
      
    */

#define NMASK ((((int64_t)1) << 61) - 1)
  /* Mask for valid edge numbers, {0..2^61-1}. */

uint64_t oct_edge_num(oct_edge_t E)
  { return E->num; }

void oct_set_edge_num(oct_edge_t E, uint64_t num)
  { demand((num & (~ NMASK)) == 0, "edge number is too big"); 
    E->num = num;
  }

uint64_t oct_arc_num(oct_arc_t *e)
  { return ((EDGE(e))->num << 3) | TUMBLECODE(e); }

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

oct_arc_t oct_walk(oct_arc_t e, int r, int n)
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
    oct_edge_t E = (oct_edge_t)notnull(malloc(ESIZE), "no mem");
    /* Get the edge's base arc: */
    oct_arc_t e = BASEARC(E);
    /* Set the fields for a map on {S^2} with {e} as the single (non-loop) edge {e}: */
    E->next[0] = NOP(e);
    E->next[1] = TOR(e);
    E->next[2] = SYM(e);
    E->next[3] = ROT(e);
    E->num = 0;
    return e;
  }

void oct_destroy_edge(oct_arc_t e)
  { /* Grab the edge record {E}: */
    oct_edge_t E = EDGE(e);
    /* Make sure that the edge is isolated: */
    oct_arc_t f = SYM(e);
    if (EDGE(E->next[0]) != E) { oct_splice(e, OPREV(e)); }
    if (EDGE(E->next[2]) != E) { oct_splice(f, OPREV(f)); }
    /* Paranoia check -- assumes that the oct-edge is well-formed: */
    int i;
    for (i = 0; i < 4; i++) { assert(EDGE(E->next[i]) == E); }
    /* OK, reclaim the edge record: */
    free(E);
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

void oct_write_arc(FILE *wr, oct_arc_t e, int width)
  { oct_edge_t E = EDGE(e);
    oct_bits_t et = TUMBLECODE(e);
    fprintf
      ( wr, ("%*" uint64_u_fmt ":%u%u%u"), width, 
        E->num, (unsigned int)((et & 4u) >> 2), (unsigned int)((et & 2u) >> 1), (unsigned int)(et & 1u)
      );
  }

oct_arc_t oct_read_arc(FILE *rd, oct_arc_vec_t *et)
  { /* Parse the edge number {num} and get the edge {E} from the edge table {ET}: */
    uint64_t num = fget_uint64(rd, 10);
    demand(num < et->ne, "invalid edge number");
    oct_edge_t E = oct_edge(et->e[num]);
    /* Parse ":" and the tumble code {t} in binary: */
    fget_match(rd, ":");
    uint64_t ut = fget_uint64(rd, 2);
    demand(ut <= 7, "invalid tumble code");
    /* Put it all together: */
    return oct_orient(E, (oct_bits_t)ut);
  }

/* MAP INPUT/OUTPUT */

#define FILE_TYPE "oct-edge"
#define FILE_VERSION "2007-01-16"

void oct_write_map(FILE *wr, oct_arc_vec_t *root, oct_arc_vec_t *et)
  { 
  
    int nE = 0;  /* Count of reachable octets (edge records). */
    
    auto bool_t renumber_edge(oct_arc_t p);
      /* Visit-proc that renumbers the base edge of {p} with {nE}
        and increments {nE}. */
    
    /* Write the header line: */
    filefmt_write_header(wr, FILE_TYPE, FILE_VERSION);
    /* Grab the number {nr} of roots and the root index width {wa}: */
    int nr = root->ne;
    /* Compute the width in digits {wr} of the root index: */
    int dr = (nr < 10 ? 1 : digits(nr-1));
    /* Make sure that {et} is non-null and points to an empty table: */
    oct_arc_vec_t etb = oct_arc_vec_new(0); /* A local edge table. */
    if (et == NULL) 
      { et = &etb; }
    else
      { oct_arc_vec_trim(et, 0); }
    /* Renumber all edge records reachable from {root[0..nr-1]}: */
    oct_enum_octets(*root, &renumber_edge, et);
    assert(et->ne == nE);
    /* Determine the width {eE} in digits of the edge number: */
    int dE = (nE < 10 ? 1 : digits(nE-1));
    /* We should have zero edges if and only if we have zero roots: */
    assert((nr == 0) == (nE == 0)); 
    /* Write the root and edge counts: */
    fprintf(wr, "roots = %d\n", nr);
    fprintf(wr, "edges = %d\n", nE);
    /* Write the roots, one per line: */
    unsigned int i;
    for (i = 0; i < nr; i++)
      { /* Write {i} and the root arc number {i}: */
        fprintf(wr, "%*u ", dr, i);
        oct_write_arc(wr, root->e[i], dE);
        fputc('\n', wr);
      }
    /* Write the edge records, one per line: */
    uint64_t num;
    for (num = 0; num < nE; num++)
      { /* Get the reachable edge number {num}: */
        oct_edge_t E = oct_edge(et->e[num]);
        /* Check whether the renumbering worked: */
        assert(E->num == num);
        /* Write the edge's number and its {onext} links: */
        fprintf(wr, ("%*" uint64_u_fmt), dE, num);
        int r;
        for (r = 0; r < 4; r++)
          { fputc(' ', wr); oct_write_arc(wr, E->next[r], dE); }
        fputc('\n', wr);
      }
    /* Write the footer line: */
    filefmt_write_footer(wr, FILE_TYPE);
    /* Reclaim the edge table if local: */
    if (et == &etb) { oct_arc_vec_trim(&etb, 0); }
    return;
    
    /* IMPLEMENTATIONS OF LOCAL PROCS */
    
    bool_t renumber_edge(oct_arc_t p)
      {
        oct_edge_t E = oct_edge(p);
        E->num = nE; nE++;
        return FALSE;
      }
  }

void oct_read_map(FILE *rd, oct_arc_vec_t *root, oct_arc_vec_t *et)
  { /* Check and consume the header line: */
    filefmt_read_header(rd, FILE_TYPE, FILE_VERSION);
    /* Parse the root count {nr} and the edge count {nE}: */
    int nr = nget_int(rd, "roots"); fget_eol(rd);
    int nE = nget_int(rd, "edges"); fget_eol(rd);
    /* Make sure that {et} is non-null and points to a table with {nE} slots: */
    oct_arc_vec_t etb = oct_arc_vec_new(0); /* A local edge table. */
    if (et == NULL) { et = &etb; }
    oct_arc_vec_trim(et, nE);
    /* Create the edge records, save their base arcs in {et->e[0..nE-1]}: */
    uint64_t num;
    for (num = 0; num < nE; num++) 
      { et->e[num] = oct_make_edge(); 
        oct_set_edge_num(oct_edge(et->e[num]), num);
      }
    /* (Re)alocate the root record and read the root arcs, one per line: */
    oct_arc_vec_trim(root, nr);
    uint64_t i;
    for (i = 0; i < nr; i++)
      { /* Parse the root index {i} and the root arc {root->e[i]}: */
        uint64_t iread = fget_uint64(rd, 10);
        demand(iread == i, "root index mismatch");
        /* Parse the root arc number {i}, save it in {root}: */
        root->e[i] = oct_read_arc(rd, et); 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    /* Read the contents of the edge records {0..nE-1}: */
    for (num = 0; num < nE; num++) 
      { /* Parse the edge number {num}: */
        uint64_t numread = fget_uint64(rd, 10);
        demand(numread == num, "edge number mismatch");
        /* Get the edge {E} from the edge table {et}: */
        oct_edge_t E = oct_edge(et->e[num]);
        /* Read its links {E->next[0..3]}: */
        int r;
        for (r = 0; r < 4; r++) { E->next[r] = oct_read_arc(rd, et); } 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    
    /* Check and consume the footer line: */
    filefmt_read_footer(rd, FILE_TYPE);
    /* Reclaim the edge table if local: */
    if (et == &etb) { oct_arc_vec_trim(&etb, 0); }
  }

vec_typeimpl(oct_arc_vec_t,oct_arc_vec,oct_arc_t);
vec_typeimpl(oct_edge_vec_t,oct_edge_vec,oct_edge_t);
