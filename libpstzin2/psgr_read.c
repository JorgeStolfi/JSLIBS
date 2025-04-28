/* See {psgr_read.h} */
/* Last edited on 2025-04-27 03:01:30 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsfile.h>
#include <filefmt.h>
#include <affirm.h>
#include <bool.h>
#include <nget.h>
#include <fget.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_unsafe.h>
#include <psgr_path.h>

#include <psgr_read.h>

#define UNMARKED psgr_UNMARKED

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE
  /* Shorter names. */  

psgr_arc_t psgr_read_arc(FILE *rd, uint32_t NE);
  /* Reads from {rd} an arc spec, either "-1" (meaning {ANONE})
    or "{ei}:{db}" where {ei} is in {0..NE-1} and {db} is 
    either 0 or 1 (meaning the arc with edge index {ei}
    and direction bit {db}). */

/* IMPLEMENTATIONS */

psgr_t* psgr_read_named(char *fname, bool_t verbose)
  { FILE *rd = open_read(fname, verbose);
    psgr_t *gr = psgr_read_file(rd, verbose);
    fclose(rd);
    return gr;
  }

psgr_t* psgr_read_file(FILE *rd, bool_t verbose)
  {
    bool_t debug = FALSE;
    
    filefmt_read_header(rd, psgr_FILE_TYPE, psgr_FILE_VERSION);
    
    uint32_t NV = nget_uint32(rd, "NV", 10); fget_eol(rd);
    uint32_t NE = nget_uint32(rd, "NE", 10); fget_eol(rd);
    demand(NV <= psgr_MAX_VERTICES, "too many vertices");
    demand(NE <= psgr_MAX_EDGES, "too many edges");
    
    uint32_t NX = nget_uint32(rd, "NX", 10); fget_eol(rd);
    uint32_t NY = nget_uint32(rd, "NY", 10); fget_eol(rd);
    demand(NX <= psgr_MAX_IMG_SIZE, "too many image cols");
    demand(NY <= psgr_MAX_IMG_SIZE, "too many image rows");
    demand((NX == 0) == (NY == 0), "inconsistent {NX.NY}");
    
    psgr_t *gr = psgr_new(NV, NE, NX, NY); 
    
    auto void read_vertex(void);
    auto void read_edge(void);
    
    psgr_arc_t *aouts = talloc(NV, psgr_arc_t); /* The out arc of each vertex. */

    if (debug) { fprintf(stderr, "reading the %d vertices ...\n", NV); }
    fget_skip_formatting_chars(rd); 
    fget_match(rd, "vertices"); fget_eol(rd);
    fget_skip_formatting_chars(rd); 
    for (uint32_t kv = 0; kv < NV ; kv++) 
      { if (debug) { fprintf(stderr, "  reading vertex %d ...\n", kv); }
        assert(psgr_num_vertices(gr) == kv);
        read_vertex();
      }
    assert(psgr_num_vertices(gr) == NV);
      
    if (debug) { fprintf(stderr, "reading the %d edges ...\n", NE); }
    fget_skip_formatting_chars(rd); 
    fget_match(rd, "edges"); fget_eol(rd);
    fget_skip_formatting_chars(rd); 
    for (uint32_t ke = 0; ke < NE ; ke++)
      { if (debug) { fprintf(stderr, "  reading edge %d ...\n", ke); }
        assert(psgr_num_edges(gr) == ke);
        read_edge();
      }
    assert(psgr_num_edges(gr) == NE);

    /* Set the {.aout}, even though the edge will be created later: */
    for (uint32_t kv = 0; kv < NV ; kv++) 
      { psgr_vertex_t v = (psgr_vertex_t){ kv };
        psgr_unsafe_vertex_set_out_arc(gr, v, aouts[kv]);
      }

    filefmt_read_footer(rd, psgr_FILE_TYPE);
    
    free(aouts);
    return gr;

    void read_vertex(void)
      { uint32_t id = fget_uint32(rd, 10);
        demand(id == psgr_num_vertices(gr), "vertex index error");
        fget_skip_spaces_and_test_char(rd, '[');
        int32_t x = fget_int32(rd);
        fget_skip_spaces_and_test_char(rd, ',');
        int32_t y = fget_int32(rd);
        fget_skip_spaces_and_test_char(rd, ']');

        fget_skip_spaces_and_match(rd, "out");
        psgr_arc_t aout = psgr_read_arc(rd, NE);
        demand((aout.ix == ANONE.ix) || (aout.ix < 2*NE), "invalid outgoing arc index");

        /* Create the disconnected vertex: */
        i2_t ip = (i2_t){{ x, y }};
        psgr_vertex_t v = psgr_add_vertex(gr, &ip);
        assert(v.ix == psgr_num_vertices(gr)-1);
        
        /* Save the outgoing arc for later, since the edge data does not exist yet: */ 
        aouts[v.ix] = aout;
        
        fget_eol(rd);
      }

    void read_edge(void)
      { uint32_t id = fget_uint32(rd, 10); 
        demand(id == psgr_num_edges(gr), "edge index error");
        psgr_vertex_t org[2];
        psgr_arc_t eprev[2];
        psgr_arc_t enext[2];
        
        for (int32_t db = 0; db <= 1; db++)
          { org[db].ix = fget_uint32(rd, 10);
            demand(org[db].ix < NV, "invalid origin vertex index");
            eprev[db] = psgr_read_arc(rd, NE);
            enext[db] = psgr_read_arc(rd, NE);
          }

        double delta = fget_double(rd);
        demand(isfinite(delta), "invalid edge delta");

        double weight = fget_double(rd);
        demand(isfinite(weight) && (weight >= 0), "invalid edge weight");

        fget_skip_spaces_and_match(rd, "path");
        psgr_path_t P = psgr_path_read(rd, NX, NY);
        
        psgr_arc_t a = psgr_unsafe_add_unlinked_edge(gr, org[0], org[1], P);
        psgr_arc_set_delta(gr, a, delta);
        psgr_arc_set_weight(gr, a, weight);
        
        psgr_edge_t e = psgr_arc_edge(a);
        assert(e.ix == psgr_num_edges(gr)-1);
        
        /* Set the edge links as read from the file */
        for (uint32_t db = 0; db <= 1; db++)
          { psgr_arc_t b = psgr_orient_edge(e, (psgr_dir_bit_t)db);
            psgr_unsafe_arc_set_onext(gr, b, enext[db]);
            psgr_unsafe_arc_set_oprev(gr, b, eprev[db]);
          }

        fget_eol(rd);
      }

  }
    
psgr_arc_t psgr_read_arc(FILE *rd, uint32_t NE)
  { int32_t e_ix = fget_int32(rd); /* Must be signed int. */
    if (e_ix < 0) { return ANONE; }
    fget_match(rd, ":");
    uint32_t db = fget_uint32(rd, 10);
    demand((db == 0) || (db == 1), "invalid direction bit");
    return psgr_orient_edge((psgr_edge_t){ (uint32_t)e_ix }, (psgr_dir_bit_t)db);
  }
