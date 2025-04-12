/* See {pst_gr_read.h} */
/* Last edited on 2025-03-14 21:55:41 by stolfi */
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

#include <pst_gr.h>
#include <pst_gr_path.h>

#include <pst_gr_read.h>

#define UNMARKED pst_gr_UNMARKED

#define NONE pst_gr_NONE

#define pst_gr_MAX_VERTICES (4096*4096)
#define pst_gr_MAX_EDGES (3*pst_gr_MAX_VERTICES + 2)

pst_gr_arc_t pst_gr_read_arc(FILE *rd, uint32_t NE);
  /* Reads from {rd} an arc spec, either "-1" (meaning {pst_gr_NONE})
    or "{ei}:{db}" where {ei} is in {0..NE-1} and {db} is 
    either 0 or 1 (meaning the arc with edge index {ei}
    and direction bit {db}). */

/* IMPLEMENTATIONS */

pst_gr_t* pst_gr_read_named(char *fname, bool_t verbose)
  { FILE *rd = open_read(fname, verbose);
    pst_gr_t *gr = pst_gr_read_file(rd, verbose);
    fclose(rd);
    return gr;
  }

pst_gr_t* pst_gr_read_file(FILE *rd, bool_t verbose)
  {
    bool_t debug = FALSE;
    
    filefmt_read_header(rd, pst_gr_FILE_TYPE, pst_gr_FILE_VERSION);
    
    uint32_t NV = nget_uint32(rd, "NV", 10); fget_eol(rd);
    uint32_t NE = nget_uint32(rd, "NE", 10); fget_eol(rd);
    demand(NV <= pst_gr_MAX_VERTICES, "too many vertices");
    demand(NE <= pst_gr_MAX_EDGES, "too many edges");
    
    uint32_t NX = nget_uint32(rd, "NX", 10); fget_eol(rd);
    uint32_t NY = nget_uint32(rd, "NY", 10); fget_eol(rd);
    demand(NX <= pst_gr_MAX_IMG_SIZE, "too many image cols");
    demand(NY <= pst_gr_MAX_IMG_SIZE, "too many image rows");
    demand((NX == 0) == (NY == 0), "inconsistent {NX.NY}");
    
    pst_gr_t *gr = pst_gr_new(NV, NE, NX, NY); 
    assert((gr->NV == 0) && (gr->NV_max == NV));
    assert((gr->NE == 0) && (gr->NE_max == NE));
    assert((gr->NX == NX));
    assert((gr->NY == NY));


    auto void read_vertex(void);
    auto void read_edge(void);

    if (debug) { fprintf(stderr, "reading the %d vertices ...\n", NV); }
    fget_skip_formatting_chars(rd); 
    fget_match(rd, "vertices"); fget_eol(rd);
    fget_skip_formatting_chars(rd); 
    for (pst_gr_vertex_t vi = 0; vi < NV ; vi++) 
      { if (debug) { fprintf(stderr, "  reading vertex %d ...\n", vi); }
        assert(gr->NV == vi);
        read_vertex();
      }
    assert(gr->NV == NV);
      
    if (debug) { fprintf(stderr, "reading the %d edges ...\n", NE); }
    fget_skip_formatting_chars(rd); 
    fget_match(rd, "edges"); fget_eol(rd);
    fget_skip_formatting_chars(rd); 
    for (pst_gr_edge_t ei = 0; ei < NE ; ei++)
      { if (debug) { fprintf(stderr, "  reading edge %d ...\n", ei); }
        assert(gr->NE == ei);
        read_edge();
      }
    assert(gr->NE == NE);

    filefmt_read_footer(rd, pst_gr_FILE_TYPE);
    return gr;
    

  void read_vertex(void)
    { uint32_t id = fget_uint32(rd, 10); demand(id == gr->NV, "vertex index error");
      fget_skip_spaces_and_test_char(rd, '[');
      int32_t x = fget_int32(rd);
      fget_skip_spaces_and_test_char(rd, ',');
      int32_t y = fget_int32(rd);
      fget_skip_spaces_and_test_char(rd, ']');

      fget_skip_spaces_and_match(rd, "at");
      r2_t p;
      fget_skip_spaces_and_test_char(rd, '(');
      p.c[0] = fget_double(rd); 
      fget_skip_spaces_and_test_char(rd, ',');
      p.c[1] = fget_double(rd); 
      fget_skip_spaces_and_test_char(rd, ')');
      demand(isfinite(p.c[0]) && isfinite(p.c[1]), "invalid vertex coordinates");

      fget_skip_spaces_and_match(rd, "out");
      pst_gr_arc_t aout = pst_gr_read_arc(rd, NE);
      demand((aout == NONE) || (aout < 2*NE), "invalid outgoing arc index");
      
      /* Create the disconnected vertex: */
      pst_gr_vertex_t vi = pst_gr_add_vertex(gr, x, y, p);
      assert(vi == gr->NV-1);
 
      /* Set the {aout} link (hopefully valid once edges are added): */
      gr->vdata[vi].aout = aout;

      fget_eol(rd);
    }

    void read_edge(void)
      { uint32_t id = fget_uint32(rd, 10); demand(id == gr->NE, "edge index error");
        pst_gr_vertex_t org[2];
        pst_gr_arc_t eprev[2];
        pst_gr_arc_t enext[2];
        
        for (int32_t db = 0; db <= 1; db++)
          { org[db] = fget_uint32(rd, 10);
            demand(org[db] < NV, "invalid origin vertex index");
            eprev[db] = pst_gr_read_arc(rd, NE);
            enext[db] = pst_gr_read_arc(rd, NE);
          }

        double delta = fget_double(rd);
        demand(isfinite(delta), "invalid edge delta");

        double weight = fget_double(rd);
        demand(isfinite(weight) && (weight >= 0), "invalid edge weight");

        fget_skip_spaces_and_match(rd, "path");
        pst_gr_path_t P = pst_gr_path_read(rd);
        
        pst_gr_arc_t ai = pst_gr_add_edge(gr, org[0], org[1], delta, weight, P, FALSE);
        pst_gr_edge_t ei = pst_gr_arc_edge(ai);
        assert(ei == gr->NE-1);
        
        /* Set the edge links as read from the file */
        pst_gr_edge_data_t *ed = &(gr->edata[ei]);
        for (int32_t db = 0; db <= 1; db++)
          { ed->org[db] = org[db];
            ed->eprev[db] = eprev[db];
            ed->enext[db] = enext[db];
          }

        fget_eol(rd);
      }

  }
    
pst_gr_arc_t pst_gr_read_arc(FILE *rd, uint32_t NE)
  { int32_t ei = fget_int32(rd);
    if (ei < 0) { return NONE; }
    fget_match(rd, ":");
    uint32_t db = fget_uint32(rd, 10);
    demand((db == 0) || (db == 1), "invalid direction bit");
    return pst_gr_orient_edge((pst_gr_edge_t)ei, (pst_gr_dir_bit_t)db);
  }
