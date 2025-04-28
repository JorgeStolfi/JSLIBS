/* See {psgr_write.h} */
/* Last edited on 2025-04-27 02:58:49 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsfile.h>
#include <filefmt.h>
#include <affirm.h>
#include <bool.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_path.h>

#include <psgr_write.h>

void psgr_write_vertex(FILE *wr, psgr_t *gr, psgr_vertex_t v);
  /* Writes the vertex {gr->vertex[v]} to {wr} in a format compatible
    with {psgr_read_vertex}. */

void psgr_write_edge(FILE* wr,psgr_t* gr, psgr_edge_t e);

/* IMPLEMENTATIONS */

void psgr_write_named(char *fname, psgr_t* gr, bool_t verbose)
  { FILE *wr = open_write(fname, verbose);
    psgr_write_file(wr, gr, verbose);
    fclose(wr);
  }

void psgr_write_file(FILE *wr, psgr_t *gr, bool_t verbose)
  {
    uint32_t NV = psgr_num_vertices(gr);
    uint32_t NE = psgr_num_edges(gr);
    i2_t sz = psgr_image_size(gr);
    filefmt_write_header(wr, psgr_FILE_TYPE, psgr_FILE_VERSION);
    fprintf(wr, "NV = %d\n", NV);
    fprintf(wr, "NE = %d\n", NE);
    fprintf(wr, "NX = %d\n", sz.c[0]);
    fprintf(wr, "NY = %d\n", sz.c[1]);
    fprintf(wr, "\n");
    fprintf(wr, "vertices\n");
    for (uint32_t kv = 0; kv < NV; kv++)
      { psgr_vertex_t v = (psgr_vertex_t){ kv };
        psgr_write_vertex(wr, gr, v);
      }
    fprintf(wr, "\n");
    fprintf(wr,"edges\n");
    for (uint32_t ke = 0; ke < NE; ke++)
      { psgr_edge_t e = (psgr_edge_t){ ke };
        psgr_write_edge(wr, gr, e); 
      }
    filefmt_write_footer(wr, psgr_FILE_TYPE);
  }

void psgr_write_vertex(FILE *wr, psgr_t *gr, psgr_vertex_t v)
  { uint32_t NV = psgr_num_vertices(gr);
    demand(v.ix <= NV, "invalid vertex index {v}");
    fprintf(wr, "%d", v.ix); 
    i2_t ip = psgr_vertex_pixel(gr, v);
    fprintf(wr, " [%d,%d]", ip.c[0], ip.c[1]);
    psgr_arc_t aout = psgr_vertex_out_arc(gr, v);
    fprintf(wr, " out "); psgr_arc_print(wr, gr, aout);
    fprintf(wr, "\n");
  }

void psgr_write_edge(FILE *wr, psgr_t *gr, psgr_edge_t e)
  { uint32_t NE = psgr_num_edges(gr);
    demand(e.ix < NE, "invalid edge {e}");
    fprintf(wr, "%d", e.ix); 
    
    for (int32_t db = 0; db <= 1; db++)
      { psgr_arc_t a = psgr_orient_edge(e, (psgr_dir_bit_t)db);
        fprintf(wr, "  ");
        fprintf(wr, " %d", psgr_arc_org(gr, a).ix);
        fprintf(wr, " ");
        psgr_arc_print(wr, gr, psgr_arc_oprev(gr, a));
        fprintf(wr, " ");
        psgr_arc_print(wr, gr, psgr_arc_onext(gr, a));
      }
        
    psgr_arc_t a0 = psgr_orient_edge(e, 0);
    fprintf(wr," %+.8f", psgr_arc_delta(gr, a0));
    fprintf(wr," %.6f", psgr_arc_weight(gr, a0));

    fprintf(wr, " path ");
    psgr_path_write(wr, psgr_arc_path(gr, a0));
    fprintf(wr,"\n");
  }
