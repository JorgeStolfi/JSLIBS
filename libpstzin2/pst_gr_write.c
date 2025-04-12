/* See {pst_gr_write.h} */
/* Last edited on 2025-03-14 19:22:23 by stolfi */
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

#include <pst_gr.h>
#include <pst_gr_path.h>

#include <pst_gr_write.h>

void pst_gr_write_vertex(FILE *wr, pst_gr_t *gr, uint32_t vi);
  /* Writes the vertex {gr->vertex[vi]} to {wr} in a format compatible
    with {pst_gr_read_vertex}. */

void pst_gr_write_edge(FILE* wr,pst_gr_t* gr, pst_gr_arc_t e);

/* IMPLEMENTATIONS */

void pst_gr_write_named(char *fname, pst_gr_t* gr, bool_t verbose)
  { FILE *wr = open_write(fname, verbose);
    pst_gr_write_file(wr, gr, verbose);
    fclose(wr);
  }

void pst_gr_write_file(FILE *wr, pst_gr_t *gr, bool_t verbose)
  {
    filefmt_write_header(wr, pst_gr_FILE_TYPE, pst_gr_FILE_VERSION);
    fprintf(wr, "NV = %d\n", gr->NV);
    fprintf(wr, "NE = %d\n", gr->NE);
    fprintf(wr, "NX = %d\n", gr->NX);
    fprintf(wr, "NY = %d\n", gr->NY);
    fprintf(wr, "\n");
    fprintf(wr, "vertices\n");
    for (pst_gr_vertex_t vi = 0; vi < gr->NV; vi++)
      { pst_gr_write_vertex(wr, gr, vi); }
    fprintf(wr, "\n");
    fprintf(wr,"edges\n");
    for (pst_gr_edge_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_write_edge(wr, gr, ei); }
    filefmt_write_footer(wr, pst_gr_FILE_TYPE);
  }

void pst_gr_write_vertex(FILE *wr, pst_gr_t *gr, pst_gr_vertex_t vi)
  { demand(vi <= gr->NV, "invalid vertex index {vi}");
    fprintf(wr, "%d", vi); 
    pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);
    fprintf(wr, " [%d,%d] at (%f,%f) out ", vd->x, vd->y, vd->coords.c[0], vd->coords.c[1]);
    pst_gr_arc_print(wr, gr, vd->aout);
    fprintf(wr, "\n");
  }

void pst_gr_write_edge(FILE *wr, pst_gr_t *gr, pst_gr_edge_t ei)
  { demand(ei < gr->NE, "invalid edge index {ei}");
    fprintf(wr, "%d", ei); 
    pst_gr_edge_data_t *ed = &(gr->edata[ei]);
    
    for (int32_t db = 0; db <= 1; db++)
      { fprintf(wr, "  ");
        fprintf(wr, " %d", ed->org[db]);
        fprintf(wr, " ");
        pst_gr_arc_print(wr, gr, ed->eprev[db]);
        fprintf(wr, " ");
        pst_gr_arc_print(wr, gr, ed->enext[db]);
      }
        
    fprintf(wr," %+.8f %.6f", ed->delta, ed->weight);

    fprintf(wr, " path ");
    pst_gr_path_write(wr, ed->path);
    fprintf(wr,"\n");
  }
