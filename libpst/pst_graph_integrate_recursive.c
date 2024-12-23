/* See {pst_graph_integrate_recursive.h}. */
/* Last edited on 2024-12-22 21:42:17 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <float_image.h>
#include <rn.h>

#include <pst_graph.h>
#include <pst_graph_integrate.h>

#include <pst_graph_integrate_recursive.h>
 
void pst_graph_estimate_from_shrunk(pst_graph_t* g, pst_graph_t* jg, float_image_t* OZ);

void pst_graph_integration_recursive
  ( pst_graph_t *g,
    float_image_t *OZ,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose,
    uint32_t level
  )
  { char *filename_grph = jsprintf("graph-%02d.txt", level);
    FILE *wr_grph = open_write(filename_grph, FALSE);
    pst_graph_write(wr_grph, g);
    fclose(wr_grph);

    if (g->num_vertex > 2)
      { fprintf(stderr, "reducing Graph level[%d] with [%d] vertices\n", level, g->num_vertex);
        pst_graph_t *jg = pst_graph_shrink(g);
        uint32_t maxIterSub = (uint32_t)ceil(((double)maxIter)*M_SQRT2);
        pst_graph_integration_recursive(jg, OZ, maxIterSub, convTol/M_SQRT2, para, szero, verbose, level+1);
        pst_graph_estimate_from_shrunk(g, jg, OZ);

        char *filename = jsprintf("guess-%02d.fni", level);
        FILE *wr_dump = open_write(filename, FALSE);
        float_image_write(wr_dump, OZ);
        fclose(wr_dump);
        pst_graph_free(jg);
      }
    else
      { fprintf(stderr, "end of Recursion - returning\n");
        float_image_fill_channel(OZ, 0, 0);
      }
    fprintf(stderr, "solving Graph level[%d] with [%d] vertices\n", level, g->num_vertex);
    pst_imgsys_t *S = pst_graph_build_integration_system(g, (uint32_t)OZ->sz[1], (uint32_t)OZ->sz[2]);

    char *filename = jsprintf("system-%02d.txt", level);
    FILE *wr_dump = open_write(filename,  FALSE);
    free(filename);
    pst_imgsys_write(wr_dump, S);
    fclose(wr_dump);

    filename = jsprintf("height-%02d.fni", level);
    wr_dump = open_write(filename, FALSE);
    free(filename);
    float_image_write(wr_dump, OZ);
    fclose(wr_dump);

    pst_graph_solve_system(g, S, OZ, maxIter, convTol, para, szero, verbose);
    pst_imgsys_free(S);
  }

void pst_graph_estimate_from_shrunk(pst_graph_t *g, pst_graph_t *jg, float_image_t *OZ)
  {
    uint32_t ind_vt = 0;
    uint32_t NX_Z = (uint32_t)OZ->sz[1];
    uint32_t NY_Z = (uint32_t)OZ->sz[2];
    for (uint32_t i = 0; i < g->num_vertex; i++)
    {
      pst_vertex_t *v = &(g->vertices[i]);
      if (((v->x + v->y)%2) != 0) { continue; }
      while (jg->vertices[ind_vt].id != v->id)
        { ind_vt++;
          assert(ind_vt < jg->num_vertex);
        }
      int32_t x, y;
      pst_graph_restore_vertex_index((int32_t)v->id, NX_Z, NY_Z, &x, &y);
      assert((x != -1) && (y != -1));
      double z = float_image_get_sample(OZ, 0, (int32_t)x, (int32_t)y);
      uint32_t has_xp = (uint32_t)(v->vtxp != -1); 
      uint32_t has_xm = (uint32_t)(v->vtxm != -1);
      uint32_t has_yp = (uint32_t)(v->vtyp != -1);
      uint32_t has_ym = (uint32_t)(v->vtym != -1);
      if (has_xp)
          { int32_t vxp;
            double dxp, wxp;
            pst_graph_vertex_get_neighbour(g, i, 1, 0, &vxp, &dxp, &wxp);
            assert(vxp != -1);
            int32_t ix, iy;
            pst_graph_restore_vertex_index((int32_t)g->vertices[vxp].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int32_t)ix, (int32_t)iy, (float)(z - dxp));
          }
        if (has_yp)
          { int32_t vyp;
            double dyp, wyp;
            pst_graph_vertex_get_neighbour(g, i, 0, 1, &vyp, &dyp, &wyp);
            assert(vyp != -1);
            int32_t ix, iy;
            pst_graph_restore_vertex_index((int32_t)g->vertices[vyp].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int32_t)ix, (int32_t)iy, (float)(z - dyp));
          }
        if (has_xm)
          { int32_t vxm;
            double dxm, wxm;
            pst_graph_vertex_get_neighbour(g, i, -1, 0, &vxm, &dxm, &wxm);
            assert(vxm != -1);
            int32_t ix, iy;
            pst_graph_restore_vertex_index((int32_t)g->vertices[vxm].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int32_t)ix, (int32_t)iy, (float)(z - dxm));
          }
        if (has_ym)
          { int32_t vym;
            double dym, wym;
            pst_graph_vertex_get_neighbour(g, i, 0, -1, &vym, &dym, &wym);
            assert(vym != -1);
            int32_t ix, iy;
            pst_graph_restore_vertex_index((int32_t)g->vertices[vym].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int32_t)ix, (int32_t)iy, (float)(z - dym));
          }
      }
  }
