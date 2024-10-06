/* See {tosl_mesh_obj_write.h} */
/* Last edited on 2024-10-06 16:51:41 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

#include <haf.h>

#include <tosl.h>
#include <tosl_mesh_obj_write.h>

void tosl_mesh_obj_write(FILE *wr, tosl_mesh_t *mesh)
  {
    /* Find coordinate ranges: */
    tosl_point_t vmin, vmax;
    tosl_mesh_coord_range_get(mesh, &vmin, &vmax);
    
    fprintf(wr, "# Created by %s\n", __FUNCTION__);
    fprintf(wr, "# Vertices\n");
    for (tosl_vert_id_t iv = 0; iv < mesh->NV; iv++)
      { tosl_point_t *v = &(mesh->Vpos[iv]);
        fprintf(wr, "v");
        for (int32_t j = 0; j < 3; j++) 
          { double cj = 2*((double)v->c[j] - vmin.c[j])/((double)vmax.c[j] - vmin.c[j] + 1) - 1.0;
            fprintf(wr, " %+8.6f", cj);
          }
        if ((mesh->Vlab != NULL) && (mesh->Vlab[iv] != NULL))
          { fprintf(wr, " # %s", mesh->Vlab[iv]); }
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
    
    fprintf(wr, "# Faces\n");
    int32_t NA = 2*mesh->NE;
    uint8_t *seen = malloc(NA*sizeof(uint8_t));
    for (tosl_arc_id_t ia = 0; ia < NA; ia++) { seen[ia] = 0; }
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { if (seen[ia] == 0)
          { fprintf(wr, "f");
            tosl_arc_id_t ka = ia;
            do 
              { assert(seen[ka] == 0); 
                seen[ka] = 1; 
                tosl_vert_id_t kv = mesh->Arc[ka].ivorg;
                fprintf(wr, " %d", kv+1);
                ka = mesh->Arc[ka].skip;
              } while (ka != ia);
            fprintf(wr, "\n");
          }
      }
    fflush(wr);
    free(seen);
  }
