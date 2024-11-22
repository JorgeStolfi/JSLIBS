/* See {tosl_mesh_make_keg.h} */
/* Last edited on 2024-11-20 05:21:37 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_mesh_make_keg.h>

tosl_mesh_t *tosl_mesh_make_keg
  ( int32_t NS,
    int32_t NR,
    int32_t NB,
    tosl_coord_t Zmax,
    int32_t debug
  )
  {
    int32_t NV = NS*(NR+1) + NS*NR*NB;
    int32_t NE = NS*(NR+1) + NS*NR*(NB + 1);
    
    tosl_mesh_t *mesh = tosl_mesh_new(NE, NV);
   
    /* Half-height of barrel {ZH} is even and a bit bigger than {Zmax}: */
    tosl_coord_t ZH = Zmax + 3;
    ZH = 2*(ZH/2);
    
    double lat_max = M_PI/6; /* Max latitude of barrel wall. */
    double Rw = ((double)ZH)/sin(lat_max); /* Approx radius of curvature of meridians. */

    auto tosl_point_t main_vert_pos(int32_t s, int32_t r);
      /* Returns the the coords of the main (not subdivision)
        vertex on meridian {s} and parallel {r}. */
    
    fprintf(stderr, "  creating %d rings of %d vertices and edges...\n", NR+1, NS);
    tosl_arc_id_t iaring[NR+1];
    for (int32_t r = 0; r <= NR; r++)
      { 
        char *pref = NULL; 
        char *pref = jsprintf("h.r%d.s", r);
        iaring[r] = tosl_mesh_add_ring(NS, pref, mesh);
        /* Set the vertex coordinates: */
        tosl_arc_id_t ia = iaring[r];
        for (int32_t s = 0; s < NS; s++)
          { tosl_point_t vsr = main_vert_pos(s, r);
            mesh->Vpos[mesh->Arc[ia].ivorg] = vsr;
            ia = mesh->Arc[ia].skip;
          }
        free(pref);
      }
    
    auto tosl_point_t interp_vert_pos(tosl_vert_id_t kv0, tosl_vert_id_t kv1, int32_t b);
      /* Returns the coordinates of the vertex with local index {b} on the line between main
        vertices {kv0} and {kv1}.  */
    
    fprintf(stderr, "  Connecting the rings with paths of %d edges...\n", NB+1);
    for (int32_t r = 0; r < NR; r++)
      { tosl_arc_id_t ia0 = iaring[r];
        tosl_arc_id_t ia1 = iaring[r+1];
        
        tosl_arc_id_t ja0 = mesh->Arc[ia0].skip;
        tosl_arc_id_t ja1 = mesh->Arc[ia1].skip;
        for (int32_t s = 0; s < NS; s++)
          { tosl_arc_id_t ka1 = tosl_sym(ja1);
            char *pref = NULL; 
            char *pref = jsprintf("v.r%d.s%d.b", r, s);
            tosl_mesh_add_path(NB, ia0, ka1, pref, mesh);
            
            /* Set the interpolated vertex coordinates: */
            tosl_vert_id_t kv0 = mesh->Arc[ja0].ivorg;
            tosl_vert_id_t kv1 = mesh->Arc[ja1].ivorg;
            tosl_arc_id_t iab = mesh->Arc[ia0].skip;
            for (int32_t b = 0; b < NB; b++)
              { tosl_point_t vsrb = interp_vert_pos(kv0, kv1, b);
                tosl_arc_id_t kab = tosl_sym(iab);
                tosl_vert_id_t kv = mesh->Arc[kab].ivorg;
                mesh->Vpos[kv] = vsrb;
                if (debug) 
                  { fprintf(stderr, "    r = %d  s = %d  b = %d", r, s, b);
                    tosl_mesh_arc_print(stderr, "  kab = ", kab, "\n", mesh);
                  }
                iab = mesh->Arc[iab].skip;
              }
            free(pref);
            ia0 = ja0; ja0 = mesh->Arc[ja0].skip;
            ia1 = ja1; ja1 = mesh->Arc[ja1].skip;
          }
     }
    return mesh;                   
    
    tosl_point_t main_vert_pos(int32_t s, int32_t r)
      {
        double fr = ((double)r)/((double)NR);
        double lat = 2*(fr - 0.5)*lat_max;
        tosl_coord_t Zv = (tosl_coord_t)floor(Rw*sin(lat) + 0.5);
        Zv = 2*(Zv/2);
        
        assert(cos(lat) > 0.65); 
        double Rr = Rw*(cos(lat) - 0.5);
        double fs = ((double)s)/((double)NS);
        double lon = fs*2*M_PI;
        tosl_coord_t Xv = (tosl_coord_t)floor(Rr*cos(lon) + 0.5);
        Xv = 2*(Xv/2);
        tosl_coord_t Yv = (tosl_coord_t)floor(Rr*sin(lon) + 0.5);
        Yv = 2*(Yv/2);
        
        return (tosl_point_t){{ Xv, Yv, Zv }};
      }
    
    tosl_point_t interp_vert_pos(tosl_vert_id_t kv0, tosl_vert_id_t kv1, int32_t b)
      {
        /* Get the coordinates of the two vertices: */
        assert((kv0 >= 0) && (kv0 < mesh->NV));
        tosl_point_t *v0 = &(mesh->Vpos[kv0]);
        assert((kv1 >= 0) && (kv1 < mesh->NV));
        tosl_point_t *v1 = &(mesh->Vpos[kv1]);
        
        double f = ((double)b+1)/((double)NB+1);
        
        tosl_coord_t Xv = (tosl_coord_t)floor((1-f)*v0->c[0] + f*v1->c[0] + 0.5);
        Xv = 2*(Xv/2);
        tosl_coord_t Yv = (tosl_coord_t)floor((1-f)*v0->c[1] + f*v1->c[1] + 0.5);
        Yv = 2*(Yv/2);
        tosl_coord_t Zv = (tosl_coord_t)floor((1-f)*v0->c[2] + f*v1->c[2] + 0.5);
        Zv = 2*(Zv/2);
        
        return (tosl_point_t){{ Xv, Yv, Zv }};
      }
  }
