/*  Last edited on 2024-10-06 16:59:45 by stolfi */
/* Test of {tosl_mesh.h} */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <tosl.h>
#include <tosl_mesh.h>

int32_t main(int32_t argc, char **argv);

void do_test_ring(int32_t NS, int8_t debug);
  /* Tests the procedure {tosl_mesh_add_ring(NS,...)} with different
    values of {NS}. */

void do_test_path(int32_t NB, int8_t debug);
  /* Tests the procedure {tosl_mesh_add_path(NB,...)}
    with different values of {NB}. */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615);
    
    do_test_ring(7, 1);
    do_test_ring(100, 0);
    
    do_test_path(0, 1);
    do_test_path(4, 1);
    do_test_path(100, 0);
    
    return 0;
  }
  
void do_test_ring(int32_t NS, int8_t debug)
  {
    fprintf(stderr, "start test of {tosl_mesh_add_ring} with NS = %d...\n", NS);
        
    /* Radius of ring: */
    tosl_coord_t R = 20*NS;

    fprintf(stderr, "creating mesh...\n");
    int32_t NE = NS; /* Expected number of edges. */
    int32_t NV = NS; /* Expected number of vertices. */
    tosl_mesh_t *mesh = tosl_mesh_new(NE, NV);

    fprintf(stderr, "adding the ring...\n");
    tosl_arc_id_t ia = tosl_mesh_add_ring(NS, "v.", mesh);
    assert(ia == 0);
    assert(mesh->NE == NE);
    assert(mesh->NV == NV);
    
    fprintf(stderr, "defining vertex coordinates...\n");
    double tilt = M_PI/6;
    for (tosl_vert_id_t s = 0; s < NV; s++)
      { double ang = ((double)s)/((double)NV)*2*M_PI;
        tosl_coord_t Xv = (tosl_coord_t)floor(R*cos(ang)*cos(tilt) + 0.5); Xv = 2*(Xv/2);
        tosl_coord_t Yv = (tosl_coord_t)floor(R*sin(ang)*cos(tilt) + 0.5); Yv = 2*(Yv/2);
        tosl_coord_t Zv = (tosl_coord_t)floor(R*sin(ang)*sin(tilt) + 0.5); Zv = 2*(Zv/2);
        mesh->Vpos[s] = (tosl_point_t){{ Xv, Yv, Zv }};
      }
    
    fprintf(stderr, "printing mesh...\n");
    if (debug) { tosl_mesh_print(stderr, mesh); } 
     
    fprintf(stderr, "checking topological consistency of mesh...\n");
    tosl_mesh_check(mesh);
  }
  
void do_test_path(int32_t NB, int8_t debug)
  {
    fprintf(stderr, "start test of {tosl_mesh_add_path} with NB = %d...\n", NB);
        
    /* Size and radius of ring: */
    int32_t NS = 2*NB + 5;
    tosl_coord_t R = 20*NS;

    fprintf(stderr, "creating mesh...\n");
    int32_t NE = NS + NB + 1; /* Expected total number of edges. */
    int32_t NV = NS + NB;     /* Expected total number of vertices. */
    tosl_mesh_t *mesh = tosl_mesh_new(NE, NV);

    fprintf(stderr, "adding a ring with %d edges...\n", NS);
    tosl_arc_id_t ia_ring = tosl_mesh_add_ring(NS, "r.", mesh);
    assert(ia_ring == 0);
    assert(mesh->NE == NS);
    assert(mesh->NV == NS);
    
    fprintf(stderr, "defining ring vertex coordinates...\n");
    for (tosl_vert_id_t s = 0; s < NS; s++)
      { double ang = ((double)s)/((double)NS)*2*M_PI;
        tosl_coord_t Xv = (tosl_coord_t)floor(R*cos(ang)+ 0.5); Xv = 2*(Xv/2);
        tosl_coord_t Yv = (tosl_coord_t)floor(R*sin(ang)+ 0.5); Yv = 2*(Yv/2);
        tosl_coord_t Zv = 0;
        mesh->Vpos[s] = (tosl_point_t){{ Xv, Yv, Zv }};
      }
    
    fprintf(stderr, "choosing places to attach the path...\n");
    tosl_arc_id_t ia0 = 2*(NS-1); /* Arc before arc 0. */
    assert(ia0 == 2*(NS-1));
    int32_t skip = NB + 2; /* How many ring edges to bypass. */
    tosl_arc_id_t ia1 = (ia0 + 2*skip) % (2*NS);
    
    fprintf(stderr, "adding the path across the ring from ");
    tosl_arc_id_print(stderr, "dst(", ia0, ")");
    tosl_arc_id_print(stderr, " to dst(", ia1, ")\n");
    
    tosl_mesh_add_path(NB, ia0, ia1, "p.", mesh);
     
    fprintf(stderr, "defining path vertex coordinates...\n");
    tosl_vert_id_t iv0 = mesh->Arc[tosl_sym(ia0)].ivorg; tosl_point_t *v0 = &(mesh->Vpos[iv0]);
    tosl_vert_id_t iv1 = mesh->Arc[tosl_sym(ia1)].ivorg; tosl_point_t *v1 = &(mesh->Vpos[iv1]);
    for (tosl_vert_id_t b = 0; b < NB; b++)
      { double f = ((double)b)/((double)NB);
        tosl_vert_id_t iv = NS + b;
        assert(iv < mesh->NV);
        tosl_point_t *v = &(mesh->Vpos[iv]);
        for (uint32_t j = 0;  j < 3; j++)
          { tosl_coord_t Cv = (tosl_coord_t)floor((1-f)*v0->c[j] + f*v1->c[j] + 0.5); Cv = 2*(Cv/2);
            v->c[j] = Cv;
          }
      }
    
    fprintf(stderr, "printing mesh...\n");
    if (debug) { tosl_mesh_print(stderr, mesh); } 

    fprintf(stderr, "checking mesh consistency...\n");
    tosl_mesh_check(mesh);
  }
