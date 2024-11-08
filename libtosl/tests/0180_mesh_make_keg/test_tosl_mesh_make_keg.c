/*  Last edited on 2024-10-08 18:56:45 by stolfi */
/* Test of {tosl_mesh_make_keg.h} */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_mesh_obj_write.h>
#include <tosl_mesh_make_keg.h>

int32_t main(int32_t argc, char **argv);

void do_test(int32_t NS, int32_t NR, int32_t NB, int8_t debug);
  /* Tests the procedure {tosl_mesh_make_keg(NS,NR,NB,...)} with different
    values of {NS,NR,NB}. */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615);
    
    do_test(4, 1, 0, 1);
    do_test(4, 1, 1, 1);
    do_test(4, 1, 2, 1);
    do_test(9, 3, 1, 0);
    do_test(9, 3, 2, 0);
    do_test(9, 4, 2, 0);
    do_test(9, 6, 3, 0);
    do_test(100, 200, 10, 0);
    
    return 0;
  }
  
void do_test(int32_t NS, int32_t NR, int32_t NB, int8_t debug)
  {
    fprintf(stderr, "start test with NS = %d NR = %d NB = %d...\n", NS, NR, NB);
        
    /* Decide the size of the barrel: */
    tosl_coord_t Zmax = 100*NR*(NB+1);
    
    /* Estimate number of edges: */
    int32_t NE_exp = NS*(NR+1) + NS*NR*(NB+1); 

    fprintf(stderr, "creating mesh...\n");
    tosl_mesh_t *mesh = tosl_mesh_make_keg(NS, NR, NB, Zmax, debug);
    assert(mesh->NE == NE_exp);
    
    if (debug) { tosl_mesh_print(stderr, mesh); } 
     
    fprintf(stderr, "checking mesh...\n");
    tosl_mesh_check(mesh);
    
    fprintf(stderr, "writing as OBJ file...\n");
    char *fname_obj= NULL;
    asprintf(&fname_obj, "out/test_ns%05d_nr%05d_nb%05d.obj", NS, NR, NB);
    FILE *wr_obj = fopen(fname_obj, "w");
    assert(wr_obj != NULL);
    tosl_mesh_obj_write(wr_obj, mesh);
    fclose(wr_obj);
    free(fname_obj);
    tosl_mesh_free(mesh);
  }
