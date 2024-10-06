/*  Last edited on 2024-10-06 16:49:44 by stolfi */
/* Test of {tosl_mesh_make_ico.h} */

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
#include <tosl_mesh_make_ico.h>

int32_t main(int32_t argc, char **argv);

void do_test(int8_t debug);
  /* Tests the procedure {tosl_mesh_make_ico(,...)}. */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615);
    do_test(1);
    return 0;
  }
  
void do_test(int8_t debug)
  {
    fprintf(stderr, "start test ...\n");
        
    /* Decide the apothem radius: */
    double R = 100;
    
    fprintf(stderr, "creating mesh...\n");
    tosl_mesh_t *mesh = tosl_mesh_make_ico(R);
    assert(mesh->NV == 12);
    assert(mesh->NE == 30);
    
    if (debug) { tosl_mesh_print(stderr, mesh); } 
     
    fprintf(stderr, "checking mesh...\n");
    tosl_mesh_check(mesh);
    
    fprintf(stderr, "writing as OBJ file...\n");
    char *fname_obj= NULL;
    asprintf(&fname_obj, "out/test.obj");
    FILE *wr_obj = fopen(fname_obj, "w");
    assert(wr_obj != NULL);
    tosl_mesh_obj_write(wr_obj, mesh);
    fclose(wr_obj);
    free(fname_obj);
    tosl_mesh_free(mesh);
  }
