#define PROG_NAME "test_mroc_basic"
#define PROG_DESC "Basic tests of the marching octahedra algorithm"
#define PROG_VERS "1.0"
/* Last edited on 2021-07-08 02:25:37 by jstolfi */

#define test_mroc_basic_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"
  
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <r3.h>
#include <i3.h>

#include <mroc.h>
#include <mroc_stl.h>

/* INTERNAL PROTOTYPES */

#define tmb_RAD (3.0)
  /* Radius of sphere. */
  
#define tmb_CTR ((r3_t){{ tmb_RAD + 3.5, tmb_RAD + 3.5, tmb_RAD + 3.5 }})
  /* Center of sphere. */

typedef double tmb_func_t(r3_t *p);
  /* Type of function that computes the signed distance to surface of some object. */

double tmb_fudge(double v);
  /* Returns {v}, but pushing it a bit away from 0.5 if it is too close to 0.5. */
  
void do_tests
  ( tmb_func_t func,
    char *fname, 
    bool_t showTetras
  );
    
void tmb_eval_voxel_centers(tmb_func_t func, int32_t iz, int32_t NEX, int32_t NEY, double fcL[]);
  /* Evaluates the function {func} at the centers of the voxels in layer {iz},
    whose {Z} coordinates are {iz+0.5}.  Assumes that the layer has {NEX} by {NEY} voxels
    spanning the rectange {[-1 _ NEX-1] x [-1 _ NEY-1]}.  The values are stored into
    {fcL[0..NEX*NEY-1]} with the {X} index varying faster. */
    
void tmb_eval_voxel_corners(tmb_func_t func, int32_t iz, int32_t NVX, int32_t NVY, double fvL[]);
  /* Evaluates the function {func} at the bottom vertices of the voxels
    in layer {iz}, awhose {Z} coordinates are {iz}. Assumes that the
    voxel layer has {NVX-1} by {NVY-1} voxels and hence the vertex layer
    has {NVX} by {NVY} vertices, spanning the rectange {[0 _ NVX-1] x [0
    _ NVY-1]}. The values are stored into {fvL[0..NVX*NVY-1]} with the
    {X} index varying faster. */

double sphere_func(r3_t *p);
  /* Returns the distance to the surface of the sphere with center 
    {(5.5,5.5,5.5)} and radius 3.0, scaled from {[-2 _ +2]} to {[0 _ 1]}
    and clipped. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    do_tests(sphere_func, "out/basic_tet.stl", TRUE);
    do_tests(sphere_func, "out/basic_iso.stl", FALSE);
    
    return 0;
  }

double sphere_func(r3_t *p)
  { r3_t ctr = tmb_CTR;
    double rad = tmb_RAD;
    double d = r3_dist(p, &ctr) - rad;
    /* Map [-2 _ +2]} to {[0 _ 1]} and clip: */
    double v = fmax(0.0, fmin(1.0, (d + 2)/4));
    return v;
  }

void do_tests
  ( tmb_func_t func,
    char *fname, 
    bool_t showTetras
  )
  {
    float eps = 0.005f;
    fprintf(stderr, "--- testing basic {mroc.h} functions\n");
    
    bool_t debug = FALSE;
    bool_t verbose = TRUE;
    
    FILE *wr = open_write(fname, TRUE);
    fprintf(wr, "solid\n");
    
    /* Number of cells (voxels) along each axis of the grid: */
    r3_t ctr = tmb_CTR;
    int32_t NCX = (int32_t)ceil(2*ctr.c[0]);
    int32_t NCY = (int32_t)ceil(2*ctr.c[1]);
    int32_t NCZ = (int32_t)ceil(2*ctr.c[2]);
 
    /* Number of cells in each cell layer, including extrapolated cells: */
    int32_t NEX = NCX + 2;
    int32_t NEY = NCY + 2;
    
    /* Number of vertices (voxel corners) in each vertex layer: */
    int32_t NVX = NCX + 1;
    int32_t NVY = NCY + 1;
    
    /* Storage for func values at cell centers of three consecutive horizontal voxel layers: */
    double *thiscL = notnull(malloc(NEX*NEY*sizeof(double)), "no mem");
    double *nextcL = notnull(malloc(NEX*NEY*sizeof(double)), "no mem");
    
    /* Storage for func values at vertices of each horizontal cell layer: */
    double *botvL = notnull(malloc(NVX*NVY*sizeof(double)), "no mem"); /* Bot vertex layer */
    double *topvL = notnull(malloc(NVX*NVY*sizeof(double)), "no mem"); /* Top vertex layer */
    
    /* Initialize the arrays {botvL,prevcL,thiscL} to process voxel layer 0: */
    
    /* Initialize {botvL} with values at layer 0 of vertex grid: */
    tmb_eval_voxel_corners(func, 0, NVX,NVY, botvL);  

    /* Initialize {thiscL} with  values from cell layer 0: */
    tmb_eval_voxel_centers(func, 0, NEX,NEY, thiscL);
    
    int32_t nt = 0; /* Number of triangles written out. */
    int32_t ne = 0; /* Number of triangles discarded. */
    int32_t lastz = -2; /* Last {z} coordinate seen. */
    
    /* Tetrahedron-processing function: */
    auto void tetra_proc(r3_t *p[], double f[]);
      /* Processes the tetrahedron  with corners {p[0..3]} and occupancy 
        values {f[0..3]}. */

    for (uint32_t iz = 0;  iz < NCZ; iz++)
      { 
        /* Processing horizontal layer {iz} of voxels.
        
          Invariant: At this point, the arrays contain the value estimated
          at the following points:
            
            {botvL}: vertices at bottom of voxel layer {iz}.
            
            {thiscL}: centers of cells of voxel layer {iz}.
        
        */
        
        /* Evaluate function at top vertices of voxel layer {iz}: */
        tmb_eval_voxel_corners(func, iz+1, NVX,NVY, topvL);

        /* Apply marching octahedra to tetrahedra inside layer {iz}: */
        mroc_2D_intra(iz, NCX,NCY, thiscL, botvL, topvL, tetra_proc);
        
        /* Fill {nextcL} to voxel values from cell layer {iz+1}: */
        tmb_eval_voxel_centers(func, iz+1, NEX,NEY, nextcL);

        /* Apply marching octahedra to tetrahedra spanning layers {iz} and {iz+1}: */
        mroc_2D_inter(iz+1, NCX,NCY, thiscL, nextcL, topvL, tetra_proc);
        
        /* Shuffle the array pointers to prepare for the next layer: */
        { double *tmp = thiscL; thiscL = nextcL; nextcL = tmp; }
        { double *tmp = botvL; botvL = topvL; topvL = tmp; }
      }          
    
    free(thiscL); free(nextcL);
    free(botvL); free(topvL);
    
    fprintf(wr, "\n");
    fprintf(wr, "endsolid\n");
    fclose(wr);

    if (verbose) { fprintf(wr, "\n"); }
    fprintf(stderr, "%d triangles eliminated\n", ne);
    fprintf(stderr, "%d triangles written\n", nt);
    
    return;
    
    /* INTERNAL IMPLEMENTATION */
    void tetra_proc(r3_t *p[], double f[])
      {
        if (debug) { fprintf(stderr, "."); }
        bool_t newz = FALSE;

        for (uint32_t k = 0;  k < 4; k++)
          { int32_t thisz = (int32_t)floor(p[k]->c[2] - 0.499999);
            if (thisz > lastz) { newz = TRUE; lastz = thisz; }
          }
        if (newz && verbose) { fprintf(stderr, " %d", lastz); }
        r3_t u[4];
        for (uint32_t k = 0;  k < 4; k++) { u[k] = *(p[k]); }
        if (showTetras)
          { mroc_stl_write_tetra_faces(wr, u, f, &nt); }
        else
          { mroc_stl_write_surface_in_tetra(wr, u, f, eps, &nt, &ne); }
      }
  }

void tmb_eval_voxel_centers(tmb_func_t func, int32_t iz, int32_t NEX, int32_t NEY, double fcL[])
  {
    double *fcR = fcL; /* Start of next row in slice. */
    r3_t p; p.c[2] = iz + 0.5;
        
    for (uint32_t ky = 0;  ky < NEY; ky++)
      { /* Evaluate on layer : */
        p.c[1] = ky - 0.5;
        for (uint32_t kx = 0;  kx < NEX; kx++)
          { p.c[0] = kx - 0.5;
            fcR[kx] = tmb_fudge(func(&p));
          }
        fcR += NEX;
      }
  }

void tmb_eval_voxel_corners(tmb_func_t func, int32_t iz, int32_t NVX, int32_t NVY, double fvL[])
  { 
    r3_t p; p.c[2] = iz;
    double *fvR = fvL; /* Start of next row in slice. */
    for (uint32_t iy = 0;  iy < NVY; iy++)
      { p.c[1] = iy;
        for (uint32_t ix = 0;  ix < NVX; ix++)
          { p.c[0] = ix;
            fvR[ix] = tmb_fudge(func(&p));
          }
        fvR += NVX;
      }
  } 

double tmb_fudge(double v)
  {
    double eps = 1.0e-4;
    if (fabs(v - 0.5) < eps)
      { /* Too close to 0.5, fudge it: */
        if (v <= 0.5)
          { v = 0.5 - eps; }
        else 
          { v = 0.5 + eps; }
      }
    return v;
  }
  
