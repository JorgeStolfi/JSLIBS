/* See mroc.h */
/* Last edited on 2016-03-10 20:47:34 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <bool.h>
#include <ppv_types.h>
#include <r3.h>
#include <affirm.h>

#include <mroc.h>

/* INTERNAL PROTOTYPES */


void mroc_1D_inter
  ( int32_t iy,
    int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *thiscL, 
    double *botvL, 
    double *topvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  );
  /* Processes the tetrahedra that span rows {iy-1} and {iy} in 
    the horizontal layer {iz} of the tomogram.  Assumes that the
    cell-centered values are stored in
    the array {thiscL} with {NCY+2} rows of {NCX+2} elements
    (including padding border elements).  Assumes that {botvL} 
    and {topvL} are arrays with {NCY+1} rows of {NCX+1} elements,
    containing the values interpolated at the 
    vertices respectively below and above that voxel layer. */

void mroc_1D_intra
  ( int32_t iy,
    int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *thiscL, 
    double *botvL, 
    double *topvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  );
  /* Processes the tetrahedra that span adjacent voxels of row {iy} in 
    the horizontal layer {iz} of the tomogram.  Assumes that the
    cell-centered values are stored in
    the array {thiscL} with {NCY+2} rows of {NCX+2} elements
    (including padding border elements).  Assumes that {botvL} 
    and {topvL} are arrays with {NCY+1} rows of {NCX+1} elements,
    containing the values interpolated at the 
    vertices respectively below and above that voxel layer. */

/* IMPLEMENTATIONS */

void mroc_2D_inter
  ( int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *antcL, 
    double *poscL, 
    double *midvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  )
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "enter %s(%d)\n", __FUNCTION__, iz); }

    /* Size of vertex arrays: */
    ppv_size_t NVX = NCX + 1;
    
    /* Size of cell arrays, including extrapolated cells: */
    ppv_size_t NEX = NCX + 2;

    int32_t ix, iy;
    
    /* Corners of an octahedron spanning a voxel face, and their values: */
    r3_t pu[4];    /* Equatorial points: vertices of a voxel face, in circular order. */
    double fu[4];  /* Values at those vertices. */
    r3_t pv[2];    /* Polar points: centers of two voxels incident to that face. */
    double  fv[2]; /* Values at those cell centers. */
    
    for (iy = 0; iy < NCY; iy++)
      { 
        /* Get the relevant rows of the cell value arrays (skipping the padding cells): */
        double *antcR = &(antcL[(iy+1)*NEX + 1]);
        double *poscR = &(poscL[(iy+1)*NEX + 1]);

        /* Get the relevant rows of the vertex value array: */
        double *antvR = &(midvL[iy*NVX]);
        double *posvR = &(midvL[(iy+1)*NVX]);
        
        for (ix = 0; ix < NCX; ix++)
          { 
            /* Get the equatorial vertices: */
            pu[0] = (r3_t){{ ix,   iy,   iz }}; fu[0] = antvR[ix];
            pu[1] = (r3_t){{ ix+1, iy,   iz }}; fu[1] = antvR[ix+1];
            pu[2] = (r3_t){{ ix+1, iy+1, iz }}; fu[2] = posvR[ix+1];
            pu[3] = (r3_t){{ ix,   iy+1, iz }}; fu[3] = posvR[ix];
            
            /* Get the polar vertices: */
            pv[0] = (r3_t){{ ix+0.5, iy+0.5, iz-0.5 }}; fv[0] = antcR[ix];
            pv[1] = (r3_t){{ ix+0.5, iy+0.5, iz+0.5 }}; fv[1] = poscR[ix];
            
            /* Process the octahedron: */
            mroc_process_octahedron(pu, fu, pv, fv, tetra_proc);
         }
     }
    if (debug) { fprintf(stderr, "\nexit %s\n", __FUNCTION__); }
  }
  
void mroc_2D_intra
  ( int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *thiscL, 
    double *botvL, 
    double *topvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  )
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "enter %s(%d)\n", __FUNCTION__, iz); }

    int32_t iy;
    for (iy = 0; iy < NCY; iy++)
      { 
        if (iy == 0) 
          { /* Process the tetrahedra that span rows {iy-1} and{iy}: */
            mroc_1D_inter(iy, iz, NCX,NCY, thiscL, botvL, topvL, tetra_proc); 
          }
          
        /* Process the tetrahedra that span adjacent voxels of row {iy}: */
        mroc_1D_intra(iy, iz, NCX,NCY, thiscL, botvL, topvL, tetra_proc);
         
        /* Process the tetrahedra that span rows {iy} and {iy+1}: */
        mroc_1D_inter(iy+1, iz, NCX,NCY, thiscL, botvL, topvL, tetra_proc);
         
      }
    if (debug) { fprintf(stderr, "\nexit %s\n", __FUNCTION__); }
  }

void mroc_1D_inter
  ( int32_t iy,
    int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *thiscL, 
    double *botvL, 
    double *topvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  )
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  enter %s(%d,%d)\n", __FUNCTION__, iy,iz); }

    /* Size of vertex arrays: */
    ppv_size_t NVX = NCX + 1;
    
    /* Size of cell arrays, including extrapolated cells: */
    ppv_size_t NEX = NCX + 2;
    
    /* Get the relevant rows {iy-1} and {iy} of the cell value arrays (skipping the padding cells): */
    double *antcR = &(thiscL[iy*NEX + 1]);
    double *poscR = &(thiscL[(iy+1)*NEX + 1]);

    /* Get the relevant rows {iy} of the vertex value arrays: */
    double *botvR = &(botvL[iy*NVX]);
    double *topvR = &(topvL[iy*NVX]);

    /* Corners of an octahedron spanning a voxel face, and their values: */
    r3_t pu[4];    /* Equatorial points: vertices of a voxel face, in circular order. */
    double fu[4];  /* Values at those vertices. */
    r3_t pv[2];    /* Polar points: centers of two voxels incident to that face. */
    double fv[2];  /* Values at those cell centers. */

    int32_t ix;
    
    for (ix = 0; ix < NCX; ix++)
      { 
        /* Process the tetrahedra that span voxels {ix} in rows {iy-1} and {iy}: */

        /* Get the equatorial vertices: */
        pu[0] = (r3_t){{ ix,   iy,   iz   }}; fu[0] = botvR[ix];
        pu[1] = (r3_t){{ ix+1, iy,   iz   }}; fu[1] = botvR[ix+1];
        pu[2] = (r3_t){{ ix+1, iy,   iz+1 }}; fu[2] = topvR[ix+1];
        pu[3] = (r3_t){{ ix,   iy,   iz+1 }}; fu[3] = topvR[ix];

        /* Get the polar vertices: */
        pv[0] = (r3_t){{ ix+0.5, iy-0.5, iz+0.5 }}; fv[0] = antcR[ix];
        pv[1] = (r3_t){{ ix+0.5, iy+0.5, iz+0.5 }}; fv[1] = poscR[ix];

        /* Process the octahedron: */
        mroc_process_octahedron(pu, fu, pv, fv, tetra_proc);
      }
    if (debug) { fprintf(stderr, "\n  exit %s\n", __FUNCTION__); }
  }

void mroc_1D_intra
  ( int32_t iy,
    int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *thiscL, 
    double *botvL, 
    double *topvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  )
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  enter %s(%d,%d)\n", __FUNCTION__, iy,iz); }
    
    /* Size of vertex arrays: */
    ppv_size_t NVX = NCX + 1;
    
    /* Size of cell arrays, including extrapolated cells: */
    ppv_size_t NEX = NCX + 2;
    
    /* Get the relevant row {iy} of the cell value array (INCLUDIG the first padding cell of row): */
    double *thiscR = &(thiscL[(iy+1)*NEX]);

    /* Get the relevant rows {iy} and {iy+1} of the vertex value arrays: */
    double *antbotvR = &(botvL[iy*NVX]);
    double *posbotvR = &(botvL[(iy+1)*NVX]);
    double *anttopvR = &(topvL[iy*NVX]);
    double *postopvR = &(topvL[(iy+1)*NVX]);

    /* Corners of an octahedron spanning a voxel face, and their values: */
    r3_t pu[4];    /* Equatorial points: vertices of a voxel face, in circular order. */
    double fu[4];  /* Values at those vertices. */
    r3_t pv[2];    /* Polar points: centers of two voxels incident to that face. */
    double  fv[2]; /* Values at those cell centers. */

    int32_t ix;
    
    for (ix = 0; ix < NVX; ix++)
      { 
        /* Process the tetrahedra that span voxels {ix-1} and {ix} in row {iy}: */

        /* Get the equatorial vertices: */
        pu[0] = (r3_t){{ ix, iy,   iz   }}; fu[0] = antbotvR[ix];
        pu[1] = (r3_t){{ ix, iy+1, iz   }}; fu[1] = posbotvR[ix];
        pu[2] = (r3_t){{ ix, iy+1, iz+1 }}; fu[2] = postopvR[ix];
        pu[3] = (r3_t){{ ix, iy,   iz+1 }}; fu[3] = anttopvR[ix];

        /* Get the polar vertices: */
        pv[0] = (r3_t){{ ix-0.5, iy+0.5, iz+0.5 }}; fv[0] = thiscR[ix];
        pv[1] = (r3_t){{ ix+0.5, iy+0.5, iz+0.5 }}; fv[1] = thiscR[ix+1];

        /* Process the octahedron: */
        mroc_process_octahedron(pu, fu, pv, fv, tetra_proc);
      }
    if (debug) { fprintf(stderr, "\n  exit %s\n", __FUNCTION__); }
 }
    
void mroc_process_octahedron
  ( r3_t pu[], 
    double fu[], 
    r3_t pv[], 
    double fv[],
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  )    
  {
    r3_t *p[4];
    double f[4];
    p[0] = &(pv[0]);
    f[0] = fv[0];
    p[1] = &(pv[1]);
    f[1] = fv[1];
    
    int32_t k0;
    for (k0 = 0; k0 < 4; k0++)
      { p[2] = &(pu[k0]);
        f[2] = fu[k0];
        int32_t k1 = (k0+1)%4;
        p[3] = &(pu[k1]);
        f[3] = fu[k1];
    
        tetra_proc(p, f);
     }
  
  }
    
