/* See mroc_ppv.h */
/* Last edited on 2021-07-08 15:54:17 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <ppv_types.h>
#include <ppv_array.h>
#include <r3.h>
#include <affirm.h>

#include <mroc.h>
#include <mroc_ppv.h>

/* INTERNAL PROTOTYPES */
    
void mroc_ppv_fill_array(ppv_size_t NX, ppv_size_t NY, double el[], double smp);
  /* Fills the array {el}, which is assumed to have {NX*NY} elements, with {smp}. */

void mroc_ppv_eval_voxel_centers(ppv_array_t *A, int32_t iz, double pad, double smpL[]);
  /* Extracts the horizontal layer {iz} from the pixel array {A},
    converts each sample to {double} in {[0 _ 1]}, and stores the
    result in {smpL}, with a frame of extra elements all around, set to
    value {pad}. If the first two dimensions of {A} are {NX,NY},
    then {smpL} is assumed to have {(NX+2)*(NY+2)} elements. */

void mroc_ppv_eval_voxel_corners(ppv_size_t NCX, ppv_size_t NCY, double botcL[], double topcL[], double midvL[]);
  /* Computes the values for one layer of the grid vertices, given the 
    cell center values {botcL,topcL} for the two adjacent cell layers.  The value at
    a vertex is the average of the center values of the 8 voxels incident to it.
    The arrays {botcL,topcL} should have {NCY2} rows of {NCX+2} elements,
    including a layer of padding all around.  The output array {midvL}
    must have {NCY+1} rows of {NCX+1} elements each. */

void mroc_ppv_process_array
  ( ppv_array_t *A,
    double pad,
    mroc_tetra_proc_t tetra_proc
  )
  {
    /* Number of cells (voxels) along each axis of the tomogram: */
    ppv_size_t NCX = A->size[2];  assert(NCX >= 1);
    ppv_size_t NCY = A->size[1];  assert(NCY >= 1);
    ppv_size_t NCZ = A->size[0];  assert(NCZ >= 1);
 
    /* Number of cells in each cell array, including extrapolated cells: */
    ppv_size_t NEX = NCX + 2;
    ppv_size_t NEY = NCY + 2;
    
    /* Values at cell centers of three consecutive horizontal voxel layers: */
    double *prevcL = notnull(malloc(NEX*NEY*sizeof(double)), "no mem");
    double *thiscL = notnull(malloc(NEX*NEY*sizeof(double)), "no mem");
    double *nextcL = notnull(malloc(NEX*NEY*sizeof(double)), "no mem");
    
    /* Number of vertices (voxel corners) along each axis: */
    ppv_size_t NVX = NCX + 1;
    ppv_size_t NVY = NCY + 1;
    
    /* Interpolated values at vertices of each horizontal voxel layer: */
    double *botvL = notnull(malloc(NVX*NVY*sizeof(double)), "no mem");
    double *topvL = notnull(malloc(NVX*NVY*sizeof(double)), "no mem");
    
    /* Initialize the arrays {botvL,prevcL,thiscL} to process voxel layer 0: */
    
    /* Initialize {prevcL} with the padding cell value: */
    mroc_ppv_fill_array(NEX,NEY, prevcL, pad);
    
    /* Initialize {thiscL} with  values from cell layer 0: */
    mroc_ppv_eval_voxel_centers(A, 0, pad, thiscL);
    
    /* Initialize {botvL} with values at layer 0 of vertex grid: */
    mroc_ppv_eval_voxel_corners(NCX,NCY, prevcL, thiscL, botvL);  
    
    /* Loop on cell layers: */
    int32_t iz;
    for (iz = 0; iz < NCZ; iz++)
      { 
        /* Processing horizontal layer {iz} of voxels.
        
          Invariant: At this point, the arrays contain the value estimated
          at the following points:
            
            {prevcL}: centers of cells of voxel layer {iz-1} (or padding if {iz=0}).
            
            {thiscL}: centers of cells of voxel layer {iz}.
        
            {botvL}: vertices at bottom of voxel layer {iz}.
            
        */
        
        /* Set {nextcL}: */
        if (iz < NCZ-1)
          { /* Set {nextcL} to voxel values from cell layer {iz+1} */
            mroc_ppv_eval_voxel_centers(A, iz+1, pad, nextcL);
          }
        else
          { /*  Fill {nextcL} with the padding cell value: */
            mroc_ppv_fill_array(NEX,NEY, nextcL, pad);
          }

        /* Set {topvL} by inerpolating cell center values of layers {iz} and {iz+1} */
        mroc_ppv_eval_voxel_corners(NCX,NCY, thiscL, nextcL, topvL);

        if (iz == 0)
          { /* Apply marching octahedra to tetrahedra spanning layers {-1} and 0: */
            mroc_2D_inter(iz, NCX,NCY, prevcL, thiscL, botvL, tetra_proc);
          }
        
        /* Apply marching octahedra to tetrahedra inside layer {iz}: */
        mroc_2D_intra(iz, NCX,NCY, thiscL, botvL, topvL, tetra_proc);
        
        /* Apply marching octahedra to tetrahedra spanning layers {iz} and {iz+1}: */
        mroc_2D_inter(iz+1, NCX,NCY, thiscL, nextcL, topvL, tetra_proc);
        
        /* Shuffle the array pointers to prepare for the next layer: */
        { double *tmp = thiscL; thiscL = nextcL; nextcL = prevcL; prevcL = tmp; }
        { double *tmp = botvL; botvL = topvL; topvL = tmp; }
      }          
    
    free(prevcL); free(thiscL); free(nextcL);
    free(botvL); free(topvL);
  }

void mroc_ppv_fill_array(ppv_size_t NX, ppv_size_t NY, double el[], double fsmp)
  { int32_t N = (int32_t)(NX*NY);
    int32_t k;
    for (k = 0; k < N; k++) { el[k] = fsmp; }
  }

void mroc_ppv_eval_voxel_centers(ppv_array_t *A, int32_t iz, double pad, double fsmpL[])
  {
    ppv_dim_t d = A->d;
    assert(d == 3);
    
   /* Number of cells (voxels) along each axis of the tomogram: */
    ppv_size_t NCX = A->size[2];
    ppv_size_t NCY = A->size[1];
    ppv_size_t NCZ = A->size[0];
 
    /* Number of cells in each cell array, including extrapolated cells: */
    ppv_size_t NEX = NCX + 2;

    demand((iz >= 0) && (iz < NCZ), "bad layer index");

    ppv_index_t ix[d];
    int32_t kx, ky;
    ix[0] = iz; 
    double *fsmpR = fsmpL; /* Start of next row in slice. */
    
    /* Store padding row before first row: */
    for (kx = 0; kx < NEX; kx++) { fsmpR[0 + kx] = pad; }
    fsmpR += NEX;
    
    for (ky = 0; ky < NCY; ky++)
      { /* Extract next row {ky} of voxels, floatize, store: */
        ix[1] = ky;
        fsmpR[0] = pad;
        for (kx = 0; kx < NCX; kx++)
          { ix[2] = kx;
            ppv_pos_t pos = ppv_sample_pos(A, ix);
            ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
            fsmpR[kx + 1] = mroc_ppv_floatize(smp, A->maxsmp);
          }
        fsmpR[NEX-1] = pad;
        fsmpR += NEX;
      }
      
    /* Store padding after last row: */
    for (kx = 0; kx < NEX; kx++) { fsmpR[0 + kx] = pad; }
  }

double mroc_ppv_floatize(ppv_sample_t smp, ppv_sample_t maxsmp) 
  { demand(smp <= maxsmp, "invalid sample value");
    double fsmp = ((double)smp)/((double)maxsmp);
    assert(fsmp != 0.5);
    return fsmp; 
  }

void mroc_ppv_eval_voxel_corners(ppv_size_t NCX, ppv_size_t NCY, double botcL[], double topcL[], double midvL[])
  { 

    /* Number of cells in each cell array, including extrapolated cells: */
    ppv_size_t NEX = NCX + 2;
    
    /* Number of vertices (voxel corners) along each axis: */
    ppv_size_t NVX = NCX + 1;
    ppv_size_t NVY = NCY + 1;
    
    int32_t ix, iy;
    
    for (iy = 0; iy < NVY; iy++)
      { /* Compute row {iy} of {midvL}: */
        double *antbotcR = &(botcL[iy*NEX]);
        double *posbotcR = &(botcL[(iy+1)*NEX]);
        double *anttopcR = &(topcL[iy*NEX]);
        double *postopcR = &(topcL[(iy+1)*NEX]);
        double *midvR = &(midvL[iy*NVX]);
        for (ix = 0; ix < NVX; ix++)
          { /* Compuet the average {avg} of the 8 surrounding voxels: */
            double bot_sum = antbotcR[ix] + antbotcR[ix+1] + posbotcR[ix] + posbotcR[ix+1];
            double top_sum = anttopcR[ix] + anttopcR[ix+1] + postopcR[ix] + postopcR[ix+1];
            double avg = (bot_sum + top_sum)/8.0;
            
            /* Make sure that the value is not 0.5: */
            double eps = 1.0e-4;
            if (fabs(avg - 0.5) < eps)
              { /* Too close to 0.5, fudge it: */
                if (avg <= 0.5)
                  { avg = 0.5 - eps; }
                else 
                  { avg = 0.5 + eps; }
              }
            
            /* Store it: */
            midvR[ix] = avg;
          }
      }
  } 

