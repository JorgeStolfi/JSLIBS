/* See {neuromat_eeg_image_basis.h}. */
/* Last edited on 2014-01-12 20:46:35 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <r3_extra.h>
#include <rn.h>
#include <rmxn.h>
#include <lsq.h>
#include <float_image.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_image.h>
#include <neuromat_eeg_image.h>
#include <neuromat_eeg_image_basis.h>
                
void neuromat_eeg_image_basis_compute_shepard_weights(int ne, r3_t pos3D[], double erad[], int order, r3_t *p, double wraw[]);
  /* Evaluates {ne} Shepard radial elements of given {order}, for
    centers {pos3D[0..ne-1]} and scaling radii {erad[0..ne-1]}, at the
    point {*p}. Assumes that all those points are on the unit sphere.
    Returns the weights in {wraw[0..ne-1]}, not normalized. */

float_image_t **neuromat_eeg_image_basis_make(int btype, float_image_t *msk, int ne, r2_t pos[], r2_t *ctr, r2_t *rad)
  {
    demand(msk->sz[0] == 1, "mask should be monochromatic");
    int NX = (int)msk->sz[1];
    int NY = (int)msk->sz[2];
    
    /* Convert 2D positions in unit disk to 3D points on the unit spehre: */
    r3_t *pos3D = notnull(malloc(ne*sizeof(r3_t)), "no mem");
    int ie;
    for (ie = 0; ie < ne; ie++) { pos3D[ie] = neuromat_eeg_geom_3D_from_2D(&(pos[ie])); }
    
    /* Allocate the basis images: */
    float_image_t **bas = notnull(malloc(ne*sizeof(float_image_t*)), "no mem");
    for (ie = 0; ie < ne; ie++) { bas[ie] = float_image_new (1, NX, NY); }
    
    switch(btype)
      { 
        case neuromat_eeg_image_basis_type_SHEPARD0:
          { neuromat_eeg_image_basis_fill_shepard_0(msk, ne, bas, ctr, rad, pos3D);
          }
          break;
        case neuromat_eeg_image_basis_type_SHEPARD1:
          { neuromat_eeg_image_basis_fill_shepard_1(msk, ne, bas, ctr, rad, pos3D);
          }
          break;
        case neuromat_eeg_image_basis_type_RADIAL:
          { double rho = 2.0;
            neuromat_eeg_image_basis_fill_radial(msk, ne, bas, ctr, rad, pos3D, rho);
          }
          break;
        case neuromat_eeg_image_basis_type_VORONOI:
          { neuromat_eeg_image_basis_fill_voronoi(msk, ne, bas, ctr, rad, pos3D);
          }
          break;
        default:
          assert(FALSE);
      }
      
    free(pos3D);
    return bas;
  }

void neuromat_eeg_image_basis_find_elem_radii(int ne, r3_t pos3D[], double erad[]);
  /* Given {ne} electrode position on the unit sphere {pos3D[0..ne-1]} Stores into {erad[0..ne-1]}
    the distance from each electrode to its nearest neighbor. If there is only one
    electrode, sets {erad[0]} to {+INF}. */

void neuromat_eeg_image_basis_find_elem_radii(int ne, r3_t pos3D[], double erad[])
  {
    int ie, je;
    for (ie = 0; ie < ne; ie++)
      { /* Find the distance squared {d2min} from {pos[ie]} to its nearest neighbor: */
        r3_t *pi = &(pos3D[ie]);
        double d2min = +INF;
        for (je = 0; je < ne; je++)
          { if (ie != je) 
              { r3_t *pj = &(pos3D[je]);
                double d2j = r3_dist_sqr(pi, pj);
                if (d2j < d2min) { d2min = d2j; }
              }
          }
        erad[ie] = sqrt(d2min);
        demand(d2min > 1.0e-5, "two electrodes are nearly coincident");
      }
  }

void neuromat_eeg_image_basis_fill_shepard_0(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[])
  {
    bool_t verbose = FALSE;

    if (ne == 0) { return; }
    
    demand(msk->sz[0] == 1, "mask must be monochromatic");
    int NX = (int)msk->sz[1]; 
    int NY = (int)msk->sz[2]; 
    
    /* Check image sizes: */
    int ie;
    for (ie = 0; ie < ne; ie++)
      { demand(bas[ie]->sz[0] == 1, "basis element must be monochromatic");
        demand(bas[ie]->sz[1] == NX, "basis element with wrong {NX}");
        demand(bas[ie]->sz[2] == NY, "basis element with wrong {NY}");
      }
    
    /* Choose the nominal radius squared {erad[i]} of each element, namely {dist(nearest_nbr)}: */
    double *erad = rn_alloc(ne);
    neuromat_eeg_image_basis_find_elem_radii(ne, pos3D, erad);
    
    /* Now loop on pixels, paint every image: */
    double *wraw = rn_alloc(ne); /* Element values at a certain point. */
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { bool_t debugpx = (ix == NX/2) & (iy == NY/2);
            /* Get mask weight of pixel: */
            double mxy = float_image_get_sample(msk, 0, ix, iy);
            demand((mxy >= 0) && (mxy <= 1.0), "invalid mask value");
            /* Set {wraw[ie]} to value of element {ie} at center of pixel {ix,iy}: */
            if (mxy <= 0)
              { for (ie = 0; ie < ne; ie++) { wraw[ie] = 0; } }
            else
              { /* Compute the position {pxyz} of the pixel center on the unit sphere: */
                r2_t qxy = (r2_t){{ ix + 0.5, iy + 0.5 }}; /* Point in image domain. */
                r2_t pxy = neuromat_eeg_geom_disk_from_ellipse(&qxy, ctr, rad); /* Corresp point in unit-disk schematic. */
                r3_t pxyz = neuromat_eeg_geom_3D_from_2D(&pxy); /* Corresp point on unit sphere. */
                if (verbose && debugpx){ r3_gen_print(stderr, &pxyz, "%+8.5f", "  pxyz = ( ", " ", " )\n"); }

                /* Set {wraw[0..ne-1]} to the normalized Shepard weights of each element at {pxyz}: */
                neuromat_eeg_image_basis_compute_shepard_weights(ne, pos3D, erad, 0, &pxyz, wraw);

                /* Now normalize so that elements add to 1: */
                double sum = rn_sum(ne, wraw);
                assert(sum > 0);
                for (ie = 0; ie < ne; ie++) 
                  { float fval = (float)(mxy*wraw[ie]);
                    if (verbose && debugpx){ fprintf(stderr, "    wraw[%3d] = %+8.5f fval = %+8.5f\n", ie, wraw[ie], fval); }
                    float_image_set_sample(bas[ie], 0, ix, iy, fval);
                  }
              }
            /* Now set pixels to {wraw[0..ne-1]}, multiplied by mask weight: */ 
            for (ie = 0; ie < ne; ie++) { float_image_set_sample(bas[ie], 0, ix, iy, (float)(mxy*wraw[ie])); }
          }
        if (verbose) { fprintf(stderr, "."); if ((iy == NY-1) || ((iy % 50) == 0)) { fprintf(stderr, "\n"); } }
      }
      
    free(erad);
    free(wraw);
  }

void neuromat_eeg_image_basis_fill_shepard_1(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[])
  {
    bool_t verbose = FALSE;

    if (ne == 0) { return; }
    
    demand(msk->sz[0] == 1, "mask must be monochromatic");
    int NX = (int)msk->sz[1]; 
    int NY = (int)msk->sz[2]; 
    
    /* Check image sizes: */
    int ie;
    for (ie = 0; ie < ne; ie++)
      { demand(bas[ie]->sz[0] == 1, "basis element must be monochromatic");
        demand(bas[ie]->sz[1] == NX, "basis element with wrong {NX}");
        demand(bas[ie]->sz[2] == NY, "basis element with wrong {NY}");
      }
    
    /* Choose the nominal radius squared {erad[i]} of each element, namely {dist(nearest_nbr)}: */
    double *erad = rn_alloc(ne);
    neuromat_eeg_image_basis_find_elem_radii(ne, pos3D, erad);
    
    /* Now loop on pixels, paint every image: */
    double *bval = rn_alloc(ne); /* Element values at a certain point. */
    double *wraw = rn_alloc(ne); /* (WORK) Raw shepard weights at certain point. */
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { bool_t debugpx = (ix == NX/2) & (iy == NY/2);
            /* Get mask weight of pixel: */
            double mxy = float_image_get_sample(msk, 0, ix, iy);
            demand((mxy >= 0) && (mxy <= 1.0), "invalid mask value");
            /* Set {wraw[ie]} to value of element {ie} at center of pixel {ix,iy}: */
            if (mxy <= 0)
              { for (ie = 0; ie < ne; ie++) { bval[ie] = 0; } }
            else
              { /* Compute the position {pxyz} of the pixel center on the unit sphere: */
                r2_t qxy = (r2_t){{ ix + 0.5, iy + 0.5 }}; /* Point in image domain. */
                r2_t pxy = neuromat_eeg_geom_disk_from_ellipse(&qxy, ctr, rad); /* Corresp point in unit-disk schematic. */
                r3_t p = neuromat_eeg_geom_3D_from_2D(&pxy); /* Corresp point on unit sphere. */
                if (verbose && debugpx){ r3_gen_print(stderr, &p, "%+8.5f", "  p = ( ", " ", " )\n"); }

                /* Set {wraw[0..ne-1]} to the unnormalized Shepard weights of each element at {pxyz}: */
                neuromat_eeg_image_basis_compute_shepard_weights(ne, pos3D, erad, 1, &p, wraw);
                
                /* 
                  Now find for each electrode {pos3D[k]} a linear combination {s} of the LSQ basis functions
                  {\phi[0..nb-1]} such that, for each {j}, {s(pos3D[j])} is as close to {(j==k)} as possible,
                  in the LSQ sense with weight {wraw[j]}.
                */

                int nb = 3; /* Number of elements in the LSQ basis. */
                
                /* Choose two arbitrary unit vectors {u[0..1]} orthogonal to {p}: */
                r3_t u[2];
                { double pm = r3_pick_ortho(&p, &(u[0])); assert(pm >= 0.1);
                  double ue = r3_dir(&(u[0]), &(u[0])); assert(fabs(ue - 1.0) < 0.001);
                  r3_cross(&p, &(u[0]), &(u[1]));
                  double ve = r3_dir(&(u[1]), &(u[1])); assert(fabs(ve - 1.0) < 0.001);
                }

                auto double phi(int j, r3_t *q);
                  /* Computes the LSQ basis function {\phi[j]} at point {q}, for {j} in {0..nb-1}. 
                    Namely, {\phi[0](q) = dot(q-p,u)}, {\phi[1](q) = dot(q-p,v)}, {\phi[2](q) = 1}. */

                double phi(int j, r3_t *q)
                  { assert((j >= 0) && (j < nb));
                    if (j <= 1)
                      { r3_t d;
                        r3_sub(q, &p, &d);
                        return r3_dot(&d, &(u[j]));
                      }
                    else if (j == 2)
                      { return 1.0; }
                    else
                      { assert(FALSE); }
                  }
                    
                auto void gen_case(int i, int nv, double vi[], int nf, double fi[], double *wiP);
                  /* Generates the LSQ basis values {vi[0..nv-1]} at electrode {i}, the 
                    corresponding function values {fi[0..nf-1]}, and the weight {*wiP}.
                    Namely, sets {vi[j] = \phi[j](pos3D[i])}, {fi[k] = (k==i)}, and {*wiP = wraw[i]}.
                    Requires {nv==nb} and {nf==ne}. */ 
                    
                void gen_case(int i, int nv, double vi[], int nf, double fi[], double *wiP)
                  { assert(nv == nb);
                    assert(nf == ne);
                    assert((i >= 0) && (i < ne));
                    int j, k;
                    for (j = 0; j < nv; j++) { vi[j] = phi(j, &(pos3D[i])); }
                    for (k = 0; k < nf; j++) { fi[j] = (k == i ? 1 : 0); }
                    (*wiP) = wraw[i];
                  }
                  
                /* Compute the linear function: */
                bool_t verbose_lsq = FALSE;
                double *U = rmxn_alloc(nb, ne);
                lsq_fit(nb, ne, ne, gen_case, U, verbose_lsq);
                
                /* Now evaluate the linear functions at {p}: */
                /* Should be just {U[nb-1,ie]} since {\phi(p) = (0,0,1)}, but let's play safe: */
                double phip[nb];
                { int j; for (j = 0; j < nb; j++) { phip[j] = phi(j, &p); } }
                rmxn_map_row(nb, ne, phip, U, bval);
                
                free(U);
              }
            /* Now set pixels: */ 
            for (ie = 0; ie < ne; ie++) 
              { float fval = (float)(mxy*bval[ie]);
                if (verbose && debugpx){ fprintf(stderr, "    bval[%3d] = %+8.5f fval = %+8.5f\n", ie, bval[ie], fval); }
                float_image_set_sample(bas[ie], 0, ix, iy, fval);
              }
          }
        if (verbose) { fprintf(stderr, "."); if ((iy == NY-1) || ((iy % 50) == 0)) { fprintf(stderr, "\n"); } }
      }
      
    free(erad);
    free(wraw);
    free(bval);
  }
    
void neuromat_eeg_image_basis_compute_shepard_weights(int ne, r3_t pos3D[], double erad[], int order, r3_t *p, double wraw[])
  {
    int ie;
    for (ie = 0; ie < ne; ie++)
      { r3_t *pi = &(pos3D[ie]);
        double Ri = 4*erad[ie];
        /* Compute {z2} = distance to electrode {ie}, relative to {4*erad[ie]}, squared: */
        double d2 = r3_dist_sqr(p, pi);
        double z2 = d2/(Ri*Ri);
        /* Weight is {1/z^2} times a Gaussian: */
        double val = ( z2 == 0 ? 1.0e+30 : ( z2 >= 90 ? 0.0 : exp(-z2/2)/z2 ) );
        /* Make sure it is not zero: */
        wraw[ie] = val + 1.0e-30;
      }
  }

void neuromat_eeg_image_basis_fill_radial(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[], double rho)
  {
    /* bool_t verbose = FALSE; */
    bool_t debug_rval = FALSE;
    
    if (ne == 0) { return; }

    demand(msk->sz[0] == 1, "mask must be monochromatic");
    int NX = (int)msk->sz[1]; 
    int NY = (int)msk->sz[2]; 
    
    /* Check image sizes: */
    int ie;
    for (ie = 0; ie < ne; ie++)
      { demand(bas[ie]->sz[0] == 1, "basis element must be monochromatic");
        demand(bas[ie]->sz[1] == NX, "basis element with wrong {NX}");
        demand(bas[ie]->sz[2] == NY, "basis element with wrong {NY}");
      }

    /* We start with a basis of radial functions on the units sphere. 
      Each element is centered at one electrode position and has a gauss 
      times sinc shape with first root at the nearest neighbor.  Then this 
      basis is converted to an interpolating one. */
    
    /* Choose the nominal radius squared {erad[i]} of each element, namely {dist(nearest_nbr)}: */
    double *erad = rn_alloc(ne);
    neuromat_eeg_image_basis_find_elem_radii(ne, pos3D, erad);
    
    auto double mexican(double d, double sigma, double tau);
      /* The mother function of the radial basis, with width parameter {sigma}
        and node spacing {tau}, evaluated at distance {d} from
        the center. */
    
    auto double radial(int j, r3_t *p);
      /* The raw basis element centered at {pos3D[j]} evaluated at {p}. */

    /* Compute the Lagrangian matrix {L}: */
    double *L = rmxn_alloc(ne, ne); /* Matrix that converts radial basis to interp basis. */
    neuromat_eeg_image_basis_lagrangian_matrix(ne, pos3D, radial, L);

    /* Now fill the images with the Lagrangian basis derived from the radial basis: */
    double *rval = rn_alloc(ne); /* Raw element values at a certain point. */
    double *ival = rn_alloc(ne); /* Interp element values at a certain point. */
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { /* Get mask weight of pixel: */
            double mxy = float_image_get_sample(msk, 0, ix, iy);
            demand((mxy >= 0) && (mxy <= 1.0), "invalid mask value");
            /* Set {bval[ie]} to value of element {ie} at center of pixel {ix,iy}: */
            if (mxy <= 0)
              { for (ie = 0; ie < ne; ie++) { rval[ie] = ival[ie] = 0; } }
            else
              { /* Compute the position {pxyz} of the pixel center on the unit sphere: */
                r2_t qxy = (r2_t){{ ix + 0.5, iy + 0.5 }}; /* Point in image domain. */
                r2_t pxy = neuromat_eeg_geom_disk_from_ellipse(&qxy, ctr, rad); /* Corresp point in unit-disk schematic. */
                r3_t pxyz = neuromat_eeg_geom_3D_from_2D(&pxy); /* Corresp point on unit sphere. */

                /* Set {bval[0..ne-1]} to the raw radial basis values of each element at {pxyz}: */
                for (ie = 0; ie < ne; ie++) { rval[ie] = radial(ie, &pxyz); }
                
                /* Convert to {rval[0..ne-1]} to interpolating basis values {ival[0..ne-1]}: */
                if (debug_rval)
                  { rn_copy(ne, rval, ival); }
                else
                  { rmxn_map_col(ne, ne, L, rval, ival); }

                /* Scale by mask weight: */
                for (ie = 0; ie < ne; ie++) { ival[ie] = ival[ie] * mxy; }
              }
            /* Now set pixels: */ 
            for (ie = 0; ie < ne; ie++) { float_image_set_sample(bas[ie], 0, ix, iy, (float)(ival[ie])); }
          }
      }
      
    free(erad);
    free(rval);
    free(ival);
    free(L);

    return;
                    
    /* Implementation of local functions: */
               
    double radial(int j, r3_t *p)
      { r3_t *cj = &(pos3D[j]);
        double d = r3_dist(cj, p);
        return mexican(d, rho*erad[j], erad[j]); 
      }
    
    double mexican(double d, double sigma, double tau)
      { if (d == 0) { return 1.0; }
        double tbell = d/sigma;
        if (tbell >= 8.5) { return 0.0; }
        double bell = exp(-tbell*tbell/2);
        double tsinc = d/tau;
        double sinc = sin(M_PI*tsinc)/(M_PI*tsinc);
        return bell*sinc;
      }
  }
    
void neuromat_eeg_image_basis_fill_voronoi(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[])
  {
    /* bool_t verbose = TRUE; */

    if (ne == 0) { return; }

    demand(msk->sz[0] == 1, "mask must be monochromatic");
    int NX = (int)msk->sz[1]; 
    int NY = (int)msk->sz[2]; 
    
    /* Check image sizes: */
    int ie;
    for (ie = 0; ie < ne; ie++)
      { demand(bas[ie]->sz[0] == 1, "basis element must be monochromatic");
        demand(bas[ie]->sz[1] == NX, "basis element with wrong {NX}");
        demand(bas[ie]->sz[2] == NY, "basis element with wrong {NY}");
      }
    
    /* Now fill the basis images: */
    int *nnear = notnull(malloc(ne*sizeof(int)), "no mem"); /* Counts subsamples of pixel in each voronoi region. */
    double *bval = rn_alloc(ne); /* Element values at a certain point. */
    int msub = 5;  /* Subsampling order. */
    int ix, iy, k;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { /* Get mask weight of pixel: */
            double mxy = float_image_get_sample(msk, 0, ix, iy);
            demand((mxy >= 0) && (mxy <= 1.0), "invalid mask value");
            /* Set {bval[ie]} to value of element {ie} at center of pixel {ix,iy}: */
            if (mxy <= 0)
              { for (ie = 0; ie < ne; ie++) { bval[ie] = 0; } }
            else
              { /* Sample {msub} by {msub} points {p} in pixel {ix,iy}, count in {nnear[k]}:  */
                for (k = 0; k < ne; k++) { nnear[k] = 0; }
                int dx, dy;
                for (dy = 0; dy < msub; dy++) 
                  { for (dx = 0; dx < msub; dx++) 
                      { /* Compute the subsampling point {pxyz} on the unit sphere: */
                        r2_t qxy = (r2_t) {{ ix + (dx + 0.5)/msub,  iy + (dy + 0.5)/msub }}; /* Pt in image domain. */
                        r2_t pxy = neuromat_eeg_geom_disk_from_ellipse(&qxy, ctr, rad); /* Corresp point in unit-disk schematic. */
                        r3_t pxyz = neuromat_eeg_geom_3D_from_2D(&pxy); /* Corresp point on unit sphere. */
                        /* Find the nearest electrode to point {p}: */
                        int jmin = -1;
                        double d2min = +INF;
                        int je;
                        for (je = 0; je < ne; je++) 
                          { r3_t *pj = &(pos3D[je]);
                            double d2j = r3_dist_sqr(pj, &pxyz);
                            if (d2j < d2min) { d2min = d2j; jmin = je; }
                          }
                        assert(jmin >= 0);
                        /* Bump counter of nearest electrode: */
                        nnear[jmin]++;
                      }
                  }
                /* Assign basis values for this pixel according to hits, times {mxy}: */
                for (ie = 0; ie < ne; ie++) { bval[ie] = mxy*((double)(nnear[ie]))/(msub*msub); }
              }
            /* Now set pixels: */ 
            for (ie = 0; ie < ne; ie++) { float_image_set_sample(bas[ie], 0, ix, iy, (float)(bval[ie])); }
          }
      }
    free(nnear);
    free(bval);
  }

void neuromat_eeg_image_basis_lagrangian_matrix(int ne, r3_t pos3D[], neuromat_image_basis_elem_t *elem, double *L)
  {
    bool_t verbose = FALSE;

    /* Build colocation matrix {A} such that {A[ie,je]} is the value of element {ie} on point {pos3d[je]}: */
    double *A = rmxn_alloc(ne, ne);
    int ie, je;
    for (je = 0; je < ne; je++)
      { /* Fill column {je} of matrix {A}:  */
        r3_t *pj = &(pos3D[je]);
        for (ie = 0; ie < ne; ie++)
          { double Aij = elem(ie, pj);
            if (verbose && (Aij != 0)) { fprintf(stderr, "  A[%3d,%3d] = %+12.7f\n", ie, je, Aij); }
            A[ie*ne + je] = Aij;
          }
      }
               
    /* Invert matrix {A} to get the matrix {L} that maps from raw elem values to interp elem values: */
    assert(ne == ne);
    (void)rmxn_inv(ne, A, L);
    /* Check inversion: */
    { double *R = rmxn_alloc(ne, ne);
      rmxn_mul(ne, ne, ne, L, A, R);
      int nerr = 0;
      for (ie = 0; ie < ne; ie++)
        { for (je = 0; je < ne; je++)
            { double Rij = R[ie*ne + je];
              double Iij = (ie == je ? 1.0 : 0.0);
              double err = fabs(Rij - Iij);
              if (err > 1.0e-7) 
                { fprintf(stderr, "  inversion error R[%3d,%3d] = %+12.8f\n", ie, je, Rij);
                  nerr++;
                }
            }
        }
      demand(nerr == 0, "aborted");
      free(R);
    }
    free(A);
  }

    
