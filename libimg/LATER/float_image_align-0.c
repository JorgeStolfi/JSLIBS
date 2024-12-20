/* See {float_image_align.h}. */
/* Last edited on 2024-12-04 23:19:45 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <float_image.h>
#include <rn.h>

#include <float_image_align.h>

#define INF MAXDOUBLE

alignment_t *fial_alloc_alignment_set(int32_t ncan);

void fial_generate_all_alignments(int32_t nims, int32_t dmax, double step, int32_t nal, alignment_t *al);
  /* Generate all relative alignments for {nims} images whose entries 
    range in {[-dmax .. +dmax]} times {step}.
    Expects {nal} to be the exact number of such alignments. 
    Returns the alignments in {al[0..*nal-1]}. */

image_rec_t *fial_reduce_images(image_rec_t *iml);
  /* Reduces all imagaes in {iml} by a linear factor of 2;
    returns a list of the results.  Also computes the coordinates
    of the centers of the reduced images. */

int32_t fial_num_rel_shifts(int32_t nims, int32_t dmax, int32_t sum);
  /* Computes the number of integer vectors in {[-dmax..+dmax]^nims}
    which add to {sum}. */
    
void fial_recenter_displacements(int32_t nims, double *dcol, double *drow);
  /* Shifts all displacements so that their sum is zero.  */

void fial_adjust_alignment
  ( int32_t nims, 
    image_rec_t *iml, 
    double dmax, 
    double epsilon, 
    double radius, 
    alignment_t *a
  );
  /* Adjusts the shifts {(al->dcol[i],al->drow[i])} for each image {i}, so
    as to minimize the average pixel variance {al->tvar}. The maximum
    adjustment is {dmax} and the nominal precision required is {epsilon}, in
    each direction. The pixel variances are weighted with a Gaussian of
    characteristic radius {radius}. A displacement of zero means to
    align the center of the image at the origin. The displacements are
    adjusted so as to add to zero.  */

void fial_average_images(int32_t nims, image_rec_t *iml, int32_t iex, Image *avg, double *dcol, double *drow);
  /* Computes average of all images, excluding image {iex}, shifting each one 
    so that its center displaced by {dcol[j],drow[j]} matches the center
    of {avg}.  */
    
interval_t fial_compute_mismatch
  ( Image *tim, 
    double dcol, 
    double drow, 
    Image *sim, 
    double scale, 
    int32_t nwt, 
    double *gwt
  );
  /* Computes the mismatch between the image {tim}, shifted by {(dcol,drow)},
    against the square image {sim}.  The vector {wt} is a unidimensional weight
    distribution.  Assumes that the length of {wt} is the same as the number
    of rows and columns of {sim}.  

    The mismatch is the sum of the squared pixel differences, each multiplied
    by the {scale} factor and the bidimensional weight distribution {wt[col]*wt[row]}.  */

void fial_recompute_pixel_variance(int32_t nims, image_rec_t *iml, alignment_t *al, Image *avg, int32_t nwt, double *wt);
  /* Recomputes the mismatches {al->mism[i]} and the total pixel variance {al->tvar},
    from scratch. */
    
int32_t fial_max_chns(image_rec_t *iml);
  /* Returns the maximum {chns} among all images in {iml}.  */
  
void fial_sort_and_prune_alignments(int32_t nims, int32_t *ncanP, alignment_t *can, double epsilon, int32_t max);

double *fial_gaussian_distr(int32_t width, double sigma);
  /* Gaussian distribution {g[i] = A*exp(-((i+0.5)-m)^2/sigma^2)}
    where {m = width/2} and {A} is such that sum is 1.  */
    
void fial_debug_displacement(char *label, int32_t i, double col, double row, interval_t s, char *tail);

double fial_compare_cands(int32_t nims, alignment_t *al, alignment_t *bl);
  /* Returns the maximum difference between the shifts of {al} and {bl}.  */
  
int32_t fial_compare_mismatches(interval_t *a, interval_t *b); 
  /* Returns +1 if {a} looks smaller than {b}, -1 if {a} looks larger than {b}, 
    and 0 if they look equal or incomparable. */

void fial_optimal_shifts
  ( int32_t nims, 
    float_image_t *iml, 
    double dxmax, 
    double dymax, 
    double epsilon.
    double rx, 
    double ry, 
    double *dx, 
    double *dy
  );

void fial_multiscale_optimal_shifts
  ( int32_t nims, 
    float_image_t *iml, 
    double dxmax, 
    double dymax, 
    double epsilon, 
    double rx, 
    double ry, 
    int32_t mult, 
    alignment_t *opt, 
    int32_t *noptP
  )
  {
    double guess_step = 0.5;   /* Initial guesses are spaced every half-pixel: */
    int32_t idmax = (int32_t)(floor(dmax/guess_step));
    double allshifts = (double)fial_num_rel_shifts(nims, idmax, 0); /* Watch oflow... */
    double allsols = allshifts*allshifts; /* Total possible alignments. */
    alignment_t *can = NULL;
    int32_t ncan;
    int32_t j, k;
    fprintf(stderr, "enter compute_multiscale_optimal_shifts\n");
    fprintf(stderr, "  allsols = %g\n", allsols);
    if (mult >= allsols)
      { /* The caller asked for all possible solutions. */
        int32_t nall = (int32_t)(floor(allsols+0.5));
        fprintf(stderr, "  computing all %d solutions...\n", nall);
        can = fial_alloc_alignment_set(nall);
        fial_generate_all_alignments(nims, idmax, guess_step, nall, can);
        ncan = nall;
      }
    else 
      { /* Generate {4*mult} solutions for coarser scale, then refine and select. */
        image_rec_t *rdl = fial_reduce_images(iml);
        /* Compute 4*mult solutions for coarse images */
        can = fial_alloc_alignment_set(4*mult);
        fprintf(stderr, "  recursion to smaller scale...\n");
        compute_multiscale_optimal_shifts(nims, rdl, 4*mult, dmax/2, guess_step, radius/2, can, &ncan);
        for (k = 0; k < ncan; k++) 
          { /* Map displacements to current scale: */
            alignment_t *al = &(can[k]);
            image_rec_t *rdli, *imli;
            int32_t i;
            for (i=0, imli=iml, rdli=rdl; i<nims; i++, imli=imli->next, rdli=rdli->next)
              { al->dcol[i] = 2*al->dcol[i] + (rdli->ccol - imli->ccol);
                al->drow[i] = 2*al->drow[i] + (rdli->crow - imli->crow);
                al->mism[i] = (interval_t){0.00, 0.25};
              }
          }
      }

    /* Refine each sol by moving at most {epsilon} from current pos: */
    fprintf(stderr, "  local refinement...\n");
    for (k = 0; k < ncan; k++) 
      { fial_adjust_alignment(nims, iml, guess_step, epsilon, radius, &(can[k])); }

    /* Select {mult} best solutions: */
    fprintf(stderr, "  sorting and pruning...\n");
    fial_sort_and_prune_alignments(nims, &ncan, can, epsilon, mult);
    (*noptP) = (ncan < mult ? ncan : mult);
    for (j = 0; j < (*noptP); j++) { opt[j] = can[j]; }
    for (j = (*noptP); j < ncan; j++) 
      { opt[j].dcol = NULL; opt[j].drow = NULL; 
        free(can[j].dcol); free(can[j].drow); free(can[j].mism);
      }
    free(can);
    fprintf(stderr, "exit compute_multiscale_optimal_shifts\n");
  }

int32_t fial_num_rel_shifts(int32_t nims, int32_t dmax, int32_t sum)
  {
    fprintf(stderr, "(%d)", nims);
    if (abs(sum) > nims*dmax)
      { return 0; }
    else if (nims <= 1) 
      { return 1; }
    else if (nims == 2) 
      { return 2*dmax + 1; }
    else if (nims == 3) 
      { return 3*dmax*(dmax + 1) + 1; }
    else if (nims == 4) 
      { return 8*dmax*(dmax*(2*dmax + 3) + 1)/3 + 1; }
    else if (nims == 5) 
      { /* There must be a better way... */
        int32_t i, s = 0;
        for (i=-dmax; i <= dmax; i++) { s += fial_num_rel_shifts(nims-1,dmax,sum-i); }
        return s;
      }
    else
      { pnm_error("too many images");
        return 0;
      }
  }

image_rec_t *fial_reduce_images(image_rec_t *iml)
  {
    fprintf(stderr, "r");
    if (iml == NULL)
      { return NULL; }
    else
      { Image *im = iml->im;
        double rcol = iml->ccol;
        double rrow = iml->crow;
        Image *rd = reduce_image(im, &rcol, &rrow);
        return cons_image(rd, rcol, rrow, fial_reduce_images(iml->next));
      }
  }

void fial_sort_and_prune_alignments(int32_t nims, int32_t *ncanP, alignment_t *can, double epsilon, int32_t max)
  {
    /* !!! must eliminate duplicates !!! */
    int32_t ncan = *ncanP;
    int32_t nok, i;
    /* Select the best {max} alignments: */
    nok = 1;
    for (i = 1; i < ncan; i++)
      { int32_t j;
        alignment_t t = can[i];
        j = nok-1;
        while ((j >= 0) && (fial_compare_cands(nims, &(can[j]), &t) > epsilon/2)) { j--; }
        if (j < 0)
          { j = nok;
            while ((j > 0) && (fial_compare_mismatches(&(can[j-1].tvar), &(t.tvar)) > 0))
              { if (j < max) { can[j] = can[j-1]; }
                j--;
              }
            can[j] = t; if (nok < max) { nok++; }
          }
      }
    (*ncanP) = nok;
  }

double fial_compare_cands(int32_t nims, alignment_t *al, alignment_t *bl)
  {
    int32_t i;
    double dmax = 0;
    for (i = 0; i < nims; i++)
      { double dc = fabs(al->dcol[i] - bl->dcol[i]);
        double dr = fabs(al->drow[i] - bl->drow[i]);
        if (dc > dmax) { dmax = dc; }
        if (dr > dmax) { dmax = dr; }
      }
    return dmax; 
  }

int32_t fial_compare_mismatches(interval_t *a, interval_t *b)
  {
    double amin = 0.75*a->lo + 0.25*a->hi;
    double bmin = 0.75*b->lo + 0.25*b->hi;
    if (amin < bmin)
      { return -1; }
    else if (amin > bmin)
      { return +1; }
    else
      { return 0; }
  }

alignment_t *fial_alloc_alignment_set(int32_t ncan)
  {
    int32_t i;
    alignment_t *als = (alignment_t *)malloc(ncan*sizeof(alignment_t));
    if (als == NULL) { pnm_error("out of memory for alignments"); }
    for (i = 0; i < ncan; i++) 
      { alignment_t *alsi = &(als[i]);
        alsi->dcol = NULL; alsi->drow = NULL; alsi->mism = NULL; 
        alsi->tvar = (interval_t){0.00, 0.25};
      }
    return als;
  }

alignment_t alloc_alignment(int32_t nims)
  {
    alignment_t al;
    int32_t i;
    al.dcol = rn_alloc(nims);
    al.drow = rn_alloc(nims);
    al.mism = (interval_t *)malloc(nims*sizeof(interval_t));
    if ((al.dcol == NULL) || (al.drow == NULL) || (al.mism == NULL))
      { pnm_error("out of memory for displacement vectors"); }
    for (i=0; i < nims; i++) 
      { al.dcol[i] = 0; al.drow[i] = 0; 
        al.mism[i] = (interval_t){0.00, 0.25};
      }
    al.tvar = (interval_t){0.00, 0.25};
    return al;
  }

int32_t imin(int32_t x, int32_t y)
  { 
    return (x < y ? x : y);
  }

int32_t imax(int32_t x, int32_t y)
  { 
    return (x > y ? x : y);
  }

void fial_generate_all_alignments(int32_t nims, int32_t dmax, double step, int32_t nal, alignment_t *al)
  {
    int32_t sc, sr;
    int32_t *dc = (int32_t *)malloc(nims*sizeof(int32_t));
    int32_t *dr = (int32_t *)malloc(nims*sizeof(int32_t));
    int32_t k;
    int32_t carry; /* boolean */
    int32_t nfound = 0;

    k = nims; sc = 0; sr = 0; carry = 0;
    while(1)
      { 
        /* Here {sc} is the sum of {dc[k..nims-1]}, ditto for {sr,dr}. */
        /* Generate the smallest solution with current {dc,dr[k..nfound-1]}: */
        int32_t j;
        while (k > 0)
          { /* Move down and set digit to minimum value: */
            k--;
            dc[k] = imax(-dmax, -(k*dmax + sc));
            dr[k] = imax(-dmax, -(k*dmax + sr));
            sc += dc[k]; sr += dr[k];
          }     

        /* Store this solution: */
        if (nfound >= nal) 
          { pnm_error("num of alignments underestimated"); }
        if ((sc != 0) || (sr != 0)) 
          { pnm_error("displacements don't add to zero"); }
        { alignment_t a = alloc_alignment(nims);
          fprintf(stderr, "\n");
          for (j = 0; j < nims; j++)
            { a.dcol[j] = step*dc[j];
              a.drow[j] = step*dr[j];
              a.mism[j] = (interval_t){0.00, 0.25};
              fial_debug_displacement("  ", j, a.dcol[j], a.drow[j], a.mism[j], "\n");
            }
          a.tvar = (interval_t){0.00, 0.25};
          al[nfound] = a; nfound++;
        }

        /* Find first position {k} that can be incremented, and increment it: */
        carry = 1;
        while (carry)
          { /* Now {sc,sr} is sum of {dc,dr[k..nims-1]} */
            if (k >= nims)
              { /* Exhausted all possibilities, stop: */
                if (nfound != nal) { pnm_error("number of alignments doesn't check"); }
                free(dc); free(dr);
                return;
              }
            sc -= dc[k]; sr -= dr[k]; 
            if (dr[k] < imin(dmax, k*dmax - sr))
              { /* Bump {dr[k]}, reset carry: */
                dr[k]++;
                carry = 0; 
              }
            else if (dc[k] < imin(dmax, k*dmax - sc))
              { /* Bump {dc[k]}, reset {dr[k]}, reset carry: */
                dc[k]++; dr[k] = imax(-dmax, -(k*dmax+sr)); 
                carry = 0; 
              }
            else 
              { /* This position can't be incremented, carry to the next one: */
                k++; 
              }
            /* Now {sc,sr} is sum of {dc,dr[k+1-carry..nims-1]} */
          }
        sc += dc[k]; sr += dr[k];
      }
  }

void fial_adjust_alignment
  ( int32_t nims, 
    image_rec_t *iml, 
    double dmax, 
    double epsilon, 
    double radius, 
    alignment_t *a
  )
  {
    if (nims >= 2)
      { 
        int32_t arad = ceil(3*radius); /* Weights are negligible beyond this. */
        int32_t nwt = 2*arad+1;
        Image *avg = allocate_image(nwt, nwt, fial_max_chns(iml));
        double *gwt = fial_gaussian_distr(nwt, radius);
        double s1 = ((double)(nims-1))/((double)nims);
        double corr = s1*s1; /* Variance correction factor. */
        image_rec_t *impi;
        double step;
        int32_t i;

        /* Loop until convergence: */
        i = 0; impi = iml;
        step = 0.5*(dmax + epsilon);
        while(1)
          { int32_t krowi, kcoli;
            interval_t prev_best; 
            /* Adjust displacements so that they add to zero: */
            fial_recenter_displacements(nims, al->dcol, al->drow);
            /* Compute the average of all images excluding image {i}: */
            pnm_message("averaging images except %d...", i);
            fial_average_images(nims, iml, i, avg, al->dcol, al->drow);
            if (avg->cols < 20) { debug_image("avg", "", avg); }
            /* Compare image {i} against average, in all possible displacements: */
            pnm_message("computing current mismatch of image %d...", i);
            al->mism[i] = fial_compute_mismatch(impi->im, al->dcol[i], al->drow[i], avg, corr, nwt, gwt);
            fial_debug_displacement("INIT", i, al->dcol[i], al->drow[i], al->mism[i], "\n");
            prev_best = al->mism[i];
            for (krowi = -1; krowi <= +1; krowi++)
              { double drowi = al->drow[i] + step*krowi;
                for (kcoli = -1; kcoli <= +1; kcoli++)
                  { double dcoli = al->dcol[i] + step*kcoli;
                    interval_t s = fial_compute_mismatch(impi->im, dcoli, drowi, avg, corr, nwt, gwt);
                    if (fial_compare_mismatches(&s, &(al->mism[i])) < 0) 
                      { fial_debug_displacement("good", i, dcoli, drowi, s, "\n");
                        al->mism[i] = s; al->dcol[i] = dcoli; al->drow[i] = drowi;
                      }
                  }
              }
            if (fial_compare_mismatches(&(al->mism[i]), &prev_best) < 0) 
              { fial_debug_displacement("BEST", i, al->dcol[i], al->drow[i], al->mism[i], "\n");
              }
            impi = impi->next; i++;
            if (i >= nims) 
              { i = 0; impi = iml; 
                if (step <= 1.00001*epsilon) 
                  { fial_recompute_pixel_variance(nims, iml, al, avg, nwt, gwt); return; } 
                step /= 2;
              }
          }
      }
  }

void fial_recompute_pixel_variance(int32_t nims, image_rec_t *iml, alignment_t *al, Image *avg, int32_t nwt, double *wt)
  {
    int32_t i = 0;
    image_rec_t *impi = iml;
    double qt = 1.0/((double)nims);
    al->tvar = (interval_t){0, 0};
    fial_average_images(nims, iml, nims, avg, al->dcol, al->drow);
    for (i = 0, impi=iml; i < nims; i++, impi = impi->next)
      { interval_t s = fial_compute_mismatch(impi->im, al->dcol[i], al->drow[i], avg, 1.0, nwt, wt);
        al->mism[i] = s;
        al->tvar.lo += qt*s.lo;
        al->tvar.hi += qt*s.hi;
      }
  }

void fial_recenter_displacements(int32_t nims, double *dcol, double *drow)
  {
    double tcol = 0, trow = 0;
    int32_t i;
    for (i = 0; i < nims; i++) { tcol += dcol[i]; trow += drow[i]; }
    tcol /= (double)nims;
    trow /= (double)nims;
    for (i = 0; i < nims; i++) { dcol[i] -= tcol; drow[i] -= trow; }
  }

double *fial_gaussian_distr(int32_t width, double sigma)
  { double m = ((double)width)/2.0;
    double *g = rn_alloc(width);
    double s = 0;
    int32_t i;
    if (g == NULL) { pnm_error("out of memory for gaussian weight"); }
    for (i=0; i < width; i++)
      { double x = (i + 0.5 - m)/sigma;
        double y = exp(-x*x); 
        g[i] = y; s += y;
      }
    for (i=0; i < width; i++)
      { g[i] /= s; }
    return g;
  }

void write_optimal_positions(int32_t nims, image_rec_t *iml, double *dcol, double *drow)
  {
    int32_t i;
    image_rec_t *impi;
    for (i = 0, impi = iml; i < nims; i++, impi = impi->next)
      { double tcol = dcol[i] + impi->ccol;
        double trow = drow[i] + impi->crow;
        printf("%03d %8.3f %8.3f\n", i, tcol, trow);
      }
    fflush(stdout);
  }

interval_t fial_compute_mismatch
  ( Image *tim, 
    double dcol, 
    double drow, 
    Image *sim, 
    double scale, 
    int32_t nwt, 
    double *gwt
  )
  {
    int32_t srow, scol;
    double wt2;
    interval_t d2;
    interval_t td2 = (interval_t){0.0, 0.0}; /* Total mismatch. */
    if ((sim->rows != nwt) || (sim->cols != nwt)) 
      { pnm_error("ref image has wrong size"); } 
    for (srow = 0; srow < sim->rows; srow++)
      { for (scol = 0; scol < sim->cols; scol++)
          { int32_t sindex = (srow * sim->cols + scol)*sim->chns;
            interval_t *sv = &(sim->el[sindex]);
            interval_t tv[MAX_CHANNELS];
            double tcol = scol + ((double)(tim->cols - sim->cols))/2.0 - dcol;
            double trow = srow + ((double)(tim->rows - sim->rows))/2.0 - drow;

            /* debug = ((scol % 56) == 23) & ((srow % 56) == 23); */
            eval_at_point(tim, tcol, trow, &(tv[0]));
            if (debug) 
              { debug_floatized_pixel("imi", tcol, trow, sim->chns, &(tv[0]), "\n"); } 
            if (tim->chns != sim->chns)
              { convert_pixel(tim->chns, &(tv[0]), sim->chns); }
            subtract_pixel(sim->chns, &(sv[0]), &(tv[0]));
            if (debug) 
              { debug_floatized_pixel("dif", scol, srow, sim->chns, &(tv[0]), "\n\n"); } 
            d2 = pixel_norm(sim->chns, &(tv[0]));
            wt2 = scale*wt[srow]*wt[scol];
            td2.lo += wt2*d2.lo;
            td2.hi += wt2*d2.hi;
          }
      }
    return td2;
  }

void fial_average_images(int32_t nims, image_rec_t *iml, int32_t iex, Image *avg, double *dcol, double *drow)
  {
    int32_t row, col;
    int32_t navg = ((iex >= 0) && (iex < nims) ? nims-1 : nims);
    double qt = 1.0/((double)navg);
    double tqt;
    for (row = 0; row < avg->rows; row++)
      { for (col = 0; col < avg->cols; col++)
          { int32_t aindex = (row * avg->cols + col)*avg->chns;
            interval_t *fva = &(avg->el[aindex]);
            int32_t j, k;
            image_rec_t *impj;
            for (k = 0; k < avg->chns; k++) { fva[k] = (interval_t){ 0.0, 0.0 }; }
            for (j = 0, impj = iml; j < nims; j++, impj = impj->next)
              { if (j != iex) 
                  { Image *imj = impj->im;
                    interval_t fvj[MAX_CHANNELS];
                    double colj = col + ((double)(imj->cols - avg->cols))/2.0 - dcol[j];
                    double rowj = row + ((double)(imj->rows - avg->rows))/2.0 - drow[j];
                    eval_at_point(imj, colj, rowj, &(fvj[0]));
                    if (imj->chns != avg->chns)
                      { convert_pixel(imj->chns, &(fvj[0]), avg->chns); }
                    accum_pixel(avg->chns, &(fvj[0]), qt, fva, &tqt); 
                  }
              }
          }
      }
  }

int32_t fial_max_chns(image_rec_t *iml)
  {
    int32_t chns = 0;
    while (iml != NULL)
      { Image *im = iml->im;
        if (im->chns > chns) { chns = im->chns; }
        iml = iml->next;
      }
    return chns;
  }

void fial_debug_displacement(char *label, int32_t i, double col, double row, interval_t s, char *tail)
  {
    fprintf(stderr, "%s im[%d] d=[%8.3f,%8.3f] s = [%15.6f _ %15.6f]", label, i, row, col, s.lo, s.hi);
    fprintf(stderr, "%s", tail);
  }

