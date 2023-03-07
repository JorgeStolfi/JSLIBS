/* See rn_classif.h. */
/* Last edited on 2023-03-07 17:15:58 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <rn_classif.h>
#include <filefmt.h>
#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <frgb_ops.h>

void rn_classif_class_count(int NS, int class[], int NC, int num[])
  {
    int cl;
    for (cl = 0; cl <= NC; cl++) { num[cl] = 0; }
    int i;
    for (i = 0; i < NS; i++)
      { cl = class[i];
        assert((cl >= 0) && (cl <= NC));
        num[cl]++;
      }
  }

void rn_classif_class_count_print(FILE *wr, int NC, int num[])
  { 
    /* Compute total samples for percentage: */
    int cl;
    int NS = 0;
    for (cl = 0; cl <= NC; cl++) { NS += num[cl]; }
    /* Print table: */
    fprintf(wr, "samples per class:\n");
    for (cl = 0; cl <= NC; cl++)
      { if ((cl > 0) || (num[cl] > 0))
          { double pct = (100.0*num[cl])/NS;
            fprintf(wr, "  %5d %9d %6.2f %%\n", cl, num[cl], pct);
          }
      }
    fprintf(wr, "  %5s %9d\n", "TOTAL", NS);
  }

rn_classif_dataset_t *rn_classif_dataset_new(int NS, int NA)
  {
    rn_classif_dataset_t *D = notnull(malloc(sizeof(rn_classif_dataset_t)), "no mem");
    double **smp = notnull(malloc(NS*sizeof(double*)), "no mem");
    int i;  for (i = 0; i < NS; i++)  { smp[i] = NULL; }
    (*D) = (rn_classif_dataset_t) { .NS = NS, .NA = NA, .smp = smp };
    return D;
  }

void rn_classif_dataset_free(rn_classif_dataset_t *D)
  { if (D != NULL)
      { if (D->smp != NULL) 
          { int i;
            for (i = 0; i < D->NS; i++) { free(D->smp[i]); }
            free(D->smp); 
          }
        free(D);
      }
  }

#define FILE_TYPE "rn_classif_dataset_t"

#define FILE_VERSION "2010-05-24"

void rn_classif_dataset_write(FILE *wr, rn_classif_dataset_t *D, int cl, int class[])
  {  
    filefmt_write_header(wr, FILE_TYPE, FILE_VERSION);
    fprintf(wr, "attributes = %d\n", D->NA);
    fprintf(wr, "samples = %d\n", D->NS);
    int ndigs = digits(D->NS-1);
    int i;
    for (i = 0; i < D->NS; i++)
      { if ((class == NULL) || (class[i] == cl))
          { fprintf(wr, "%*d", ndigs, i);
            double *ati = D->smp[i];
            rn_classif_sample_print(wr, " ", D->NA, ati, " ", "\n");
          }
      }
    filefmt_write_footer(wr, FILE_TYPE);
    fflush(wr);
  }  

void rn_classif_sample_print(FILE *wr, char *pre, int NA, double p[], char *sep, char *suf)
  {  
    if (pre != NULL) { fputs(pre, wr); }
    int t;
    for (t = 0; t < NA; t++) 
      { if ((sep != NULL) && (t > 0)) { fputs(sep, wr); }
        fprintf(wr, "%+10.7f", p[t]);
      }
    if (suf != NULL) { fputs(suf, wr); }
  }  

void rn_classif_dataset_label(rn_classif_dataset_t *D, int NAL, int NCL, rn_classif_labeler_t *lab, int class[])
  {
    int NAD = D->NA;
    int NSD = D->NS;
    demand(NAL <= NAD, "not enough attributes in dataset");
    
    int i;
    for (i = 0; i < NSD; i++)
      { /* Get attributes {ati} and nominal class {old} of sample {i}: */
        double *ati = D->smp[i];
        int cl = lab(NAL, NCL, ati);
        assert((cl >= 0) && (cl <= NCL));
        class[i] = cl;
      }
  }

void rn_classif_compare(int NS, int cold[], int cnew[], int *nsP, int *nfP)
  {
    int ns = 0; /* Number of class matches. */
    int nf = 0; /* Number of class discrepancies. */
    int i;
    for (i = 0; i < NS; i++)
      { int old = cold[i];
        int new = cnew[i];
        if (old != new) { nf++; } else { ns++; }
      }
    (*nsP) = ns;
    (*nfP) = nf;
  }

void rn_classif_cross_matrix_build(int NS, int cold[], int NCold, int cnew[], int NCnew, int **nccP)
  {
    int cols = NCnew+1; /* Number of columns. */
    int rows = NCold+1; /* Number of rows. */
     
    int NTOT = rows*cols; /* Total elements in matrix. */
    int *ncc = notnull(malloc(NTOT*sizeof(int)), "no mem");
    
    int k; 
    for (k = 0; k < NTOT; k++) { ncc[k] = 0; }
    int i;
    for (i = 0; i < NS; i++)
      { int old = cold[i]; assert((old >= 0) && (old <= NCold));
        int new = cnew[i]; assert((new >= 0) && (new <= NCnew));
        ncc[old*cols + new]++;
      }
    (*nccP) = ncc;
  }
  
void rn_classif_cross_matrix_print(FILE *wr, int NCold, int NCnew, int *ncc, bool_t score)
  { 
    int cols = NCnew + 1; /* Number of columns in array. */
    
    int old, new;
    fprintf(wr, "row is old class, column is new class:\n");
    int ntot = 0; /* Sum of all matrix entries. */
    fprintf(wr, "%4s", "");
    for (new = 0; new <= NCnew; new++) { fprintf(wr, " %8d", new); }
    fprintf(wr, "\n");
    fprintf(wr, "%4s", "----");
    for (new = 0; new <= NCnew; new++) { fprintf(wr, " %8s", "--------"); }
    fprintf(wr, "\n");
    int ns = 0, nf = 0;
    for (old = 0; old <= NCold; old++)
      { fprintf(wr, "%4d", old);
        for (new = 0; new <= NCnew; new++) 
          { int count = ncc[old*cols + new];
            fprintf(wr, " %8d", count);
            ntot += count;
            if (old == new) { ns += count; } else { nf += count; }
          }
        fprintf(wr, "\n");
      }
    if (score)
      { double ps = (100.0*ns)/ntot;
        double pf = (100.0*nf)/ntot;
        fprintf(wr, "\n");
        fprintf(wr, "total %7d samples tested\n", ns + nf);
        fprintf(wr, "      %7d (%5.2f%%) successes\n", ns, ps);
        fprintf(wr, "      %7d (%5.2f%%) failures\n", nf, pf);
      }
  }

void rn_classif_find_nearest_neighbor_index(int NS, int NA, double p[], rn_classif_pi_dist_t *dist, int *iminP, double *dminP)
  {
    int imin = -1; double dmin = +INF;
    int i;
    for (i = 0; i < NS; i++)
      { double dpi = dist(NA, p, i);
        /* fprintf(stderr, "i = %d  dpi = %6.4f\n", i, dpi ); */
        if (dpi < dmin) { imin = i; dmin = dpi; }
      }
    (*iminP) = imin;
    (*dminP) = dmin;
  }

int rn_classif_find_nearest_in_dataset(rn_classif_dataset_t *M, double HM[], int NA, double p[], rn_classif_pq_dist_t *dist)
  {
    auto double pi_dist(int NA, double p[], int i);
      /* Distance from {p} to {M->smp[i]} with handicap {HM[i]} if given. */
      
    /* Proc body: */
    int imin; double dmin;
    rn_classif_find_nearest_neighbor_index(M->NS, NA, p, pi_dist, &imin, &dmin);
    return imin;
      
    /* Local procedure implementations: */
    
    double pi_dist(int NA, double p[], int i)
      { double *q = M->smp[i];
        double Hq = (HM == NULL ? 0 : HM[i]); 
        return fmax(dist(NA,p,q), Hq);
      }
  }

void rn_classif_nn_label_dataset
  ( rn_classif_dataset_t *M,    /* Model samples. */
    double HM[],                 /* {HM[i]} is the handicap of sample {M.smp[i]}. */
    int classM[],               /* {classM[i]} is the class of sample {M.smp[i]}. */
    rn_classif_dataset_t *C,    /* Sanples to be classified. */
    rn_classif_pq_dist_t *dist, /* Distance metric. */
    int classC[]                /* (OUT) {classC{j]} is the class assigned to {C.smp[j]}. */
  )
  {
    demand(C->NA <= M->NA, "incompatible domain dimensions");
    int j;
    for (j = 0; j < C->NS; j++)
      { double *p = C->smp[j];
        int i = rn_classif_find_nearest_in_dataset(M, HM, C->NA, p, dist);
        classC[j] = classM[i];
      }
  }

void rn_classif_nn_improve
  ( rn_classif_dataset_t *M,    /* (IN/OUT) Model samples. */
    double HM[],                 /* (IN/OUT) {HM[i]} is the handicap of sample {M.smp[i]}. */
    int classM[],               /* (IN/OUT) {classM[i]} is the given class of sample {M.smp[i]}. */
    rn_classif_dataset_t *R,    /* (IN/OUT) Training samples. */
    int classR[],               /* (IN/OUT) {classR[i]} is the given class of sample {R.smp[i]}. */
    int classX[],               /* (OUT) {classX[i]} is the class of {R.smp[i]} assigned by {M}. */
    rn_classif_pq_dist_t *dist, /* Sample distance function. */
    rn_classif_handicap_proc_t *hproc, /* Procedure that computes handicaps for {M}. */
    int maxIters,               /* Maximum iterations. */
    bool_t sameClass,           /* TRUE only swaps samples of the same class. */
    bool_t verbose              /* TRUE to show the progress at each iteration. */
  )
  {
    auto int evalM(void);
      /* Evaluates the current {M} by recomputing the handicaps {HM}, 
        using {M,HM,classM} to label {R}, and comparing the resulting
        labeling {classX} with the given labeling {classR}. R
        eturns the number of mismatches. */
        
    auto void swap_samples(int i, int j);
       /* Swaps {M.smp[i]} with {R.smp[j]} and {classM[i]} with {classM[j]}. */
        
    if (verbose) { fprintf(stderr, "classifying R with M ..."); }
    int nf; /* Number of errors in classifying {R} with {M} */
    nf = evalM();
    int candM = 0;  /* Next sample of {M} that is candidate to be moved out. */
    int candR = 0;  /* Next sample of {R} that is a candidate to be moved in. */
    
    int try;
    for (try = 0; (try < maxIters) && (nf > 0); try++)
      { if (verbose) { fputc('\n', stderr); fprintf(stderr, "iteration %5d ", try); }
        /* Choose a misclassified point {candR} in {R}: */
        while (classR[candR] == classX[candR]) { candR = (candR + 1) % R->NS; }
        /* Pick a candidate for removal {candM} in {M}: */
        /* !!! Should pick {candM} in a smarter way !!! */
        if (sameClass)
          { /* Tries to find a point in {M} starting at {candM} with class {classR[candR]}. */
            /* If there If there is no such point, returns candM. */
            int nt = M->NS; 
            while ((nt > 0) && (classM[candM] != classR[candR]))
              { candM = (candM + 1) % M->NS; nt--; }
          }
        /* Swap sample {candM} of {M} by {candR}: */
        /* !!! Should save the {HM} vector and restore it when unswapping. !!! */
        if (verbose) { fprintf(stderr, "swapping M[%5d] and R[%5d] ...  ", candM, candR ); }
        swap_samples(candM, candR);
        int nfold = nf;
        nf = evalM();
        if (nf > nfold)
          { /* Replacement was a bad idea, restore the removed representative: */
            if (verbose) { fprintf(stderr, "unswapping M[%5d] and R[%5d] ...", candM, candR ); }
            swap_samples(candM, candR);
            nf = evalM();
            affirm(nf == nfold, "why?");
          }
        /* Round-robin of replacement candidates: */
        candM = (candM + 1) % M->NS;
      }
    
    /* Local procedure implementations */
    
    int evalM(void)
      { if (hproc == NULL)
          { int i; for (i = 0; i < M->NS; i++) { HM[i] = 0; } } 
        else if (HM != NULL) 
          { hproc(M, classM, HM); }
        rn_classif_nn_label_dataset(M, HM, classM, R, dist, classX);
        int ns, nf;
        rn_classif_compare(R->NS, classR, classX, &ns, &nf);
        if (verbose) { fprintf(stderr, " nf = %6d (%5.1f)\n", nf, (100.0*nf)/R->NS ); }
        return nf;
      }

    void swap_samples(int i, int j)
      {
        { double *tsmp = M->smp[i]; M->smp[i] = R->smp[j]; R->smp[j] = tsmp; }
        { int tcl = classM[i]; classM[i] = classR[j]; classR[j] = tcl; }
      }
  }

frgb_t *rn_classif_pick_class_colors(int NC)
  {
    frgb_t *cmap = notnull(malloc((NC+1)*sizeof(frgb_t)), "no mem"); 
    cmap[0] = (frgb_t){{ 1,1,1 }}; /* The `none' class is black. */
    double YMin = 0.299; /* Min luminance. */
    double YMax = 0.800; /* Max luminance. */
    int cl;
    for (cl = 1; cl <= NC; cl++)
      { /* Map the class from {1..NC} to a real {r} in {[0_1]}: */
        double r = (cl - 1.0)/NC;
        /* Compute the luminance {Y} from {r}: */
        double Y = exp(log(YMin)*(1-r) + log(YMax)*r);
        /* Pick hue so that it is fairly varied when {NC} is large: */
        double H = (1 + NC/8)*r;
        /* Use the maximum saturation possible in the RGB cube: */
        double T = 1;
        frgb_t val = (frgb_t){{ (float)H, (float)T, (float)Y }};
        frgb_from_HTY(&val);
        /* fprintf(stderr, "cmap[%d] = ( %6.4f %6.4f %6.4f )\n", cl, val.c[0], val.c[1], val.c[2]); */
        cmap[cl] = val;
      }
    return cmap;
  }
