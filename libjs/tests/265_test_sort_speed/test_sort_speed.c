#define PROG_NAME "test_sort_speed"
#define PROG_DESC "tests speed of quick-sort, heap-sort, binsertion-sort and merge-sort"
#define PROG_VERS "1.1"

/* Last edited on 2013-10-25 18:48:00 by stolfilocal */

#define test_sort_speed_COPYRIGHT \
  "Copyright © 2004  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    ??? \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

#define PROG_INFO \
  "???\n" \
  "AUTHORS\n" \
  "  Created aug/2004 by Jorge Stolfi, IC-UNICAMP"

/* We must set _GNU_SOURCE to get {asprintf}. */
/* Why can't we just include "gnuio.h"??? */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <intsort.h>
#include <intsort_extra.h>
#include <affirm.h>
#include <jsfile.h>
#include <jstime.h>

static int ncmp;
static long int *data; /* Random values to be compared. */
static int ndata; /* Size of {data} vector. */

/* Max number of runs for each algorithm and each array size: */
#define MAXRUNS 50

/* Maximum array size to try: */
/* #define MAXSIZE (4*256*256) */
/* #define MAXSIZE (128*256) */
#define MAXSIZE (4*256)

typedef int int_cmp_t(int a, int b); 
  /* A signed comparison predicate for integers (indices, etc.). */

typedef void tss_sorter_t(int *h, int n, int_cmp_t cmp, int sgn);
  /* A procedure that sorts {h} so that  {sgn*cmp(h[i-1],h[i]) <= 0}
    for all {i} */
  
typedef double tss_timer_t(int n);
  /* A procedure that estimates how many comparisons (or equivalent
    ops) an algorithm will take for an {n}-element array. 
    
    The actual value returned is not important; the test procs use
    only the ratios of the estimated times among various algorithms,
    for each given {n}. */
  
typedef struct tss_alg_t 
  { tss_sorter_t *srt;  /* The sorting procedure. */
    char *name;         /* The procedure's name. */
    char *descr;        /* A short description. */
    int thr;            /* Max {n} included in the base case. */
    tss_timer_t *tmr;   /* Estimates the number of operations. */
  } tss_alg_t;

typedef enum {TP_NCMP = 0, TP_TIME = 1} stype_t;
  /* Type of statistics to print (comparisons or running times). */
  
typedef struct tss_stats_t 
  { int ncmp;    /* Number of comparisons */
    double time; /* Time in microseconds. */
  } tss_stats_t;

/* PROTOS */

void tss_all_tests(char *outname);
  /* Tests all algorithms for various array sizes, prints
    statistics to the file "{outname}-{tpname}.tex" where 
    {tpname} is "ncmp" or "time". */

tss_stats_t tss_test_sorter(int *h, int n, tss_sorter_t srt, int_cmp_t cmp, int sgn);
  /* Tests sorter {srt} on {h} with comparator {sgn*cmp}; checks
    order of result, returns statistics. */
    
void tss_print_stats(FILE *wr, char *sname, int n, tss_stats_t *st, int ntests, stype_t tp);
  /* Writes to {wr} the summarized statistics of type {tp} for sorter {sname}
    acting on arrays of {n} elements, given the statitistics {st[0..ntests-1]} 
    obtained in {ntests} test runs. */

char *tss_protect(char *name);
  /* Returns a copy of {name} wrapped in '{}'s
    with every '_' changed into '\_'. */

void tss_fill_algorithms_table(int *nalgsP, tss_alg_t *alg[]);
  /* Stores into {alg[0..nalgs-1]} a set of algorithms to test.
    Also sets {*nalgsP} the number {nalgs} of algorithms 
    defined.  Assumes that {alg} has at least {MAXALGS} entries. */

void tss_write_tex_algorithms_table(char *outname, int nalgs, tss_alg_t *alg[]);
  /* Writes to "{outname}-algs.tex" a TeX table with procedure name,
    algorithm description, and recursion threshold {SMALL}. */

void tss_begin_tex_performance_table(FILE *wr, char *outname, char *tpname, int n);
  /* Writes into {wr} the preamble of a TeX table with performance data for various
    algorithms on tables with {n} entries. */
    
void tss_finish_tex_performance_table(FILE *wr, char *tpname, int n);
  /* The file {wr} must have been opened with
    {tss_begin_tex_performance_table}. Writes into {wr} the
    postamble of the table */

int tss_cmp(int a, int b);
  /* Compares {data[a]} with {data[b]}. */

tss_alg_t *tss_new_alg(tss_sorter_t *srt, char *name, char *descr, int thr, tss_timer_t *tmr);
  /* Creates a new algorithm data record. */
  
double tss_n_logn(int n);
double tss_n2_2(int n);
double tss_n2_4(int n);
  /* Useful time estimators. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { tss_all_tests(PROG_NAME);
    return 0;
  }

tss_alg_t *tss_new_alg(tss_sorter_t *srt, char *name, char *descr, int thr, tss_timer_t *tmr)
  { tss_alg_t *pa = (tss_alg_t *)notnull(malloc(sizeof(tss_alg_t)), "no mem");
    pa->srt = srt; 
    pa->name = name;
    pa->descr = descr;
    pa->thr = thr; 
    pa->tmr = tmr;
    return pa;
  }

/* Max number of algorithms to test: */
#define MAXALGS 10

void tss_all_tests(char *outname)
  {
    char *tpname[2];
    tpname[TP_NCMP] = "ncmp"; 
    tpname[TP_TIME] = "time";
    stype_t tp;
    
    int nalgs = 0;
    tss_alg_t *alg[MAXALGS];
    
    /* Assemble a table of algorithms and parameters: */
    tss_fill_algorithms_table(&nalgs, alg);

    /* Write table describing algorithms: */
    tss_write_tex_algorithms_table(outname, nalgs, alg);

    /* Open separate files for comparison counts and times: */
    fprintf(stderr, "  %-45s %6s  %8s .. %8s  %10s ± %10s (%4s)\n\n",
      "Algorithm", "n", "min", "max", "avg", "std", "runs"); 

    /* Open performance data files: */
    FILE *wr[2];
    for (tp = TP_NCMP; tp <= TP_TIME; tp++)
      { char *fname = NULL;
        asprintf(&fname, "%s-%s.tex", outname, tpname[tp]);
        wr[tp] = open_write(fname, TRUE);
        free(fname);
      }

    int n = 0;
    while (n <= MAXSIZE) 
      { /* Allocate data array */
        long int dt[n]; data = dt; ndata = n;
        int h[n];

        /* Performance statistics for the various runs of each algorithm: */
        tss_stats_t st[nalgs*MAXRUNS];
        
        /* Define the max number of runs for each algorithm: */
        int max_runs[nalgs];
        int ialg;
        for (ialg = 0; ialg < nalgs; ialg++)
          { double relcost = alg[ialg]->tmr(n)/tss_n_logn(n);
            max_runs[ialg] = (int)(MAXRUNS/relcost);
          }

        int sgn = +1;

        fprintf(stderr, "\n");
        int it;
        for (it = 0; it < MAXRUNS; it++)
          { fprintf(stderr, "=");
            /* Generate "random" numbers and test the two sorters on them: */
            int i;
            srandom(it+1);
            for (i = 0; i < n; i++) { dt[i] = random(); }
            
            for (ialg = 0; ialg < nalgs; ialg++)
              { if (it < max_runs[ialg]) 
                  { tss_alg_t *pa = alg[ialg];
                    tss_stats_t *pst = &(st[ialg*MAXRUNS + it]);
                   (*pst) = tss_test_sorter(h, n, pa->srt, tss_cmp, sgn);
                  }
              }
          }
        fprintf(stderr, "\n");
          
        for (tp = TP_NCMP; tp <= TP_TIME; tp++)
          { tss_begin_tex_performance_table(wr[tp], outname, tpname[tp], n);
            for (ialg = 0; ialg < nalgs; ialg++)
              { tss_alg_t *pa = alg[ialg];
                tss_stats_t *pst = &(st[ialg*MAXRUNS]);
                tss_print_stats(wr[tp], pa->descr, n, pst, MAXRUNS, tp);
              }
            tss_finish_tex_performance_table(wr[tp], tpname[tp], n);
          }


        n = (n == 0 ? 1 : 2*n); 
      }
    
    for (tp = TP_NCMP; tp <= TP_TIME; tp++) { fclose(wr[tp]); }

  }
    
void tss_fill_algorithms_table(int *nalgsP, tss_alg_t *alg[])
  {
    int nalgs = 0;
    
    #define defalg(DESCR,SRT,THR,TMR) \
      do { \
        assert(nalgs < MAXALGS); \
        alg[nalgs] = tss_new_alg(&SRT, #SRT, DESCR, THR, &TMR); \
        nalgs++; \
      } while (0)
    
    defalg(
      "Heapsort (classic)",
      isrt_heapsort_classic, 
      0,
      tss_n_logn   
    );
    defalg(
      "Heapsort (vacsink)",
      isrt_heapsort_vacsink,     
      0,
      tss_n_logn   
    );
    defalg(
      "Heapsort (libjs)",
      isrt_heapsort,       
      0,
      tss_n_logn   
    );
    defalg(
      "Quick-sort (middle elem)",
      isrt_quicksort_middle,
      isrt_quicksort_middle_SMALL,
      tss_n_logn   
    );
    defalg(
      "Quick-sort (median-of-three)",
      isrt_quicksort_median3,    
      isrt_quicksort_median3_SMALL,
      tss_n_logn   
    );
    defalg(
      "Mergesort (in-place, piv)",
      isrt_mergesort_pivot,    
      isrt_mergesort_pivot_SMALL,
      tss_n_logn   
    );
    defalg(
      "Mergesort (in-place, sym)",
      isrt_mergesort_symsplit,  
      isrt_mergesort_symsplit_SMALL,
      tss_n_logn   
    );
    defalg(
      "Mergesort (libjs)",
      isrt_mergesort,      
      0,
      tss_n_logn   
    );
    defalg(
      "Insertion sort (libjs)",
      isrt_binssort,
      0,
      tss_n2_4   
    );
    defalg(
      "Binary insertion sort (libjs)",
      isrt_binssort,
      0,
      tss_n2_2   
    );

    (*nalgsP) = nalgs;
  }

void tss_begin_tex_performance_table(FILE *wr, char *outname, char *tpname, int n)
  {
    fprintf(wr, "\\advance\\endlinechar by -256\n");
    fprintf(wr, "\\jstab\n");
    fprintf(wr, "  {t.%s.%010d}\n", tpname, n);
    fprintf(wr, "  {\n");
    fprintf(wr, "     %% max runs for each n = %d\n", MAXRUNS);
    fprintf(wr, "    \\begin{tabular}{|l|r|r|r|r|r|r|}\n");
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "      %-45s & %6s & %8s & %8s & %10s & %10s & %4s\\\\\n",
      "Algorithm", "n", "min", "max", "avg", "std", "runs"); 
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "      \\hline\n"); 
  }

void tss_finish_tex_performance_table(FILE *wr, char *tpname, int n)
  {
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "    \\end{tabular}\n"); 
    fprintf(wr, "  }\n");
    fprintf(wr, "  {\n");
    fprintf(wr, "    Performance data (\\texttt{%s})", tpname);
    fprintf(wr, " for $n = %d$ entries.\n", n);
    fprintf(wr, "  }\n");
    fprintf(wr, "\\advance\\endlinechar by 256\n");
  }
  
void tss_write_tex_algorithms_table(char *outname, int nalgs, tss_alg_t *alg[])
  {
    char *fname = NULL;
    asprintf(&fname, "%s-algs.tex", outname);
    FILE *wr = open_write(fname, TRUE);
    free(fname);

    fprintf(wr, "\\advance\\endlinechar by -256\n");
    fprintf(wr, "\\jstab\n");
    fprintf(wr, "  {t.algs}\n");
    fprintf(wr, "  {\n");
    fprintf(wr, "    \\begin{tabular}{|l|l|r|}\n");
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "      %-35s & %-40s & %5s\\\\\n", "Procedure", "Algorithm", "SMALL");
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "      \\hline\n"); 
    int ialg;
    for (ialg = 0; ialg < nalgs; ialg++)
      { tss_alg_t *pa = alg[ialg]; 
        fprintf(wr, "      \\pn%-35s & %-40s & %5d \\\\\n", 
          tss_protect(pa->name), pa->descr, pa->thr
        );
        fprintf(wr, "      \\hline\n"); 
      }
    fprintf(wr, "    \\end{tabular}\n");
    fprintf(wr, "  }\n");
    fprintf(wr, "  {\n");
    fprintf(wr, "    Algorithms and their parameters.\n");
    fprintf(wr, "  }\n");
    fprintf(wr, "\\advance\\endlinechar by 256\n");
  } 

double tss_n_logn(int n)
  { double x = (double)(n + 1);
    return x*log(x)/log(2.0);
  }
  
double tss_n2_2(int n)
  { double x = (double)(n + 1);
    return x*x/2.0;
  }
  
double tss_n2_4(int n)
  { double x = (double)(n + 1);
    return x*x/4.0;
  }
  
/* One second in the units used by {now}: */
#define SECOND 1.0e6

char *tss_protect(char *name)
  {
    /* Compute length of new string: */
    int nnew = 0; /* Length of protected string. */
    char *p = name;
    while ((*p) != 0) { if ((*p) == '_') { nnew++; } p++; nnew++; }
    char *new = notnull(malloc((nnew+3)*sizeof(char)), "no mem");
    /* Copy {name} into {new}, wrapping and protecting all '_'s: */
    p = name;
    char *q = new;
    (*q) = '{'; q++;
    while ((*p) != 0) 
      { if ((*p) == '_') { (*q) = '\\'; q++; }
        (*q) = (*p); 
        q++; p++;
      }
    (*q) = '}'; q++;
    /* Schwarzeneggerize the new string: */
    (*q) = 0;
    return new;      
  }

void tss_print_stats(FILE *wr, char *sname, int n, tss_stats_t *st, int ntests, stype_t tp)
  {
    /* Extract the relevant statistics: */
    double x[ntests];
    int nok = 0;
    int it;
    for (it = 0; it < ntests; it++)
      { tss_stats_t *sti = &(st[it]);
        double stx;
        switch(tp)
          { case TP_NCMP:
              stx = (double)sti->ncmp;
              break;
            case TP_TIME:
              stx = sti->time;
              break;
            default:
              affirm(FALSE, "bad tp");
          }
        if (stx >= 0.0) { x[nok] = stx; nok++; }
      }
    
    fprintf(wr, "      \\pn%-35s & %6d", tss_protect(sname), n);
    fprintf(stderr, "  %-45s n = %6d", sname, n);
    
    /* Compute min and max: */
    double xmin, xmax;
    xmin = xmax = x[0];
    for (it = 1; it < nok; it++)
      { double xi = x[it]; 
        if (xi < xmin) 
          { xmin = xi; }
        else if
          (xi > xmax)
        { xmax = xi; }
      }
    
    /* Print min and max: */
    switch(tp)
      { case TP_NCMP:
          { int nmin = (int)(xmin+0.5), nmax = (int)(xmax+0.5);
            fprintf(wr, " & %8d & %8d", nmin, nmax); 
            fprintf(stderr, "  %8d .. %8d", nmin, nmax);
          }
          break;
        case TP_TIME:
          { fprintf(wr, " & %8.3f & %8.3f", xmin/SECOND, xmax/SECOND);
            fprintf(stderr, " %8.3f .. %8.3f", xmin/SECOND, xmax/SECOND);
          }
          break;
        default:
          affirm(FALSE, "bad tp");
      }
    
    /* Compute average: */
    double avg;
    if (nok >= 1)
      { double tc = 0;
        for (it = 0; it < nok; it++)
          { double xi = x[it]; 
            tc += xi;
          }
        avg = tc/nok;
      }
    else
      { avg = 0.0; }
    
    /* Compute standard deviation: */
    double var, std;
    if (nok >= 2)
      { double tcc = 0;
        for (it = 0; it < nok; it++)
          { double dc = x[it] - avg; 
            tcc += dc*dc;
          }
        var = tcc/((double)nok-1);
        std = (var < 0.0 ? 0.0 : sqrt(var));
      }
    else
      { var = std = 0.0; }
    
    /* Print average and deviation out: */
    switch(tp)
      { case TP_NCMP:
          fprintf(wr, " & %10.1f & %10.1f ", avg, std);
          fprintf(stderr, "  %10.1f ± %10.1f", avg, std);
          break;
        case TP_TIME:
          fprintf(wr, " & %10.3f & %10.3f", avg/SECOND, std/SECOND);
          fprintf(stderr, "  %10.3f ± %10.3f", avg/SECOND, std/SECOND);
          break;
        default:
          affirm(FALSE, "bad tp");
      }
    fprintf(wr, " & %4d ", nok);
    fprintf(stderr, " (%4d)", nok);
    fprintf(wr, " \\\\\n");
    fprintf(stderr, "\n");
    fflush(wr);
  }

tss_stats_t tss_test_sorter(int *h, int n, tss_sorter_t srt, int_cmp_t cmp, int sgn)
  { 
    tss_stats_t st;
    /* Start with trivial perm */
    int i;
    for (i = 0; i < n; i++) { h[i] = i; }
    /* Sort the data: */
    ncmp = 0;
    double start = user_cpu_time_usec();
    srt(h, n, cmp, sgn);
    st.time = user_cpu_time_usec() - start;
    st.ncmp = ncmp;
    /* Check range: */
    for (i = 0; i < n; i++)
      { int hi = h[i]; 
        if ((hi < 0) || (hi >= n)) {  affirm(FALSE, "not 0..n-1"); }
      }
    /* Check ordering: */
    for (i = 1; i < n; i++)
      { if (sgn*cmp(h[i-1],h[i]) > 0) 
          { affirm(FALSE, "out of order"); }
      }
    /* Check permutation: */
    int seen[n]; 
    for (i = 0; i < n; i++) { seen[i] = 0; }
    for (i = 0; i < n; i++)
      { int hi = h[i]; 
        if (seen[hi]) {  affirm(FALSE, "not perm"); }
        seen[hi] = 1;
      }
    return st;
  }

int tss_cmp(int a, int b)
  {
    affirm((a >= 0) && (a < ndata), "tss_cmp: bad a");
    affirm((b >= 0) && (b < ndata), "tss_cmp: bad b");
    ncmp++;
    long int da = data[a], db = data[b];
    if (da < db) 
      { return -1; }
    else if (da > db)
      { return +1; }
    else
      { return 0; }
  }
