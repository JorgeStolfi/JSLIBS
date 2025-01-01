#define PROG_NAME "test_sort_speed"
#define PROG_DESC "tests speed of quick-sort, heap-sort, binsertion-sort and merge-sort"
#define PROG_VERS "1.1"

/* Last edited on 2024-12-30 18:25:51 by stolfi */

#define test_sort_speed_COPYRIGHT \
  "Copyright © 2004  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME ""

#define PROG_INFO \
  "SYNOPSIS\n" \
  "" PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Tests various sorting algorithms and writes to" \
  " directory \"out\" tables with comparison counts" \
  " and running times.\n" \
  "\n" \
  "AUTHORS\n" \
  "  Created aug/2004 by Jorge Stolfi, IC-UNICAMP"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <intsort.h>
#include <intsort_extra.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jstime.h>
#include <jsmath.h>
#include <jsrandom.h>
  
/* One second in the units used by {now}: */
#define SECOND (1.0e6)
#define MILLISECOND (1.0e3)
#define MICROSECOND (1.0)

/* Max number of runs for each algorithm and each array size: */
#define MAXRUNS 50

/* Number of different sizes to test with: */
#define NUMSIZES 11

/* Maximum total time allowed for all test runs of each algoritm (usec): */
#define MAXALGTIME (300*SECOND)

typedef int32_t int32_t_cmp_t(int32_t a, int32_t b); 
  /* A signed comparison predicate for integers (indices, etc.). */

typedef void tss_sorter_t(int32_t ix[], uint32_t n, int32_t_cmp_t cmp, int32_t sgn);
  /* A procedure that sorts {ix} so that  {sgn*cmp(ix[i-1],ix[i]) <= 0}
    for all {i} */
  
typedef double tss_cost_pred_t(uint32_t n);
  /* A procedure that predicts the relative number of operations an
    algorithm will perform when sorting an {n}-element array.
    This may be the number of comparisons, or, for some
    algorithms, the number of element swaps or shifts.
    
    The actual value returned is not important, only the ratio between
    values for different {n}. In any case, the result must be 1.0 for
    {n=0}, and must be increasing with {n}. */
  
typedef enum {TP_NOPS = 0, TP_TIME = 1} stype_t;
  /* Type of performance metric for an algorithm. Metric {TP_TIME} is
    CPU time in usec, and metric {TP_NOPS} is number of operations
    performed. */
#define TP_FIRST TP_NOPS
#define TP_LAST TP_TIME

typedef struct tss_alg_t 
  { tss_sorter_t *srt;        /* The sorting procedure. */
    char *name;               /* The procedure's name. */
    char *descr;              /* A short description. */
    uint32_t n_base;          /* Max {n} included in the base case. */
    tss_cost_pred_t *asymp;   /* Estimator for relative cost as function of {n}. */
    double cfactor[2];        /* Constant factor to multiply into {asymp}. */
    bool_t shifter;           /* Count estimated shifts as operations. */
  } tss_alg_t;
  /* The array {cfactor} is indexed with {stype_t}:
    {cfactor[tp]*asymp(n)} is the theoretical predicted average of the
    {tp} metric (CPU time in usec if {tp=TP_TIME}, operation count if
    {tp=TP_NOPS}) for a run of the algorithm on {n} random inputs.
    
    If the {shifter} flag is false, the {TP_NOPS} metric for a run of
    the arlgorithm will be the number of comparisons actually made in
    the run.
    
    The {shifter} flag is true for those algorithms, like insertion sort
    with binary search, that performs relatively few ({\O(n \log n)}
    comparisons but many ({\O(n^2)}) element moves. It says that the
    {TP_NOPS} metric should be set to {n*(n-1)/4} instead of the actual
    number of comparisons. */
  
typedef struct tss_run_stats_t 
  { double data[2]; } tss_run_stats_t;
  /* Performance metrics for a single test run of a specific algorithm. The {data}
    array is indexed by {stype_t}: {data[TP_NOPS]} is the numer of
    operations done, and {data[TP_TIME]} is the CPU time in usec. */
  
typedef struct tss_alg_stats_t 
  { uint32_t nruns;
    double asymp;     /* Asymptotic formula. */
    double min[2];    /* Minimum. */
    double max[2];    /* Maximum. */
    double avg[2];    /* Average. */
    double dev[2];    /* Standard deviation. */
    double exp[2];    /* Expected value. */
  } tss_alg_stats_t;
  /* Condensed performance metrics for zero or more runs of 
    an algorithm for a specific value of {n}.
    
    The {nruns} field is the number of runs that were actually performed.
    
    The arrays {min,max,avg,dev,exp} are indexed by {stype_t} as in {tss_run_stats_t}.
    They summarize the metric {tp} over all those runs:
      
      * {min[tp]} and {max[tp]} are the min and max the metric,
        or {NAN} if {nruns} is zero.  
      
      * {avg[tp]} is the average of the metrics,
        or {NAN} if {nruns} is zero.  
    
      * {dev[tp]} is the standard deviation of that metric,
        or {NAN} if {nruns<2}.  
        
      * {exp[tp]} is the value of that metric predicted by the 
        theoretical formula.
        
    The CPU time data are in usec. */

/* PROTOS */

void tss_do_all_tests(char *outdir);
  /* Tests all algorithms for various array sizes, prints
    statistics to stderr and to files "{outdir}/{N}-{tpname}.tex" where 
    {tpname} is "nops" or "time" and {N} is the number of items 
    in the test array, formatted as "%010d". */
    
tss_run_stats_t tss_test_sorter(uint32_t n, tss_sorter_t srt, int32_t sgn);
  /* Runs sorter {srt} on an array of {ix[0..n-1]} integers, then checks
    whether the output is ordered and is a permutation of 
    the input.  Returns statistics.
    
    The integers {ix[0..n-1]} to be sorted are the integers {0..n-1},
    interpreted as indices into an array {data[0..n-1]} of distinct
    random 64-bit integers, generated internally. That is, the array
    {ix} is considered properly sorted if {ix[0..n-1]} is a permutation
    of {0..n-1}, and {data[ix[i]] < data[ix[j]]} for all {i,j} in
    {0..n-1} with {i<j}.
    
    The caller should initialize the random generator calling
    {srandom(seed)}, where the {seed} is reproducible but distinct for
    each call of {tss_test_sorter} with same {n} and {srt}. */
     
int64_t *tss_random_data(uint32_t n);
  /* Returns a vector of {n} random 64-bit integers, many negative,
    all distinct. */

tss_alg_stats_t tss_summarize_run_stats
  ( uint32_t nruns,
    double asymp,
    tss_run_stats_t st[],
    double nops_pred,
    double time_pred
  );
  /* Returns a single {tss_alg_stats_t} record {ast} that summarizes the
    performance records {st[0..nruns-1]} for {nruns} runs of some
    algorithm with random input arrays of some fixed size {n}.
    The {asymp} should be the asymptotic relative prediction of the 
    cost (without the constant factor). */
    
void tss_print_stats_line(char *sname, uint32_t n, tss_alg_stats_t *ast, stype_t tp);
  /* Assumes that {ast} is the summarized statistics for sorter {sname}
    acting on arrays of {n} elements. Writes to stderr a line with the 
    statistics of type {tp}. */
     
void tss_write_tex_table_line(FILE *wr, char *sname, uint32_t n, tss_alg_stats_t *ast, stype_t tp);
  /*  Assumes that {ast} is the summarized statistics for sorter {sname}
    acting on arrays of {n} elements.  Writes to {wr} a TeX table line with the 
    statistics of type {tp}. */
   
void tss_print_stderr_col_headers(void);
  /* Prints to {stderr} the column headers for the {stderr} output of {tss_print_stats}. */

char *tss_protect(char *name);
  /* Returns a copy of {name} wrapped in '{}'s
    with every '_' changed into '\_'. */

void tss_fill_algorithms_table(uint32_t *nalgsP, tss_alg_t *alg[]);
  /* Stores into {alg[0..nalgs-1]} a set of algorithms to test.
    Also sets {*nalgsP} the number {nalgs} of algorithms 
    defined.  Assumes that {alg} has at least {MAXALGS} entries. */

void tss_write_tex_algorithms_table(char *outdir, uint32_t nalgs, tss_alg_t *alg[]);
  /* Writes to "{outdir}/algs.tex" a TeX table with procedure name,
    algorithm description, and recursion threshold (when applicable). */

void tss_start_tex_performance_table(FILE *wr, char *tpname, char *tptitle, uint32_t n);
  /* Writes into {wr} the preamble of a TeX table with performance data for various
    algorithms on tables with {n} entries.  The {tpname} is "time" or "nops", and title 
    is a longer readable version thereof. */
    
void tss_finish_tex_performance_table(FILE *wr, char *tpname, char *tptitle, uint32_t n);
  /* The file {wr} must have been opened with
    {tss_start_tex_performance_table(wr,tpname,tptitle,n)}. Writes into
    {wr} the postamble of the table. */

tss_alg_t *tss_new_alg
  ( tss_sorter_t *srt,
    char *name,
    char *descr,
    uint32_t n_base,
    tss_cost_pred_t *asymp,
    double nops_factor,
    double time_factor,
    bool_t shifter
  );
  /* Creates a new algorithm data record.
    The {time_factor} should be in usec. */
  
double tss_n_logn(uint32_t n);
  /* Returns approximately {n log2(n)}, more precisely {(n+1)log2(n+2)}
    so that it is 1 for {n=0}. */
  
double tss_n2(uint32_t n);
  /* Returns approximately {n^2}, more precisely {(n+1)^2} so that it is
    1 for {n=0}. */
  
int32_t main (int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { tss_do_all_tests("out");
    return 0;
  }

/* Max number of algorithms to test: */
#define MAXALGS 10

void tss_do_all_tests(char *outdir)
  {
    char *tpname[2];
    tpname[TP_NOPS] = "nops"; 
    tpname[TP_TIME] = "time";
    
    char *tptitle[2];
    tptitle[TP_NOPS] = "operation counts"; 
    tptitle[TP_TIME] = "running time (msec)";
    
    uint32_t nalgs = 0;
    tss_alg_t *alg[MAXALGS];
    
    /* Assemble a table of algorithms and parameters: */
    tss_fill_algorithms_table(&nalgs, alg);

    /* Write table describing algorithms: */
    tss_write_tex_algorithms_table(outdir, nalgs, alg);
    
    /* The values of {n} to use in the tests: */
    uint32_t size[NUMSIZES] = { 0, 1, 2, 4, 64, 96, 128, 192, 256, 64*64, 64*64*64 }; 

    for (uint32_t kn = 0; kn < NUMSIZES; kn++) 
      { /* Pick the number {n} of items to sort in each run: */
        uint32_t n = size[kn];
        fprintf(stderr, "testing all algorithms with n = %d\n", n);
      
        /* Theoretical predictions for each alg with random {n} inputs: */
        double nops_pred[nalgs]; /* Average num operations. */
        double time_pred[nalgs]; /* Average running time (usec). */
        double *pred[2] = { [TP_NOPS] = nops_pred, [TP_TIME] = time_pred };
        
        tss_alg_stats_t ast[nalgs]; /* Summary stats for all runs of each algorithm. */

        for (uint32_t ialg = 0; ialg < nalgs; ialg++)
          { /* Compute theoretical performance: */
            double asymp = alg[ialg]->asymp(n);
            for (stype_t tp = TP_FIRST; tp <= TP_LAST; tp++)
              { pred[tp][ialg] = asymp*alg[ialg]->cfactor[tp];
                fprintf(stderr, "  predicted %s = %20.6f for %s\n", tpname[tp], pred[tp][ialg], alg[ialg]->name); 
              }
            
            /* Define the max number of runs for this algorithm: */
            uint32_t nruns = (uint32_t)(MAXALGTIME/time_pred[ialg]);
            if (nruns > MAXRUNS) { nruns = MAXRUNS; }
            fprintf(stderr, "  doing %d runs ...\n", nruns);
            
            /* Performance statistics for the various runs of this algorithm: */
            tss_run_stats_t rst[nruns]; /* Run {r} of alg {ialg} is {st[ialg*MAXRUNS + r]} */
        
            int32_t sgn = +1;
            for (uint32_t it = 0; it < nruns; it++)
              { tss_alg_t *pa = alg[ialg];
                srandom(666*418+2*it+1);
                rst[it] = tss_test_sorter(n, pa->srt, sgn);
                if (pa->shifter) 
                  { /* Use the estimated number of shifts for insertion sorters. */
                    double dn = (double)n;
                    rst[it].data[TP_NOPS] = (dn*(dn-1))/4.0;
                  }
                sgn = -sgn;
              }
            fprintf(stderr, "  finished the %d runs\n", nruns);
            ast[ialg] = tss_summarize_run_stats(nruns, asymp, rst, nops_pred[ialg], time_pred[ialg]);
          }
        

        /* Write the tables for this {n}: */
        for (stype_t tp = TP_FIRST; tp <= TP_LAST; tp++)
          { 
            /* Open performance data file for {n} and {tp}: */
            FILE *wr;
            char *fname = jsprintf("%s/%010d-%s.tex", outdir, n, tpname[tp]);
            wr = open_write(fname, TRUE);
            free(fname);

            tss_start_tex_performance_table(wr, tpname[tp], tptitle[tp], n);
            fprintf(stderr, "  -- %s --\n", tptitle[tp]);
            tss_print_stderr_col_headers();
            
            for (uint32_t ialg = 0;  ialg < nalgs; ialg++)
              { tss_alg_t *pa = alg[ialg];
                tss_print_stats_line(pa->descr, n, &(ast[ialg]), tp);
                tss_write_tex_table_line(wr, pa->descr, n, &(ast[ialg]), tp);
              }
              
            fprintf(stderr, "\n");
            tss_finish_tex_performance_table(wr, tpname[tp], tptitle[tp], n);
            fclose(wr);
          }
        fprintf(stderr, "\n");
      }
    
    return;
  }

tss_run_stats_t tss_test_sorter
  ( uint32_t n,
    tss_sorter_t srt,
    int32_t sgn
  )
  { 
    /* Allocate and fill the data array: */
    int64_t *data = tss_random_data(n);

    /* Comparison counter: */
    uint64_t ncmp;

    auto int32_t compare_data(int32_t a, int32_t b);
      /* Compares {data[a]} with {data[b]}, and increments {ncmp}. */
    
    /* Allocate index array */
    int32_t *ix = talloc(n, int32_t);

    /* Start {ix} with trivial perm of indices {0..n-1}. */
    for (int32_t i = 0;  i < n; i++) { ix[i] = i; }

    /* Sort {ix} according to {data}: */
    tss_run_stats_t st;
    double start = user_cpu_time_usec();
    ncmp = 0;
    srt(ix, n, compare_data, sgn);
    st.data[TP_TIME] = user_cpu_time_usec() - start;
    st.data[TP_NOPS] = (double)ncmp;
    
    if (n > 0)
      { /* Check that output array is a permutation of {0..n-1}: */
        bool_t seen[n]; 
        for (uint32_t i = 0;  i < n; i++) { seen[i] = FALSE; }
        for (uint32_t i = 0;  i < n; i++)
          { int32_t hi = ix[i]; 
            affirm((hi >= 0) && (hi < n), "sorted array elemnet is not in {0..n-1}");
            affirm (! seen[hi], "sorted array is not a permutation");
            seen[hi] = TRUE;
          }

        /* Check ordering of output array: */
        ncmp = 0;
        for (uint32_t i = 1;  i < n; i++)
          { if (sgn*compare_data(ix[i-1],ix[i]) > 0) { affirm(FALSE, "out of order"); } }
        assert(ncmp == (n == 0 ? 0 : n-1)); /* Check on {compare}. */
     }

    free(data);
    free(ix);
    /* Return performance data: */
    return st;

    int32_t compare_data(int32_t a, int32_t b)
      {
        affirm((a >= 0) && (a < n), "tss_cmp: bad a");
        affirm((b >= 0) && (b < n), "tss_cmp: bad b");
        ncmp++;
        int64_t da = data[a], db = data[b];
        if (da < db) 
          { return -1; }
        else if (da > db)
          { return +1; }
        else
          { return 0; }
      }
  }
  
void tss_print_stderr_col_headers(void)
  { char *fmt = "  %-45s %6s  %14s %14s  %14s %14s  %14s %12s %14s  %4s\n";
    fprintf
      ( stderr, fmt, 
        "Algorithm", "n", 
        "min", "max", 
        "avg", "std", 
        "pred", "avg/pred", "avg/asymp", 
        "runs"
      ); 
    fprintf
      ( stderr, fmt, 
        "---------------------------------------------", "------", 
        "--------------", "--------------", 
        "--------------", "--------------", 
        "--------------", "------------", "--------------", 
        "----"
      ); 
  }

tss_alg_t *tss_new_alg
  ( tss_sorter_t *srt,
    char *name,
    char *descr,
    uint32_t n_base,
    tss_cost_pred_t *asymp,
    double nops_factor,
    double time_factor,
    bool_t shifter
  )
  { tss_alg_t *pa = (tss_alg_t *)notnull(malloc(sizeof(tss_alg_t)), "no mem");
    pa->srt = srt; 
    pa->name = name;
    pa->descr = descr;
    pa->n_base = n_base; 
    pa->asymp = asymp;
    pa->cfactor[TP_NOPS] = nops_factor;
    pa->cfactor[TP_TIME] = time_factor;
    pa->shifter = shifter;
    return pa;
  }
    
void tss_fill_algorithms_table(uint32_t *nalgsP, tss_alg_t *alg[])
  {
    uint32_t nalgs = 0;
    
    #define defalg(DESCR,SRT,N_BASE,ASYMP,NOPSF,TIMEF,SHIFTER) \
      do { \
        assert(nalgs < MAXALGS); \
        alg[nalgs] = tss_new_alg(&SRT, #SRT, DESCR, N_BASE, &ASYMP, NOPSF, TIMEF, SHIFTER); \
        nalgs++; \
      } while (0)
    
    defalg(
      "Heapsort (classic)",
      isrt_heapsort_classic, 
      0,
      tss_n_logn, 1.8507, 0.19556094, FALSE   
    );
    defalg(
      "Heapsort (vacsink)",
      isrt_heapsort_vacsink,     
      0,
      tss_n_logn, 1.0547, 0.14718884, FALSE   
    );
    defalg(
      "Heapsort (libjs)",
      isrt_heapsort,       
      0,
      tss_n_logn, 1.0547, 0.16295728, FALSE 
    );
    defalg(
      "Quicksort (middle elem)",
      isrt_quicksort_middle,
      isrt_quicksort_middle_SMALL,
      tss_n_logn, 1.2355, 0.17148143, FALSE  
    );
    defalg(
      "Quicksort (median-of-three)",
      isrt_quicksort_median3,    
      isrt_quicksort_median3_SMALL,
      tss_n_logn, 1.1393, 0.15575583, FALSE  
    );
    defalg(
      "Mergesort (in-place, piv)",
      isrt_mergesort_pivot,    
      isrt_mergesort_pivot_SMALL,
      tss_n_logn, 1.0038, 0.17781383, FALSE   
    );
    defalg(
      "Mergesort (in-place, sym)",
      isrt_mergesort_symsplit,  
      isrt_mergesort_symsplit_SMALL,
      tss_n_logn, 2.4507, 0.38175505, FALSE   
    );
    defalg(
      "Mergesort (libjs)",
      isrt_mergesort,      
      0,
      tss_n_logn, 1.0038, 0.18447458, FALSE 
    );
    defalg(
      "InsertSort (libjs)",
      isrt_binssort,
      0,
      tss_n2, 0.2500, 0.00051229, TRUE   
    );
    defalg(
      "BinaryInsertsort (libjs)",
      isrt_binssort,
      0,
      tss_n2, 0.2500, 0.00054182, TRUE
    );

    (*nalgsP) = nalgs;
  }

tss_alg_stats_t tss_summarize_run_stats
  ( uint32_t nruns,
    double asymp,
    tss_run_stats_t st[],
    double nops_pred,
    double time_pred
  )
  {
    tss_alg_stats_t ast;
    ast.nruns = nruns;
    ast.asymp = asymp;
    
    ast.exp[TP_NOPS] = nops_pred;
    ast.exp[TP_TIME] = time_pred;
    
    for (stype_t tp = TP_FIRST; tp <= TP_LAST; tp++)
      { ast.min[tp] = NAN;
        ast.max[tp] = NAN;
        ast.avg[tp] = NAN;
        ast.dev[tp] = NAN;
        
        if (nruns >= 1)
          { /* Compute {min}, {max}, {avg}: */
            double dmin = +INF, dmax = -INF;
            double sum_d = 0;
            for (uint32_t i = 0;  i < nruns; i++)
              { double di = st[i].data[tp];
                assert(isfinite(di) && (di >= 0));
                if (di < dmin) { dmin = di; }
                if (di > dmax) { dmax = di; }
                sum_d += di;
              }
            ast.min[tp] = dmin;
            ast.max[tp] = dmax;
            ast.avg[tp] = sum_d/nruns;
            
            if (nruns >= 2)
              { /* Compute {dev}: */
                double sum_d2 = 0;
                for (uint32_t i = 0;  i < nruns; i++)
                  { double di = st[i].data[tp] - ast.avg[tp];
                    sum_d2 += di*di;
                  }
                ast.dev[tp] = sqrt(sum_d2/(nruns-1));
              }
          }
      }
    return ast;
  }

void tss_start_tex_performance_table(FILE *wr, char *tpname, char *tptitle, uint32_t n)
  {
    fprintf(wr, "\\advance\\endlinechar by -256\n");
    fprintf(wr, "\\jstab\n");
    fprintf(wr, "  {t.%s.%010d}\n", tpname, n);
    fprintf(wr, "  {\n");
    fprintf(wr, "    \\begin{tabular}{|l|r|r|r|r|r|r|r|r|r|}\n");
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "      %-45s & %6s & %14s & %14s & %14s & %14s & %14s & %12s & %14s & %4s\\\\\n",
      "Algorithm", "n", "min", "max", "avg", "std", "pred", "avg/pred", "avg/asymp", "runs"); 
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "      \\hline\n"); 
  }

void tss_print_stats_line
  ( char *sname,
    uint32_t n,
    tss_alg_stats_t *ast,
    stype_t tp
  )
  {
    fprintf(stderr, "  %-45s %6d", sname, n);
    double min = ast->min[tp];
    double max = ast->max[tp];
    double avg = ast->avg[tp];
    double dev = ast->dev[tp];
    double exp = ast->exp[tp];
    switch(tp)
      { case TP_NOPS:
          { if (! isnan(min))
              { uint64_t nmin = (uint64_t)(min+0.5), nmax = (uint64_t)(max+0.5);
                fprintf(stderr, "  %14lu %14lu", nmin, nmax);
              }
            else
              { fprintf(stderr, "  %14s %14s", "", ""); }

            if (! isnan(avg))
              { fprintf(stderr, "  %14.1f", avg); }
            else
              { fprintf(stderr, "  %14s", ""); }

            if (! isnan(dev))
              { fprintf(stderr, " %14.1f", dev); }
            else
              { fprintf(stderr, " %14s", ""); }
              
            assert(! isnan(exp));
            fprintf(stderr, "  %14.1f", exp);
            
            if (! isnan(avg))
              { fprintf(stderr, " %12.4f %14.4f", avg/exp, avg/ast->asymp); }
            else
              { fprintf(stderr, " %12s %14s", "", ""); }
          }
          break;
        case TP_TIME:
          { if (! isnan(min))
              { fprintf(stderr, "  %14.3f %14.3f", min/MILLISECOND, max/MILLISECOND); }
            else
              { fprintf(stderr, "  %14s %14s", "", ""); }

            if (! isnan(avg))
              { fprintf(stderr, "  %14.4f", avg/MILLISECOND); }
            else
              { fprintf(stderr, "  %14s", ""); }

            if (! isnan(dev))
              { fprintf(stderr, " %14.4f", dev/MILLISECOND); }
            else
              { fprintf(stderr, " %14s", ""); }
              
            assert(! isnan(exp));
            fprintf(stderr, "  %14.4f", exp/MILLISECOND);
            
            if (! isnan(avg))
              { fprintf(stderr, " %12.6f %14.8f", avg/exp, avg/ast->asymp/MICROSECOND); }
            else
              { fprintf(stderr, " %12s %14s", "", ""); }
          }
          break;
        default:
          affirm(FALSE, "bad tp");
      }
    fprintf(stderr, "  %4d", ast->nruns);
    fprintf(stderr, "\n");
  }

void tss_write_tex_table_line
  ( FILE *wr,
    char *sname,
    uint32_t n,
    tss_alg_stats_t *ast,
    stype_t tp
  )
  {
    fprintf(wr, "      \\pn%-35s & %6d", tss_protect(sname), n);
    double min = ast->min[tp];
    double max = ast->max[tp];
    double avg = ast->avg[tp];
    double dev = ast->dev[tp];
    double exp = ast->exp[tp];
    switch(tp)
      { case TP_NOPS:
          { if (! isnan(min))
              { uint64_t nmin = (uint64_t)(min+0.5), nmax = (uint64_t)(max+0.5);
                fprintf(wr, " & %14lu & %14lu", nmin, nmax);
              }
            else
              { fprintf(wr, " & %14s & %14s", "\\none", "\\none"); }

            if (! isnan(avg))
              { fprintf(wr, " & %14.1f", avg); }
            else
              { fprintf(wr, " & %14s", "\\none"); }


            if (! isnan(dev))
              { fprintf(wr, " & %14.1f", dev); }
            else
              { fprintf(wr, " & %14s", "\\none"); }
              
            assert(! isnan(exp));
            fprintf(wr, " & %14.1f", exp);
            
            if (! isnan(avg))
              { fprintf(wr, " & %12.4f & %14.8f", avg/exp, avg/ast->asymp); }
            else
              { fprintf(wr, " & %12s & %14s", "\\none", "\\none"); }
          }
          break;
        case TP_TIME:
          { if (! isnan(min))
              { fprintf(wr, " & %14.3f & %14.3f", min/MILLISECOND, max/MILLISECOND);
              }
            else
              { fprintf(wr, " & %14s & %14s", "\\none", "\\none"); }

            if (! isnan(avg))
              { fprintf(wr, " & %14.4f", avg/MILLISECOND); }
            else
              { fprintf(wr, " & %14s", "\\none"); }


            if (! isnan(dev))
              { fprintf(wr, " & %14.4f", dev/MILLISECOND); }
            else
              { fprintf(wr, " & %14s", "\\none"); }
              
            assert(! isnan(exp));
            fprintf(wr, " & %14.4f", exp/MILLISECOND);
            
            if (! isnan(avg))
              { fprintf(wr, " & %12.4f & %14.8f", avg/exp, avg/ast->asymp/MICROSECOND); }
            else
              { fprintf(wr, " & %12s & %14s", "\\none", "\\none"); }
          }
          break;
        default:
          affirm(FALSE, "bad tp");
      }
    fprintf(wr, " & %4d ", ast->nruns);
    fprintf(wr, " \\\\\n");
    fflush(wr);
  }

void tss_finish_tex_performance_table(FILE *wr, char *tpname, char *tptitle, uint32_t n)
  {
    fprintf(wr, "      \\hline\n"); 
    fprintf(wr, "    \\end{tabular}\n"); 
    fprintf(wr, "  }\n");
    fprintf(wr, "  {\n");
    fprintf(wr, "    Average %s (\\texttt{%s})", tptitle, tpname);
    fprintf(wr, " for sorting $n = %d$ items.\n", n);
    fprintf(wr, "  }\n");
    fprintf(wr, "\\advance\\endlinechar by 256\n");
  }
  
void tss_write_tex_algorithms_table(char *outdir, uint32_t nalgs, tss_alg_t *alg[])
  {
    char *fname = jsprintf("%s/algs.tex", outdir);
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

    for (uint32_t ialg = 0;  ialg < nalgs; ialg++)
      { tss_alg_t *pa = alg[ialg]; 
        fprintf(wr, "      \\pn%-35s & %-40s & %5d \\\\\n", 
          tss_protect(pa->name), pa->descr, pa->n_base
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

double tss_n_logn(uint32_t n)
  { double x = (double)(n + 1);
    return x*log(x + 1)/log(2.0);
  }
  
double tss_n2(uint32_t n)
  { double x = (double)(n + 1);
    return x*x;
  }

char *tss_protect(char *name)
  {
    /* Compute length of new string: */
    uint32_t nnew = 0; /* Length of protected string. */
    char *p = name;
    while ((*p) != 0) { if ((*p) == '_') { nnew++; } p++; nnew++; }
    char *new = talloc(nnew+3, char);
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
    
int64_t *tss_random_data(uint32_t n)
  { demand(n >= 0, "invalid {n}");
    int64_t *data = talloc(n, int64_t);
    if (n > 0) 
      { /* Fill data with a null value, for paranois: */
        for (uint32_t i = 0;  i < n; i++) { data[i] = INT64_MIN; }

        /* Now set {data[i]} to {±q*n + p(i)} where {q} is random and
          {p} is some permutation of {0..n-1}. This way the elements
          {data[0..n-1]} will be in random order but all distinct. */

        /* Find range of {q}: */
        int64_t dmax = (1L << 60) - 1;     /* Max value of {data[i]}. */
        int64_t qmax = dmax/((int64_t)n);  /* Max value of {q}, */
        assert(qmax > 0); /* Paranoia. */

        /* Find a {stride} close to {sqrt(n)} relatively prime to it: */
        uint32_t stride = (uint32_t)floor(sqrt(n));
        while ((stride > 1) && (gcd(stride, n) != 1)) { stride--; }
        assert(gcd(stride, n) == 1);

        /* Fill {q} in steps of {stride} to improve randomness: */
        uint32_t k = 0;
        int64_t sgn = +1;
        for (uint32_t i = 0;  i < n; i++)
          { int64_t q = int64_abrandom(0,qmax);
            int64_t dk = q*((int64_t)n) + (int64_t)i;
            assert(dk > INT64_MIN);
            data[k] = dk;
            k = (k + stride) % n;
            sgn = -sgn;
          }
      }
    return data;
  }
