#define PROG_NAME "test_sRGB_bias"
#define PROG_DESC "finds the {bias} that best fits the sRGB encoding"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-20 14:41:59 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_sRGB_bias_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -expo {EXPO} ] \\\n" \
  "    [ -bias {BIAS} ] \\\n" \
  "    [ -mean | -max ] \\\n" \
  "    [ -weight { lows | uniform | highs } ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2007-08-01 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2007-08-01 Created  by J. Stolfi.\n" \
  "  2017-06-15 Added \"-midrange\" weighting.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_sRGB_bias_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program writes a table comparing the" \
  " functions {sample_conv_encode_sRGB} and {sample_conv_decode_sRGB} with" \
  " {sample_conv_gamma} with bias {BIAS} and exponents" \
  " {EXPO} and {1/EXPO}, respectively.\n" \
  "\n" \
  "  The program may be asked" \
  " to find the value of {BIAS} or {EXPO} that" \
  " minimizes the discrepancy between the two models.  The" \
  " criterion is the discrepancy" \
  " (absolute or relative, root-mean-squared or" \
  " maximum) between the two models, in both directions.\n" \
  "\n" \
  "  The program also writes a file \"out/optim.txt\" with the" \
  " result(s) of parameter optimization.  The file has one line" \
  " with \"{PARAM} = {OPT_VALUE}\" where {PARAM} is either \"bias\" or \"expo\".  The" \
  " line is blank if there was no optimization."

#define PROG_INFO_OPTS \
  "  -expo {EXPO}\n" \
  "    Value of {expo} to use.  If not specified, the program" \
  " tries to find the optimum expo.  At least one of \"-bias\"" \
  " and \"-expo\" must be specified.\n" \
  "\n" \
  "  -bias {BIAS}\n" \
  "    Value of {bias} to use.  If not specified, the program" \
  " tries to find the optimum bias.  At least one of \"-bias\"" \
  " and \"-expo\" must be specified.\n" \
  "\n" \
  "  -mean\n" \
  "  -max\n" \
  "    These options specify whether the programs should" \
  " consider the root-mean-squared discrepancy (\"-mean\") or" \
  " the maximum discrepancy (\"-max\") when choosing" \
  " the bias.  The default is \"-mean\".\n" \
  "\n" \
  "  -weight { lows | uniform | highs }\n" \
  "    These options specify whether the weight {wt} of each" \
  " squared error term in the mismatch should be {1/v} (\"lows\")" \
  " 1 (\"uniform\"), or {1-(1-v)^4} (\"highs\".  The default is \"-weight uniform\".\n"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <float_image.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <bool.h>
#include <minu_brent.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>

typedef struct options_t
  { double expo;     /* User-given {expo}, or NAN if not given. */
    double bias;     /* User-given {bias}, or NAN if not given. */
    bool_t mean;     /* RMS or maximum mismatch? */
    int weight;      /* Weighting option: 0 = lows, 1 = uniform, 2 = highs. */
  } options_t;

options_t *tb_parse_options(int argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int main(int argc, char **argv);

void tb_apply_both
  ( float v, 
    double expo, 
    double bias, 
    bool_t encode, 
    float *y_sr, 
    float *y_sg
  );
  /* If {encode} is true, evaluates {sample_conv_encode_sRGB(v)} and
    {sample_conv_gamma(v, expo, bias)}. If {encode} is false, evaluates
    the inverse functions {sample_conv_decode_sRGB(v)} and
    {sample_conv_gamma(v, 1/expo, bias)}. In either case, stores the
    results in {y_sr} and {y_sg}, respectively. */

double tb_mismatch(int nv, double expo, double bias, bool_t mean, int weight);
  /* Computes a quadratic discrepancy function between the sRGB 
    gamma encoding and decoding and the generic
    {sample_conv_gamma} with the given {bias} and
    exponents {expo} and {1/expo}, respectively.
    
    The discrepancy is computed at {nv} equally spaced values of the
    argument {v} in {[0 _ 1]}, avoiding and 1.
    
    If {weight=0} each squared error term {e^2} is multiplied
    by the weight {wt = 1/v}, thus giving more importance to low 
    values. If {weight=1}, it is left unweighted ({wt = 1}).
    If {weight=2}, it is multiplied by the weight {wt = 1-(1-v)^4},
    thus giving more importance to mid and high values.
    
    If {mean} is true, returns the (weighted) average of the squared
    error terms. If {mean} is false, returns the largest (weighted)
    squared error term.  Note that {mean=FALSE} and {weight!=1}
    mean minimizing the maximum value of {fabs(e)*sqrt(wt)}. */

void tb_print_comparison(FILE *wr, int nv, double expo, double bias, bool_t encode);
  /* Analyzes the difference between the two models in the {encode}
    direction, and the round-trip error incurred when by applying one
    model in the encode direction followed by the other model in the
    opposite direction. Prints the summary statistics to {wr}. */

void tb_write_table
  ( char *name, 
    char *suff, 
    int nv, 
    double expo,
    double bias,
    bool_t encode, 
    options_t *o
  );
  /* Writes to file "{name}{suff}" a list of {2*nv} lines 
    in the format
      "{Y} {Y_SR} {Y_SG} {Y_DIFF} {DY_SR} {DY_SG}"
    where {Y} is a linear-scale intensity value, 
    {Y_SR} is the corresponding sRGB value, {Y_SG}
    is the corresponding value computed by {sample_conv_gamma},
    {Y_DIFF} is {Y_SG - Y_SR}, and
    {DY_SR} and {DY_SG} are the numerical derivatives
    of {Y_SR} and {Y_SG}.
    
    The values {Y} span {[0 _ 1]} with {nv+1} equal steps.
    
    Also prints to {stderr} some error statistics. */

double tb_find_bias(int nv, double expo, bool_t mean, int weight);
  /* Computes a {bias} value that minimizes the mismatch
    {tb_mismatch(nv, bias, mean, weight)}. */

double tb_find_expo(int nv, double bias, bool_t mean, int weight);
  /* Computes a {expo} value that minimizes the mismatch
    {tb_mismatch(nv, expo, bias, mean, weight)}. */

void tb_write_opt(char *pname, double val);
  /* Writes the file with result {val} of optimizing parameter {pname}. */

int main(int argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    int nv_eval = 4096;
    int nv_plot = 4096;
    
    double bias = o->bias;
    double expo = o->expo;
    if (isnan(bias))  
      { bias  = tb_find_bias (nv_eval, expo, o->mean, o->weight); tb_write_opt("bias", bias); }
    else if (isnan(expo)) 
      { expo = tb_find_expo(nv_eval, bias,  o->mean, o->weight); tb_write_opt("expo", expo); }
    else
      { tb_write_opt(NULL, NAN); }
    
    int idir; /* Map direction: {-1} = decoding, {+1} = encoding. */
    for (idir = -1; idir < 2; idir += 2)
      { bool_t encode = (idir > 0);
      
        /* Print a sumamry of the differences: */
        tb_print_comparison(stderr, nv_eval, expo, bias, encode);
        fprintf(stderr, "\n");
        
        /* Write a comparison table for plotting: */
        char *oname = (encode ? "out/enc" : "out/dec"); /* Table filename. */
        tb_write_table(oname, ".txt", nv_plot, expo, bias, encode, o);
      }
    /* Cleanup: */
    free(o); o = NULL;
    return 0;
  }
  
double tb_find_bias(int nv, double expo, bool_t mean, int weight)
  {
    demand(! isnan(expo), "{expo} is undefined");
    
    auto bool_t eval(void *parms, double x, double *fx, double *dfx);
      /* Evaluates the goal function at {x}, result in {fx}. Sets {dfx} to NAN. */
    
    int nEvals = 0;
    int maxEvals = 200;
    
    bool_t eval(void *parms, double x, double *fx, double *dfx)
      { if (nEvals >= maxEvals) { return TRUE; }
        double x_bias = x;
        (*fx) = tb_mismatch(nv, expo, x_bias, mean, weight);
        (*dfx) = NAN;
        nEvals++;
        return FALSE;
      }

    bool_t debug = FALSE;
    
    /* Input/output variables: */
    double bias, fbias, dfbias;
    double a, b;
    
    /* Search parameters: */
    double tol = 0.0000001;
    double dist = 0.5;
    
    /* Initial range and guess: */
    a = 0.001; b = 0.500;
    bias = (a+b)/2; 
    fbias = tb_mismatch(nv, expo, bias, mean, weight); 
    dfbias = NAN;

    minu_brent_minimize
      ( /*parms*/ NULL,
        /*eval*/  &eval,    
        /*x*/     &bias,        
        /*fx*/    &fbias,       
        /*dfx*/   &dfbias,      
        /*tol*/   tol,       
        /*dist*/  dist,      
        /*a*/     &a,        
        /*b*/     &b,        
        /*check*/ NULL,
        /*debug*/ debug       
      );

    fprintf(stderr, "evaluations =  %12d\n", nEvals);
    fprintf(stderr, "fitted bias =  %12.7f\n", bias);
    fprintf(stderr, "mismatch =     %12.7f\n", sqrt(fbias));
    return bias;
  }

double tb_find_expo(int nv, double bias, bool_t mean, int weight)
  {
    demand(! isnan(bias), "{bias} is undefined");
    
    auto bool_t eval(void *parms, double x, double *fx, double *dfx);
      /* Evaluates the goal function at {x}, result in {fx}. Sets {dfx} to NAN. */
    
    int nEvals = 0;
    int maxEvals = 200;
    
    bool_t eval(void *parms, double x, double *fx, double *dfx)
      { if (nEvals >= maxEvals) { return TRUE; }
        double x_expo = x;
        (*fx) = tb_mismatch(nv, x_expo, bias, mean, weight);
        (*dfx) = NAN;
        nEvals++;
        return FALSE;
      }

    bool_t debug = FALSE;
    
    /* Input/output variables: */
    double expo, fexpo, dfexpo;
    double a, b;
    
    /* Search parameters: */
    double tol = 0.000001;
    double dist = 0.5;
    
    /* Initial range and guess: */
    a = 0.200; b = 0.700;
    expo = (a+b)/2; 
    fexpo = tb_mismatch(nv, expo, bias, mean, weight); 
    dfexpo = NAN;

    minu_brent_minimize
      ( /*parms*/ NULL,
        /*eval*/  &eval,    
        /*x*/     &expo,        
        /*fx*/    &fexpo,       
        /*dfx*/   &dfexpo,      
        /*tol*/   tol,       
        /*dist*/  dist,      
        /*a*/     &a,        
        /*b*/     &b,        
        /*check*/ NULL,
        /*debug*/ debug       
      );

    fprintf(stderr, "evaluations =  %12d\n", nEvals);
    fprintf(stderr, "fitted expo = %12.7f\n", expo);
    fprintf(stderr, "mismatch =     %12.7f\n", sqrt(fexpo));
    return expo;
  }

double tb_mismatch(int nv, double expo, double bias, bool_t mean, int weight)
  {
    demand(! isnan(expo), "{expo} is undefined");
    demand(expo > 0, "invalid {expo}");
    demand(! isnan(bias), "{bias} is undefined");
    demand(bias >= 0, "invalid {bias}");

    double sum_wt_dy2 = 0.0;
    double sum_wt = 0.0;
    double max_wt_dy2 = 0.0;
    int iv;
    for (iv = 0; iv < nv; iv++)
      { float v = (float)((iv + 0.5)/nv);
        int idir; /* Map direction: {-1} = decoding, {+1} = encoding. */
        for (idir = -1; idir < 2; idir += 2)
          { bool_t encode = (idir > 0);
        
            float y_sr, y_sg;
            tb_apply_both(v, expo, bias, encode, &y_sr, &y_sg);
            assert(y_sr > 0);
            assert(y_sg > 0);
        
            double dy = ((double)y_sr) - ((double)y_sg);
            double wt; /* Importance weight. */
            if (weight == 0)
              { wt = 1.0/v; /* Prioritize lows. */ }
            else if (weight == 1)
              { wt = 1.0;  /* Uniform weights. */ }
            else if (weight == 2)
              { /* Prioritize mids and highs */
                double tt = 1 - v; tt = tt*tt; tt = tt*tt;
                wt = 1-tt; 
              }
            else
              { assert(FALSE); }
            sum_wt_dy2 += wt*dy*dy;
            sum_wt += wt;
            max_wt_dy2 = fmax(wt*dy*dy, max_wt_dy2);
          }
      }
    return (mean ? sum_wt_dy2/sum_wt : max_wt_dy2);
  }
  
void tb_print_comparison(FILE *wr, int nv, double expo, double bias, bool_t encode)
  {
    demand(! isnan(expo), "{expo} is undefined");
    demand(expo > 0, "invalid {expo}");
    demand(! isnan(bias), "{bias} is undefined");
    demand(bias >= 0, "invalid {bias}");
    
    /* Differences in value {y} between the two models: */
    double sum_dy2 = 0.0;     /* Sum of unweighted squared diffs. */
    double sum_wt_dy2 = 0.0;  /* Sum of weighted squared diffs. */
    double max_ry = 0.0;      /* Largest relative diff. */
    double max_dy = 0.0;      /* Largest diff. */
    
    /* Differences in arg {v} for round-trip across both models: */
    double sum_dv2 = 0.0;     /* Sum of unweighted squared round-trip diffs. */
    double sum_wt_dv2 = 0.0;  /* Sum of weighted squared round-trip diffs. */
    double max_rv = 0;        /* Largest relative round-trip diff. */
    double max_dv = 0;        /* Largest round-trip diff. */

    /* Weight sums: */
    double sum_ = 0.0;        /* Number of terms. */
    double sum_wt = 0.0;      /* Sum of weights. */

    int iv;
    for (iv = 0; iv < nv; iv++)
      { float v = (float)((iv + 0.5)/nv);
        double wt = 1.0/v; /* Importance weight. */
        
        /* Apply both models to {v} in the {encode} direction: */
        float y_sr, y_sg;
        tb_apply_both(v, expo, bias, encode, &y_sr, &y_sg);
        assert(y_sr > 0);
        assert(y_sg > 0);

        /* Compute the absolute and relative diff between {y} results: */
        double dy = ((double)y_sr) - ((double)y_sg);
        double ry = 2 * dy/((double)y_sr) + ((double)y_sg);
        
        /* Accumulate the {y} differences: */
        sum_dy2 += dy*dy;
        sum_wt_dy2 += wt*dy*dy;
        if (fabs(dy) > fabs(max_dy)) { max_dy = dy; }
        if (fabs(ry) > fabs(max_ry)) { max_ry = ry; }
        
        /* Compute the absolute and relative round-trip error in {v}: */
        float v_sr_sr, v_sr_sg, v_sg_sr, v_sg_sg;
        tb_apply_both(y_sr, expo, bias, !encode, &v_sr_sr, &v_sr_sg);
        assert(fabs(v_sr_sr - v) < 1.0e-6*v);
        tb_apply_both(y_sg, expo, bias, !encode, &v_sg_sr, &v_sg_sg);
        assert(fabs(v_sg_sg - v) < 1.0e-6*v);
        double dv_sr_sg = ((double)v_sr_sg) - ((double)v);
        double dv_sg_sr = ((double)v_sg_sr) - ((double)v);
        double dv = (fabs(dv_sr_sg) > fabs(dv_sg_sr) ? dv_sr_sg : dv_sg_sr);
        double rv = dv/v;
        
        /* Accumulate the round-trip errors: */
        sum_dv2 += dv*dv;
        sum_wt_dv2 += wt*dv*dv;
        if (fabs(dv) > fabs(max_dv)) { max_dv = dv; }
        if (fabs(rv) > fabs(max_rv)) { max_rv = rv; }

        sum_ += 1;
        sum_wt += wt;
      }

    /* Names for {encode} and {!encode} directions: */
    char *xdir = (char*[2]){ "decode", "encode" }[encode];
    char *xinv = (char*[2]){ "decode", "encode" }[! encode];
    
    /* Compute and print the statistics: */
    double rms_wt_dy = sqrt(sum_wt_dy2/sum_wt);
    double rms_dy = sqrt(sum_dy2/sum_);
    fprintf(wr, "Assuming expo = %9.7f bias = %9.7f\n", expo, bias);
    fprintf(wr, "Differences in the %s direction:\n", xdir);
    fprintf(wr, "  rms weight =     %12.7f\n", rms_wt_dy);
    fprintf(wr, "  rms unweight =   %12.7f\n", rms_dy);
    fprintf(wr, "  maximum absolute = %12.7f\n", max_dy);
    fprintf(wr, "  maximum relative = %12.7f\n", max_ry);
    
    double rms_wt_dv = sqrt(sum_wt_dv2/sum_wt);
    double rms_dv = sqrt(sum_dv2/sum_);
    fprintf(wr, "Round-trip (%s -- %s) errors:\n", xdir, xinv);
    fprintf(wr, "  rms weighted =     %12.7f\n", rms_wt_dv);
    fprintf(wr, "  rms unweighted =   %12.7f\n", rms_dv);
    fprintf(wr, "  maximum absolute = %12.7f\n", max_dv);
    fprintf(wr, "  maximum relative = %12.7f\n", max_rv);
  }
  
void tb_write_opt(char *pname, double val)
  { FILE *wr = open_write("out/optim.txt", TRUE);
    if (! isnan(val)) { fprintf(wr, "%s = %.7f\n", pname, val); }
    fclose(wr);
  }

void tb_apply_both
  ( float v, 
    double expo, 
    double bias, 
    bool_t encode, 
    float *y_sr, 
    float *y_sg
  )
  { if (encode)
      { (*y_sr) = sample_conv_sRGB_encode(v);
        (*y_sg) = sample_conv_gamma(v, expo, bias);
      }
    else
      { (*y_sr) = sample_conv_sRGB_decode(v);
        (*y_sg) = sample_conv_gamma(v, 1/expo, bias);
      }
  }

void tb_write_table
  ( char *name, 
    char *suff, 
    int nv, 
    double expo,
    double bias,
    bool_t encode, 
    options_t *o
  )
  { 
    demand(! isnan(expo), "{expo} is undefined");
    demand(expo > 0, "invalid {expo}");
    demand(! isnan(bias), "{bias} is undefined");
    demand(bias >= 0, "invalid {bias}");

    char *fname = jsprintf("%s%s", name, suff);
    FILE *wr = open_write(fname, TRUE);
    double dx = 0.5/((double)nv); /* Step for finite differencing. */
    int ix;
    for (ix = 0; ix <= nv; ix++)
      { float x = (float)(((double)ix)/((double)nv));
        /* Evaluate the functions: */
        float y_sr, y_sg;
        tb_apply_both(x, expo, bias, encode, &y_sr, &y_sg);
        /* Evaluate the numerical derivatives: */
        float x0 = (float)(x - dx), y_sr0, y_sg0;
        tb_apply_both(x0, expo, bias, encode, &y_sr0, &y_sg0);
        float x1 = (float)(x + dx), y_sr1, y_sg1;
        tb_apply_both(x1, expo, bias, encode, &y_sr1, &y_sg1);
        double dy_sr = ((double)(y_sr1 - y_sr0))/dx/2;
        double dy_sg = ((double)(y_sg1 - y_sg0))/dx/2;
        fprintf(wr, "  %9.7f", x);
        fprintf(wr, "  %9.7f %9.7f", y_sr, y_sg);
        fprintf(wr, "  %+10.7f", y_sg - y_sr);
        fprintf(wr, "  %+10.7f %+10.7f", dy_sr, dy_sg);
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
  }

options_t *tb_parse_options(int argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* PARSE KEYWORD ARGUMENTS: */
    
    if (argparser_keyword_present(pp, "-expo"))
      { o->expo = argparser_get_next_double(pp, 1.0e-6, 1.0e+6); }
    else
      { o->expo = NAN; }
    
    if (argparser_keyword_present(pp, "-bias"))
      { o->bias = argparser_get_next_double(pp, 0.0, 0.999999); }
    else
      { o->bias = NAN; }
      
    if (isnan(o->expo) && isnan(o->bias))
      { argparser_error(pp, "must specify at least one of \"-expo\" or \"-bias\""); }
    
    if (argparser_keyword_present(pp, "-weight"))
      { char *wt = argparser_get_next(pp);
        if (strcmp(wt, "lows") == 0)
          { o->weight = 0; }
        else if (strcmp(wt, "uniform") == 0)
          { o->weight = 1; }
        else if (strcmp(wt, "highs") == 0)
          { o->weight = 2; }
        else
          { argparser_error(pp, "invalid \"-weight\" option"); }
      }
    else 
      { o->weight = 1; /* Uniform weights. */ }
      
    if (argparser_keyword_present(pp, "-mean"))
      { o->mean = TRUE; }
    else if (argparser_keyword_present(pp, "-max"))
      { o->mean = FALSE; }
    else 
      { o->mean = TRUE; }

    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
