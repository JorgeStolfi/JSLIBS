#ifndef btc_bubble_parms_read_H
#define btc_bubble_parms_read_H

/* Reading BTC price bubble parameters. */
/* Last edited on 2015-04-29 23:43:18 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_parms_read(char* fName, int nd, char* dt[], int* nbP, btc_bubble_t** bpP);
  /* Read the bubble parameters {bp[0..nb-1]} from file with name
    {fName}, in the format described in the string {btc_bubble_parms_read_INFO} below. The
    number of bubbles {nd} and the array {bp} (allocated by the
    procedure) are returned in {*nbP} and {*bpP}. Uses the array
    {dt[0..nd-1]} to convert the bubble peak dates into day indices
    {id_fin_up,id_ini_dn} in {0..nd-1}. */
    
#define btc_bubble_parms_read_INFO \
  "  The model for the price series is the linear combination of one or more" \
  " /bubble functions/.  In general, each bubble function consists of three" \
  " consecutive parts: a non-decreasing exponential (/rally/) that ends" \
  " with value 1.0 on some date {DT_FIN_UP}, a constant section with value 1.0," \
  " and a non-increasing exponential (/decay/) that starts on some later" \
  " date {DT_INI_DN} with value 1.0.  The constant section may be absent, in" \
  " which case the two exponentials are joined, with value 1 on a single" \
  " date {DT_FIN_UP = DT_INI_DN}.\n" \
  "\n" \
  "  A file that describes a model must contain one data line for each bubble in the" \
  " model. Comments start with '#' and contine to the end of the" \
  " line.  Blank lines, and lines that contain only a" \
  " comment, are ignored.  Each data line line must " \
  "have the format\n" \
  "\n" \
  "    \"{COEF} {DT_INI_SG} {R_UP} {DT_FIN_UP} {DT_INI_DN} {R_DN} {DT_FIN_SG} {BTAG} {COLOR}\"\n" \
  "\n" \
  "  The {COEF} field is the coefficient of the bubble in the linear combination, i.e. the" \
  "\n" \
  " amplitude of the bubble at its peak. The fields {R_UP} and {R_DN} are the daily rates of change during" \
  " the rally phase and the decay phase, respectively.  The rate {R_UP} must" \
  " begreater than or equal to 1; if {R_UP} is 1, the bubble function is" \
  " constant (1.0) up to {DT_FIN_UP}.  The rate {R_DN} must be less" \
  " than or equal to 1; if {R_DN} is 1, the bubble function is" \
  " constant (1.0) after {DT_INI_DN}.  When the rate {R_UP} (or {R_DN}) is 1.0, the" \
  " dates {DT_FIN_UP} and {DT_INI_DN} must be the same.  If both rates are 1.0, the" \
  " dates are immaterial (but must still be equal).\n" \
  "\n" \
  "  The fields {DT_INI_SG} and {DT_FIN_SG} are dates that approximately bracket" \
  " the range of dates where the bubble is dominant or most relevant.  They do not affect the model, but" \
  " may be used by plotting or other scripts.  The" \
  " field {BTAG} is an alphanumeric identifier for the bubble, used in" \
  " debugging, headers, and in the name of some output files.  The field {COLOR} is" \
  " an RGB color value encoded as six hexadecimal digits, used by some plotting scripts."

#endif
