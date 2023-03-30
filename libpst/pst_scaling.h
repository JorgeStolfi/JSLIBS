#ifndef pst_scaling_H
#define pst_scaling_H

/* pst_scaling.h -- tools for scaling sample values. */
/* Created 2006-04-20 by Jorge Stolfi, IC-UNICAMP. */
/* Last edited on 2023-03-19 15:27:58 by stolfi */

#include <bool.h>
#include <vec.h>
#include <argparser.h>
#include <float_image.h>

#include <pst_basic.h>

/* PARSING SCALING PARAMETERS FROM THE COMMAND LINE */
  
double_vec_t pst_scaling_parse_range_option(argparser_t *pp, char *key, int *NC);
  /* If the keyword {key} is present, marks it as parsed, then 
    parses the following arguments as a tuple of one more numbers,
    using {pst_double_vec_parse}.
    
    If {NC} is is NULL, or only one numeric argument is present (with
    optional denominator), ignores {NC}. Otherwise, if {*NC} is
    negative, sets {*NC} to the number of elements read. Otherwise
    demands and parses exactly {*NC} numeric arguments (with an
    optional denominator). */

#define pst_scaling_option_list_INFO \
  "\"-min\", \"-max\", \"-center\", and \"-width\""

/* Low endpoint */

#define pst_scaling_parse_min_one_HELP \
  "-min {VMIN} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_min_RGB_HELP \
  "-min {VMIN_R} {VMIN_G} {VMIN_B} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_min_any_HELP \
  "-min {VMIN}.. " pst_double_vec_spec_den_HELP

#define pst_scaling_parse_min_INFO \
  "Specifies the lower bound {VMIN}"

/* High endpoint */

#define pst_scaling_parse_max_one_HELP \
  "-max {VMAX} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_max_RGB_HELP \
  "-max {VMAX_R} {VMAX_G} {VMAX_H} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_max_any_HELP \
  "-max {VMAX}.. " pst_double_vec_spec_den_HELP

#define pst_scaling_parse_max_INFO \
  "Specifies the upper bound {VMAX}"

/* Range width */

#define pst_scaling_parse_center_one_HELP \
  "-center {VCTR} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_center_RGB_HELP \
  "-center {VCTR_R} {VCTR_G} {VCTR_H} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_center_any_HELP \
  "-center {VCTR}.. " pst_double_vec_spec_den_HELP

#define pst_scaling_parse_center_INFO \
  "Specifies the center {VCTR = (VMAX+VMIN)/2}" \
  

/* Range center */

#define pst_scaling_parse_width_one_HELP \
  "-width {VWID} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_width_RGB_HELP \
  "-width {VWID_R} {VWID_G} {VWID_H} " pst_double_vec_spec_den_HELP
#define pst_scaling_parse_width_any_HELP \
  "-width {VWID}.. " pst_double_vec_spec_den_HELP

#define pst_scaling_parse_width_INFO \
  "Specifies the width {VWID = VMAX-VMIN}"

bool_t pst_scaling_parse_uniform(argparser_t *pp, bool_t next);
  /* Parses the \"-uniform\" switch; returns TRUE if 
     present, FALSE otherwise.  
     
     If {next} is TRUE, the procedure looks for the keyword only at
     the next argument; otherwise it looks among all arguments that
     are still unparsed.  */

#define pst_scaling_parse_uniform_HELP \
  "-uniform"

#define pst_scaling_parse_uniform_INFO \
  "Specifies that the same sample conversion must be used" \
  " in all color channels. If this option is given, any scaling" \
  " parameters (" pst_scaling_option_list_INFO ") that are" \
  " specified should have either a single value, or the" \
  " same value for all channels."

#define pst_scaling_parse_uniform_HELP_INFO \
  "  " pst_scaling_parse_uniform_HELP "\n" \
  "    " pst_scaling_parse_uniform_INFO 

/* UNIFORMIZATION AND COMPLETION OF SCALING PARAMETERS */

void pst_scaling_fix_channels(int NC, int32_vec_t *channel);
  /* Adjusts a list of channel indices to have exactly {NC} elements.
    
    If the {*channel} vector has exactly {NC} channels, the procedure
    does nothing. If {NC} is positive but the given {*channel} vector
    has zero elements, the procedure extends it to {NC} elements,
    which are set to consecutive indices {0,1,...NC-1}. If {NC} is 2
    or more but {*channel} has a single element, the procedure expands
    it to {NC} elements, all equal to element 0. Otherwise the
    procedure fails. */
    
void pst_scaling_use_actual_range
  ( float smin,
    float smax,
    double *min, 
    double *max, 
    double *ctr, 
    double *wid
  );  
  /* Makes sure that no more than two of the scaling range parameters
    {*min,*max,*ctr,*wid} are unspecified ({±INF}). The needed
    defaults are obtained from the given range [{smin} _ {smax}].
    
    If exactly two of the arguments are already defined, no
    adjustments are made. Fails if three or more of the arguments are
    already defined. The other cases are described in
    {pst_scaling_set_defaults_from_actual_range_INFO}. */

#define pst_scaling_use_actual_range_INFO \
  "If none of the four arguments is present, {VMIN} and {VMAX} are equal" \
  " to {SMIN} and {SMAX}, respectively the minimum and maximum" \
  " sample values in the input.  If only {VCTR} is given," \
  " then {VWID} is the smallest" \
  " width such that {SMIN} and {SMAX} lie in the range" \
  " [{VCTR-VWID/2} _ {VCTR+VWID/2}}.  If only {VMIN}" \
  " is given, then {VMAX} is {max(VMIN,SMAX)}.  If only" \
  " {VMAX} is given, then {VMIN} is {min(SMIN,VMAX)}.  If only" \
  " {VWID} is given, then {VCTR} is {(SMIN+SMAX)/2}."
  
#define pst_scaling_use_actual_range_with_uniform_INFO \
  "If the \"-uniform\" option is specified, the values of {SMIN} and {SMAX}" \
  " are computed over all relevant channels of the image. Otherwise," \
  " they are computed separately for each channel."

void pst_scaling_fix_params
  ( int NC,
    bool_t uniform,
    double_vec_t *min,
    double_vec_t *max,
    double_vec_t *ctr,
    double_vec_t *wid,
    float_image_t *fim, 
    int32_vec_t *channel
  );
  /* Adjusts a bunch of channel selection and scaling parameters for
    {NC} channels, as may be given in the command line.
    
    Each of the input parameter vectors {min,max,ctr,wid}
    must have either zero elements, or one element, or {NC} elements.
    
    If a parameter is unspecified (represented by the value {INF}),
    the procedure either sets it to a suitable default value, or
    computes it from other known parameters. If a single value was
    given for all channels, the procedure replicates it for all
    channels.
    
    After the procedure returns with success, the parameter vectors
    {channel,min,max,ctr,wid} will have exactly {NC} elements, no
    elements will be {INF}, and the parameters of each channel will be
    consistent with each other. 
    
    If there are less than two known parameters for each channel, and
    {fim} is not NULL, the procedure tries to use the minimum and
    maximum sample values {smin,smax} from the image {*fim} to provide
    suitable defaults, as explained for
    {pst_scaling_use_actual_range}. In that case, if
    {channel} is not NULL, the procedure assumes that the channels to
    analyze are {channel.e[0..NC-1]}. If {channel} is NULL, channels
    {0,1,...NC-1} are analyzed.
    
    In any case, once there are exactly two known parameter values for
    each channel, the other two are computed with
    {pst_scaling_complete_params}. The procedure fails if, for any
    channel, more than two parameters are defined.
    
    If {uniform} is TRUE, forces the scaling parameters to be the same
    across all channels. In that case, all input vectors which are not
    empty or trivial must contain {NC} equal values; and, if the
    acrual extreme samples {smin,smax} are needed, they are computed
    once over all relevant channels of {fim}, rather than separately
    for each channel. */ 

#define pst_scaling_num_values_INFO \
  "If only one value is specified, it applies to all" \
  " channels; otherwise there should be one value for each channel."

#define pst_scaling_use_actual_range_without_uniform_INFO \
  "The values of {SMIN} and {SMAX} are computed separately for" \
  " each channel of the image." 

void pst_scaling_complete_params
  ( double *min, 
    double *max, 
    double *ctr, 
    double *wid
  );  
  /* Provides values for any scaling range arguments
    {*min,*max,*ctr,*wid} that are unspecified ({±INF}),
    by using the fundamental equations.  Exactly two of those
    four values must have been specified; the procedure fails 
    otherwise. */

#define pst_scaling_complete_params_INFO \
  "The missing values for the arguments" \
  " {VMIN}, {VMAX}, {VCTR} and {VWID} are chosen so as to" \
  " satisfy the equations {VWID = VMAX - VMIN} and" \
  " {VCTR = (VMIN + VMAX)/2}.  Thus, if any two of those" \
  " four arguments are given, the other two are determined" \
  " from these equations."

#endif
