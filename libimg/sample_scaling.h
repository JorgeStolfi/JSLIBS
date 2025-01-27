#ifndef sample_scaling_H
#define sample_scaling_H

/* sample_scaling.h -- tools for scaling sample values. */
/* Created 2006-04-20 by Jorge Stolfi, IC-UNICAMP. */
/* Last edited on 2025-01-21 19:31:47 by stolfi */

#include <bool.h>
#include <vec.h>
#include <argparser.h>
#include <float_image.h>

#include <pst_basic.h>

/* PARSING SCALING PARAMETERS FROM THE COMMAND LINE */
  
typedef struct sample_scaling_options_t
  { double_vec_t min;      /* Low endpoint of scaling range; empty vec if not given. */
    double_vec_t max;      /* High endpoint of scaling range; empty vec if not given. */
    double_vec_t ctr;      /* Center of scaling range; empty vec if not given. */
    double_vec_t wid;      /* Width of scaling range; empty vec if not given. */
    bool_t uniform;        /* Obtains default scaling args from whole image rather than single channel. */
  } sample_scaling_options_t;
  /* Sample scaling options as read from the command line. */

#define sample_scaling_options_parse_HELP \
  "    [ -min {VMIN}.. " argparser_double_vec_den_HELP " ] \\\n" \
  "    [ -max {VMAX}.. " argparser_double_vec_den_HELP " ] \\\n" \
  "    [ -center {VCTR}.. " argparser_double_vec_den_HELP " ] \\\n" \
  "    [ -width {VWID}.. " argparser_double_vec_den_HELP " ] \\\n" \
  "    [ -uniform ]"
  
#define sample_scaling_options_parse_HELP_INFO \
  "  -min {VMIN}.. " argparser_double_vec_den_HELP "\n" \
  "  -max {VMAX}.. " argparser_double_vec_den_HELP "\n" \
  "  -center {VCTR}.. " argparser_double_vec_den_HELP "\n" \
  "  -width {VWID}.. " argparser_double_vec_den_HELP "\n" \
  "    These options specify a range [{VMIN} _ {VMAX}] used to scale samples between the {float} values and discrete pixels.\n" \
  "\n" \
  "    If none of these four arguments is present, {VMIN} and {VMAX} will be set" \
  " to {SMIN} and {SMAX}, respectively the minimum and maximum" \
  " sample values in the input.  If only {VCTR} is given," \
  " then {VWID} is the smallest" \
  " width such that {SMIN} and {SMAX} lie in the range" \
  " [{VCTR-VWID/2} _ {VCTR+VWID/2}}.  If only {VMIN}" \
  " is given, then {VMAX} is {max(VMIN,SMAX)}.  If only" \
  " {VMAX} is given, then {VMIN} is {min(SMIN,VMAX)}.  If only" \
  " {VWID} is given, then {VCTR} is {(SMIN+SMAX)/2}\n" \
  "\n" \
  "    In any case, at most two of these four parameters may be specified..\n" \
  "\n" \
  "  -uniform\n" \
  "    This option, if present, specifies that the same sample scaling" \
  " range {[VMIN _ VMAX]} must be used" \
  " in all color channels. If this option is given, any scaling" \
  " parameters (\"-min\", \"-max\", \"-center\", \"-width\") that are" \
  " specified should have either a single value, or the" \
  " same value for all channels. Also, the values of {SMIN} and {SMAX}" \
  " are computed over all relevant channels of the image. Otherwise," \
  " they are computed separately for each channel."

sample_scaling_options_t sample_scaling_parse_options(argparser_t *pp, int32_t *NC_P); 
  /* Parses the options "-min", "-max", "-center", "-width", and "-uniform" from 
    the command line.
    
    If {NC_P} is is NULL, or only one numeric argument is present (with
    optional denominator), ignores {NC_P}. Otherwise, if {*NC_P} is
    negative, sets {*NC_P} to the number of elements read. Otherwise
    demands and parses exactly {*NC_P} numeric arguments (with an
    optional denominator). */

/* UNIFORMIZATION AND COMPLETION OF SCALING PARAMETERS */

void sample_scaling_fix_channels(int32_t NC, int32_vec_t *channel);
  /* Adjusts a list of channel indices to have exactly {NC} elements.
    
    If the {*channel} vector has exactly {NC} channels, the procedure
    does nothing. If {NC} is positive but the given {*channel} vector
    has zero elements, the procedure extends it to {NC} elements,
    which are set to consecutive indices {0,1,...NC-1}. If {NC} is 2
    or more but {*channel} has a single element, the procedure expands
    it to {NC} elements, all equal to element 0. Otherwise the
    procedure fails. */

void sample_scaling_fix_params
  ( sample_scaling_options_t *sop, 
    int32_vec_t *channel, 
    int32_t NC,
    float_image_t *fim
  );
  /* Adjusts the channel selection {channel.e[0..channel.ne-1]} and scaling 
    parameters {sop}, as may be given in the command line, for {NC} channels.
    
    Each of the input parameter vectors {min,max,ctr,wid}
    must have either zero elements, or one element, or {NC} elements.
    
    If a parameter is unspecified (represented by the value {±INF} of
    {NAN}), the procedure either sets it to a suitable default value, or
    computes it from other known parameters. If a single value was given
    for all channels, the procedure replicates it for all channels.
    
    After the procedure returns with success, the parameter vectors
    {channel} and {sop.{min,max,ctr,wid}} will have exactly {NC} elements, no
    elements will be {±INF} or {NAN}, and the parameters of each channel will be
    consistent with each other. 
    
    If there are less than two known parameters for each channel, and
    {fim} is not NULL, the procedure tries to use the minimum and
    maximum sample values {smin,smax} from the image {*fim} to provide
    suitable defaults, as explained in
    {sample_scaling_options_parse_HELP_INFO}. In that case, if
    {channel} is not NULL, the procedure assumes that the channels to
    analyze are {channel.e[0..NC-1]}. If {channel} is NULL, channels
    {0,1,...NC-1} are analyzed.
    
    In any case, once there are exactly two known parameter values for
    each channel, the other two are computed with
    {sample_scaling_complete_params}. The procedure fails if, for any
    channel, more than two parameters are defined.
    
    If {sop.uniform} is TRUE, forces the scaling parameters to be the same
    across all channels. In that case, all input vectors which are not
    empty or trivial must contain {NC} equal values; and, if the
    actual extreme samples {smin,smax} are needed, they are computed
    once over all relevant channels of {fim}, rather than separately
    for each channel. */ 

#endif
