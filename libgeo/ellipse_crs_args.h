#ifndef ellipse_crs_args_H
#define ellipse_crs_args_H

/* ellipse_crs_args.h -- tools to parse ellipse parameters from the command line */
/* Last edited on 2021-06-09 20:27:34 by jstolfi */

#define _GNU_SOURCE
#include <r2.h>
#include <ellipse_crs.h>
#include <argparser.h>

/* PARSING SPHERE OUTLINE SPECS FROM THE COMMAND LINE */

void ellipse_crs_args_parse
  ( argparser_t *pp,   /* The command line parsing state. */
    ellipse_crs_t *EP, /* (OUT) The ellipse parameters. */
    double *ctrAdj,    /* (OUT) The max adjustment to {EP.ctr} coords, or NULL. */
    double *radAdj,    /* (OUT) The max adjustment to {EP.rad}, or NULL. */ 
    double *strAdj     /* (OUT) The max adjustment to {EP.str} coords, or NULL. */
  );
  /* Parses from the command line the geometric parameters
    (center, radius and stretch vector) of an ellipse
    and packs them as an {ellipse_crs_t}.
    
    The syntax is described by {ellipse_crs_args_XXX_HELP} and
    {ellipse_crs_args_XXX_INFO}, where {XXX} is {center}, {radius}, or
    {stretch}. All the parameters of the same ellipse must appear
    together in the command line. See {argparser.h} for an explanation
    of the {pp} parameter.
    
    If any of the keywords "center", "radius", or "stretch" is not specified in
    the command line, the corresponding field of {E} is set to {NAN}.
    This applies also if the keyword appears but is followed immediately
    by "adjust" rather than a numeric value.  
    
    If {ctrAdj} is not NULL, looks for the keyword "adjust" after the
    center coordinates (or after the "center" keyword, if there are
    no coordinates).  If "adjust" is present, the following
    float value {CTR_AJUST} is stored in {*ctrAdj}. If that
    keyword is not present, stores 0 in {ctrAdj}.  If {ctrAdj} 
    is NULL, the modifiers "adjust {AMOUNT}" are not accepted. 
    
    The {radAdj} and {strAdj} provide the same alternative for the
    "radius" and "stretch" keywords, respectively. If these
    keyword-value pairs are given, they must be consecutive unparsed
    arguments, starting at the next argument of {pp}. */
    
void ellipse_crs_args_print(FILE *wr, ellipse_crs_t *E, char *fmt);
  /* Prints the parameters of {E} (without adjustment amounts)
    in a format compatible with {ellipse_crs_args_parse}. */

void ellipse_crs_args_adjust_print
  ( FILE *wr, 
    ellipse_crs_t *E, 
    double ctrAdj, 
    double radAdj, 
    double strAdj,
    char *fmt
  );
  /* Prints the parameters of {E} (and their adjustments, if nonzero)
    in a format compatible with {ellipse_crs_args_parse}. */

/* ---------------------------------------------------------------------- */
/* Documentation for "center": */

#define ellipse_crs_args_center_HELP \
  "center [ {CTRX} {CTRY} ]"

#define ellipse_crs_args_center_adjust_HELP \
  ellipse_crs_args_center_HELP " [ adjust {CTR_ADJUST} ]"
  
#define ellipse_crs_args_center_INFO \
  "Specifies the center of the ellipse on the image."

#define ellipse_crs_args_center_adjust_INFO \
  ellipse_crs_args_center_INFO \
  "  If \"adjust\" is present, {CTR_ADJUST} is the" \
  " maximum amount by which" \
  " {CTRX} and/or {CTRY} may be adjusted."

#define ellipse_crs_args_center_HELP_INFO \
  "      " ellipse_crs_args_center_HELP "\n" \
  "        " ellipse_crs_args_center_INFO 

#define ellipse_crs_args_center_adjust_HELP_INFO \
  "      " ellipse_crs_args_center_adjust_HELP "\n" \
  "        " ellipse_crs_args_center_adjust_INFO

/* ---------------------------------------------------------------------- */
/* Documentation for "radius": */

#define ellipse_crs_args_radius_HELP \
  "radius [ {RAD} ]"
  
#define ellipse_crs_args_radius_adjust_HELP \
  ellipse_crs_args_radius_HELP " [ adjust {RAD_ADJUST} ]"

#define ellipse_crs_args_radius_INFO \
  "Specifies the transverse radius (minor semidiameter)" \
  " of the ellipse."

#define ellipse_crs_args_radius_adjust_INFO \
  ellipse_crs_args_radius_INFO \
  "  If \"adjust\" is present, {RAD_ADJUST} is the maximum amount by which" \
  " the minor radius can be adjusted."

#define ellipse_crs_args_radius_HELP_INFO \
  "      " ellipse_crs_args_radius_HELP "\n" \
  "        " ellipse_crs_args_radius_INFO 
  
#define ellipse_crs_args_radius_adjust_HELP_INFO \
  "      " ellipse_crs_args_radius_adjust_HELP "\n" \
  "        " ellipse_crs_args_radius_adjust_INFO


/* ---------------------------------------------------------------------- */
/* Documentation for "stretch": */

#define ellipse_crs_args_stretch_HELP \
  "stretch [ {STRX} {STRY} ]"

#define ellipse_crs_args_stretch_adjust_HELP \
  ellipse_crs_args_stretch_HELP " [ adjust stretch {STR_ADJUST} ]"

#define ellipse_crs_args_stretch_INFO \
  "Specifies the direction and and amount of stretching" \
  " of the ellipse.  The longest radius (the major" \
  " semidiameter) of the ellipse will be parallel to the vector" \
  " {(STRX,STRY)}, and its length will be the minor" \
  " semidiameter {RAD} plus the length of that" \
  " vector."

#define ellipse_crs_args_stretch_adjust_INFO \
  ellipse_crs_args_stretch_INFO \
  "  If \"adjust\" is present, {STR_ADJUST} is the maximum amount by which" \
  " {STRX} and/or {STRY} can be adjusted."

#define ellipse_crs_args_stretch_HELP_INFO \
  "      " ellipse_crs_args_stretch_HELP "\n" \
  "        " ellipse_crs_args_stretch_INFO 
  
#define ellipse_crs_args_stretch_adjust_HELP_INFO \
  "      " ellipse_crs_args_stretch_adjust_HELP "\n" \
  "        " ellipse_crs_args_stretch_adjust_INFO

#endif
 
