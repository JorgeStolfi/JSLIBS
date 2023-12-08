/* Panpage for {test_r2_opt.c}. */
/* Last edited on 2023-11-25 17:14:28 by stolfi */

#ifndef test_r2_opt_info_H
#define test_r2_opt_info_H

#define _GNU_SOURCE

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {FUNC_NAME} {BIAS} {MONOSCALE} {QUADOPT} \\\n" \
  "  {NI}  {ADJ_RAD_X} {ADJ_RAD_Y}  {ADJ_STEP_X} {ADJ_STEP_Y}  \\\n" \
  "  {CMP_SCALE_X} {CMP_SCALE_Y}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program tests {r2_opt} on an artificial target function selected by the string {FUNC_NAME}.  The argument to be optimized is a list {p} of {NI} 2-vectors.  The target function {F(p)} is constructed so that its optimum is at some argument {popt} chosen by the program.  The optimization starts at some point {pini} near {popt}.\n" \
  "\n" \
  "  The \"indiff\" target function is constant (independent of {p}).\n" \
  "\n" \
  "  The \"optdst\" target function is the square of the Euclidean distance between {p} and {popt}.\n" \
  "\n" \
  "  The \"imgmis\" target function is the squared discrepancy between {NI} images, where each image {J[k]} is evaluated in the neighborhood of {p[k]}.  More precisely, it is the sum of {wx[ix]*wy[iy]\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  ??\n" \
  "\n" \
  "OPTIONS\n" \
  "  ??\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  dm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created sometime in 2020 maybe by J. Stolfi.\n" \
  "MODIFICATION HISTORY\n" \
  "  ??.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " dm_test_180_discrim_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#endif
