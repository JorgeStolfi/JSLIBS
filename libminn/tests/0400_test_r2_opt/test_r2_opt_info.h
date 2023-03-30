/* Panpage for {test_r2_opt.c}. */
/* Last edited on 2017-06-05 22:57:55 by stolfilocal */

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
  "  The program then trims each of the candidates to various" \
  " lengths, obtaining a set of nested sub-candidates.  Each of these" \
  " sub-candidates is optionally refined (once), and its length and score" \
  " are written out to disk.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "  The main output files are called \"{OUT_NAME}-{LL}-{K}.dat\" where" \
  " {LL} is a two digit filtering level, and {K} indicates the kind" \
  " of pairing (\"t\" for true pairings, \"f\" for false).  In these files," \
  " each line contains the following fields:\n" \
  "\n" \
  "    {NX} number of positions in the X sequence that are spanned by the sub-candidate.\n" \
  "    {NY} number of positions in the Y sequence that are spanned by the sub-candidate.\n" \
  "    {RSCORE} the score of the sub-candidate, with refining if requested.\n" \
  "    {USCORE} the score of the sub-candidate before refining.\n" \
  "\n" \
  "  The data for all sub-candidates of the same candidate appear" \
  " in consecutive lines.  For the sake of {gnuplot}, blank lines separate" \
  " data from different candidates.\n" \
  "\n" \
  "  The program also writes out the candidates for each" \
  " level {LL} and eack kind {K}, to a file called \"{OUT_NAME}-{LL}-{K}.cdv\"," \
  " using the procedure {msm_cand_vec_write}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -seqA {ID_A} {NAME_A} {FILENAME_A}\n" \
  "  -seqB {ID_B} {NAME_B} {FILENAME_B}\n" \
  "    These  mandatory arguments specify the" \
  " file names where the sequences are to be read from.  The files must be in the" \
  " format used by {dnae_seq_read_from_nucleic_file}: one line of header, starting" \
  " with \">\", and then one or more lines wit the nucleotide sequence, using" \
  " letters [ATCG].  The sequences will have internal names {NAME_A.NAME_B} and" \
  " internal ID numbers {ID_A,ID_B}; these are used only for debugging.\n" \
  "\n" \
  "  -candsTrue {CDV_FILE_TRUE}\n" \
  "  -candsFalse {CDV_FILE_FALSE}\n" \
   "    These  mandatory arguments specify the" \
  " names of the files with the `true' and `false' pairings.  The files must be in the format" \
  " written by {msm_cand_vec_write}.\n" \
  "\n" \
  "  -minLevel {MIN_LEVEL}\n" \
  "  -maxLevel {MAX_LEVEL}\n" \
  "    These  mandatory arguments specify the" \
  " minimum and maximum level of filtering to be applied.  Level 0 is the" \
  " unfiltered sequence, and each subsequent level is sampled at half" \
  " the frequency of the previous one.\n" \
  "\n" \
  "  -maxCands {MAX_CANDS}\n" \
  "    This optional argument, if given, causes the input candidate files to be" \
  " truncated to a maximum of {MAX_CANDS} candidates in each.  If not" \
  " specified, or zero, all candidates are used.\n" \
  "\n" \
  "  -minSamples {MIN_SAMPLES}\n" \
  "    This mandatory argument specifies the minimum number of actual" \
  " samples (not just sequence positions) in each sequence and also" \
  " the minimum number of rungs in each candidate, at any scale, which" \
  " are needed for meaningful analysis.  The maximum filtering" \
  " level at which the test is performed is such that at" \
  " least one candidate with {MIN_SAMPLES} rungs remains.\n" \
  "\n" \
  "  -initFilter " wt_table_args_HELP " \n" \
  "  -incrFilter " wt_table_args_HELP " \n" \
  "    These mandatory arguments specify the" \
  " weights of the filter to be used at the" \
  " first filtering step and at subsequent" \
  " filtering steps.  " wt_table_args_norm_sum_INFO "\n" \
  "\n" \
  "  -refineFull\n" \
  "    This optional argument specifies that each candidates must" \
  " be refined (once) before extracting the sub-candidates.\n" \
  "\n" \
  "  -refineSub\n" \
  "    This optional argument specifies that each sub-candidate must be refined.\n" \
  "\n" \
  "  -scores " dm_score_args_HELP "\n" \
  "    This mandatory argument specifies the" \
  " scoring of steps in a pairing.  " dm_score_args_INFO "\n" \
  "\n" \
  "  -delta {DELTA}\n" \
  "    This optional argument specifies the" \
  " amount of adjustment allowed for the X and Y coordinates of" \
  " internal rungs of each pairing.  The default is 3.\n" \
  "\n" \
  "  -kappa {KAPPA}\n" \
  "    This optional argument specifies the" \
  " amount of X and Y extension allowed" \
  " at either end of each pairing.  The default is 6.\n" \
  "\n" \
  "  -maxUnp {MAX_UNP}\n" \
  "    This optional argument specifies the" \
  " maximum unpaired samples at each step used for candidate refining.\n" \
  "\n" \
  "  -expand {EXPAND}\n" \
  "    This optional argument specifies by how much the {R}-range" \
  " of a refined candidate may extend beyond the" \
  " original candidate's {R}-range, in both directions.  If omitted, or" \
  " if {EXPAND} is zero, the refined" \
  " {R}-range will be a subset of the orginal range.\n" \
  "\n" \
  "  -shrink {SHRINK}\n" \
  "    This optional argument specifies by how much the {R}-range of" \
  " the refined candidate may shrink into the original" \
  " candidate's {R}-range, in each direction.  If omitted, or" \
  " if {SHRINK} is zero, the refined {R}-range" \
  " will be a superset of the original range.\n" \
  "\n" \
  "  -lenStep {LENGTH_STEP}\n" \
  "    This mandatory parameter specifies the lenght increment of the analyzed" \
  " candidate pieces at level 0.  The pieces at higher levels will be shorter in" \
  " proportion to the sequence lengths.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  dm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 21/dec/2006 by J. Stolfi.\n" \
  "MODIFICATION HISTORY\n" \
  "  2015-11-02 J. Stolfi:\n" \
  "    * Truncates to {maxCands} candidates only after checking size.\n" \
  "    * Option \"-stepSamples\" renamed \"-lenStep\" and scaled according to level.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " dm_test_180_discrim_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#endif
