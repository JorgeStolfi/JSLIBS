#define PROG_NAME "pnmift"
#define PROG_DESC "compute the Image Foresting Transform (IFT) of a PBM/PGM/PPM image"
#define PROG_VERS "1.0"

#define pnmift_C_COPYRIGHT "Copyright © 2001 by the State University of Campinas (UNICAMP)"

/* Last edited on 2007-12-26 20:31:59 by stolfi */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  -radius RADIUS -costFunction COSTFN \\\n" \
  "  [ -seedImage PGMFILE ] [ -seedPixels TXTFILE ] \\\n" \
  "  { -output  PREFIX | \\\n" \
  "    {  -costs | -preds | -labels | -extract LABEL } OUTPGM \\\n" \
  "  } \\\n" \
  "  [ -roots  OUTPGM ] [ -spread OUTPNM ] \\\n" \
  "  [ -boxes OUTTXT [ -margin MRG ]] \\\n" \
  "  [ -plot OUTEPS ]  \\\n" \
  "  [ -bgColor {VALUE} {VALUE} {VALUE} | -bgGray {VALUE} ] \\\n" \
  "  [ -maxCostValue MAXCOST ] [ -maxCostPixel MAXPIXEL ] \\\n" \
  "  [ -lifo | -fifo ]  [ -reverse ] \\\n" \
  "  [ -verbose ] [ -listFunctions ] \\\n" \
  "  [ PNMFILE ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPT1 "\n" \
  "\n" \
  PROG_INFO_OPT2 "\n" \
  "\n" \
  PROG_INFO_OPT3 "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgm(5), ppm(5), pnm(5), pgmselect(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in sep/2001 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  IFT concept and algorithm by A.X.Falcão, R.Lotufo and J.Stolfi." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmift_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program Reads a portable anymap {PNMFILE} and computes" \
  " its Image Foresting Transform (IFT), for a given set of seed" \
  " pixels and a specified path cost function.\n" \
  "\n" \
  "  The seed pixels can be specified in two ways: as a greyscale" \
  " image, or through a text file.  If neither method is used," \
  " every pixel is assumed to be a seed.\n" \
  "\n" \
  "  The IFT is written as three separate graymaps," \
  " containing the label, predecessor, and cost of each pixel," \
  " respectively.  Each image can be requested and named separately" \
  " with the \"-cost\", \"-pred\", and \"-label\" directives.  Alternatively," \
  " one may specify a common filename prefix as the \"-output {PREFIX}\" parameter;" \
  " the three images will be written out as \"{PREFIX}-costs.pgm\"," \
  " \"{PREFIX}-preds.pgm\", and \"{PREFIX}-labels.pgm\".\n" \
  "\n" \
  "  Optionally, the program may also writes a binary image showing" \
  " the root pixels (see the \"-roots\" directive); and/or a copy" \
  " of the original image where each root pixel has been replicated" \
  " throughout its influence zone (see the \"-spread\" directive); and/or" \
  " a text file with the bounding box of each region (see the \"-boxes\" directive)."

#define PROG_INFO_OPT1 \
  "  -seedImage {PGMFILE}\n" \
  "    where {PGMFILE} is the name of a a grayscale image file in" \
  " the \"pgm\" format, where the non-zero pixels are the seeds," \
  " and their values are the respective labels.\n" \
  "\n" \
  "  -seedPixels {TEXTFILE}\n" \
  "    Specifies the name of an ascii file where each line contains" \
  " the H and V coordinates of a seed pixel, and the corresponding" \
  " label (an integer in the range [0..2^16-1]).  These labels are" \
  " set after the seed image, if any.  If both \"-seedImage\" and" \
  " \"-seedPixels\" are specified, the {PGMFILE} is read first, then" \
  " the {TEXTFILE} is applied over it.  If neither \"-seedImage\"" \
  " nor \"-seedPixels\"is specified, every pixel is considered" \
  " to be a seed, with label 1.\n" \
  "\n" \
  "  -radius {NUM}\n" \
  "    Specifies the radius of the neighborhood to use in the" \
  " propagation algorithm.  Use 1.0 for 4-neighbor topology, 1.5 for" \
  " 8-neighbor.  This parameter is required, and must be within 1.0" \
  " and 127.0, inclusive."

#define PROG_INFO_OPT2 \
  "  -costFunction {NAME}\n" \
  "    The name of the path-cost function to use.  The valid names are:\n" \
  ift_functions_INFO "\n" \
  "\n" \
  "  -output {PREFIX}\n" \
  "    The common prefix for output file names.  Use either this" \
  " option, or the options \"-costs\", \"-labels\", \"-preds\".\n" \
  "\n" \
  "  -costs {PGMFILE}\n" \
  "    Specifies the name of the output image file that contains" \
  " the pixel costs.\n" \
  "\n" \
  "  -preds {PGMFILE}\n" \
  "    Specifies the name of the output image file that contains" \
  " the predecessors of the pixels in the minimum-path forest.  Each" \
  " predecessor is given as two relative displacements (dh,dv), each" \
  " a signed integer in [-127..+127], packed into a 16-bit pixel" \
  " value 256*(dv+128) + (dh+128). \n" \
  "\n" \
  "  -labels {PGMFILE}\n" \
  "    Specifies the name of the output image file that contains" \
  " the propagated pixel labels.  Each label is stored as" \
  " a 16-bit pixel value.  A zero label marks those pixels" \
  " that cannot be reached from any seed at finite path cost.\n" \
  "\n" \
  "  -roots {PGMFILE}\n" \
  "    Causes the program to write the named binary output" \
  " image image (a graymap with maxval=1), showing which pixels" \
  " are roots (value 1) or non-roots (value 0).\n" \
  "\n" \
  "  -spread {PNMFILE}\n" \
  "    Causes the program to write an image of the same type and" \
  " depth as the original, where every pixel has been set to the" \
  " color of its corresponding root pixel. \n" \
  "\n" \
  "  -extract {LABEL} {PNMFILE}\n" \
  "    Causes the program to write an image of the same type and" \
  " depth as the original, where every pixel of class {LABEL} is copied" \
  " from the original image, while all other pixels are set to the" \
  " background pixel value. \n" \
  "\n" \
  "  -boxes {TXTFILE}\n" \
  "    Causes the program to write a text file containg one" \
  " line for each possible label, with the bounding box" \
  " of the corresponding region, in the format\n" \
  "\n" \
  "      {LABEL} {HMIN} {VMIN} {HSIZE} {VSIZE}\n" \
  "\n" \
  "    If the region is empty, the last four fields are zero."

#define PROG_INFO_OPT3 \
  "  -margin {MARGIN}\n" \
  "    Used in conjuntion with \"-boxes\".  Expands the bounding" \
  " box of each non-empty region by the given number of pixels, in" \
  " all four directions, then clips the result against the image bounds.\n" \
  "\n" \
  "  -plot {EPSFILE}\n" \
  "    Causes the program to write an Encapsulated Postscript version" \
  " of the image, overlaid with a drawing of the minimum-path forest.\n" \
  "\n" \
  "  -bgColor {VALUE} {VALUE} {VALUE}\n" \
  "  -bgGray {VALUE}\n" \
  "    These options specify the background pixel value to use for" \
  " the \"-extract\" option.  Only one of the two should be" \
  " specified.  If \"-bgColor\" is used for gray" \
  " image input, only the first {VALUE} is used.  The default" \
  " is \"-bgColor 0 0 0\" = \"-bgGray 0\".\n" \
  "\n" \
  "  -maxCostPixel {MAXPIXEL}\n" \
  "    Specifies the maximum pixel value to be used when writing the" \
  " cost layer of the IFT.  Defaults to 2^16-1 = 65535.\n" \
  "\n" \
  "  -maxCostValue {MAXCOST}\n" \
  "    Specifies the maximum nominal path cost, for the purpose of" \
  " output scaling.  The units are the same used internally by the" \
  " IFT procedure (see {ift.h}). Defaults to the maximum path cost" \
  " observed in the IFT.  When writing the cost layer, infinte path" \
  " costs will be mapped to MAXPIXEL.  Finite path costs in" \
  " [0..MAXCOST] will be mapped linearly to [0..{MAXPIXEL-1} ].   Finite" \
  " costs higher than {MAXCOST} will be mapped to {MAXPIXEL-1}.\n" \
  "\n" \
  "  -lifo\n" \
  "    Specifies LIFO precedence for selecting between equal-cost paths.\n" \
  "\n" \
  "  -fifo\n" \
  "    Specifies FIFO precedence for selecting between equal-cost paths. \n" \
  "\n" \
  "  -reverse\n" \
  "    Reverses the scanning order of seeds and outgoing arcs.\n" \
  "\n" \
  "  -verbose\n" \
  "    Produces diagnostic output."

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>

#include "ift.h"
#include "ift_functions.h"
#include "ift_plot.h"
#include "ift_image.h"

#include <bool.h>
#include <affirm.h>
#include <argparser.h>
#include <jspnm_image.h>
#include <jspnm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <pswr.h>

/* PROTOTYPES */

typedef struct options_t
  { char *input_img_name;
    char *seed_img_name;
    char *seed_txt_name;
    double radius;
    PathCost max_cost_value;
    pnm_sample_t max_cost_pixel;
    pnm_sample_t bgColor[3];
    PathCostFn *pf;
    ift_tie_breaking tbreak;
    ift_scan_order order;
    char *cost_img_name;
    char *pred_img_name;
    char *label_img_name;
    char *root_img_name;
    char *spread_img_name;
    char *extract_img_name;
    SeedLabel extract_label;
    char *region_boxes_name;
    int box_margin;
    char *plot_eps_name;
    int verbose;
  } options_t;
  /* Command line arguments. */

int main(int argc, char* argv[]);

options_t *pnmift_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

ImageGraph *pnmift_read_image_as_graph(options_t *o, pnm_sample_t *maxvalp);
  /* Reads the PNM file "{o->input_img_name}" and transforms it into
    an IFT graph. Also updates {*input_maxvalp} to be at least as big
    as the nominal maximum value of any sample in the image (the PNM
    {maxval} parameter). */

void pnmift_set_seeds_from_image_file(options_t *o, ImageGraph *G, pnm_sample_t *maxvalp);
  /* Sets the seeds from the PNM image file named
    "{o->seed_img_name}". Also updates {*maxvalp} to to be at least as
    big as the nominal maximum value of any sample in the image (the
    PNM {maxval} parameter). */

void pnmift_set_seeds_from_text_file(options_t *o, ImageGraph *G, pnm_sample_t *maxlabelp);
  /* Sets the seeds from the text file named "{o->seed_txt_name}".  Also
    sets {*maxlabelp} to the maximum label value read from the file. */

void pnmift_write_cost_image(options_t *o, ImageGraph *G, PathCost actual_max_cost);
void pnmift_write_label_image(options_t *o, ImageGraph *G, SeedLabel maxlabel);
void pnmift_write_pred_image(options_t *o, ImageGraph *G);
void pnmift_write_root_image(options_t *o, ImageGraph *G);
void pnmift_write_spread_image(options_t *o, ImageGraph *G, pnm_sample_t maxval);
void pnmift_write_extract_image(options_t *o, ImageGraph *G, pnm_sample_t maxval);
  /* These procedure construct various output images (by the corresponding procedures
     in {ift.h}) and write them to files called "{o->cost_img_name}",
     "{label_img_name}", etc.. */

void pnmift_write_region_bboxes(options_t *o, ImageGraph *G, SeedLabel maxlabel);
  /* Writes a text file containing a rectangular bounding box for each
    region of {G} in {0..maxlabel}.  See {ift_write_boxes} in {ift.h}. */

void pnmift_write_postscript_plot(options_t *o, ImageGraph *G);
  /* Writes an Encapsulated Postscript file with a plot of the
    spanning forest contained in {G}. */

void skip_rest_of_line(FILE *f);

/* REAL CODE */

int main(int argc, char** argv)
  {
    /* Parse command line arguments. */
    options_t *o = pnmift_parse_options(argc, argv);

    pnm_sample_t input_maxval;       /* Max sample value of input image. */
    pnm_sample_t seed_maxlabel = 0;  /* Maximum seed label, for output. */
    PathCost actual_max_cost;        /* Actual maximum cost observed in the IFT. */

    /* Get input image, and turn it into a graph {G}: */
    ImageGraph *G = pnmift_read_image_as_graph(o, &input_maxval); /* Graph of input image. */

    /* Mark the seed pixels given as a PNM image: */
    ift_set_all_seed_labels(G, NON_SEED_LABEL);
    if (o->seed_img_name != NULL)
      { pnmift_set_seeds_from_image_file(o, G, &(seed_maxlabel)); }

    /* Mark the seed pixels given as a text file: */
    if (o->seed_txt_name != NULL)
      { if (o->verbose) { fprintf(stderr, "reading seed file %s...\n", o->seed_txt_name); }
        pnmift_set_seeds_from_text_file(o, G, &(seed_maxlabel));
      }

    /* If no seed pixels are specified, mark all pixels as seeds: */
    if ((o->seed_img_name == NULL) && (o->seed_txt_name == NULL))
      { if (o->verbose) { fprintf(stderr, "setting all seed labels to %d...\n", DEFAULT_SEED_LABEL); }
        ift_set_all_seed_labels(G, DEFAULT_SEED_LABEL);
      }

    /* Compute the IFT: */
    if (o->verbose) { fprintf(stderr, "computing the IFT...\n"); }
    ift_compute_forest(G, o->pf, o->tbreak, o->order, &actual_max_cost);

    /* Write the output images: */

    /* Cost layer: */
    if (o->cost_img_name != NULL)
      { pnmift_write_cost_image(o, G, actual_max_cost); }

    /* Label layer: */
    if (o->label_img_name != NULL)
      { pnmift_write_label_image(o, G, seed_maxlabel); }

    /* Predecessor layer: */
    if (o->pred_img_name != NULL)
      { pnmift_write_pred_image(o, G); }

    /* Roots indicator file: */
    if (o->root_img_name != NULL)
      { pnmift_write_root_image(o, G); }

    /* Root-spread image: */
    if (o->spread_img_name != NULL)
      { pnmift_write_spread_image(o, G, input_maxval); }

    /* Single-label image: */
    if (o->extract_img_name != NULL)
      { pnmift_write_extract_image(o, G, input_maxval); }

    /* Region bounding boxes: */
    if (o->region_boxes_name != NULL)
      { pnmift_write_region_bboxes(o, G, seed_maxlabel); }

    /* Postscript plot: */
    if (o->plot_eps_name != NULL)
      { pnmift_write_postscript_plot(o, G); }

    /* All done: */
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    exit(0);
  }

ImageGraph *pnmift_read_image_as_graph(options_t *o, pnm_sample_t *input_maxvalp)
  { pnm_image_t *input_img;
    if (o->verbose) { fprintf(stderr, "reading file %s...\n", o->input_img_name); }
    input_img = pnm_image_read(o->input_img_name);
    (*input_maxvalp) = input_img->maxval;
    if (o->verbose) { fprintf(stderr, "building image graph...\n"); }
    ImageGraph *G = ift_make_graph(input_img->cols, input_img->rows, o->radius);
    if (o->verbose) { fprintf(stderr, "cols = %d rows = %d\n", G->cols, G->rows); }
    if (o->verbose) { fprintf(stderr, "channels = %d\n", G->channels); }
    if (o->verbose) { fprintf(stderr, "nodes = %ld arcs/node = %d\n", G->nodes, G->arcs); }
    ift_set_values_from_image(input_img, G);
    pnm_image_free(input_img);
    return G;
  }

void pnmift_set_seeds_from_image_file(options_t *o, ImageGraph *G, pnm_sample_t *maxlabelp)
  { pnm_image_t *seed_img;
    if (o->verbose) { fprintf(stderr, "reading seed image %s...\n", o->seed_img_name); }
    seed_img = pnm_image_read(o->seed_img_name);
    if (seed_img->chns != 1) { pnm_error("seed image must be monochromatic"); }
    if ((*maxlabelp) < seed_img->maxval) { (*maxlabelp) = seed_img->maxval; }
    if (o->verbose) { fprintf(stderr, "marking seeds in graph...\n"); }
    ift_set_seeds_from_image(seed_img, G);
    pnm_image_free(seed_img);
  }

void pnmift_set_seeds_from_text_file(options_t *o, ImageGraph *G, pnm_sample_t *maxlabelp)
  {
    FILE *f = open_read(o->seed_txt_name, o->verbose);
    int maxlabel, retcode, row, col, label, ch;
    maxlabel = 1;
    ch = fgetc(f);
    while (ch != EOF)
      { if ((ch == ' ') || (ch == '\n') || (ch == '\011') || (ch == '\015'))
          { /* Ignore. */ }
        else if (ch == '#')
          { while ((ch != '\n') && (ch != EOF)) { ch = fgetc(f); } }
        else if (ch == '(')
          { retcode = fscanf(f, " %d %d ) = %d", &col, &row, &label);
            if (retcode != 3) { IFT_ERROR("format error in seed file"); }
            if ((row < 0) || (row >= G->rows)) { IFT_ERROR("bad row in seed file"); }
            if ((col < 0) || (col >= G->cols)) { IFT_ERROR("bad col in seed file"); }
            if ((label < 0) || (label > MAX_LABEL)) { IFT_ERROR("bad label in seed file"); }
            ift_set_seed_label(G, col, row, label);
            if (label > maxlabel) { maxlabel = label; }
            skip_rest_of_line(f);
          }
        else
          { fprintf(stderr, "ch = «%c»\n", ch);
            IFT_ERROR("unexpected char in seed file");
          }
        if (ch != EOF) { ch = fgetc(f); }
      }
    fclose(f);
    if ((*maxlabelp) < maxlabel) { (*maxlabelp) = maxlabel; }
  }

void pnmift_write_cost_image(options_t *o, ImageGraph *G, PathCost actual_max_cost)
  { pnm_image_t *cost_img;
    if (o->verbose) { fprintf(stderr, "extracting the cost image...\n"); }
    fprintf(stderr, "maximum observed path cost = %.0f (%g)\n",
      actual_max_cost, actual_max_cost/((double)MAX_FINITE_ARC_COST)
    );
    if (o->max_cost_value == 0) { o->max_cost_value = actual_max_cost; }
    if (o->verbose)
      { fprintf
         ( stderr, "scaling [0..%.0f] ([0 _ %g]) to [0..%d], +oo to %d\n",
           o->max_cost_value,
           o->max_cost_value/((double)MAX_FINITE_ARC_COST),
           o->max_cost_pixel-1,
           o->max_cost_pixel
         );
      }
    cost_img = ift_get_cost_image(G, o->max_cost_value, o->max_cost_pixel);
    if (o->verbose) { fprintf(stderr, "writing the cost image %s...\n", o->cost_img_name); }
    pnm_image_write(o->cost_img_name, cost_img, 0);
    pnm_image_free(cost_img);
  }

void pnmift_write_label_image(options_t *o, ImageGraph *G, SeedLabel maxlabel)
  { pnm_image_t *label_img;
    if (o->verbose) { fprintf(stderr, "extracting the label image...\n"); }
    label_img = ift_get_label_image(G, maxlabel);
    if (o->verbose) { fprintf(stderr, "writing the label image %s...\n", o->label_img_name); }
    pnm_image_write(o->label_img_name, label_img, 0);
    pnm_image_free(label_img);
  }

void pnmift_write_pred_image(options_t *o, ImageGraph *G)
  {
    if ((int)PNM_FILE_MAX_MAXVAL < 65535)
      { fprintf(stderr, "PNM_FILE_MAX_MAXVAL = %d too small - pred image not generated\n", PNM_FILE_MAX_MAXVAL); }
    else
      {
        pnm_image_t *pred_img;
        if (o->verbose) { fprintf(stderr, "extracting the pred image...\n"); }
        pred_img = ift_get_pred_image(G);
        if (o->verbose) { fprintf(stderr, "writing the pred image %s...\n", o->pred_img_name); }
        pnm_image_write(o->pred_img_name, pred_img, 0);
        pnm_image_free(pred_img);
      }
  }

void pnmift_write_root_image(options_t *o, ImageGraph *G)
  { pnm_image_t *root_img;
    if (o->verbose) { fprintf(stderr, "generating the root mask...\n"); }
    root_img = ift_get_root_image(G);
    if (o->verbose) { fprintf(stderr, "writing the root mask %s...\n", o->root_img_name); }
    pnm_image_write(o->root_img_name, root_img, 0);
    pnm_image_free(root_img);
  }

void pnmift_write_spread_image(options_t *o, ImageGraph *G, pnm_sample_t maxval)
  { pnm_image_t *spread_img;
    if (o->verbose) { fprintf(stderr, "extracting the spread image (maxval = %d)...\n", maxval); }
    spread_img = ift_get_spread_image(G, maxval);
    if (o->verbose) { fprintf(stderr, "writing the spread image %s...\n", o->spread_img_name); }
    pnm_image_write(o->spread_img_name, spread_img, 0);
    pnm_image_free(spread_img);
  }

void pnmift_write_extract_image(options_t *o, ImageGraph *G, pnm_sample_t maxval)
  { pnm_image_t *extract_img;
    if (o->verbose) { fprintf(stderr, "extracting the single_label image (maxval = %d)...\n", maxval); }
    extract_img = ift_get_single_label_image(G, o->extract_label, o->bgColor, maxval);
    if (o->verbose) { fprintf(stderr, "writing the single_label image %s...\n", o->extract_img_name); }
    pnm_image_write(o->extract_img_name, extract_img, 0);
    pnm_image_free(extract_img);
  }

void pnmift_write_region_bboxes(options_t *o, ImageGraph *G, SeedLabel maxlabel)
  { FILE *wr = open_write(o->region_boxes_name, o->verbose);
    if (o->verbose) { fprintf(stderr, "writing the bounding boxes...\n"); }
    ift_write_boxes(wr, G, maxlabel, o->box_margin);
    fclose(wr);
  }

void pnmift_write_postscript_plot(options_t *o, ImageGraph *G)
  { double hsize = 12.0*(double)G->cols;
    double vsize = 12.0*(double)G->rows;
    if (o->verbose) { fprintf(stderr, "generating the tree plot %s...\n", o->plot_eps_name); }
    FILE *psfile = open_write(o->plot_eps_name, o->verbose);
    PSStream *ps = pswr_new_stream(TRUE, NULL, psfile, "letter", hsize + 4.0, vsize + 4.0);
    pswr_new_canvas(ps, "trees");
    pswr_set_window
      ( ps,
        0.0, (double)G->cols, 0.0, (double)G->rows,
        2.0, 2.0 + hsize, 2.0, 2.0 + vsize
      );
    pswr_set_grid(ps, G->cols, G->rows);
    ift_plot_pixel_values(ps, G, 0.5);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.5, 0.0, 0.0);
    ift_plot_forest_edges(ps, G);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
    ift_plot_forest_nodes(ps, G, 1, 1.0,  1.0,1.0,1.0);
    ift_plot_forest_nodes(ps, G, 0, 0.5,  0.0,0.0,0.0);
    pswr_close_stream(ps);
  }

options_t *pnmift_parse_options(int argc, char **argv)
  {
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");

    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    if (argparser_keyword_present(pp, "-seedImage"))
      { o->seed_img_name = argparser_get_next(pp); }
    else
      { o->seed_img_name = NULL; }

    if (argparser_keyword_present(pp, "-seedPixels"))
      { o->seed_txt_name = argparser_get_next(pp); }
    else
      { o->seed_txt_name = NULL; }

    if (argparser_keyword_present(pp, "-radius"))
      { o->radius = argparser_get_next_double(pp, 1.0, 127.0); }
    else
      { argparser_error(pp, "must specify a radius"); }

    if (argparser_keyword_present(pp, "-costFunction"))
      { char *pf_str = argparser_get_next(pp);
        o->pf = get_cost_function(pf_str);
      }
    else
      { argparser_error(pp, "must specify a cost function"); }

    if (argparser_keyword_present(pp, "-output"))
      { /* Define matched filenames for costs, labels, and preds images: */
        char *output_str = argparser_get_next(pp);
        o->cost_img_name = txtcat(output_str, "-costs.pgm");
        o->label_img_name = txtcat(output_str, "-labels.pgm");
        o->pred_img_name = txtcat(output_str, "-preds.pgm");
      }
    else
      { /* Allow individual filenames for costs, labels, and preds images: */
        if (argparser_keyword_present(pp, "-costs"))
          { o->cost_img_name = argparser_get_next(pp); }
        else
          { o->cost_img_name = NULL; }

        if (argparser_keyword_present(pp, "-preds"))
          { o->pred_img_name = argparser_get_next(pp); }
        else
          { o->pred_img_name = NULL; }

        if (argparser_keyword_present(pp, "-labels"))
          { o->label_img_name = argparser_get_next(pp); }
        else
          { o->label_img_name = NULL; }
      }

    if (argparser_keyword_present(pp, "-roots"))
      { o->root_img_name = argparser_get_next(pp); }
    else
      { o->root_img_name = NULL; }

    if (argparser_keyword_present(pp, "-spread"))
      { o->spread_img_name = argparser_get_next(pp); }
    else
      { o->spread_img_name = NULL; }

    if (argparser_keyword_present(pp, "-extract"))
      { o->extract_label = argparser_get_next_int(pp, 0, INT_MAX);
        o->extract_img_name = argparser_get_next(pp);
      }
    else
      { o->extract_label = 0;
        o->extract_img_name = NULL; 
      }

    if (argparser_keyword_present(pp, "-boxes"))
      { o->region_boxes_name = argparser_get_next(pp); }
    else
      { o->region_boxes_name = NULL; }

    if
      ( (o->cost_img_name == NULL) &&
        (o->pred_img_name == NULL) &&
        (o->label_img_name == NULL) &&
        (o->root_img_name == NULL) &&
        (o->spread_img_name == NULL) &&
        (o->region_boxes_name == NULL)
      )
      { fprintf(stderr, "what, no output files?\n"); }

    if (argparser_keyword_present(pp, "-margin"))
      { o->box_margin = argparser_get_next_int(pp, 0, INT_MAX); }
    else
      { o->box_margin = 0; }

    if (argparser_keyword_present(pp, "-plot"))
      { o->plot_eps_name = argparser_get_next(pp); }
    else
      { o->plot_eps_name = NULL; }

    if (argparser_keyword_present(pp, "-maxCostValue"))
      { o->max_cost_value = argparser_get_next_int(pp, 0, INT_MAX); }
    else
      { o->max_cost_value = 0; /* 0 means actual max */ }

    if (argparser_keyword_present(pp, "-bgColor"))
      { int ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = argparser_get_next_int(pp, 0, INT_MAX); }
      }
    else if (argparser_keyword_present(pp, "-bgGray"))
      { pnm_sample_t val = argparser_get_next_int(pp, 0, INT_MAX);
        int ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = val; }
      }
    else
      { int ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = 0; }
      }

    if (argparser_keyword_present(pp, "-maxCostPixel"))
      { o->max_cost_pixel = argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL); }
    else
      { o->max_cost_pixel = PNM_FILE_MAX_MAXVAL; }

    if (argparser_keyword_present(pp, "-lifo"))
      { o->tbreak = ift_tbreak_LIFO; }
    else if (argparser_keyword_present(pp, "-fifo"))
      { o->tbreak = ift_tbreak_FIFO; }
    else
      { o->tbreak = ift_tbreak_FIFO; }

    if (argparser_keyword_present(pp, "-reverse"))
      { o->order = ift_order_DN; }
    else
      { o->order = ift_order_UP; }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Get optional input file name: */
    argparser_skip_parsed(pp);
    if (argparser_next(pp) != NULL)
      { o->input_img_name = argparser_get_next(pp); }
    else
      { o->input_img_name = "-"; }

    /* Check for extraneous arguments: */
    argparser_finish(pp);

    return o;
  }

void skip_rest_of_line(FILE *f)
  {
    /* Skip blanks/comments until end of line, complain if anything else: */
    int ch;
    do
      { ch = fgetc(f);
      if (ch == EOF) { IFT_ERROR("missing newline in seed file"); }
      if (ch == '#')
  	{
  	  do { ch = fgetc(f);
  	  if (ch == EOF) { IFT_ERROR("missing newline in seed file"); }
  	  } while (ch != '\n');
  	}
      if ((ch != ' ') && (ch != '\n') && (ch != '\011') && (ch != '\015'))
  	{ IFT_ERROR("format error in seed file"); }
      }
    while (ch != '\n');
  }
