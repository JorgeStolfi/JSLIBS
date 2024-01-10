#ifndef bz_plot_2D_H
#define bz_plot_2D_H

/* Plotting of a 2D Bézier patch with fixed resolution (for debugging) */
/* Last edited on 2009-08-30 19:08:17 by stolfi */

#include <pswr.h>
#include <pswr_plot_2D.h>
#include <bool.h>
#include <bz_patch.h>

PSStream *bz_plot_2D_init
  ( char *prefix,
    bool_t epsformat,
    char *paperSize,
    char *caption,
    interval_t bbox[]
  );
  /* Creates a Postscript stream.
  
    If {epsformat} is FALSE, the output will be a stand-alone
    Postscript document called "{prefix}doc.ps", with page dimensions
    defined by the {papersize} string ("letter", "a4", etc.). The
    {caption} will be printed under the figure.
    
    If {epsformat} is TRUE, the ouptut will be an EPS file called
    "{prefix}{NNNNNN}.eps", where {NNNNNN} is a six-digit sequential
    number starting from 0. Each figure will have a nominal size of 
    {8 × 6} inches.
    
    In either case, the scale and clipping is set so that the visible
    range of client plotting coordinates is the rectangle {bbox[0] ×
    bbox[1]}, which is mapped to fill most of the available space,
    with equal scales. */

void bz_plot_2D_finish ( PSStream *ps );
  /* Terminates a Postscript stream that was started with {bz_plot_2D_init}. */

void bz_plot_2D_patch_and_func
  ( PSStream *ps, 
    bz_patch_t *bz,
    bz_patch_rdim_t nf,
    pswr_plot_2D_func_t *func,
    int minRank,
    int maxRank,
    double tol,
    interval_t box[],
    int nx,
    int ny,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  );
  /* Plots a Bézier patch {bz} concatenated with a function {func}
    over the rectangle {B[0] × B[1]}.
    
    The Bézier patch {bz}, if not NULL, must have domain dimension
    {d==bz.m==2} and rank dimension {nb == bz.n >= 2}. So, it is a
    function that, for any argument {x[0..d-1]}, returns a real vector
    {bz(x) = b[0..nb-1]}. If {bz} is NULL, it is assumed to be a
    degenerate patch with {d == bz.m == 2} and {nb = bz.n == 0}.
    
    The function {func} is assumed to take an argument {x[0..1]}
    and return a result vector {func(x) = f[0..nf-1]} with {nf}
    reals. If {func} is NULL, then {nf} must be 0.
    
    The function {T} to be plotted is the concatenation of the values
    {bz(x)} with the values of {func(x)}; namely, takes {x[0..d-1]} to a
    vector {T(x) = t[0..nt-1]} where {nt == nb + nf}, {t[0..nb-1] ==
    b[0..nb-1]}, and {t[nb..nt-1] == f[0..nf-1]}
    
    In any case, the first 2 components {t[0..1]} of {T(x)} are
    interpreted as the plotting coordinates {P(x)} of the point {x}.
    The remaining coordinates {t[2..nt-1]} are plotted using colors
    and/or isolines.
    
    The domain rectangle {B[0..d-1]} is subdivided recursively in half
    along alternating axes, at least {minRank} times. The subdivision
    continues until {maxRank} splits, or until the plot coordinates
    {t[0..d-1]} (which come from the Bézier patch) deviate less than
    {tol} from the multiaffine patch with the same corners. At that
    point the patch is divided into an array of {nx × ny} tiles, each
    tile is divided into 4 triangles with a shared center vertex, and
    the function {T} is plotted as if its values (shape and colors)
    were bilinear in that square.  For the meaning of the parameters
    {st,fill,draw}, see {pswr_plot_2D.h} */

#endif
