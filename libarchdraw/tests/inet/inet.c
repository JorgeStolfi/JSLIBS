/* Diagram of the internet */
/* Last edited on 2012-12-07 21:02:18 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>

#include <archdraw.h>

typedef enum 
  { 
    inet_block_type_BKB, /* Backbone (polygon). */
    inet_block_type_BKC, /* Backbone connection (dot). */
    inet_block_type_DNS, /* DNS server (polygon). */
    inet_block_type_GOO, /* Google server (polygon). */
    inet_block_type_BAN, /* Bank server (polygon). */
    inet_block_type_ETC  /* Etc. (must be last). */
  } inet_block_type_t;
  /* Code for sub-diagram type. */
  
#define inet_block_type_MAX (inet_block_type_ETC)
  /* The largest {block_type_t} value. */

adrw_unit_style_t **inet_define_type_styles(void);
  /* Returns an array, indexed by {inet_block_type_t}, with the style to use
    when plotting each block type. */

char **inet_define_type_tags(void);

adrw_point_vec_t inet_define_points(void);
  /* Returns the reference points of the drawing. */

void inet_define_points_in_diagram(adrw_point_vec_t *P, int *np);
  /* Appends to {P} the main reference points of the diagram.
    Assumes {P->e[0..np-1]} are already defined and updates {*np}.
    Points along a vertical edge of a module are numbered {N+0..N+ncy}
    where {N} is the number of the module's lower corner. Reference point {N+ncy} 
    is the upper corner and coincides with point {N+MMY+0}, the lower corner
    of the module above it. */

void inet_plot_all(epswr_figure_t *epsf, adrw_building_t *B, int nx, int ny, bool_t show_dots);

void inet_append_backbone(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void inet_append_services(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void inet_append_connections(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);

int main (int argc, char **argv)
  {
    /* Define block types, names and styles: */
    adrw_unit_style_t **style = inet_define_type_styles();
    /* char **type_tag = inet_define_type_tags(); */
    
    /* Build diagram: */
    adrw_point_vec_t P = inet_define_points();
    adrw_building_t *B = adrw_make_building();
    
    inet_append_backbone(B, &P, style);
    inet_append_services(B, &P, style);
    inet_append_connections(B, &P, style);
    
    /* Print the point table: */
    { FILE *wr = open_write("out/inet-p.txt", TRUE);
      adrw_print_points(wr, &P);
      fclose(wr);
    }

    /* Plot the diagram: */
    { bool_t eps_fmt = TRUE;
      epswr_figure_t *epsf = pswr_new_stream("out/inet-", NULL, eps_fmt, "d", "letter", FALSE, 788.0, 508.0);
      inet_plot_all(epsf, B, 1, 1, FALSE);
      pswr_close_stream(epsf);
    }
    
    return 0;
  }

adrw_point_vec_t inet_define_points(void)
  {
    adrw_point_vec_t P = adrw_point_vec_new(100);
    int np = 0;
    
    fprintf(stderr, "=== DEFINING POINTS ===\n");
    
    inet_define_points_in_diagram(&P, &np);
    
    adrw_point_vec_trim(&P, np);
    
    return P;
  }
  
/* DIAGRAM POINTS */

#define PTBASE 1000
  /* Implicilty added to all point indices. */

#define SZX 20
#define SZY 20
  /* Size of a block module. */

#define NBX 2
#define NBY 8
  /* Size of backbone block in modules. */

#define NSX 2
#define NSY 1
  /* Size of service blocks in modules. */

#define NCY 8
  /* Number of reference points per module. */

#define DIAG_SZX ((2*NSX+NBX)*SZX)
#define DIAG_SZY (NBY*SZY)
  /* Size of entire diagram. */

#define MMX 100 /* Increment in point index per module column. */
#define MMY 10  /* Increment in point index per module row. */

#define RRND 1.0
  /* Corner rounding radius (diagram units). */

void inet_define_points_in_diagram(adrw_point_vec_t *P, int *np)
  {
    /* All dimensions in millimeters. */

    double Xor = 0; /* X of diagram's lower left corner. */
    double Yor = 0; /* Y of diagram's lower left corner. */
    double Zor = 0; /* Z of diagram (arbitrary). */

    auto void s(int dj, double dX, double dY, int di);
      /* Defines {P[dj]+(dX,dY,Zor)} as the coordinates of
        point {P[di]}.  If {P[di]} was defined
        previously, checks whether the definitions agree.
        The point numbers {di,dj} are added to {PTBASE}. */
       
    void s(int dj, double dX, double dY, int di)
      { char *lab = NULL;
        int i = di + (PTBASE);
        int j = (dj < 0 ? dj : dj + (PTBASE));
        double dZ = Zor;
        asprintf(&lab, "P%04d", i);
        adrw_append_point(lab, i, j, j, j, dX, dY, dZ, P, np); 
      } 
    
    fprintf(stderr, "--- module corners and connection points ---\n");
    int ndx = NBX + 2*NSX; /* Modules in diagram (X). */
    int ndy = NBY;  /* Modules in diagram (Y). */
    int ncy = NCY;  /* Reference points per module (Y). */
    assert(MMY >= (ncy+1));
    assert(MMX > MMY*(ndy+1));
    int ix, iy, cy;
    for (ix = 0; ix <= ndx; ix++)
      for (iy = 0; iy <= ndy; iy++)
        for (cy = 0; cy <= ncy; cy++)
          { int ip = MMX*(ix+1) + MMY*(iy+1) + cy;
            double xp = Xor+SZX*ix;
            double yp = Yor+SZY*(iy + ((double)cy)/ncy);
            s( -1, xp, yp, ip);
          }
  }

/* DIAGRAM OBJECTS */

#define ADDPOLY(DIAG,PTS,LABEL,DESC,TYPE,ROUND,...) \
  do \
    { fprintf(stderr, " [%s]", LABEL); \
      int v[] = { __VA_ARGS__ , -1 }; \
      int i = 0; while (v[i] >= 0) { v[i] += (PTBASE); i++; } \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),1,(PTS),v,(ROUND),0,(TYPE),style[(TYPE)]); \
      adrw_append_unit((DIAG), rm); \
    } \
  while(0)

#define ADDBOX(DIAG,PTS,LABEL,DESC,TYPE,ROUND,CTR,WDX,WDY) \
  do \
    { fprintf(stderr, " [%s]", LABEL); \
      adrw_unit_t *rm = adrw_make_box((LABEL),(DESC),(PTS),(CTR)+(PTBASE),0.0,0.0,0.0,(WDX),(WDY),(ROUND),0,(TYPE),style[(TYPE)]); \
      adrw_append_unit((DIAG), rm); \
    } \
  while(0)

#define ADDDOT(DIAG,PTS,PTIX,TYPE) \
  do \
    { fprintf(stderr, " [P%04d]", (PTIX)+(PTBASE)); \
      adrw_unit_t *rm = adrw_make_dot(NULL,NULL,(PTS),(PTIX)+(PTBASE), 0.0,0.0,0.0, 0,(TYPE),style[(TYPE)]); \
      adrw_append_unit((DIAG), rm); \
    } \
  while(0)

void inet_append_backbone(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- backbone modules ---\n");
    ADDPOLY
      ( B, P, "BKB1", "Backbone 1", inet_block_type_BKB, RRND,
        310, 410, 420, 320, 310
      );
    ADDPOLY
      ( B, P, "BKB2", "Backbone 2", inet_block_type_BKB, RRND,
        410, 510, 530, 430, 440, 340, 320, 420, 410
      );
    ADDPOLY
      ( B, P, "BKB3", "Backbone 3", inet_block_type_BKB, RRND,
        530, 550, 350, 340, 440, 430, 530
      );
    ADDPOLY
      ( B, P, "BKB4", "Backbone 4", inet_block_type_BKB, RRND,
        350, 550, 590, 390, 350
      );
    fprintf(stderr, "\n");
  }

void inet_append_services(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- DNS servers ---\n");
    ADDPOLY
      ( B, P, "DNS1", "DNS server 1", inet_block_type_DNS, RRND,
        211, 311, 317, 217, 211
      );
    ADDPOLY
      ( B, P, "DNS2", "DNS server 2", inet_block_type_DNS, RRND,
        551, 651, 657, 557, 551
      );
    fprintf(stderr, "--- Google ---\n");
    ADDPOLY
      ( B, P, "GOO", "Google", inet_block_type_GOO, RRND,
        121, 321, 327, 127, 121
      );
    fprintf(stderr, "--- Banco ---\n");
    ADDPOLY
      ( B, P, "BAN", "Banco", inet_block_type_BAN, RRND,
        131, 331, 337, 137, 131
      );
    fprintf(stderr, "\n");
  }

void inet_append_connections(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- backbone connection dots ---\n");
    int nsx = NSX; /* Modules in service blocks (X). */
    int nbx = NBX; /* Modules in backbone (X). */
    int nby = NBY; /* Modules in backbone (Y). */
    int ix, iy, kc;
    for (ix = nsx; ix <= nsx+nbx; ix += nbx)
      for (iy = 0; iy < nby; iy++)
        { int ib = MMX*(ix+1) + MMY*(iy+1);  /* Lower corner. */
          for (kc = 2; kc <= 6; kc += 2)
            { int ip = ib + kc;
              ADDDOT(B, P, ip, inet_block_type_BKC);
            }
        }
    fprintf(stderr, "\n");
  }


void inet_plot_all(epswr_figure_t *epsf, adrw_building_t *B, int nx, int ny, bool_t show_dots)
  {
    double xwid = DIAG_SZX;
    double ywid = DIAG_SZY;
    double xmrg =  5;
    double ymrg =  5;
    double xmin = 00 - xmrg, xmax = xwid + xmrg;
    double ymin = 00 - ymrg, ymax = ywid + ymrg;

    int ox, oy;
    for (ox = 0; ox < nx; ox++)
      for (oy = 0; oy < ny; oy++)
        {
          fprintf(stderr, "=== PLOTTING PAGE [%d,%d] OF [%d,%d] ===\n", ox,oy,nx,ny);
          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Situacao atual");
          pswr_set_label_font(epsf, "Courier", 8.0);
          adrw_plot_building(epsf, B, show_dots);
        }
  }

adrw_unit_style_t **inet_define_type_styles(void)
  {
    int ntypes = inet_block_type_MAX + 1;
    adrw_unit_style_t **style = (adrw_unit_style_t **)notnull(malloc(ntypes*sizeof(adrw_unit_style_t *)), "no mem");
    frgb_t fill_rgb[ntypes];
    frgb_t draw_rgb[ntypes];
    double pen_width[ntypes];
    frgb_t dots_rgb[ntypes];
    double dot_radius[ntypes];
    
    /* Defaults: */
    int i;
    for (i = 0; i < ntypes; i++)
      { fill_rgb[i] = (frgb_t){{ -1,-1,-1 }};
        draw_rgb[i] = (frgb_t){{ 0.000, 0.000, 0.000 }}; /* Y = 0.0000 */
        pen_width[i] = 0.10;
        dots_rgb[i] = (frgb_t){{ -1,-1,-1 }};
        dot_radius[i] = -1.0;
      }

    /* Polygon fill colors: */
    fill_rgb[inet_block_type_BKB] = (frgb_t){{ 0.000, 0.600, 0.550 }}; /* Y = 0.4150 */
    fill_rgb[inet_block_type_DNS] = (frgb_t){{ 0.650, 0.350, 1.000 }}; /* Y = 0.5050 */
    
    fill_rgb[inet_block_type_GOO] = (frgb_t){{ 0.790, 0.910, 0.890 }}; /* Y = 0.8720 */
    fill_rgb[inet_block_type_BAN] = (frgb_t){{ 0.900, 0.850, 0.750 }}; /* Y = 0.8550 */
    
    
    // color[inet_block_type_DOC] = (frgb_t){{ 0.000, 0.520, 1.000 }}; /* Y = 0.4120 */
    // color[inet_block_type_INF] = (frgb_t){{ 0.150, 0.680, 0.000 }}; /* Y = 0.4530 */
    // 
    // color[inet_block_type_REU] = (frgb_t){{ 0.820, 0.820, 0.820 }}; /* Y = 0.8200 */
    // color[inet_block_type_DEP] = (frgb_t){{ 0.930, 0.780, 0.530 }}; /* Y = 0.8000 */
    // color[inet_block_type_POS] = (frgb_t){{ 0.850, 0.850, 0.650 }}; /* Y = 0.8300 */
    // color[inet_block_type_SRV] = (frgb_t){{ 1.000, 0.700, 0.850 }}; /* Y = 0.8050 */
    // 
    // color[inet_block_type_ETC] = (frgb_t){{ 0.900, 0.900, 0.900 }}; /* Y = 0.9000 */
    // color[inet_block_type_ARR] = (frgb_t){{ 1.000, 0.300, 0.000 }}; /* Y = 0.4800 */
    // 
    // color[inet_block_type_NEX] = (frgb_t){{ 0.200, 0.200, 0.200 }}; /* Y = 0.2000 */
    // color[inet_block_type_PIL] = (frgb_t){{ 0.000, 0.300, 0.100 }}; /* Y = 0.1900 */
    // color[inet_block_type_SIT] = (frgb_t){{ 0.500, 0.100, 0.700 }}; /* Y = 0.2800 */

    /* Polygon (and dot) stroke colors: */
    draw_rgb[inet_block_type_BKC] = (frgb_t){{ -1, -1, -1 }}; /* Y = 0.0000 */
    
    /* Polygon (and dot) stroke linewidths: */
    pen_width[inet_block_type_BKC] = -1.0;
    
    /* Dot fill colors: */
    dots_rgb[inet_block_type_BKC] = (frgb_t){{ 0.000, 0.000, 0.000 }}; /* Y = 0.0000 */

    /* Dot radius: */
    dot_radius[inet_block_type_BKC] = 1.0;
    
    for (i = 0; i <= inet_block_type_MAX; i++)
      { style[i] = adrw_make_unit_style
          ( &(fill_rgb[i]), &(draw_rgb[i]), pen_width[i], &(dots_rgb[i]), dot_radius[i] );
      }
    return style;
  }

char **inet_define_type_tags(void)
  { int ntypes = inet_block_type_MAX + 1;
    char **tag = (char **)notnull(malloc(ntypes*sizeof(char *)), "no mem");
    tag[inet_block_type_BKB] = "BKB";
    tag[inet_block_type_DNS] = "DNS";
    
    return tag;
  } 

