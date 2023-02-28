/* See {archdraw.h}. */
/* Last edited on 2023-02-21 04:43:04 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <values.h>
#include <string.h>
#include <stdint.h>

#include <epswr.h>
#include <r2.h>
#include <r3.h>
#include <bool.h>
#include <affirm.h>
#include <jsfile.h>

#include <archdraw.h>

/* INTERNAL PROTOTYPES */

void check_defined(r3_t *p, int32_t ip);
  /* Bombs out if {p} is not defined. The {ip} is used in the message. */

void r3_shift(r3_t *p, double dx, double dy, double dz, r3_t *r);
  /* Same as {r3_add(p, &d, r)} where {d = (r3_p){{dx,dy,dz}}}. */

char *adrw_make_derived_label(char *main_label, int32_t k);
  /* Creates a label for a point derived from another point.
    Namely returns "{main_label}.{k}". */

void adrw_print_unit_data
  ( FILE *wr,
    bool_t TeX, 
    char *label, 
    char *descr, 
    char *type, 
    double modules, 
    double area, 
    double length
  );
  /* Prints a line of a table with data for a unit (or total of various units). 
   If any {modules,area,length} is NAN, that field (and its separator) is not printed. */

void adrw_print_unit_table_sep_line
  ( FILE *wr,
    bool_t TeX, 
    bool_t print_modules, 
    bool_t print_area, 
    bool_t print_length
  );
  /* Prints a separator line for a table with data for building units. */

adrw_point_t adrw_get_point(adrw_point_vec_t *P, int32_t ip);
  /* If {ip} is positive, returns a copy of point number {ip}
    from the point list {P}, checks whether it is defined,
    and marks it as used.
    
    If {ip} is negative, returns a point with position {(0,0,0)},
    {.lab = NULL}, and {.krot =0}. */

/* IMPLEMENTATIONS */

bool_t is_defined(r3_t *p)
  { return (! (isnan(p->c[0]) || isnan(p->c[1]) || isnan(p->c[2]))); }
  
void check_defined(r3_t *p, int32_t ip)
  { if (! is_defined(p)) 
      { fprintf(stderr, "** point %04d is undefined\n", ip);
        exit(1);
      }
  }

void adrw_append_point
  ( char *lab, 
    int32_t ip, 
    int32_t jpx, 
    int32_t jpy, 
    int32_t jpz, 
    double dX, 
    double dY, 
    double dZ, 
    adrw_point_vec_t *P, 
    int32_t *nP
  )
  {
    /* Print the data: */
    fprintf(stderr, "%-8s P[%04d]", lab, ip); 
    fprintf(stderr, " = P[ ");
    r3_t r;
    for (int32_t i = 0; i < 3; i++) 
      { int32_t jp = (int32_t[]){ jpx, jpy, jpz }[i];
        adrw_point_t apj = adrw_get_point(P, jp);
        if (i > 0) { fprintf(stderr, ", "); }
        if (jp < 0)
          { r.c[i] = 0.0;
            fprintf(stderr, "%4s", "ORIG");
          }
        else
          { r.c[i] = apj.p.c[i];
            fprintf(stderr, "%04d", jp);
          }
      }
    fprintf(stderr, " ]");
    fprintf(stderr, " + (%5.0f %5.0f %5.0f)", dX, dY, dZ); 
    
    if (! is_defined(&r))
      { fprintf(stderr, "\n");
        demand(FALSE, "undefined ref point"); 
      }
    
    double X = r.c[0] + dX;
    double Y = r.c[1] + dY;
    double Z = r.c[2] + dZ;
    
    /* Make sure that {P.e[ip]} exists: */
    int32_t n = (*nP);
    adrw_point_vec_expand(P, ip); 

    /* If {ip >= n}, mark all points {n..ip} as undefined, bump {n}: */
    while(n <= ip) 
      { r3_t u = (r3_t){{NAN, NAN, NAN}};
        P->e[n] = (adrw_point_t){ .p = u, .lab = NULL, .krot = (60*n) % 360, .used = FALSE };
        n++;
      }
    /* Get the point's record: */
    adrw_point_t *api = &(P->e[ip]);
      
    /* Define the coordinates and label of point {ip}: */
    fprintf(stderr, " = (%5.0f %5.0f %5.0f)", X, Y, Z); 
    fprintf(stderr, " rot = %3d", api->krot); 
    fprintf(stderr, "\n"); 
    r3_t *ppi = &(api->p);
    if(is_defined(ppi))
      { if ((ppi->c[0] != X) || (ppi->c[1] != Y) || (ppi->c[2] != Z))
          { fprintf(stderr, "** error: P[%04d] (%s) previously defined", ip, lab);
            fprintf(stderr, " as (%10.5f %10.5f %10.5f)", ppi->c[0], ppi->c[1], ppi->c[2]);
            fprintf(stderr, " error = (%10.5f %10.5f %10.5f)\n", ppi->c[0]-X, ppi->c[1]-Y, ppi->c[2]-Z);
          }
      }
    else
      { api->p = (r3_t){{X, Y, Z}}; 
        api->lab = lab;
      }
    (*nP) = n;
  }
  
void show_holes(adrw_point_vec_t *P, int32_t ini, int32_t fin)
  { int32_t ip;
    demand(ini >= 0, "invalid ini");
    int32_t hole_ini = -999;  /* Beginning of most recent gap. */
    int32_t hole_fin = -998;  /* End of most recent gap. */
    for (ip = ini; ip <= fin+1; ip++)
      { bool_t defined;
        if (ip >= fin) 
          { defined = TRUE; /* Hack to flush out last gap. */ }
        else if (ip >= P->ne) 
          { defined = FALSE; }
        else 
          { adrw_point_t *api = &(P->e[ip]);
            r3_t *p = &(api->p);
            defined = is_defined(p);
          }
        if ((hole_ini <= hole_fin) && (hole_fin == ip-1))
          { /* Item follows a gap: */
            if(defined)
              { /* Gap ended: */
                fprintf(stderr, "P[%04d] .. P[%04d] undefined\n", hole_ini, hole_fin);
              }
            else
              { /* Extend gap: */
                hole_fin = ip;
              }
          }
        else
          { if (! defined)
              { /* Starts a new gap: */
                hole_ini = hole_fin = ip;
              }
          }
      }
  }
  
adrw_unit_style_t *adrw_make_unit_style
  ( frgb_t *fill_rgb,
    frgb_t *draw_rgb,
    double pen_width,
    frgb_t *dots_rgb,
    double dot_radius
  )
  { adrw_unit_style_t *st = notnull(malloc(sizeof(adrw_unit_style_t)), "no mem");
    (*st) = (adrw_unit_style_t)
      { 
        .fill_rgb = (fill_rgb != NULL ? (*fill_rgb) : (frgb_t){{ -1,-1,-1 }}),
        .draw_rgb = (draw_rgb != NULL ? (*draw_rgb) : (frgb_t){{ -1,-1,-1 }}),
        .pen_width = pen_width,
        .dots_rgb = (dots_rgb != NULL ? (*dots_rgb) : (frgb_t){{ -1,-1,-1 }}),
        .dot_radius = dot_radius
      };
    return st;
  }

bool_t adrw_polygon_is_closed(adrw_point_vec_t *P)
  { 
    int32_t n = P->ne;
    if (n < 3) { return FALSE; }
    double x0 = P->e[0].p.c[0], xn = P->e[n-1].p.c[0];
    double y0 = P->e[0].p.c[1], yn = P->e[n-1].p.c[1];
    return ((x0 == xn) && (y0 == yn));
  }

void adrw_compute_area_and_length(adrw_point_vec_t *P, double *areaP, double *lengthP)
  {
    int32_t n = P->ne;
    double length = 0.0;          /* Length (or perimeter) in 3D. */
    double area = 0.0;            /* Area of XY projection. */
    int32_t ip;
    bool_t closed = adrw_polygon_is_closed(P);
    for (ip = 1; ip < n; ip++)
      { adrw_point_t *ap0 = &(P->e[ip-1]);
        adrw_point_t *ap1 = &(P->e[ip]);
        r3_t *p0 = &(ap0->p);
        r3_t *p1 = &(ap1->p);
        length += r3_dist(p0, p1);
        if (closed)
          { double x0 = p0->c[0], y0 = p0->c[1];
            double x1 = p1->c[0], y1 = p1->c[1];
            area += 0.5*(x0*y1 - x1*y0);
          }
      }
    (*areaP) = fabs(area);
    (*lengthP) = length;
  }

adrw_unit_t *adrw_make_poly
  ( char *label, 
    char *descr, 
    double modules,
    adrw_point_vec_t *P, 
    int32_t v[], 
    double round,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  )
  { 
    /* Count corners: */
    int32_t m = 0;  while(v[m] >= 0) { m++; }
    demand(m > 0, "no points in polygon");
    /* Convert corner indices to points and pack as {adrw_unit_t}: */
    adrw_unit_t *rm = notnull(malloc(sizeof(adrw_unit_t)), "no mem");
    (*rm) = (adrw_unit_t)
      { .label = label,
        .descr = descr,
        .pt = adrw_point_vec_new(m),
        .round = round,
        .kfloor = kfloor,
        .type = type,
        .style = style,
        .modules = modules,
        .area = NAN,
        .length = NAN
      };
    int32_t k;
    for (k = 0; k < m; k++) { rm->pt.e[k] = adrw_get_point(P, v[k]); }
    adrw_compute_area_and_length(&(rm->pt), &(rm->area), &(rm->length));
    return rm;
  }

void r3_shift(r3_t *p, double dx, double dy, double dz, r3_t *r)
  { 
    r->c[0] = p->c[0] + dx;
    r->c[1] = p->c[1] + dy;
    r->c[2] = p->c[2] + dz;
  }
  
adrw_point_t adrw_make_derived_point(adrw_point_t *ap, int32_t k, double dx, double dy, double dz, int32_t krot)
  { 
    adrw_point_t dp;
    r3_shift(&(ap->p), dx, dy, dz, &(dp.p)); 
    if ((dx == 0) && (dy == 0) && (dz == 0))
      { dp.lab = ap->lab; 
        dp.krot = ap->krot;
      }
    else
      { dp.lab = adrw_make_derived_label(ap->lab, k); 
        dp.krot = (krot % 360);
      }
    ap->used = TRUE;
    dp.used = FALSE;
    return dp;
  }

char *adrw_make_derived_label(char *main_label, int32_t k)
  { 
    if (main_label == NULL)
      { return NULL; }
    else
      { char *label = NULL;
        asprintf(&label, "%s.%d", main_label, k);
        return label;
      }
  }

adrw_point_t adrw_get_point(adrw_point_vec_t *P, int32_t ip)
  {
    adrw_point_t dp;
    if (ip >= 0)
      { demand((ip > 0) && (ip < P->ne), "invalid reference index");
        dp = P->e[ip];
        check_defined(&(dp.p), ip);
        P->e[ip].used = TRUE;
      }
    else
      { dp = (adrw_point_t){ .p = (r3_t){{ 0,0,0 }}, .lab = NULL, .krot = 0, .used = FALSE }; }
    return dp;
  }

adrw_unit_t *adrw_make_box
  ( char *label, 
    char *descr, 
    adrw_point_vec_t *P, 
    int32_t ip, 
    double ctrx,
    double ctry,
    double ctrz,
    double wdx,
    double wdy,
    double round,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  )
  { 
    /* Get the center reference point: */
    adrw_point_t dp = adrw_get_point(P, ip);
    
    /* Create an an {adrw_unit_t} and compute its corners: */
    adrw_unit_t *rm = notnull(malloc(sizeof(adrw_unit_t)), "no mem");
    (*rm) = (adrw_unit_t)
      { .label = label,
        .descr = descr,
        .pt = adrw_point_vec_new(5),
        .round = round,
        .kfloor = kfloor,
        .type = type,
        .style = style,
        .modules = 0.0,
        .area = 0.0,
        .length = 0.0
      };
    adrw_point_t *rp = &(rm->pt.e[0]);
    rp[0] = adrw_make_derived_point(&dp, 0, ctrx-0.5*wdx, ctry-0.5*wdy, ctrz,  225);
    rp[1] = adrw_make_derived_point(&dp, 1, ctrx+0.5*wdx, ctry-0.5*wdy, ctrz,  315);
    rp[2] = adrw_make_derived_point(&dp, 2, ctrx+0.5*wdx, ctry+0.5*wdy, ctrz,   45);
    rp[3] = adrw_make_derived_point(&dp, 3, ctrx-0.5*wdx, ctry+0.5*wdy, ctrz,  135);
    rp[4] = rp[0];
    return rm;
  }

adrw_unit_t *adrw_make_dot
  ( char *label, 
    char *descr, 
    adrw_point_vec_t *P, 
    int32_t ip, 
    double ctrx,
    double ctry,
    double ctrz,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  )
  {
    /* Get the center reference point: */
    adrw_point_t ap = adrw_get_point(P, ip);
    
    /* Create an an {adrw_unit_t} for the dot: */
    adrw_unit_t *rm = notnull(malloc(sizeof(adrw_unit_t)), "no mem");
    (*rm) = (adrw_unit_t)
      { .label = label,
        .descr = descr,
        .pt = adrw_point_vec_new(1),
        .round = 0.0,
        .kfloor = kfloor,
        .type = type,
        .style = style,
        .modules = 0.0,
        .area = 0.0,
        .length = 0.0
      };
    adrw_point_t *rp = &(rm->pt.e[0]);
    rp[0] = adrw_make_derived_point(&ap, 0, ctrx, ctry, ctrz, (ap.krot + 90) % 360);
    return rm;
  }

adrw_building_t *adrw_make_building(void)
  {
    adrw_building_t *B = notnull(malloc(sizeof(adrw_building_t)), "no mem");
    (*B) = (adrw_building_t)
      { .NU = 0,
        .unit = adrw_unit_vec_new(100)
      };
    return B;
  }

void adrw_append_unit(adrw_building_t *B, adrw_unit_t *rm)
  {
    int32_t k = B->NU;
    adrw_unit_vec_expand(&(B->unit),k);
    B->unit.e[k] = rm;
    B->NU++;
  }

void adrw_append_seats
  ( adrw_building_t *B, 
    adrw_point_vec_t *P, 
    int32_t v00,
    int32_t v11,
    double szx,
    double szy,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  )
  {
    /* Get corner coordinates: */
    demand((v00 >= 0) && (v00 < P->ne), "invalid point index");
    adrw_point_t *ap00 = &(P->e[v00]);
    r3_t *p00 = &(ap00->p); double x00 = p00->c[0], y00 = p00->c[1], z00 = p00->c[2];
    check_defined(p00, v00);
    demand((v11 >= 0) && (v11 < P->ne), "invalid point index");
    adrw_point_t *ap11 = &(P->e[v11]);
    r3_t *p11 = &(ap11->p); double x11 = p11->c[0], y11 = p11->c[1], z11 = p11->c[2];
    check_defined(p11, v11);
    
    /* Ignore signs of sizes: */
    szx = fabs(szx);
    szy = fabs(szy);
    
    /* Compute number of seats: */
    int32_t nx = (int32_t)floor(fabs(x11 - x00)/szx + 0.001);
    int32_t ny = (int32_t)floor(fabs(y11 - y00)/szy + 0.001);
    
    /* Compute actual spacing of seats (signed): */
    double dx = (x11 - x00)/nx;
    double dy = (y11 - y00)/ny;
    double dz = (z11 - z00)/ny; /* Sic: {ny} not {nz}. */
    
    /* Compute signd of displacements: */
    int32_t signx = (dx < 0 ? -1 : +1);
    int32_t signy = (dy < 0 ? -1 : +1);
    
    /* Determine signed size {wdx,wdy} of seat proper, and signed margins to leave around it: */
    double wdx, wdy;
    double xmlo, xmhi, ymlo, ymhi;
    if (szx < szy)
      { /* Seats are turned towards {±Y} */
        wdx = szx*signx;
        wdy = 0.80*szx*signy;
        xmlo = xmhi = 0.5*(dx - wdx);
        ymlo = 0;
        ymhi = (dy - wdy);
      }
    else
      { /* Seats are turned towards {±X} */
        wdy = szy*signy;
        wdx = 0.80*szy*signx;
        ymlo = ymhi = 0.5*(dy - wdy);
        xmlo = 0;
        xmhi = (dx - wdx);
      }
   
    /* Append the seats as little boxes: */
    int32_t ix, iy;
    for (ix = 0; ix < nx; ix++)
      { for (iy = 0; iy < ny; iy++)
          { double xc = ix*dx + xmlo + 0.5*wdx;
            double yc = iy*dy + ymlo + 0.5*wdy;
            double zc = iy*dz;
            adrw_unit_t *rm = adrw_make_box(NULL, NULL, P, v00, xc, yc, zc, wdx, wdy, 0.5, kfloor, type, style);
            adrw_append_unit(B, rm);
          }
      }
  }

epswr_figure_t *adrw_new_figure
  ( char *dir,
    char *prefix,
    char *suffix,
    double xmin, 
    double xmax, 
    double ymin, 
    double ymax, 
    int32_t ox, 
    int32_t nx, 
    int32_t oy, 
    int32_t ny, 
    char *title
  )
  {
    bool_t segmented = ((nx > 1) || (ny > 1)); /* Only part of a multipage plot. */
    
    /* Compute the desired sub-rectangle {[xlo_xhi]×[ylo_yhi]}: */
    double dx = (xmax - xmin)/nx, mx = (nx == 0 ? 0.00 : 0.02*dx);
    double dy = (ymax - ymin)/ny, my = (ny == 0 ? 0.00 : 0.02*dy);
    
    double xlo = xmin + ox*dx - mx, xhi = xmin + (ox+1)*dx + mx;
    double ylo = ymin + oy*dy - my, yhi = ymin + (oy+1)*dy + my;

    double maxHSize = 190*epswr_pt_per_mm;
    double maxVSize = 210*epswr_pt_per_mm;
    int32_t capLines = (title != NULL ? 1 : 0) + (segmented ? 2 : 0);
    double textFontHeight = 10.0;
    bool_t eps_verbose = FALSE;
    int32_t npage = (segmented ? 10000 + oy*nx + ox : 0);
    epswr_figure_t *epsf = epswr_new_captioned_figure
      ( dir, prefix, NULL, npage, suffix,
        xlo,xhi, ylo,yhi, maxHSize, maxVSize,
        capLines, textFontHeight, eps_verbose
      );
      
    /* Draw a light gray frame: */
    epswr_set_pen(epsf,  0.800f,0.800f,0.800f,  0.10, 0,0);
    epswr_frame(epsf);
    
    /* Write the caption text: */
    epswr_set_text_font(epsf, "Courier", textFontHeight);
    epswr_text(epsf, title, FALSE, 0.5, TRUE,FALSE);
    if (segmented)
      { char *title1 = NULL;
        asprintf(&title1, "page (%d,%d) of (%d,%d)", ox,oy,nx,ny);
        epswr_text(epsf, title1, FALSE, 0.5, TRUE,FALSE);
        free(title1);
        char *title2 = NULL;
        asprintf(&title2, "[%.1f _ %.1f] x [%.1f _ %.1f]", xlo,xhi,ylo,yhi);
        epswr_text(epsf, title2, FALSE, 0.5, TRUE,FALSE);
        free(title2);
      }

    double labelFontHeight = 10.0;
    epswr_set_label_font(epsf, "Courier", labelFontHeight);

    /* epswr_set_pen(epsf,  0.000,0.000,0.000,  0.10, 0,0); */
    /* epswr_coord_line(epsf, epswr_axis_HOR, 0.0); */
    /* epswr_coord_line(epsf, epswr_axis_VER, 0.0); */
    epswr_set_pen(epsf,  0.000,0.000,0.000,  0.10, 0,0);
    epswr_set_fill_color(epsf,  0.000,0.000,0.000);
    
    return epsf;
  }

void adrw_print_points(FILE *wr, adrw_point_vec_t *P)
  { int32_t ip;
    for (ip = 0; ip < P->ne; ip++)
      { fprintf(wr, "P%04d = ", ip);
        adrw_point_t *api = &(P->e[ip]);
        r3_t *p = &(api->p);
        if(! is_defined(p))
          { fprintf(wr, "undefined"); }
        else
          { fprintf(wr, "(%5.0f %5.0f %5.0f)", p->c[0], p->c[1], p->c[2]);
            fprintf(wr, " %-10s", api->lab);
            fprintf(wr, " rot = %3d", api->krot);
            if (! api->used) { fprintf(wr, " UNUSED"); }
          }
        fprintf(wr, "\n");
      } 
  }

void adrw_plot_point
  ( epswr_figure_t *epsf, 
    char *lab, 
    int32_t ip, 
    double x, 
    double y, 
    double rot, 
    double hAlign, 
    double vAlign
  )
  { char *ptlab = NULL;
    if (lab == NULL)
      { asprintf(&ptlab, "P%04d", ip); }
    else
      { ptlab = lab; }
    epswr_dot(epsf, x, y, 0.25, TRUE, FALSE);
    epswr_set_label_font(epsf, "Courier", 4.0);
    epswr_label(epsf, ptlab, "P", x, y, rot, FALSE, hAlign, vAlign, TRUE, FALSE);
    if (lab != ptlab) { free(ptlab); }
  }

void adrw_plot_unit
  ( epswr_figure_t *epsf,
    adrw_unit_t *rm,
    bool_t show_dots
  )
  {
    bool_t closed = adrw_polygon_is_closed(&(rm->pt));
    
    /* Plot the unit: */
    adrw_unit_style_t *st = rm->style;
    frgb_t *frgb = &(st->fill_rgb);
    bool_t frgb_ok = (!isnan(frgb->c[0])) && (frgb->c[0] >= 0) && (frgb->c[0] < +INF);
    
    frgb_t *drgb = &(st->draw_rgb);
    bool_t drgb_ok = (!isnan(drgb->c[0])) && (drgb->c[0] >= 0) && (drgb->c[0] < +INF);
    
    double penw = st->pen_width;
    bool_t penw_ok = (!isnan(penw)) && (penw > 0) && (penw < +INF);
    
    frgb_t *vrgb = &(st->dots_rgb);
    bool_t vrgb_ok = (!isnan(vrgb->c[0])) && (vrgb->c[0] >= 0) && (vrgb->c[0] < +INF);
    
    double drad = st->dot_radius;
    bool_t drad_ok = (!isnan(drad)) && (drad > 0) && (drad < +INF);

    bool_t fill = FALSE, draw = FALSE;

    if (closed && frgb_ok)
      { epswr_set_fill_color(epsf, frgb->c[0],frgb->c[1],frgb->c[2]);
        fill = TRUE;
      }
    if (drgb_ok && penw_ok)
      { epswr_set_pen(epsf,  drgb->c[0],drgb->c[1],drgb->c[2],  penw, 0,0);
        draw = TRUE; 
      }
    
    /* Extract corner coordinates {x[0..n-1], y[0..n-1]} minus closer: */
    int32_t n = (closed ? rm->pt.ne - 1 : rm->pt.ne);
    double x[n], y[n];
    char *lab[n];
    int32_t krot[n];
    int32_t k;
    for (k = 0; k < n; k++)
      { adrw_point_t *api = &(rm->pt.e[k]);
        r3_t *p = &(api->p);
        x[k] = p->c[0]; 
        y[k] = p->c[1];
        lab[k] = api->lab;
        krot[k] = api->krot;
      }

    /* Draw polygon/polyline: */
    if (rm->round == 0)
      { epswr_polygon(epsf, closed, x, y, n, fill, draw, TRUE); }
    else
      { epswr_rounded_polygon(epsf, closed, x, y, n, rm->round, fill, draw, TRUE); }
    
    if (show_dots)
      { /* Plot vertices as black dots with standard size, and labels: */
        epswr_set_fill_color(epsf, 0,0,0);
        epswr_set_pen(epsf,  0.000,0.000,0.000,  0.10, 0,0);
        for (k = 0; k < n; k++) {           /* Reduce label direction angle to {-90..269}: */
          int32_t krotk = krot[k] % 360; 
          if (krotk < 0) { krotk += 360; }
          if (krotk >= 270) { krotk -= 360; }
          double rot, hAlign;
          if (krotk < 90)
            { rot = (double)krotk; hAlign = -0.10; }
          else
            { rot = (double)(krotk - 180); hAlign = 1.10; }
          double vAlign = 0.5;
          adrw_plot_point(epsf, lab[k], 0, x[k], y[k], rot, hAlign, vAlign);
        }
      }
    else if (drad_ok && (vrgb_ok || draw))
      { /* Plot vertices with specified dot size, no labels: */
        if (vrgb_ok) { epswr_set_fill_color(epsf, vrgb->c[0],vrgb->c[1],vrgb->c[2]); }
        /* Assumes the pen has already been set to {drgb,penw}. */
        for (k = 0; k < n; k++) { epswr_dot(epsf, x[k], y[k], drad, vrgb_ok, draw); }
      }
  }

void adrw_plot_building
  ( epswr_figure_t *epsf,
    adrw_building_t *B,
    bool_t show_dots
  )
  { 
    int32_t iu;
    for (iu = 0; iu < B->NU; iu++)
      { adrw_unit_t *rm = B->unit.e[iu];
        adrw_plot_unit(epsf, rm, show_dots);
      }
  }

void adrw_print_unit
  ( FILE *wr,
    adrw_unit_t *rm,
    char *type_tag[],
    bool_t print_modules,
    bool_t print_area,
    bool_t print_length,
    bool_t TeX
  )
  {
    /* Get number of modules, length, area: */
    double modules = rm->modules;   /* Number of standard office modules. */
    double length = rm->length;   /* Length (or perimeter) in 3D. */
    double area = rm->area;       /* Area of XY projection. */
    
    /* Supress unwanted fields: */
    if (! print_modules) { modules = NAN; }
    if (! print_area)    { area = NAN; }
    if (! print_length)  { length = NAN; }
    
    /* Print it: */
    adrw_print_unit_data(wr, TeX, rm->label, rm->descr, type_tag[rm->type], modules, area, length);
  }
  
void adrw_print_unit_data
  ( FILE *wr,
    bool_t TeX, 
    char *label, 
    char *descr, 
    char *type, 
    double modules, 
    double area, 
    double length
  )
  {
    auto void prnt(char *fmt, double val);
    /* Prints {val} to {wr} in format {fmt}, but changes '.' to ','. */

    void prnt(char *fmt, double val) 
      { char *s = NULL;
        asprintf(&s, fmt, val);
        int32_t i = 0; while (s[i] != 0) { if (s[i] == '.') { s[i] = ','; } i++; }
        fputs(s, wr);
        free(s);
      }
    
    if (TeX)
      { fprintf(wr, "  %7s", (label == NULL ? "~" : label));
        if (! isnan(modules))  { fprintf(wr, " & "); prnt("%4.1f", modules);  }
        if (! isnan(area))   { fprintf(wr, " & "); prnt("%6.1f", area/10000);   }
        if (! isnan(length)) { fprintf(wr, " & "); prnt("%6.1f", length/100); }
        fprintf(wr, " & %s", (type == NULL ? "~" : type));
        fprintf(wr, " & %-24s", (descr == NULL ? "~" : descr));
        fprintf(wr, " \\\\\n");
      }
    else
      {
        fprintf(wr, " | %7s", (label == NULL ? "" : label));
        if (! isnan(modules))  { fprintf(wr, " | %4.1f", modules); }
        if (! isnan(area))   { fprintf(wr, " | %6.1f", area/10000); }
        if (! isnan(length)) { fprintf(wr, " | %6.1f", length/100); }
        fprintf(wr, " | %s", (type == NULL ? "" : type));
        fprintf(wr, " | %-24s", (descr == NULL ? "" : descr));
        fprintf(wr, " |\n");
      }
  }

void adrw_print_unit_table_sep_line
  ( FILE *wr,
    bool_t TeX, 
    bool_t print_modules, 
    bool_t print_area, 
    bool_t print_length
  )
  {
    if (TeX)
      { fprintf(wr, "  \\hline\n"); }
    else
      {
        fprintf(wr, " +--------");
        if (print_modules)  { fprintf(wr, "-+-------"); }
        if (print_area)   { fprintf(wr, "-+-------"); }
        if (print_length) { fprintf(wr, "-+-------"); }
        fprintf(wr, "-+----");
        fprintf(wr, "-+---------------");
        fprintf(wr, "-+\n");
      }
  }

void adrw_print_building
  ( FILE *wr,
    adrw_building_t *B,
    bool_t select[],
    char *type_tag[],
    bool_t print_modules,
    bool_t print_area,
    bool_t print_length,
    bool_t TeX
  )
  { 
    /* Sort building units by type, label: */
    
    auto int32_t ucomp(const void *a, const void *b);
    
    int32_t ucomp(const void *a, const void *b)
      { adrw_unit_t *rma = *((adrw_unit_t **)a);
        adrw_unit_t *rmb = *((adrw_unit_t **)b);
        if (rma->type < rmb->type)
          { return -1; }
        else if (rma->type > rmb->type)
          { return +1; }
        else
          { char *laba = (rma->label == NULL ? "" : rma->label);
            char *labb = (rmb->label == NULL ? "" : rmb->label);
            return strcmp(laba, labb);
          }
      }
      
    qsort(B->unit.e, B->NU, sizeof(adrw_unit_t *), &ucomp);
    
    /* Now print the units with subtotal and total lines per type: */
    
    auto void print_sub_tot(char *descr, char *tag, double modules, double area, double length);
      /* Prints a total or subtotal line. */
      
    void print_sub_tot(char *descr, char *tag, double modules, double area, double length)
      { /* Supress unwanted subtotal fields: */
        if (! print_modules) { modules = NAN; }
        if (! print_area)    { area = NAN; }
        if (! print_length)  { length = NAN; }
        /* Print subtotal line: */
        adrw_print_unit_table_sep_line(wr, TeX, print_modules, print_area, print_length);
        adrw_print_unit_data
          ( wr, TeX, NULL, descr,  tag, modules, area, length );
        adrw_print_unit_table_sep_line(wr, TeX, print_modules, print_area, print_length);
      }
    
    adrw_print_unit_table_sep_line(wr, TeX, print_modules, print_area, print_length);
    int32_t iu;
    int32_t old_type = -1;
    double tot_modules = 0, tot_area = 0, tot_length = 0;
    double sub_modules = 0, sub_area = 0, sub_length = 0;
    for (iu = 0; iu < B->NU; iu++)
      { adrw_unit_t *rm = B->unit.e[iu];
        if ((select != NULL) && (! select[rm->type])) { continue; }
        if (rm->type != old_type)
          { if (old_type >= 0)
              { print_sub_tot("Subtotal",  type_tag[old_type], sub_modules, sub_area, sub_length); }
            sub_modules = sub_area = sub_length = 0;
            old_type = rm->type;
          }
        adrw_print_unit(wr, rm, type_tag, print_modules, print_area, print_length, TeX);
        sub_modules += rm->modules;  tot_modules += rm->modules;
        sub_area += rm->area;        tot_area += rm->area;      
        sub_length += rm->length;    tot_length += rm->length;  
      }
    if (old_type >= 0)
      { print_sub_tot("Subtotal",  type_tag[old_type], sub_modules, sub_area, sub_length); }
    print_sub_tot("Total",  "TOT", tot_modules, tot_area, tot_length);
    fflush(wr);
  }

void adrw_compute_building_stats
  ( adrw_building_t *B,
    int32_t ntypes,
    double units[],
    double modules[],
    double area[],
    double length[]
  )
  {
    int32_t type;
    for (type = 0; type < ntypes; type++)
      { if (units != NULL) units[type] = 0.0;
        if (modules != NULL) modules[type] = 0.0;
        if (area != NULL) area[type] = 0.0; 
        if (length != NULL) length[type] = 0.0;
      }
    int32_t iu;
    for (iu = 0; iu < B->NU; iu++)
      { adrw_unit_t *rm = B->unit.e[iu];
        type = rm->type;
        demand((type >= 0) && (type < ntypes), "invalid unit type");
        if (units != NULL) units[type] += 1;
        if (modules != NULL) modules[type] += rm->modules;
        if (area != NULL) area[type] += rm->area; 
        if (length != NULL) length[type] += rm->length;
      }
  }
  
void adrw_plot_type_legend
  ( char *dir,
    char *prefix,
    char *suffix,
    int32_t key_type[], 
    int32_t ntypes,
    int32_t ncols, 
    char *type_tag[], 
    adrw_unit_style_t *style[]
  )
  {
    /* Count the keys to be plotted: */
    int32_t nkeys = 0;  /* Number of keys actually plotted. */
    while (key_type[nkeys] >= 0) 
      { demand((key_type[nkeys] >= 0) && (key_type[nkeys] < ntypes), "invalid type");
        nkeys++;
      }
    
    /* Compute the size of the rectangle in mm: */
    double wdy = 4.0;
    double wdx = wdy*(1 + sqrt(5))/2;
    
    /* Row and column displacements in mm: */
    double dy = 2.0*wdy;
    double dx = 6.0*wdy;

    /* Compute the number of rows: */
    int32_t nrows = (nkeys + ncols - 1)/ncols;
    
    /* Compute the Client sub-rectangle {[xlo_xhi]×[ylo_yhi]} (mm): */
    double xlo = 0.0, xhi = ncols*dx;
    double ylo = 0.0, yhi = nrows*dy;
    
    /* Compute the size of the plot area in pt: */
    double hSize = (xhi - xlo)*epswr_pt_per_mm;
    double vSize = (yhi - ylo)*epswr_pt_per_mm;
        
    /* Choose margins and set up the plot window */
    double hvMarg = 4.0;
    bool_t eps_verbose = TRUE;
    epswr_figure_t *epsf = epswr_new_named_figure
      ( dir, prefix, NULL, -1, suffix, 
        hSize, vSize, hvMarg, hvMarg, hvMarg, hvMarg, 
        eps_verbose
      );
    epswr_set_client_window(epsf, xlo, xhi, ylo, yhi);
    epswr_set_label_font(epsf, "Courier", 14.0);
    epswr_set_pen(epsf,  0.000,0.000,0.000,  0.10, 0,0);

    int32_t slot; 
    for (slot = 0; slot < nkeys; slot++)
      { int32_t type = key_type[slot];
        int32_t row = nrows - 1 - (slot / ncols);
        int32_t col = slot % ncols;
        double rxlo = col*dx + 2.0, rxhi = rxlo + wdx;
        double rylo = row*dy + (dy - wdy)/2, ryhi = rylo + wdy;
        adrw_unit_style_t *st = style[type];
        frgb_t *rgb = &(st->fill_rgb);
        epswr_set_fill_color(epsf, rgb->c[0],rgb->c[1],rgb->c[2]);
        epswr_rectangle(epsf, rxlo, rxhi, rylo, ryhi, TRUE, TRUE);
        epswr_set_fill_color(epsf, 0,0,0);
        epswr_label(epsf, type_tag[type], "P", rxhi + 1.5, rylo + 0.5, 0.0, FALSE, 0.0,0.0, TRUE, FALSE);
      }
    epswr_end_figure(epsf);
  }
  
void adrw_plot_histogram_bar
  ( char *dir,
    char *prefix,
    char *name,
    char *suffix,
    adrw_space_type_t type,
    adrw_unit_style_t *st,
    double val,
    double vmax
  )
  {
    /* Compute the height and width of the bar in mm, including stroked border: */
    double wdx = 160.0;
    double wdy = 4.0;

    /* Compute the size of the plot area in pt: */
    double hSize = wdx * epswr_pt_per_mm;
    double vSize = wdy * epswr_pt_per_mm;
        
    /* Choose margins and set up the plot window */
    bool_t eps_verbose = TRUE;
    epswr_figure_t *epsf = epswr_new_named_figure
      ( dir, prefix, name, -1, suffix, 
        hSize, vSize, 0.0, 0.0, 0.0, 0.0, 
        eps_verbose
      );
    
    /* Pen width in mm: */
    double pwd = 0.20;
    epswr_set_pen(epsf,  0.000,0.000,0.000,  pwd, 0,0);

    /* Set Client window to plot window but in mm: */
    epswr_set_client_window(epsf, 0, wdx, 0, wdy);
    
    /* The bar rectangle (excluding half of the stroked border) in mm: */
    double rwdx = fmax(0.0, wdx*val/vmax - 1.1*pwd);               
    double rwdy = fmax(0.0, wdy - 1.1*pwd);               
    double rxlo = 1.1*pwd/2, rxhi = rxlo + rwdx;
    double rylo = 1.1*pwd/2, ryhi = rylo + rwdy;
    frgb_t *rgb = &(st->fill_rgb);
    epswr_set_fill_color(epsf, rgb->c[0],rgb->c[1],rgb->c[2]);
    epswr_rectangle(epsf, rxlo, rxhi, rylo, ryhi, TRUE, TRUE);
    epswr_end_figure(epsf);
  }

vec_typeimpl(adrw_point_vec_t, adrw_point_vec, adrw_point_t); 

vec_typeimpl(adrw_unit_vec_t, adrw_unit_vec, adrw_unit_t *);
