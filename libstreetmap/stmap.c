/* See stmap.h */
/* Last edited on 2023-02-21 21:34:10 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include <stmap.h>
#include <stheap.h>

#include <epswr.h>
#include <affirm.h>
#include <r2.h>
#include <fget.h>
#include <nget.h>
#include <jsfile.h>
#include <jsmath.h>
#include <sign.h>
#include <sign_get.h>

/* INTERNAL PROTOTYPES */

void st_map_file_error(int32_t nlin, char *msg);
void st_map_add_edge(Map *m, int32_t org, int32_t dst, int32_t ei, float c0, float c1);
bool_t st_map_insert_edge_in_ring(quad_arc_t e, quad_arc_t r);
int32_t  st_map_ccw(Point p, Point q, Point r);
bool_t st_map_vertex_is_visible(Point p, Interval xr, Interval yr);
bool_t st_map_edge_is_visible(Point p, Point q, Interval xr, Interval yr);

/* IMPLEMENTATIONS */

Map *st_map_read(FILE *f)
  { Map *m = (Map *)malloc(sizeof(Map));
    int32_t vi;
    int32_t ei, ai;
    int32_t nlin = 0;
    int32_t res;

    affirm(m != NULL, "out of memory"); 
    res = fscanf(f, "p street-graph %d %d\n", &(m->nv), &(m->ne));
    nlin++;
    fprintf(stderr, "res = %d\n", res);
    if (res != 2) { st_map_file_error(nlin, "bad header line"); }
    m->out = (quad_arc_t *)malloc(m->nv * sizeof(quad_arc_t));
    m->along = (quad_arc_t *)malloc(2 * m->ne * sizeof(quad_arc_t));
    m->vd = (VertexData **)malloc(m->nv * sizeof(VertexData *));
    m->ed = (EdgeData **)malloc(m->ne * sizeof(EdgeData *));
    
    for (vi = 0; vi < m->nv; vi++) { m->vd[vi] = NULL; m->out[vi] = quad_arc_NULL; }
    for (ei = 0; ei < m->ne; ei++) { m->ed[ei] = NULL; }
    for (ai = 0; ai < 2*m->ne; ai++) { m->along[ai] = quad_arc_NULL; }

    for (vi = 0; vi < m->nv; vi++) 
      { VertexData *vd = (VertexData *)malloc(sizeof(VertexData));
        int32_t vot, vin, id;
        res = fscanf(f, "v %d %lf %lf %d %d\n", &id, &(vd->p.c[0]), &(vd->p.c[1]), &vin, &vot);
        nlin++;
        if (res != 5) { st_map_file_error(nlin, "bad vertex line"); }
        affirm(id == vi, "vertices out of order");
        vd->id = id; vd->deg = vin+vot;
        m->vd[id] = vd;
      }
    for (ei = 0; ei < m->ne; ei++)
      { int32_t org, dst, blocked;
        float c0, c1;
        res = fscanf(f, "a %d %d %d\n", &org, &dst, &blocked);
        nlin++;
        if (res != 3) { st_map_file_error(nlin, "bad edge line"); }
        if (blocked == 1) 
          { c0 = c1 = +INF; }
        else
          { double d = r2_dist(&(m->vd[org]->p), &(m->vd[dst]->p));
            c0 = c1 = (float)d;
          }
        st_map_add_edge(m, org, dst, ei, c0, c1);
      }
    return m;
  }
  
#define DIRID(a) ((((EdgeData *)quad_ldata(a))->id << 1) | quad_sym_bit(a))

void st_map_init_costs(Map *m, float *d, quad_arc_t *e, float *c)
  { int32_t vi, ai;
    for(vi = 0; vi < m->nv; vi++) { d[vi] = +INF; e[vi] = quad_arc_NULL; }
    for(ai = 0; ai < 2*m->ne; ai++) { c[ai] = +INF; }
  }

void st_map_reset_costs(Map *m, int32_t *r, int32_t nr, float *d, quad_arc_t *e, float *c)
  { int32_t k;
    for (k = 0; k < nr; k++)
      { int32_t vi = r[k]; 
        d[vi] = +INF; e[vi] = quad_arc_NULL;
        quad_arc_t b = m->out[vi], a = b;
        do
          { EdgeData *ad = (EdgeData *)quad_ldata(a);
            int32_t s = quad_lon_bit(a);
            int32_t ai = 2 * ad->id + s;
            c[ai] = +INF;
            a = quad_onext(a);
          }
        while (a != b);
      }
  }

void st_map_compute_costs
  ( Map *m, 
    int32_t u, 
    float dMax,
    int32_t *r,
    int32_t *nr, 
    float *d, 
    quad_arc_t *e,
    float *c
  )
  {
    int32_t vi;
    st_Heap *h = st_heap_new(2*m->ne);
    float cprev = 0.0, dprev = 0.0;

    /* Let DIRID(a) be the ID number of the *directed* edge {a},
      namely {DIRID(a) = 2*ei+s} where {ei} is the ID of the
      undirected edge, and {s = quad_lon_bit(a)}. The value of {c[DIRID(a)]}
      is the cost of the optimal path from {u} that ends with the
      quad_arc_t {a}; or is {+INF} if that cost is still unknown. */
    
    affirm(dMax >= 0.0, "invalid {dMax}");

    (*nr) = 0;
    d[u] = 0.0;
    vi = u; 
    do
      { /* The next reached vertex in order of increasing cost is {vi}. */
        /* Consistency check: */
        affirm(d[vi] >= dprev, "vertices popped out of order");
        dprev = d[vi];
        /* Store {vi} in list of reached vertices: */
        r[*nr] = vi; (*nr)++;
        /* Insert in heap all fresh edges out of {vi}, truncating at {dMax}: */
        { quad_arc_t b = m->out[vi], a = b;
          do
            { EdgeData *ad = (EdgeData *)quad_ldata(a);
              int32_t s = quad_lon_bit(a);
              int32_t ai = 2 * ad->id + s;
              float cnew = d[vi] + ad->cost[s];
              affirm(cnew >= d[vi], "negative arc cost");
              if (cnew <= dMax)
                { VertexData *vd = (VertexData *)quad_ddata(a);
                  int32_t wi = vd->id;
                  affirm(c[ai] == +INF, "total arc cost c[ai] set twice");
                  c[ai] = cnew; 
                  if (cnew < d[wi]) { st_heap_insert(h, ai, c); }
                }
              a = quad_onext(a);
            }
          while (a != b);
        }
        /* Get next vertex {vi} in cost order, or -1 if none: */
        vi = -1;
        while ((h->n > 0) && (vi < 0))
          { int32_t ai = st_heap_pop(h, c);
            quad_arc_t a = m->along[ai];
            VertexData *vd = (VertexData *)quad_ddata(a);
            int32_t wi = vd->id;
            affirm(c[ai] >= cprev, "total arc costs popped out of order");
            cprev = c[ai];
            if (c[ai] < d[wi]) { d[wi] = c[ai]; e[wi] = a; vi = wi; }
          }
      }
    while (vi >= 0);
    st_heap_discard(h);
  }
     
void st_compute_coverage
  ( Map *m, 
    int32_t *u, 
    float *dMax, 
    int32_t n, 
    int32_t *vcover,
    int32_t *ecover
  )
  {
    int32_t vi, ei;
    for (vi = 0; vi < m->nv; vi++) { vcover[vi] = 0; }
    for (ei = 0; ei < m->ne; ei++) { ecover[ei] = 0; }
    if (n > 0) 
      { /* Work files for {st_map_compute_costs}: */
        float d[m->nv];
        float c[2 * m->ne];
        quad_arc_t e[m->nv];
        int32_t r[m->nv];
        int32_t nr;
        /* Count sites that cover each vertex and each edge: */
        st_map_init_costs(m, d, e, c);
        int32_t i;
        for(i = 0; i < n; i++)
          { int32_t ui = u[i];
            float dmi = dMax[i];
            /* Find vertices within cost {dMax} from {ui}: */
            st_map_compute_costs(m, ui, dmi, r, &nr, d, e, c);
            /* Tally vertex and edge coverage: */
            st_increment_coverage(m, dmi, r, nr, d, c, vcover, ecover);
            /* Reset costs for next {st_map_compute_costs}: */
            st_map_reset_costs(m, r, nr, d, e, c);
          }
      }
  }
  
void st_increment_coverage
  ( Map* m,
    float dMax, 
    int32_t *r, 
    int32_t nr, 
    float *d, 
    float *c,
    int32_t *vcover, 
    int32_t *ecover
  )
  { int32_t k;
    for (k = 0; k < nr; k++)
      { int32_t vi = r[k];
        /* Tally vertex coverage: */
        affirm(d[vi] <= dMax, "inconsistent cost");
        vcover[vi]++;
        /* Tally edge coverage: */
        quad_arc_t b = m->out[vi], a = b;
        do
          { EdgeData *ad = (EdgeData *)quad_ldata(a);
            int32_t ei = ad->id;
            int32_t ai = 2 * ei, aj = ai + 1;
            /* Check if edge {ei} is {dMax}-reachable: */
            if ((c[ai] <= dMax) || (c[aj] <= dMax))
              { /* Take care not to increment {ecover[ei]} twice: */
                VertexData *wd = (VertexData *)quad_odata(quad_sym(a));
                int32_t wi = wd->id;
                if (vi < wi) { ecover[ei]++; }
              }
            a = quad_onext(a);
          }
        while (a != b);
      }
  }

void st_map_plot(
    epswr_figure_t *eps, 
    Map *m, 
    Interval xr, Interval yr, 
    float *vwidth, RGBColor *vcolor,
    float *ewidth, RGBColor *ecolor
  )
  {
    int32_t vi, ei;
    bool_t clipping = FALSE;
    
    if ( 
        ( (xr.lo <= xr.hi) && (yr.lo <= yr.hi) ) &&
        ( (xr.lo != -INF) || (xr.hi != +INF) ||
          (yr.lo != -INF) || (yr.hi != +INF)
        )
      )
      { clipping = TRUE; }
    else
      { clipping = FALSE; }

    /* Plot edges: */
    for (ei = 0; ei < m->ne; ei++)
      { quad_arc_t a0 = m->along[2*ei];
        if (a0 != quad_arc_NULL)
          { VertexData *od = (VertexData *)quad_odata(a0);
            VertexData *dd = (VertexData *)quad_ddata(a0);
            if ((! clipping) || st_map_edge_is_visible(od->p, dd->p, xr, yr))
              { float pr,pg,pb,pwd;
                if (ecolor != NULL) 
                  { RGBColor *cp = &(ecolor[ei]); pr = cp->R; pg = cp->G; pb = cp->B; }
                else
                  { pr = 1.0; pg = 0.0; pb = 0.0; }
                pwd = (ewidth != NULL ? ewidth[ei] : 0.250f);
                if (pwd > 0.0)
                  { epswr_set_pen(eps, pr,pg,pb, pwd, 0.0, 0.0);
                    epswr_segment(eps, od->p.c[0], od->p.c[1], dd->p.c[0], dd->p.c[1]);
                  }
              }
          }
      }
      
    /* Plot vertices: */
    for (vi = 0; vi < m->nv; vi++)
      { VertexData *vd = m->vd[vi];
        if ((! clipping) || st_map_vertex_is_visible(vd->p, xr, yr))
          { float pr,pg,pb,pwd;
            if (vcolor != NULL) 
              { RGBColor *cp = &(vcolor[vi]); pr = cp->R; pg = cp->G; pb = cp->B; }
            else
              { pr = 1.000f; pg = 0.000f; pb = 0.000f; }
            pwd = (vwidth != NULL ? vwidth[vi] : (vd->deg == 2 ? 0.000f : 0.500f));
            /* Plot dot: */
            if (pwd > 0.0)
              { epswr_set_pen(eps, pr,pg,pb, pwd, 0.0, 0.0);
                epswr_segment(eps, vd->p.c[0], vd->p.c[1], vd->p.c[0], vd->p.c[1]);
              }
          }
      }
  }
  
void st_map_get_bbox(Map *m, Interval *xr, Interval *yr)
  {
    double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0;
    int32_t i;
    for (i = 0; i < m->nv; i++)
      { VertexData *vd = m->vd[i];
        double x = vd->p.c[0], y = vd->p.c[1];
        if (i == 0)
          { xlo = xhi = x; ylo = yhi = y; }
        else
          { if (x < xlo) { xlo = x; } else if (x > xhi) { xhi = x; }
            if (y < ylo) { ylo = y; } else if (y > yhi) { yhi = y; }
          }
      }
    xr->lo = xlo; xr->hi = xhi;
    yr->lo = ylo; yr->hi = yhi;
  }

int32_t st_map_nearest_vertex(Map *m, Point p)
  {
    int32_t vi;
    int32_t imin = -1;
    double d2min = +INF;
    for (vi = 0; vi < m->nv; vi++)
      { VertexData *vd = m->vd[vi]; 
        double dx = (p.c[0] - vd->p.c[0]);
        double dy = (p.c[1] - vd->p.c[1]);
        double d2 = dx*dx + dy*dy;
        if (d2 < d2min) { imin = vi; d2min = d2; }
      }
    return imin;
  }

/* INTERNAL PROCEDURES */
  
void st_map_file_error(int32_t nlin, char *msg)
   /*
     Prints {msg} and {nlin}, and exits with error status.
     Handy while parsing a map file. */
  {
    fprintf(stderr, "map format error on line %d: %s\n", nlin, msg); 
    exit(1);
  }

void st_map_add_edge(Map *m, int32_t org, int32_t dst, int32_t ei, float c0, float c1)
   /*
     Adds a new edge with number {ei} to the map's topology
     structure, between vertices {org} and {dst}, with costs {c0} and
     {c1} for forward and backwards traversal. Also stores the
     corresponding arcs in {m->out[org]}, {m->out[dst]} and
     {m->along[ei]}, and stores the edge's data record in {m->ed[ei]}.
     Assumes that the vertex data records {m->vd[...]} are already set. */
  {
    EdgeData *ed;
        
    quad_arc_t a0 = quad_make_edge();
    quad_arc_t a1 = quad_sym(a0);
    int32_t ia0 = 2*ei + quad_lon_bit(a0);
    int32_t ia1 = 2*ei + quad_lon_bit(a1);
    
    affirm(org != dst, "loop edge");

    /* Set origin and destination data pointers: */
    quad_set_odata(a0, (void*)m->vd[org]);
    quad_set_odata(a1, (void*)m->vd[dst]);
    
    /* Create edge data record and set data pointers to it: */
    ed = (EdgeData *)malloc(sizeof(EdgeData));
    ed->id = ei;
    ed->cost[0] = c0;
    ed->cost[1] = c1;
    m->ed[ei] = ed;
    quad_set_ldata(a0, (void *)ed);
    quad_set_ldata(a1, (void *)ed);

    /* Connect to adjacent edges, and store arc in topology tables: */
    { bool_t ok0 = TRUE, ok1 = TRUE;
      if (m->out[org] != quad_arc_NULL) 
        { ok0 = st_map_insert_edge_in_ring(a0, m->out[org]); }
      if (m->out[dst] != quad_arc_NULL)
        { ok1 = st_map_insert_edge_in_ring(a1, m->out[dst]); }
      if (!(ok0 & ok1))
        { /* Detach this edge (a desperate attempt to keep going): */
          fprintf(stderr, "ring insertion failed for edge %d = (%d,%d)\n", ei,org,dst);
          quad_splice(a0, quad_oprev(a0)); quad_splice(a1, quad_oprev(a1));
        }
      else
        { m->out[org] = a0; m->along[ia0] = a0;
          m->out[dst] = a1; m->along[ia1] = a1;
        }
    }
  }

bool_t st_map_insert_edge_in_ring(quad_arc_t e, quad_arc_t r)
  /* 
    Inserts {e} in the proper place among the ring of edges out 
    of some vertex {v}, of which {r} is a representative.
    Returns {FALSE} if the insertion failed for some reason. */
  {
    VertexData *od = (VertexData *)quad_odata(e);
    VertexData *ed = (VertexData *)quad_ddata(e);
    quad_arc_t a = r;
    bool_t here;
    
    affirm(quad_onext(e) == e, "edge ONEXT");
    do { 
      quad_arc_t b = quad_onext(a);
      if (b == a)
        { here = TRUE; }
      else
        { VertexData *ad = (VertexData *)quad_ddata(a);
          VertexData *bd = (VertexData *)quad_ddata(b);
          int32_t sae, seb, sba;

          affirm(quad_odata(a) == quad_odata(e), "edge origin");
          sae = st_map_ccw(od->p, ad->p, ed->p);
          seb = st_map_ccw(od->p, ed->p, bd->p);
          sba = st_map_ccw(od->p, bd->p, ad->p);
          here = (sae + seb + sba > 0);
        }
      if (here) { quad_splice(a, e); return TRUE; }
      a = b;
    } while (a != r);
    return FALSE;
  }
  
int32_t st_map_ccw(Point p, Point q, Point r)
  /*
   Sign of triangle {p,q,r} (-1, 0, or +1). */
  {
    double ux = q.c[0] - p.c[0];
    double uy = q.c[1] - p.c[1];
    double vx = r.c[0] - p.c[0];
    double vy = r.c[1] - p.c[1];
    double det = ux*vy - uy*vx;
    return sign_double(det);
  }

bool_t st_map_vertex_is_visible(Point p, Interval xr, Interval yr) 
   /* 
     Returns {TRUE} if the point {p} is within the box {xr × yr}. */
  { 
    return (p.c[0] >= xr.lo) & (p.c[0] <= xr.hi) & (p.c[1] >= yr.lo) & (p.c[1] <= yr.hi);
  }
  
double_vec_t st_read_double_vec_t(char *name)
  { FILE *rd = open_read(name, TRUE);
    double_vec_t d = st_fread_double_vec_t(rd);
    if (rd != stdin) { fclose(rd); }
    return d;
  }
  
double_vec_t st_fread_double_vec_t(FILE *rd)
  { int32_t ne = nget_int32(rd, "vertices"); fget_eol(rd);
    double_vec_t d = double_vec_new(ne);
    int32_t i;
    for (i = 0; i < ne; i++) { d.e[i] = fget_double(rd); fget_eol(rd); }
    return d;
  }
  
void st_write_double_vec_t(char *name, char *fmt, double_vec_t *d)
  { FILE *wr = open_write(name, TRUE);
    st_fwrite_double_vec_t(wr, fmt, d);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
  }

void st_fwrite_double_vec_t(FILE *wr, char *fmt, double_vec_t *d)
  { int32_t i;
    fprintf(wr, "vertices = %d\n", d->ne);
    for (i = 0; i < d->ne; i++) { fprintf(wr, fmt, d->e[i]); fputc('\n', wr); }
    fflush(wr);
  }

bool_t st_map_edge_is_visible(Point p, Point q, Interval xr, Interval yr) 
   /* 
     Returns {TRUE} if the bounding box of the segment {p--q}
     intercepts the box {xr × yr}. */
  { 
    if (p.c[0] < q.c[0])
      { if ((q.c[0] < xr.lo) || (p.c[0] > xr.hi)) { return FALSE; } }
    else
      { if ((p.c[0] < xr.lo) || (q.c[0] > xr.hi)) { return FALSE; } }
    if (p.c[1] < q.c[1])
      { if ((q.c[1] < yr.lo) || (p.c[1] > yr.hi)) { return FALSE; } }
    else
      { if ((p.c[1] < yr.lo) || (q.c[1] > yr.hi)) { return FALSE; } }
    return TRUE;
  }

void st_interval_widen(Interval *r, double margin)
  { r->lo = r->lo - margin;
    r->hi = r->hi + margin;
  }


void st_adjust_rect_shape(Interval *xr, Interval *yr, double tx, double ty)
  { double dx = xr->hi - xr->lo;
    double dy = yr->hi - yr->lo;
    if (ty*dx > tx*dy)
      { double ey = (ty*dx/tx - dy)/2.0;
        yr->lo = yr->lo - ey;
        yr->hi = yr->hi + ey;
      }
    else if (ty*dx < tx*dy)
      { double ex = (tx*dy/ty - dx)/2.0;
        xr->lo = xr->lo - ex;
        xr->hi = xr->hi + ex;
      }
  }
