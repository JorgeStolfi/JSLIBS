/* Planta basica da casa na Boaretto da Silva 113 */
/* Last edited on 2012-12-07 20:54:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <r2.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>
#include <archdraw.h>

#include <boas_define_points.h>
#include <boas_draw.h>

void plot_all(char *prefix, adrw_point_vec_t *P, int nx, int ny, bool_t show_dots);

int main (int argc, char **argv)
  {

    adrw_point_vec_t P = boas_define_points();
    
    /* Plot the floorplan: */
    plot_all("out/boas_A", &P, 1, 3, TRUE);
    plot_all("out/boas_B", &P, 1, 3, FALSE);
    plot_all("out/boas_C", &P, 1, 1, FALSE);
    
    /* Print the point table: */
    FILE *wr = open_write("out/boas_p.txt", TRUE);
    adrw_print_points(wr, &P);
    fclose(wr);

    return 0;
  }

void plot_all(char *prefix, adrw_point_vec_t *P, int nx, int ny, bool_t show_dots)
  {
    double hSize = 460;
    double vSize = 640;

    double xmin = -50, xmax = +1550;
    double ymin = -350, ymax = +3250;

    for (int ox = 0; ox < nx; ox++)
      for (int oy = ny-1; oy >= 0; oy--)
        {
          char *fname = NULL;
          asprintf(&fname, "%s_%03d_%03d.eps", prefix, ox, oy);
          FILE *wr = open_write(fname, TRUE);
          bool_t verbose = TRUE;
          double hvMarg = 4.0
          epswr_figure_t *epsf = epswr_new_figure(wr, hSize, vSize, hvMarg, hvMarg, hvMarg, hvMarg, verbose);

          /* Removal: */
          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Remocao de paredes e ladrilhos");
          boas_draw_land_plot(epsf, P, show_dots);
          boas_draw_trimmed_pavement(epsf, P, show_dots);
          boas_draw_cuts_pavement(epsf, P, show_dots);
          boas_draw_drain_pipes(epsf, P, show_dots);
          boas_draw_drain_boxes(epsf, P, show_dots);
          boas_draw_trimmed_walls(epsf, P, show_dots);
          boas_draw_cuts_walls(epsf, P, show_dots);
          boas_draw_building(epsf, P, show_dots);
 
          /* Build: */
          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Novas paredes e ladrilhos");
          boas_draw_land_plot(epsf, P, show_dots);
          boas_draw_trimmed_pavement(epsf, P, show_dots);
          boas_draw_new_pavement(epsf, P, show_dots);
          boas_draw_drain_pipes(epsf, P, show_dots);
          boas_draw_drain_boxes(epsf, P, show_dots);
          boas_draw_trimmed_walls(epsf, P, show_dots);
          boas_draw_new_walls(epsf, P, show_dots);
          boas_draw_building(epsf, P, show_dots);

          epswr_end_figure(epsf);
          free(fname);
        }
  }

