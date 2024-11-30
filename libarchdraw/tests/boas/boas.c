/* Planta basica da casa na Boaretto da Silva 113 */
/* Last edited on 2023-02-20 17:43:33 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include <r2.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>

#include <archdraw.h>

#include <boas_define_points.h>
#include <boas_draw.h>

#define OUT_PREFIX "boas"

void plot_all(char *prefix, adrw_point_vec_t *P, int32_t nx, int32_t ny, bool_t show_dots);

int32_t main (int32_t argc, char **argv)
  {

    adrw_point_vec_t P = boas_define_points();
    
    /* Print the point table: */
    FILE *wr = open_write("out/" OUT_PREFIX "_P.txt", TRUE);
    adrw_print_points(wr, &P);
    fclose(wr);
    
    /* Plot the floorplan: */
    plot_all(OUT_PREFIX "_A", &P, 1, 3, TRUE);
    plot_all(OUT_PREFIX "_B", &P, 1, 3, FALSE);
    plot_all(OUT_PREFIX "_C", &P, 1, 1, FALSE);

    return 0;
  }

void plot_all(char *prefix, adrw_point_vec_t *P, int32_t nx, int32_t ny, bool_t show_dots)
  {
    double xmin = -50, xmax = +1550;
    double ymin = -350, ymax = +3250;
    
    epswr_figure_t *epsf = NULL; /* Defined later. */

    for (uint32_t ox = 0;  ox < nx; ox++)
      { for (int32_t oy = ny-1; oy >= 0; oy--)
          { epsf = adrw_new_figure
              ( "out", prefix, "1R", xmin, xmax, ymin, ymax,
                ox, nx, oy, ny, "Remocao de paredes e ladrilhos"
              );

            /* Removal: */
            boas_draw_land_plot(epsf, P, show_dots);
            boas_draw_trimmed_pavement(epsf, P, show_dots);
            boas_draw_cuts_pavement(epsf, P, show_dots);
            boas_draw_drain_pipes(epsf, P, show_dots);
            boas_draw_drain_boxes(epsf, P, show_dots);
            boas_draw_trimmed_walls(epsf, P, show_dots);
            boas_draw_cuts_walls(epsf, P, show_dots);
            boas_draw_building(epsf, P, show_dots);
            
            epswr_end_figure(epsf);

            /* Build: */
            epsf = adrw_new_figure
              ( "out", prefix, "2N", xmin, xmax, ymin, ymax,
                ox, nx, oy, ny, "Novas paredes e ladrilhos"
              );
            boas_draw_land_plot(epsf, P, show_dots);
            boas_draw_trimmed_pavement(epsf, P, show_dots);
            boas_draw_new_pavement(epsf, P, show_dots);
            boas_draw_drain_pipes(epsf, P, show_dots);
            boas_draw_drain_boxes(epsf, P, show_dots);
            boas_draw_trimmed_walls(epsf, P, show_dots);
            boas_draw_new_walls(epsf, P, show_dots);
            boas_draw_building(epsf, P, show_dots);

            epswr_end_figure(epsf);
          }
      }
  }

