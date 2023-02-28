/* Draws various units of the house on Boaretto da Silva 113. */
/* Last edited on 2023-02-20 17:25:16 by stolfi */

#ifndef boas_draw_H
#define boas_draw_H

#define _GNU_SOURCE

#include <bool.h>
#include <archdraw.h>

void boas_draw_land_plot(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

void boas_draw_old_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void boas_draw_old_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

void boas_draw_cuts_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void boas_draw_cuts_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

void boas_draw_trimmed_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void boas_draw_trimmed_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

void boas_draw_new_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void boas_draw_new_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

void boas_draw_drain_boxes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void boas_draw_drain_pipes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

void boas_draw_building(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

#endif
