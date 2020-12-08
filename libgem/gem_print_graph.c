/* See {gem_print_graph} */
/* Last edited on 2014-07-23 19:49:37 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <jsfile.h>

#include <gem.h>

#include <gem_print_graph.h>

#define gem_print_graph_BLACK 0
#define gem_print_graph_GRAY 7

#define gem_print_graph_color_code(i) (1 << i)
//#define gem_print_graph_color_code(i) (i)

char *gem_print_graph_color_name[17] = 
  { "BLACK", "BLUE",  "GREEN", "BLACK",  "RED", 
    "BLACK", "BLACK", "BLACK", "ORANGE", "BLACK",
    "BLACK", "BLACK", "BLACK", "BLACK",  "BLACK",
    "BLACK", "YELLOW"
  };

void gem_print_graph_named(char *filename, gem_ref_t root, int d)
  { FILE *wr = open_write(filename, TRUE);
    gem_print_graph(wr, root, d);
    if (wr != stdout) { fclose(wr); }
  }
  
#define gem_print_graph_INI_NODES 200
  /* Initial size of node stacks. */

void gem_print_graph(FILE *wr, gem_ref_t root, int d)
  {
    gem_ref_vec_t nodes = gem_ref_vec_new(gem_print_graph_INI_NODES);
    int i, j, nn = 0;
    gem_traverse(root, d, &nodes, &nn);

    fprintf(wr, "graph gr { graph [fontsize=18];\n");
    fprintf(wr, "node [shape=circle,style=bold,height=\"0.25\",width=\"0.25\"];\n");
    fprintf(wr, "edge [style=bold,labelfontsize=3];\n");
    for (i = 0; i < nn; i++)
      { char *color = gem_print_graph_color_name[gem_print_graph_BLACK];
        fprintf(wr, "  %d [color=%s];\n", gem_get_label(nodes.e[i]), color);
      }
    for (i = 0; i < nn; i++)
      { for (j = 0; j < d; j++)
          { gem_ref_t q = gem_step(nodes.e[i], j);
            int ind = gem_get_label(q);
            if (ind <= i)
              { char *color = gem_print_graph_color_name[gem_print_graph_color_code(j)];
                fprintf(wr, "  %d -- %d [color=%s];\n", i, ind, color);
              }
          }
      }
    fprintf(wr, "}\n");
    fflush(wr);
    free(nodes.e);
  }
