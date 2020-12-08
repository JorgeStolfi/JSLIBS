/* See {gem_print_face.h} */
/* Last edited on 2014-07-23 19:57:02 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>

#include <jsfile.h>

#include <gem.h>

#include <gem_print_face.h>
  
void gem_print_face_2d_vertices_named(char *filename, gem_ref_t root, int d)
  { FILE *wr = open_write(filename, TRUE);
    gem_print_face_2d_vertices(wr, root, d);
    if (wr != stdout) { fclose(wr); }
  }
  
#define gem_print_face_INI_NODES 200
  /* Initial size of node stacks. */

void gem_print_face_2d_vertices(FILE *wr, gem_ref_t root, int d)
  { int i, j;

    /* Collecting a node {faces[i]} in each 2-face: */
    int R[d]; /* Wall-colors to use. */
    int nr = 0; /* Number of wall-colors to use. */
    for (i = 0; i < d; i++) { if (i != 2) { R[nr] = i; nr++; } }
    assert(nr == d-1);
    
    int S[1]; /* The wall-color 2 only. */
    S[0] = 2;
    int ns = 1;

    gem_ref_vec_t nodes = gem_ref_vec_new(gem_print_face_INI_NODES);
    gem_ref_vec_t faces = gem_ref_vec_new(gem_print_face_INI_NODES);
    int nn = 0; /* Number of nodes. */
    int nf = 0; /* Number of 2-faces. */
    gem_domains_enum(root, nr, R, 1, S, &nodes, &nn, &faces, &nf);

    /* Print all vertices of each 2-face: */
    gem_ref_vec_t vertices = gem_ref_vec_new(gem_print_face_INI_NODES);
    nr = 0; /* Number of colors except colors 0 and 2: */
    for (i = 0; i < d; i++) { if ((i != 0) && (i != 2)) { R[nr] = i; nr++; } }
    assert(nr == d-2);
    S[0] = 0; /* The wall-color 0 only. */
    for (i = 0; i < nf; i++)
      { int nv = 0; /* Number of vertices in face. */
        int nnf = 0; /* Nodes in face. */
        gem_domains_enum(faces.e[i], nr, R, ns, S, &nodes, &nnf, &vertices, &nv);
        for (j = 0; j < nv; j++) 
          { if (j > 0) { fprintf(wr, " "); }
            fprintf(wr, "%d", gem_get_data(vertices.e[j]));
          }
        fprintf(wr,"\n");
      }
    fflush(wr);
    free(nodes.e);
    free(faces.e);
    free(vertices.e);
  }
