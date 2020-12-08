/* See {stmesh_io.h} */
/* Last edited on 2016-05-02 16:27:35 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <jsfile.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <bool.h>
#include <affirm.h>

#include <stmesh.h>
#include <stmesh_rep.h>
#include <stmesh_STM.h>

/* INTERNAL PROTOTPES */

/* IMPLEMENTATIONS */
    
void stmesh_STM_write(FILE *wr, stmesh_t mesh)
  {
    /* Write header: */
    filefmt_write_header(wr, stmesh_STM_FILE_TYPE, stmesh_STM_FILE_VERSION);
    
    /* Write unit of measure: */
    fprintf(wr, "eps = %25.15e\n", stmesh_get_eps(mesh));
    
    /* Write element counts: */
    uint32_t nv = stmesh_vert_count(mesh);
    uint32_t ne = stmesh_edge_count(mesh);
    uint32_t nf = stmesh_face_count(mesh);
    fprintf(wr, "nv = %ud\n", nv);
    fprintf(wr, "ne = %ud\n", ne);
    fprintf(wr, "nf = %ud\n", nf);
    
    /* Write the vertex data: */
    { fprintf(wr, "vertices\n");
      stmesh_vert_unx_t uxv;
      for (uxv = 0; uxv < nv; uxv++) 
        { fprintf(wr, "%ud", uxv);
          stmesh_vert_t v = stmesh_get_vert(mesh, uxv);
          i3_t pos = stmesh_vert_get_pos(v);
          for (int j = 0; j < 3; j++) { fprintf(wr, " %d", pos.c[j]); }
          fprintf(wr, "\n");
        }
    }

    /* Write the edge data: */
    { fprintf(wr, "edges\n");
      stmesh_edge_unx_t uxe;
      for (uxe = 0; uxe < ne; uxe++) 
        { fprintf(wr, "%ud", uxe);
          stmesh_edge_t e = stmesh_get_edge(mesh, uxe);
          /* Get the vertex indices {uxv[0..1]} of the endpoints: */
          stmesh_vert_t v[2];
          stmesh_edge_get_endpoints(e, v);
          stmesh_vert_unx_t uxv[2];
          for (int j = 0; j < 2; j++) { uxv[j] = stmesh_vert_get_unx(mesh, v[j]); } 
          /* Ensure their order: */
          if (uxv[0] > uxv[1]) { stmesh_vert_unx_t t = uxv[0]; uxv[0] = uxv[1]; uxv[1] = t; }
          assert(uxv[0] < uxv[1]);
          /* Write the vertx numbers: */
          for (int j = 0; j < 2; j++) { fprintf(wr, " %ud", uxv[j]); }
          fprintf(wr, "\n");
        }
    }
    
    /* Write the face data: */
    { fprintf(wr, "faces\n");
      stmesh_face_unx_t uxf;
      for (uxf = 0; uxf < nf; uxf++) 
        { fprintf(wr, "%ud", uxf);
          stmesh_face_t f = stmesh_get_face(mesh, uxf);
          /* Get the edge indices {uxe[0..2]} of the sides: */
          stmesh_edge_t e[3];
          stmesh_face_get_sides(f, e);
          stmesh_edge_unx_t uxe[3];
          for (int j = 0; j < 3; j++) { uxe[j] = stmesh_edge_get_unx(mesh, e[j]); } 
          /* Ensure their order: */
          if (uxe[0] > uxe[1]) { stmesh_edge_unx_t t = uxe[0]; uxe[0] = uxe[1]; uxe[1] = t; }
          if (uxe[0] > uxe[2]) { stmesh_edge_unx_t t = uxe[0]; uxe[0] = uxe[2]; uxe[2] = t; }
          if (uxe[1] > uxe[2]) { stmesh_edge_unx_t t = uxe[1]; uxe[1] = uxe[2]; uxe[2] = t; }
          assert((uxe[0] < uxe[1]) && (uxe[1] < uxe[2]));
          for (int j = 0; j < 3; j++)  { fprintf(wr, " %ud", uxe[j]); }
          fprintf(wr, "\n");
        }
    }
    
    /* Write footer: */
    filefmt_write_footer(wr, stmesh_STM_FILE_TYPE);

    fflush(wr);
  }

void stmesh_STM_write_named (char *fileName, stmesh_t mesh, bool_t verbose)
  { FILE *wr = open_write(fileName, verbose);
    stmesh_STM_write(wr, mesh);
    fclose(wr);
  }
  
#define stmesh_STM_read_COORD_MAX ((int64_t)INT32_MAX)
  /* Max absolute quantized coordinate. */
  
stmesh_t stmesh_STM_read(FILE *rd)
  {
    bool_t validate = TRUE; /* True to perform validation checks. */
    
    /* Read header: */
    filefmt_read_header(rd, stmesh_STM_FILE_TYPE, stmesh_STM_FILE_VERSION);
    
    /* Read unit of masure: */
    float eps = (float)nget_double(rd, "eps");
    demand((! isnan(eps)) && (eps > 0.0), "invalid unit of measure");
    fget_eol(rd);
    
    /* Read element counts: */
    uint64_t nv = nget_uint64(rd, "nv", 10);
    demand(nv <= stmesh_nv_MAX, "too many vertices");
    fget_eol(rd); 

    uint64_t ne = nget_uint64(rd, "ne", 10);
    demand(ne <= stmesh_ne_MAX, "too many edges");
    fget_eol(rd); 

    uint64_t nf = nget_uint64(rd, "nf", 10);
    demand(nf <= stmesh_nf_MAX, "too many faces");
    fget_eol(rd); 

    /* Allocate the structure: */
    stmesh_t mesh = stmesh_new_desc(eps, (uint32_t)nv, (uint32_t)ne, (uint32_t)nf);
    
    /* Read the vertices: */
    { fget_skip_and_match(rd, "vertices"); fget_eol(rd);
      stmesh_vert_unx_t uxv;
      for (uxv = 0; uxv < nv; uxv++) 
        { /* Parse the vertex index: */
          uint64_t uxv_rd = fget_uint64(rd, 10);
          demand(uxv == uxv_rd, "vertex index mismatch");
          /* Parse the vertex coordinates: */
          i3_t pos;
          for (int j = 0; j < 3; j++) 
            { int64_t cj = fget_int64(rd);
              demand((cj < 0 ? -cj : +cj) <= stmesh_STM_read_COORD_MAX, "vertex coord too large");
              pos.c[j] = (int32_t)cj;
            }
          fget_eol(rd);
          /* Add the vertex: */
          stmesh_vert_unx_t uxv_ad = stmesh_add_vert(mesh, &pos);
          assert(uxv_ad == uxv);
          if (validate)
            { /* Paranoia checks: */
              stmesh_vert_t v = stmesh_get_vert(mesh, uxv);
              stmesh_vert_unx_t uxv_ck = stmesh_vert_get_unx(mesh, v);
              assert(uxv_ck == uxv);
              i3_t pos_ck = stmesh_vert_get_pos(v);
              for (int j = 0; j < 3; j++) { assert(pos.c[j] == pos_ck.c[j]); }
            }
        }
    }

    /* Read the edges: */
    { fget_skip_and_match(rd, "edges"); fget_eol(rd);
      stmesh_edge_unx_t uxe;
      for (uxe = 0; uxe < ne; uxe++) 
        { /* Parse the edge index: */
          uint64_t uxe_rd = fget_uint64(rd, 10);
          demand(uxe == uxe_rd, "edge index mismatch");
          /* Parse the edge endpoints: */
          stmesh_vert_unx_t uxv[2];
          for (int j = 0; j < 2; j++) 
            { int64_t vj = fget_int64(rd);
              demand((vj >= 0) && (vj < nv), "invalid edge endpoint");
              uxv[j] = (stmesh_edge_unx_t)vj;
            }
          fget_eol(rd);
          demand(uxv[0] < uxv[1], "edge endpoints out of order");
          /* Add the edge: */
          stmesh_edge_unx_t uxe_ad = stmesh_add_edge(mesh, uxv);
          assert(uxe_ad == uxe);
          if (validate)
            { /* Paranoia checks: */
              stmesh_edge_t e = stmesh_get_edge(mesh, uxe);
              stmesh_edge_unx_t uxe_ck = stmesh_edge_get_unx(mesh, e);
              assert(uxe_ck == uxe);
              stmesh_vert_t v_ck[2];
              stmesh_edge_get_endpoints(e, v_ck);
              stmesh_vert_unx_t uxv_ck[2];
              for (int j = 0; j < 2; j++) { uxv_ck[j] = stmesh_vert_get_unx(mesh, v_ck[j]); }
              assert
                ( ((uxv[0] == uxv_ck[0]) && (uxv[1] == uxv_ck[1])) ||
                  ((uxv[0] == uxv_ck[1]) && (uxv[1] == uxv_ck[0])) );
            }
        }
    }

    /* Read the faces: */
    { fget_skip_and_match(rd, "faces"); fget_eol(rd);
      stmesh_face_unx_t uxf;
      for (uxf = 0; uxf < nf; uxf++) 
        { /* Parse the face index: */
          uint64_t uxf_rd = fget_uint64(rd, 10);
          demand(uxf == uxf_rd, "face index mismatch");
          /* Parse the face sides: */
          stmesh_edge_unx_t uxe[2];
          for (int j = 0; j < 3; j++) 
            { int64_t ej = fget_int64(rd); /* Unoriented edge number. */
              demand((ej >= 0) && (ej < ne), "invalid face side");
              uxe[j] = (stmesh_edge_unx_t)ej;
            }
          fget_eol(rd);
          demand((uxe[0] < uxe[1]) && (uxe[1] < uxe[2]), "face sides out of order");
          /* Add the face: */
          stmesh_face_unx_t uxf_ad = stmesh_add_face(mesh, uxe);
          assert(uxf_ad == uxf);
          if (validate)
            { /* Paranoia checks: */
              stmesh_face_t f = stmesh_get_face(mesh, uxf);
              stmesh_face_unx_t uxf_ck = stmesh_face_get_unx(mesh, f);
              assert(uxf_ck == uxf);
              stmesh_edge_t e_ck[3];
              stmesh_face_get_sides(f, e_ck);
              stmesh_edge_unx_t uxe_ck[3];
              for (int j = 0; j < 3; j++) { uxe_ck[j] = stmesh_edge_get_unx(mesh, e_ck[j]); }
              if (uxe_ck[0] > uxe_ck[1]) { stmesh_edge_unx_t t = uxe_ck[0]; uxe_ck[0] = uxe_ck[1]; uxe_ck[1] = t; }
              if (uxe_ck[0] > uxe_ck[2]) { stmesh_edge_unx_t t = uxe_ck[0]; uxe_ck[0] = uxe_ck[2]; uxe_ck[2] = t; }
              if (uxe_ck[1] > uxe_ck[2]) { stmesh_edge_unx_t t = uxe_ck[1]; uxe_ck[1] = uxe_ck[2]; uxe_ck[2] = t; }
              for (int j = 0; j < 3; j++) { assert(uxe[j] == uxe_ck[j]); }
            }
        }
    }

    /* Read footer: */
    filefmt_read_footer(rd, stmesh_STM_FILE_TYPE);
    return mesh;
  }
