#ifndef tabela_H
#define tabela_H

#include <float_image.h>
#include <r2.h>
#include <r3.h>
#include <r3x3.h>
#include <approx_system.h>
#include <ellipse_crs.h>
#include <argparser.h>

#define MAX_NC 3
  /* Maximum number of channels in images. */

typedef struct LinhaVetor LinhaVetor;
typedef struct LinhaTabela LinhaTabela;
typedef struct Tabela Tabela;

Tabela* aloca_tabela_vazia(int num_luzes, int num_linhas, r3_t view_dir);
/*Cria tabela sem inicializar dados internos*/

Tabela* cria_tabela
  ( float_image_t  *G[], 
    float_image_t  *M,
    double albedo, 
    int num_luzes, 
    int canal, 
    r2_t ponto[], 
    r3_t normal[], 
    int num_pontos,
    r3_t view_dir,
    bool_t interpolate_pixels
  );
/* Creates a table {tab} of observation vectors from a list of gauge
  images {G[0..num_luzes-1]} and a list of sampling points
  {ponto[0..num_pontos-1]} with their respective surface normals
  {normal[0..num_pontos-1]}.  Each observation vector is a vector
  {GO[0..num_luzes-1]} of float values, where {GO[i]} is the intensity
  in channel {canal} of image {G[i]}, in the range [0_1], interpolated
  at the sample point. Assumes that the gauge has uniform albedo {albedo}
  and is defined only where the (grayscale) mask image {M} is nonzero. */

Tabela* cria_tabela_virtual
  ( r3_t  directions[],
    double radius,r2_t centro,
    double albedo, 
    int num_luzes, 
    int canal, 
    r2_t ponto[], 
    r3_t normal[], 
    int num_pontos,
    r3_t view_dir
  );
/* Creates a table {tab} of observation vectors from a list of virtual gauge
  with a pontual light source with direction {direction} with their respective surface normals
  {normal[0..num_pontos-1]}.  Each observation vector is a vector
  {GO[0..num_luzes-1]} of float values, where {GO[i]} is the intensity
  in channel {canal} of image {G[i]}, in the range [0_1], interpolated
  at the sample point. Assumes that the gauge has uniform albedo {albedo}
 */
void gera_normais_na_esfera_uniforme(int resolucao, r3_t** normalP, int* num_pontosP, r3_t view_dir);
void gera_normais_no_gabarito(int resolucao,double thetaMax, r3_t** normalP, int* num_pontosP, r3_t view_dir);
void gera_pontos_no_gabarito_eliptico(int resolucao,double thetaMax, r2_t centro, double raio, r2_t estica, r2_t **pontoP, r3_t **normalP, int *num_pontosP,r3_t view_dir);
/*Generates a set of sample points on a eliptical image of a light gauge.

 The procedures assimes that the gaugeś projection has
center{centro} and radius {raio} and a stretch factor{estica}.
The center and estica are in the graphics coordinate system with y down.
The parameter {resolucao} controls the
 number of points generated.  The procedure returns in {*pontoP} the
 address of a vector of sample points (in the graphics coordinate system) , in {*normalP} 
 the address of a vector containing  the surface normals at those points, and in {*num_pontosP} the number of
 points (which is approximately {resolucao^2}.
 The normals are in matheatical coordinate system (y up) */

void gera_pontos_com_mascara(float_image_t* mask,r2_t centro, double raio, r2_t estica, r2_t **pontoP, r3_t **normalP, int *num_pontosP,r3_t view_dir);
void gera_pontos_no_gabarito_old(int resolucao, r2_t centro, double raio, r2_t **pontoP, r3_t **normalP, int *num_pontosP);
void gera_pontos_no_gabarito(int resolucao, r2_t centro, double raio, r2_t **pontoP, r3_t **normalP, int *num_pontosP);
/* Similar to above but with stretch {0,0}. */


float_image_t  *gera_mascara_do_gabarito(int nx, int ny, r2_t centro, double raio,r2_t stretch);
/* Generates a grayscale image with size {nx,ny} with stretch factor {stretch} where pixels that are inside
  the gauge outline have value 1, and pixels outside have value 0. */

const double *get_intdir(Tabela* tab, int linha);
/* Returns a pointer to the normalized observation vector in line {linha} of the 
  table {tab}. */

double get_intmag(Tabela* tab, int linha);
/* Returns a the L2 norm of the observation vector in line {linha} of the 
  table {tab}. */

void set_intdir(Tabela* tab, int linha,double g[]);
/* Copy the content of a pointer to the normalized observation vector in line {linha} of the 
  table {tab}. */

void set_intmag(Tabela* tab, int linha,double Gmag);
/* Sets a the L2 norm of the observation vector in line {linha} of the 
  table {tab}. */



r3_t get_view_dir(Tabela* tab);
/*Returns the associated view direction from table*/

r3_t compute_view_dir_from_stretch(r2_t gauge_stretch, double radius);
/*Computes camera view direction as seen from gauge center.
  It is important that the stretch vector {gauge_stretch} points towards
  the optical axis. Assumes that the stretch is in the graphics coordinate system (y down)
  and view_dir is in the mathematical coordinate system (y up)
*/



r3_t get_normal(Tabela* tab, int linha); 
/* Normal associada à linha da tabela. */

void set_normal(Tabela* tab, int linha, r3_t normal); 
/* Normal associada à linha da tabela. */

int get_num_linhas(Tabela* tab);
/* Returns the number of entries in table {tab}.  Was {get_tamanho_tabela}. */

void set_num_linhas(Tabela* tab,int num_linhas);
/* Sets the number of entries in table {tab}.  Was {get_tamanho_tabela}. */

int get_num_luzes(Tabela* tab);
/* Returns the number of lights in table {tab}. */

void print_linha(Tabela* tab, int linha);

int localiza_linha_por_normal(Tabela* tab, r3_t *sn);

void DumpTabela(Tabela* tab, char* fileName);

void ShowTableData(char* prefix, Tabela* tab,r2_t ponto[], r3_t normal[], int num_pontos);
/* Saves signature table {tab} information in file {prefix}_TableData.txt with
the associated set of points {ponto} and normals {normal}, where {num_pontos}
is the number of points and normals. Notice that the number of entries in a table
dont need to be equal to {num_pontos}. Only normals that are present in the table will be dumped.*/

void loadTableData(char* filename, Tabela** table,r2_t** pontos_P, r3_t** normals_P );
/* Loads a signature table {table},their associated points {pontos_P} and normals {normals_P} in gauge image from
 file {filename}. */

void SaveTable(char* filename,Tabela* tab, bool_t check);
/*Saves a signature table {tab} inside file {filename}. It is human-readable.
 If {check{ is true, checks that all normals are in the visible hemisphere, that
 all signature coordinates are non-negative, and that all magnitudes are 
 positive. */

void LoadTable(char* filename, Tabela** tab);
/*Loads a signature table {tab} from file {filename}.  */

Tabela* criaSubTabela(Tabela* tab,int subsetSize,int* subSets, int** index);
/* Returns a subtable using main table {tab} with {subsetSize} lights specified bt {subSets}
  the output equivalent index in the main table of the subtable is stored in {index}
*/

r3_t calcula_normal_esfera(r2_t *uv);
  /* Computes the outwards normal vector to the sphere's surface at
    the point that projects onto the point {uv}. Assumes a {U,V,W}
    coordinate system where the sphere has radius 1 and center at the
    origin, and {W} points towards the camera. */
void PrintTableStats(Tabela* tab);

void LiberaTabela(Tabela* tab);

Tabela* FixTableData(Tabela* tab);
/* Given a table, returns an allocated table wich contains no anomalies.*/

Tabela* cria_tabela_from_model(int num_lights,int resolution,r3_t view_dir,double thetaMax, approx_model_t* lm,void** l_data);
/*Creates a virtual-gauge table based on a lighting model*/

/*Gauge options */

typedef struct gauge_data_t
  { char *image;          /* Filename of image containing the gauge. */
    double gamma; 
    ellipse_crs_t E;      /* Ellipse that is the gauge's projection. */
    r3_t view;            /* Viewing direction at the gauge's center. */
    double albedo[MAX_NC]; /* Albedo of gauge in each color channel. */
    /* Only for the output image: */
    int magnify;          /* Magnification factor. */
    int margin;           /* Extra margin in pixels. */
    double trim;
    
  } gauge_data_t;
  /* Describes the geometry of a gauge's projection on an image. */

gauge_data_t *parse_gauge_args(argparser_t *pp, bool_t input);
  /* Parses the specification of a gauge.  
  
    If {input} is true, assumes an input gauge image spec: requires
    the "center" and "radius" keywords, allows "stretch" and
    "view", forbids "magnify" and "margin".
    
    If {input} is false, assumes an output gauge image spec:
    forbids "center", "radius", "stretch" and "view", 
    allows "magnify" and "margin".
    
    The "albedo" keyword is allowed in both cases 
    (defaults to 1.0). The "image" keyword is required in both cases. */

    r3_t gauge_normal_at_point(r2_t *q, ellipse_crs_t *E);
  /* Returns the gauge's view-relative normal at the point {q}.  In particular,
    returns {(0,0,1)} for the center of {E}. */

  
r3_t gauge_normal_in_pixel(int x, int y, ellipse_crs_t *E);
  /* Returns the mean view-relative normal of the gauge's surface
    in the pixel whose lower left corner is {(x,y)}.  In particular,
    returns {(0,0,1)} for a pixel whose center is the center of {E}.  */
  
  int pixel_position_in_gauge_image(int x, int y, ellipse_crs_t *E);
  /* Position of the pixel whose lower left corner is {(x,y)} relative to the gauge's projection.
    Returns {+1} if totally inside, {-1} if totally outside, and 0 if the pixel straddles
    the projection's outline (or too close to tell). */
  
  double gauge_coverage_of_pixel(int x, int y, ellipse_crs_t *E);
  /* Returns 1.0 if the pixel with lower left corer {(x,y)} is entirely inside the gauge's
    projection; 0.0 if it is entirely outside; and a number strictly between 0 and 1 
    if the pixel straddles the projection's boundary. */
  
  
void extract_data_from_gauge_image
  ( float_image_t *img,
    float_image_t *xtr,
    float_image_t* mask,
    ellipse_crs_t *E, 
    double cutRadius,
    r3_t *view, 
    r3_t **XP,
    double ***FP,
    int *NPP
  );
  /* Enumerates the pixels of image {img} that lie on a spherical
    gauge described by {*ga}. Returns the corresponding surface
    normals {(*XP)[0..NP-1] and image intensities
    {(*FP)[0..NC-1][0..NP-1]}, where {NC} is the number of channels in
    the image {img} and {NP} is the number of said pixels, which is
    returned in {*NPP}. 
    
    The normals are adjusted so that the center of {E} has normal
    {view}. The data arrays are allocated by the procedure.
    
    If cutRadius is bigger than 0, the samples with distance less than ${cutRadius}
    will have cover zero assigned. 
    
    If {xtr} is not NULL, it should be a single channel image. In this
    case the procedure stores into it a coverage mask, where each
    sample is the area fraction {cov} (in {[0_1]}) of the corresponding
    pixel of {img} that is covered by the gauge's projection. If {cov}
    is neither 0 nor 1, it is only an approximation computed by
    sampling. */
  
  float_image_t *read_gauge_image(char *filename, double gamma);
  /* Reads a gauge image from file {filename}, which must have
    extensions ".fni", ".ppm", or ".pgm".  The {gamma} is used
    for PNM images. */

  
    void synthetic_gauge_color_at_point
  ( r2_t *p,
    ellipse_crs_t *E, 
    r3x3_t *VM,
    int NC,
    void** l_data,
    evaluate_function* shading,
    double value[]  /* (OUT) */
  );
  /* For each channel {c} in {0..NC-1}, computes the apparent color
    (per-channel radiance) {value[c]} of the synthetic gauge or
    background at point {p} of the synthetic image. Uses a given shading function {shading}
    for a given lighting model ${l_data}*/
  
  
  void synthetic_gauge_color_in_pixel
  ( int x,
    int y,
    ellipse_crs_t *E, 
    r3x3_t *VM,
    int NC,
    void** l_data,
    evaluate_function* shading,
    double value[]  /* (OUT) */
  );
  /* Computes the mean apparent color (per-channel radiance)
    of the synthetic image in the pixel with bottom left corner {(x,y)}. */
  
  void paint_synthetic_gauge_image
  ( float_image_t *img,
    ellipse_crs_t *E, 
    r3_t *view, 
   void** l_data,
    evaluate_function* shading
  );
  /* Paints into each channel {c} of image {img} a painting of a
    virtual gauge, given its projected {center} and {radius} in
    pixels, the viewing direction {view}, albedo {D[c]},
    and light flow {L[c]}. */


#endif
