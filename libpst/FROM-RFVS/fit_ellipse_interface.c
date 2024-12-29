/* This program uses the y axis orientation provided by grpahics cards, it means, y values grow down. The program
assums that the image has not been cropped, so that the optical center is the image center. The signs of the stretch parameter 
are adjusted so that the stretch points towards the optical center.*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <float_image.h>
#include <float_image_paint.h> 
#include <float_image_gradient.h>
#include <float_pnm_image_io.h>
#include <ellipse_ouv.h>
#include <ellipse_crs.h>
#include <ellipse_crs_args.h>
#include <math.h>
#include <jsfile.h>
#include <pst_fit_ellipse.h>
#include <argparser.h>

#include <SDL.h>

#include "SDL/SDL_gfxPrimitives.h"

#define PROG_NAME "fit_ellipse_interface"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -nImages {NUM} \\\n" \
  "    -prefix {FILE_PREFIX} \\\n" \
  "    -noise {NUM} \\\n" \
  "    -images {FILENAME}...  \\\n" \
  "    -cameraDist {CAMDIST_PX} \\\n"  \
  "    [ -width {WIDTH} ] \\\n" \
  "    [ -height {HEIGHT} ] \\\n" \
  "    [ -bpp {BPP} ] \\\n" \
  "    [ -warp ] \\\n" \
  "    [ -hw ] \\\n" \
  "    [ -autoSize ] \\\n" \
  "    [ -fullScreen ] \\\n" \
  "    [ -yflip ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  Etc. etc..\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "OPTIONS" \
  "  Etc. etc.."



char debug[500];

struct Ellipse_t{
	r2_t center;	
	ellipse_ouv_t el;
	int status; //0 - center not defined, 1 - center only defined, 2 - one axis defined, 3- both axes (complete)
	int type; //1 - ellipse, 2 - circle , 3 - sphere
};

typedef struct Ellipse_t ellipse_t;


struct Button_t{
	char* name;
	SDL_Surface* figure;
	SDL_Rect rect;
	int status; //0 normal, 1 higlighted, 2 pressed
};

typedef struct Button_t button_t;


 struct SDL_Options_t{
	int w, h;
	int desired_bpp;
	Uint32 video_flags;
	int screen_auto_size;
 };

typedef struct SDL_Options_t SDL_options_t;

struct Workspace_Options_t{
	char* outPrefix;
	char** float_image_filenames;
	double* image_noise;
	double gamma;
        double cameraDist;
	int image_count; //total of images loaded (REPEATED)
	bool_t yflip;
};

typedef struct Workspace_Options_t workspace_options_t;


struct WorkSpace_t{
	SDL_Surface* image_surface; //surface of current image (performance purposes)
	int current_image; //index of image being worked
	int image_count; //total of images loaded
	
	SDL_Rect image_rect; //defines area where image will be drawed
	SDL_Rect display_rect; //defines area where image will be displayed on screen
	SDL_Rect board_rect; //static are where user can clik in image
	SDL_Rect panel_rect; //defines control panel area
	int pan_x,pan_y;
	float_image_t** float_image_list; //list of float images
	float_image_t** gradient_list; //list of gradient of the images above
	int drag_disable;

	button_t** button_list;
	int button_count;

	/*interface options*/
	int selected_action;
	int status_action;
	int show_gradient;

	ellipse_t** ellipse_list; //first index (image), second (elipses)
	int *ellipse_count; //first index (image)

	workspace_options_t* opts;
	
};

typedef struct WorkSpace_t workspace_t;

/*Ellipse Functions*/

void DrawEllipse(SDL_Surface* screen, ellipse_t* el, int dx,int dy);
/* Draws in a SDL_Surface{screen} a ellipse with translation {dx}{dy}
   If ellipse have only U vector defined only a line will be draw.
*/

void InsertEllipse(ellipse_t** ellipse_list,int* list_size,ellipse_t elem);
/*Inserts a ellipse{elem} in the last position of a array of ellipses{ellipse_list}
with size {listsize}. The array is reallocated to fit to the new size
*/

void RemoveEllipse(ellipse_t** ellipse_list,int* list_size,int index);
/* Remove a ellipse in position{index} of a list {ellipse_list}. The array is
reallocated to fit the new size. If the array reaches 0 size it receives NULL value
*/

/*Button Functions*/

button_t* CreateButton(char* name, SDL_Surface* figure, int x, int y);
/* Create a button structure with name {name}, associating a figure {figure}
   to the button. A SDL_Rect is created internaly with position at {x,y}. The
   width and heigth of the rec is determined by the size of the name
*/
void DrawButton(SDL_Surface* screen, button_t* bt);
/* Draws a button in a SDL_Surface {screen} according  button attributes and status
*/

void BlitFloatImage(SDL_Surface* screen,float_image_t* img,int dx,int dy);
/*  Blits a float image {img} in a SDL_Surface {screen} with translation {dx,dy}
    To avoid compatibility problems with pixel format the function SDL_FillRect is used for
    each pixel of the image. It means... it is SLOW and should ot be used
    to draw at real-time execution.
*/

/*Workspace Functions*/

void SetWorkspaceImage(workspace_t* wspc, int index);

workspace_t* CreateWorkspace(int screen_w, int screen_h, workspace_options_t* opts,float_image_t** image_list);

void SetWorkspacePan(workspace_t *wspc, int dx, int dy);

void DrawWorkspace(SDL_Surface* screen, workspace_t* wspc);

void ChangeSelectedAction(workspace_t* wspc, int num_action);

/*Events Management*/

void HandleEvent(workspace_t* wspc);

void HandleKeyDown_Pan(SDL_Event* event, workspace_t* wspc);

void HandleKeyDown_ChangeImage(SDL_Event* event, workspace_t* wspc);

void HandleMouseMotion_Pan(SDL_Event* event, workspace_t* wspc);

void HandleMouseButtonDown_Panel(SDL_Event* event, workspace_t* wspc);

void HandleMouseButtonDown_Ellipse(SDL_Event* event, workspace_t* wspc);

void HandleMouseButtonDown_Circle(SDL_Event* event, workspace_t* wspc);

void HandleMouseButtonDown_Sphere(SDL_Event* event, workspace_t* wspc);

void HandleMouseMotion_Ellipse(SDL_Event* event, workspace_t* wspc);

void HandleMouseMotion_Circle(SDL_Event* event, workspace_t* wspc);

void HandleMouseMotion_Sphere(SDL_Event* event, workspace_t* wspc);



/*Ellipse fitting functions*/

void AdjustEllipseToFloatImage(workspace_t* wspc);

void fep_fit_ellipse  ( float_image_t *IMG,     ellipse_crs_t *E,     double noise  , int type);

void fep_write_ellipse_params(char *prefix, char *tag, r2_t gab_center,double gab_radius,  r2_t stretch, r3_t view_dir);

void fep_write_image_as_pnm(char *name, float_image_t *A);

float_image_t *fep_compute_msk_image(ellipse_crs_t *E, int NX, int NY);

void fep_write_sub_image_as_pnm
  ( char *prefix, 
    char *tag, 
    float_image_t *A,
    int xLo,
    int xHi,
    int yLo, 
    int yHi
  );


float_image_t *fep_compute_msk_image(ellipse_crs_t *E, int NX, int NY);
  /* Computes a monochromatic image {MSK} with size {NX} by {NY} where
    the value of each pixel is the fraction of its area that is
    covered by the ellipse described in {E}. */

void printCutParams(char *prefix, char* tag, int W, int H, int X, int Y);

void computeViewDirAndStretch(r2_t gab_center,double radius,r2_t cam_center, double cam_dist,r3_t* view_dirP, r2_t* stretchP);

/*Misc Functions*/
void printGaugeParams(FILE* arq,r2_t gab_center,double gab_radius,  r2_t stretch, r3_t view_dir);


float_image_t** ComputeGradients(float_image_t* image_list[],int image_count,double image_noise[]);

void GetOptions(int argc, char** argv,SDL_options_t* sdl_opts, workspace_options_t* wspc_opts);

void ClearScreen(SDL_Surface *screen);

void MainLoop(SDL_Surface* screen,workspace_options_t* wspc_opts,float_image_t** image_list);

SDL_Surface* InitializeSDL(SDL_options_t* opts);

/*Implementation*/


int main ( int argc, char *argv[] )
{

	SDL_options_t sdl_opts;
	workspace_options_t wspc_opts;

	SDL_Surface *screen;
	/* Title */
	fprintf (stderr,"Fit Ellipse\n");
	GetOptions(argc, argv,&sdl_opts,&wspc_opts);
	float_image_t** image_list;
	image_list = float_pnm_image_list_read(wspc_opts.image_count,wspc_opts.float_image_filenames,FALSE,wspc_opts.gamma,VIEW_BIAS,FALSE,TRUE,TRUE);
	
	if(sdl_opts.screen_auto_size == 1){
		int maxW,maxH;
		int i;
		maxH = -1;
		maxW = -1;
		for(i = 0; i < wspc_opts.image_count; i++){
			float_image_t* img = image_list[i];
			int w = img->sz[1];
			int h = img->sz[2];
			if(w > maxW){
				maxW = w;
			}
			if(h > maxH){
				maxH = h;
			}
		}
		maxH = maxH + 50;
		fprintf(stderr,"Screen resized to %d %d\n",maxW,maxH);
		sdl_opts.w = maxW;
		sdl_opts.h = maxH;
	}
	screen = InitializeSDL(&sdl_opts);
	MainLoop(screen,&wspc_opts,image_list);
	
	return(0);
}

SDL_Surface* InitializeSDL(SDL_options_t* opts){
	SDL_Surface* screen = NULL;
	
	if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
		fprintf(stderr,
			"Couldn't initialize SDL: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	/* Initialize the display */
	screen = SDL_SetVideoMode(opts->w, opts->h, opts->desired_bpp, opts->video_flags);
	if ( screen == NULL ) {
		fprintf(stderr, "Couldn't set %dx%dx%d video mode: %s\n",opts->w, opts->h, opts->desired_bpp, SDL_GetError());
		return NULL;
	}

	fprintf(stderr,"Set %dx%dx%d mode\n",screen->w, screen->h, screen->format->BitsPerPixel);
	fprintf(stderr,"Video surface located in %s memory.\n",	(screen->flags&SDL_HWSURFACE) ? "video" : "system");
	
	/* Check for double buffering */
	if ( screen->flags & SDL_DOUBLEBUF ) {
		fprintf(stderr,"Double-buffering enabled - good!\n");
	}

	/* Set the window manager title bar */
	SDL_WM_SetCaption("Fit Ellipse to Photo", "Fit Ellipse to Photo");

	return screen;
}

void MainLoop(SDL_Surface* screen,workspace_options_t* wspc_opts,float_image_t** image_list){
	/* Black screen */
	ClearScreen(screen);
	stringRGBA(screen,0,0,"loading...",255,255,255,255);
	SDL_Flip(screen);
	workspace_t* w = CreateWorkspace(screen->w, screen->h, wspc_opts,image_list);
	while(1) {
		HandleEvent(w);
		DrawWorkspace(screen, w);
		stringRGBA(screen,0,0,debug,255,0,0,255);
		SDL_Flip(screen);
		SDL_Delay(10);
	}
}

workspace_t* CreateWorkspace(int screen_w, int screen_h, workspace_options_t* opts,float_image_t** image_list){
	workspace_t* novo;
	
	if(opts->image_count < 1){
		fprintf(stderr,"No images to be worked!\n");
		return NULL;
	}
	if( (screen_h < 50) || (screen_w < 320) ){
		fprintf(stderr,"No Enough Space to draw a workspace !\n");
		return NULL;
	}
	
	novo = (workspace_t*)malloc(sizeof(workspace_t));
	
	novo->panel_rect.x = 0;
	novo->panel_rect.y = screen_h - 50;
	novo->panel_rect.w = screen_w;
	novo->panel_rect.h = 50;
	
	novo->display_rect.x = 0;
	novo->display_rect.y = 0;
	novo->display_rect.w = screen_w;
	novo->display_rect.h = screen_h - 50;

	novo->board_rect = novo->display_rect;
	novo->opts = opts;
	novo->image_count = opts->image_count;
	novo->image_surface = NULL;
	
	fprintf(stderr,"Reading Image List of %d images\n",opts->image_count);
	int i;
	//novo->float_image_list = float_pnm_image_list_read(img_count,img_filelist,1.0,0.0,FALSE,TRUE,TRUE);
	novo->float_image_list = image_list;
	novo->gradient_list = ComputeGradients(image_list,opts->image_count,opts->image_noise);
	fprintf(stderr,"Seting Workspace Image\n");
	novo->show_gradient = 0;
	SetWorkspaceImage(novo,0);
	
	novo->drag_disable = 0;
	fprintf(stderr,"Seting Buttons\n");
	char *buttons_names[5] = {
 	"Pan",
 	"Draw Ellipse",
 	"Draw Circle",
 	"Draw Sphere",
 	"Fit Form"
	};
	novo->button_count = 5;
	novo->button_list = (button_t**)malloc(sizeof(button_t*)*(novo->button_count));
	fprintf(stderr,"List Created\n");
	int advanced = 10;
	for(i = 0; i < novo->button_count; i++){
		int x = advanced;
		int y =  novo->panel_rect.y + 10;
		novo->button_list[i] = CreateButton(buttons_names[i], NULL,x,y);	
		advanced += novo->button_list[i]->rect.w +10;
	}
	
	novo->selected_action = 0;
	novo->status_action = 0;

	
	novo->ellipse_count = (int*)malloc(sizeof(int)*opts->image_count);
	novo->ellipse_list = (ellipse_t**)malloc(sizeof(ellipse_t*)*opts->image_count);
	for(i = 0; i < opts->image_count;i++){
		novo->ellipse_count[i] = 0;	
		novo->ellipse_list[i] = NULL;
	}
	debug[0] = '\0';
	
	return novo;
}

void HandleEvent(workspace_t* wspc)
{
	SDL_Event event; 
	int mx,my;

	/* Check for events */
        while ( SDL_PollEvent(&event) ) {
                        switch (event.type) {
			 case SDL_KEYDOWN:
				HandleKeyDown_ChangeImage(&event,wspc);
				if(wspc->selected_action == 0){
					HandleKeyDown_Pan(&event,wspc);
				}
				break;
			case SDL_MOUSEBUTTONDOWN:
				
 				mx = event.button.x;
				my = event.button.y;
				if(  (mx < wspc->board_rect.x) || (mx > (wspc->board_rect.x + wspc->board_rect.w))
				   || (my < wspc->board_rect.y) || (my > (wspc->board_rect.y + wspc->board_rect.h))
				){
					wspc->drag_disable = 1;
					//fprintf(stderr,"Out of board area\n");
					HandleMouseButtonDown_Panel(&event,wspc);
				}
				else{
					wspc->drag_disable = 0;
				}
				//puts("here");
				if(wspc->drag_disable == 0){
 					if(wspc->selected_action == 1){
						HandleMouseButtonDown_Ellipse(&event,wspc);
					}
					if(wspc->selected_action == 2){
						HandleMouseButtonDown_Circle(&event,wspc);
					}
					if(wspc->selected_action == 3){
						HandleMouseButtonDown_Sphere(&event,wspc);
					}	
					if(wspc->selected_action == 4){
						AdjustEllipseToFloatImage(wspc);
					}
				}
				//puts("here");
			break;
			case SDL_MOUSEMOTION:
				if(wspc->selected_action == 0){
					HandleMouseMotion_Pan(&event,wspc);
				}
				if(wspc->selected_action == 1){
					HandleMouseMotion_Ellipse(&event,wspc);
				}
				if(wspc->selected_action == 2){
					HandleMouseMotion_Circle(&event,wspc);
				}
				if(wspc->selected_action == 3){
					HandleMouseMotion_Sphere(&event,wspc);
				}
				break;
			 case SDL_QUIT:
                                        exit(0);
                                        break;
			}
	}
}


void DrawEllipse(SDL_Surface* screen, ellipse_t* el, int dx,int dy){
	SDL_Rect temp_rect;
	Uint32 color;
	sprintf(debug,"Center:[%f %f]\n U:[%f %f] A:%f V:[%f %f] B:%f",el->center.c[0],el->center.c[1],el->el.u.c[0], el->el.u.c[1],el->el.a, el->el.v.c[0], el->el.v.c[1], el->el.b);
	if(el->status > 0){
		//draw center
		color = SDL_MapRGB (screen->format, 255,0,0);
		temp_rect.x = el->center.c[0] + dx;
		temp_rect.y = el->center.c[1] + dy;
		temp_rect.w = 1;
		temp_rect.h = 1;
		SDL_FillRect(screen,&(temp_rect),color);	
	}
	//if status is 1 draw only a line towards center
	if((el->status >= 1) &&  (el->status <= 2)){
		//puts("status 1");
		int x1,x2,y1,y2;
		double a = el->el.a;
		x1 = el->center.c[0] +dx;
		y1 = el->center.c[1] +dy;
		x2 = (el->el.u.c[0]*a) + x1;
		y2 = (el->el.u.c[1]*a) + y1;
		
		lineRGBA(screen, x1,y1,x2,y2, 0,255, 0, 255);
	}
	if(el->status >= 2){
		//draw whole thing
		//puts("status 2");
		double a = el->el.a;
		double b = el->el.b;
		r2_t u = el->el.u;
		r2_t v = el->el.v;
		double passo = 0.75/a;
		double dt;
		
		if(passo <= 0.001) return;
		int x1,x2,y1,y2;
		x1 = el->center.c[0] +dx;
		y1 = el->center.c[1] +dy;
		x2 = (el->el.v.c[0]*b) + x1;
		y2 = (el->el.v.c[1]*b) + y1;
		if(el->status != 3) lineRGBA(screen, x1,y1,x2,y2, 255,0, 0, 255);
		int times = 0;
		//puts("Lined");
		color = SDL_MapRGB (screen->format, 255,0,0);
		for(dt = 0; dt < (2.0*M_PI); dt = dt + passo){
			temp_rect.x = el->center.c[0] + dx + (a*cos(dt)*u.c[0]) + (b*sin(dt)*v.c[0]);
			temp_rect.y = el->center.c[1] + dy + (a*cos(dt)*u.c[1]) + (b*sin(dt)*v.c[1]);
			temp_rect.w = 1;
			temp_rect.h = 1;
			SDL_FillRect(screen,&(temp_rect),color);
			times++;
			//printf("dt-%d\n",times);
		}
		//puts("end loop");
		
	}
	
}


void InsertEllipse(ellipse_t** ellipse_list,int* list_size,ellipse_t elem){
	fprintf(stderr,"Inserting ellipse %d\n",*list_size);
	if(*ellipse_list == NULL){
		*ellipse_list = (ellipse_t*)malloc(sizeof(ellipse_t));
		(*ellipse_list)[0] = elem;
		
	}else{
		int tam = *list_size;
		*ellipse_list = (ellipse_t*)realloc(*ellipse_list,sizeof(ellipse_t)*(tam+1));
		(*ellipse_list)[*list_size] = elem;
	}
	*list_size = *list_size+1;
}

void RemoveEllipse(ellipse_t** ellipse_list,int* list_size,int index){
	fprintf(stderr,"Removing ellipse %d\n",index);
	if(index < 0) return;
	if(index >= (*list_size)) return;
	ellipse_t* new_list;
	if(*list_size > 1){
		int i;
		int t = *list_size -1;
		new_list = (ellipse_t*)malloc(sizeof(ellipse_t)*t);
		int count = 0;
		for(i =0 ; i < *list_size; i++){
			if(i != index){
				new_list[count] = (*ellipse_list)[i];
				count++;
			}
		}
		free(*ellipse_list);
		*ellipse_list = new_list;
		*list_size = t;
	}else{
		free(*ellipse_list);
		*ellipse_list = NULL;
		*list_size = 0;
	}
	
}




button_t* CreateButton(char* name, SDL_Surface* figure, int x, int y){
	button_t* novo;
	novo = (button_t*) malloc(sizeof(button_t));
	novo->name = (char*)malloc(sizeof(char)*(strlen(name) +1));
	//puts("Aqui");
	strcpy(novo->name,name);
	//fprintf(stderr,"Button Name %s\n",name);
	novo->figure = figure;
	novo->rect.x = x;
	novo->rect.y = y;
	novo->rect.w = (strlen(name)*8) + 20;
	novo->rect.h = (8) + 20;
	novo->status = 0;
	return novo;
}

void DrawButton(SDL_Surface* screen, button_t* bt){
	
	SDL_Rect temp_rect;
	Uint32 color;
	color = SDL_MapRGB (screen->format, 0,0,0);
	SDL_FillRect(screen,&(bt->rect),color);
	temp_rect.x = bt->rect.x + 1;
	temp_rect.y = bt->rect.y + 1;
	temp_rect.h = bt->rect.h - 2;
	temp_rect.w = bt->rect.w - 2;
	if(bt->status == 0){
		color = SDL_MapRGB (screen->format, 50,50,50);
		
	}else{
		color = SDL_MapRGB (screen->format, 110,110,110);
	}
	int R,G,B;
	R = B = G = 255;
	if(bt->status == 2){
		R = G = B = 0;
	}
	SDL_FillRect(screen,&(temp_rect),color);
	stringRGBA(screen,bt->rect.x +10, bt->rect.y +10,bt->name,R,G,B,255);
	
}






void SetWorkspaceImage(workspace_t* wspc, int index){

	assert(index >= 0);
	assert(index < wspc->image_count);
	wspc->pan_x = 0;
	wspc->pan_y = 0;
	float_image_t* fimg;
	if(wspc->show_gradient == 0){
		fimg= wspc->float_image_list[index];
	}else{
		fimg= wspc->gradient_list[index];
	}	

	if(wspc->image_surface != NULL){
		SDL_FreeSurface(wspc->image_surface);
	}
	int nx = fimg->sz[1];
	int ny = fimg->sz[2];

	wspc->image_surface = SDL_CreateRGBSurface(SDL_HWSURFACE,nx,ny,32,0,0,0,0);
	BlitFloatImage(wspc->image_surface,fimg,0,0);
	
	wspc->image_rect.x = 0;
	wspc->image_rect.y = 0;
	//wspc->image_rect.w = wspc->display_rect.w;
	//wspc->image_rect.h = wspc->display_rect.h;
	wspc->image_rect.w = nx;
	wspc->image_rect.h = ny;
	
	wspc->current_image = index;
	
	SetWorkspacePan(wspc,-wspc->pan_x,-wspc->pan_y);

}

float_image_t** ComputeGradients(float_image_t* image_list[],int image_count,double image_noise[]){
	int i;
	float_image_t** gradient_list = (float_image_t **) malloc(sizeof(float_image_t*)*image_count);
	for(i = 0; i < image_count;i++){
		gradient_list[i] = float_image_gradient_sqr_relative(image_list[i], image_noise[i], TRUE);
		assert(gradient_list[i]->sz[0] == 1);
		float vMax = 1.0e-100; /* Just in case {REF} is all zeros. */
        	float_image_update_sample_range(gradient_list[i], 0, NULL, &vMax);
		float_image_rescale_samples(gradient_list[i], 0, 0.0, vMax, 0.0, 0.80);
	}
	
	return gradient_list;
    
	
}



void SetWorkspacePan(workspace_t *wspc, int dx, int dy){
	wspc->pan_x += dx;
	wspc->pan_y += dy;
	
// 	if(wspc->pan_x > 0){
// 		wspc->image_rect.x = wspc->pan_x;
// 		wspc->display_rect.x = 0;
// 	}
// 	else{
// 		wspc->display_rect.x = - wspc->pan_x ;
// 		wspc->image_rect.x = 0;
// 	}
// 	
// 	if(wspc->pan_y > 0){
// 		wspc->image_rect.y = wspc->pan_y;
// 		wspc->display_rect.y = 0;
// 	}
// 	else{
// 		wspc->display_rect.y = -wspc->pan_y;
// 		wspc->image_rect.y = 0;
// 	}
	if(wspc->pan_x >= 0){
		wspc->display_rect.x = wspc->pan_x ;
		wspc->image_rect.x = 0;
	}
	else{
		wspc->image_rect.x = -wspc->pan_x ;
		wspc->display_rect.x = 0;
	}
	if(wspc->pan_y >= 0){
		wspc->display_rect.y = wspc->pan_y ;
		wspc->image_rect.y = 0;
	}
	else{
		wspc->image_rect.y = -wspc->pan_y ;
		wspc->display_rect.y = 0;
	}

	
}

void DrawWorkspace(SDL_Surface* screen, workspace_t* wspc){
	//fill it with black tint
	Uint32 color;
	color = SDL_MapRGB (screen->format, 0,0, 0);
	SDL_FillRect(screen,NULL,color);
	//draw image being worked
	//SDL_BlitSurface(wspc->image_surface,&(wspc->image_rect),screen,&(wspc->display_rect));
	SDL_BlitSurface(wspc->image_surface,&(wspc->image_rect),screen,&(wspc->display_rect));
	
	int i;
	for(i = 0; i < wspc->ellipse_count[wspc->current_image]; i++){
	//	if((wspc->selected_action == 1) ||(wspc->selected_action == 2)){
			//printf("ploting %d circle\n",i);
			ellipse_t* el = &(wspc->ellipse_list[wspc->current_image][i]);
			DrawEllipse(screen,el,wspc->pan_x,wspc->pan_y);
			//puts("ok");
	//	}
	}
	//we need draw the menu !
	color = SDL_MapRGB (screen->format, 50,50,50);
	SDL_FillRect(screen,&(wspc->panel_rect),color);
	//buttons
	
	for(i = 0; i < wspc->button_count; i++){
		if(wspc->selected_action == i){
			wspc->button_list[i]->status=2;
		}else{
			wspc->button_list[i]->status=0;
		}
		DrawButton(screen,wspc->button_list[i]);
	}
}


void BlitFloatImage(SDL_Surface* screen,float_image_t* img,int dx,int dy){
	int i,j;
	SDL_Rect rec;
	Uint32 color;
	rec.x = 0;
	rec.y = 0;
	rec.w = 1;
	rec.h = 1;
	color = SDL_MapRGB (screen->format, 0,0, 0);
	SDL_FillRect (screen, NULL, color);
	for(i = 0; i < img->sz[1]; i++){
		for(j = 0; j < img->sz[2]; j++){
			rec.x = i + dx;
			rec.y = j + dy;
			if( (rec.x >= 0) && (rec.x < screen->w) && (rec.y >= 0) && (rec.y < screen->h) ){
			//	printf("%d %d\n",rec.x,rec.y);
				double R,G,B;
				R = float_image_get_sample(img, 0,i,j)*255;
				if( img->sz[0] == 3){
					G = float_image_get_sample(img, 1,i,j)*255;
					B = float_image_get_sample(img, 2,i,j)*255;
				}else{
					G = B = R;
				}
				color = SDL_MapRGB (screen->format, R,G, B);
				SDL_FillRect (screen, &rec, color);
			}else{
					//printf("%d %d\n",rec.x,rec.y);
			}
		}	
	}
}


void HandleKeyDown_Pan(SDL_Event* event, workspace_t* wspc){
	if(event->key.keysym.sym == SDLK_DOWN){
		SetWorkspacePan(wspc,0,1);
	}
	if(event->key.keysym.sym == SDLK_UP){
		SetWorkspacePan(wspc,0,-1);
	}
	if(event->key.keysym.sym == SDLK_RIGHT){
		SetWorkspacePan(wspc,1,0);
	}
	if(event->key.keysym.sym == SDLK_LEFT){
		SetWorkspacePan(wspc,-1,0);
	}
	if(event->key.keysym.sym == SDLK_g){
		wspc->show_gradient = 1 - wspc->show_gradient;
		SetWorkspaceImage(wspc,wspc->current_image);
	}
}

void HandleKeyDown_ChangeImage(SDL_Event* event, workspace_t* wspc){
	int next;
	if(event->key.keysym.sym == SDLK_PAGEDOWN){
		next = wspc->current_image + 1;
		if(next >= wspc->image_count) next = 0;
		if(next == wspc->current_image) return;
		SetWorkspaceImage(wspc,next);
		
	}
	else if(event->key.keysym.sym == SDLK_PAGEUP){
	        next = wspc->current_image - 1;
		if(next < 0) next = wspc->image_count -1;
		if(next == wspc->current_image) return;
		SetWorkspaceImage(wspc,next);
	}
	
}

void HandleMouseMotion_Pan(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	
	if(event->motion.state&&SDL_BUTTON(1)){
		if(wspc->drag_disable == 0){
			SetWorkspacePan(wspc,event->motion.xrel,event->motion.yrel);
		}
	}
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	sprintf(debug,"Pan %d %d Mouse [%d %d] [%d %d]",wspc->pan_x,wspc->pan_y,mx,my,x,y);
}

void ChangeSelectedAction(workspace_t* wspc, int num_action){
	wspc->selected_action = num_action;
	wspc->status_action = 0;
	sprintf(debug," ");

}


void HandleMouseButtonDown_Panel(SDL_Event* event, workspace_t* wspc){
	int i;
	int mx = event->button.x;
	int my = event->button.y;
	for(i = 0; i < wspc->button_count; i++){
		SDL_Rect rect = wspc->button_list[i]->rect;
		if( (mx >= rect.x) && (mx <= (rect.x + rect.w)) &&
		    (my >= rect.y) && (my <= (rect.y + rect.h)) ){
			 ChangeSelectedAction(wspc,i);
		}
	}
}


void HandleMouseButtonDown_Ellipse(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	int current_el = wspc->ellipse_count[wspc->current_image]-1;
	fprintf(stderr,"Event: HandleMouseButtonDown_Ellipse\n");
	if(event->button.button == SDL_BUTTON_LEFT){
		switch(wspc->status_action ){
			case 0:
				{
					ellipse_t el;
					el.center.c[0] = x;
					el.center.c[1] = y;
					el.status= 1;
					el.type = 1;
					InsertEllipse(&(wspc->ellipse_list[wspc->current_image]),&(wspc->ellipse_count[wspc->current_image]),el);
					wspc->status_action = 1;
					
				}
				break;
			case 1:
				//direction of long axis is where mouse points
				wspc->ellipse_list[wspc->current_image][wspc->ellipse_count[wspc->current_image]-1].status = 2;
				wspc->status_action = 2;
				//puts("2");
				break;
			case 2:
				//direction of long axis is where mouse points
				{
					ellipse_ouv_t* el = &(wspc->ellipse_list[wspc->current_image][current_el].el);
					if(el->a < el->b){
						double temp = el->a;
						el->a = el->b;
						el->b = temp;
						r2_t vect = el->u;
						el->u = el->v;
						el->v = vect;
					}
				}
				wspc->ellipse_list[wspc->current_image][wspc->ellipse_count[wspc->current_image]-1].status = 3;
				wspc->status_action = 0;
				//puts("0");
				break;
		}
	}else if( event->button.button == SDL_BUTTON_RIGHT){
		RemoveEllipse(&(wspc->ellipse_list[wspc->current_image]),&(wspc->ellipse_count[wspc->current_image]),wspc->ellipse_count[wspc->current_image]-1);
		wspc->status_action = 0;
	}

}

void HandleMouseButtonDown_Circle(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	ellipse_t el;
	fprintf(stderr,"Event: HandleMouseButtonDown_Circle\n");
	if(event->button.button == SDL_BUTTON_LEFT){
		switch(wspc->status_action ){
			case 0:
// 				if(wspc->ellipse_list == NULL){
// 					wspc->ellipse_list = (ellipse_t*)malloc(sizeof(ellipse_t));
// 				}
// 				wspc->ellipse_list[0].center.c[0] = x;
// 				wspc->ellipse_list[0].center.c[1] = y;
// 				wspc->status_action = 2;
// 				wspc->ellipse_count = 1;
// 				fprintf(stderr,"Center %d %d \n",x,y);
				el.center.c[0] = x;
				el.center.c[1] = y;
				el.status= 2;
				el.type = 2;
				InsertEllipse(&(wspc->ellipse_list[wspc->current_image]),&(wspc->ellipse_count[wspc->current_image]),el);
				wspc->status_action = 2;
				//puts("2");
				break;
			case 2:
				//direction of long axis is where mouse points
				wspc->ellipse_list[wspc->current_image][wspc->ellipse_count[wspc->current_image]-1].status = 3;
				wspc->status_action = 0;
			//	puts("0");
				break;
		}
	}else if( event->button.button == SDL_BUTTON_RIGHT){
		RemoveEllipse(&(wspc->ellipse_list[wspc->current_image]),&(wspc->ellipse_count[wspc->current_image]),wspc->ellipse_count[wspc->current_image]-1);
		wspc->status_action = 0;
	}
}


void HandleMouseButtonDown_Sphere(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	fprintf(stderr,"Event: HandleMouseButtonDown_Ellipse\n");
	if(event->button.button == SDL_BUTTON_LEFT){
		switch(wspc->status_action ){
			case 0:
				{
					ellipse_t el;
					el.center.c[0] = x;
					el.center.c[1] = y;
					el.status=2;
					el.type = 3;
					InsertEllipse(&(wspc->ellipse_list[wspc->current_image]),&(wspc->ellipse_count[wspc->current_image]),el);
					wspc->status_action = 2;
				}
				break;
			case 2:
				wspc->ellipse_list[wspc->current_image][wspc->ellipse_count[wspc->current_image]-1].status = 3;
				wspc->status_action = 0;
				//puts("0");
				break;
		}
	}else if( event->button.button == SDL_BUTTON_RIGHT){
		RemoveEllipse(&(wspc->ellipse_list[wspc->current_image]),&(wspc->ellipse_count[wspc->current_image]),wspc->ellipse_count[wspc->current_image]-1);
		wspc->status_action = 0;
	}
}
 
void HandleMouseMotion_Ellipse(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	if(wspc->status_action  == 0) return;
	if(wspc->ellipse_count[wspc->current_image] <= 0) return ;
	int current_el = wspc->ellipse_count[wspc->current_image]-1;
	ellipse_ouv_t* el = &(wspc->ellipse_list[wspc->current_image][current_el].el);
	double dx = x - wspc->ellipse_list[wspc->current_image][current_el].center.c[0];
	double dy = y - wspc->ellipse_list[wspc->current_image][current_el].center.c[1];
	double l;
	
	switch(wspc->status_action ){
		case 1:
			el->u.c[0] = dx;
			el->u.c[1] = dy;
			l = r2_dir (&(el->u), &(el->u));
			el->a = l;
			break;
		case 2:
			//sprintf(debug,"Type 2 [%f,%f] center %f %f",dx,dy,wspc->ellipse_list[0].center.c[0], wspc->ellipse_list[0].center.c[1]);
			el->v.c[0] = dx;
			el->v.c[1] = dy;
			l = r2_dir (&(el->v), &(el->v));
			el->b = l;
			//redefining u to be orthogonal to v
			el->u.c[0] = el->v.c[1];
			el->u.c[1] = -el->v.c[0];
			
		break;
	}
}

void HandleMouseMotion_Circle(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	if(wspc->status_action  == 0) return;
	if(wspc->ellipse_count[wspc->current_image] <= 0) return ;
	int current_el = wspc->ellipse_count[wspc->current_image]-1;
	ellipse_ouv_t* el = &(wspc->ellipse_list[wspc->current_image][current_el].el);
	double dx = x - wspc->ellipse_list[wspc->current_image][current_el].center.c[0];
	double dy = y - wspc->ellipse_list[wspc->current_image][current_el].center.c[1];
	//puts("here?");
	switch(wspc->status_action ){
		case 2:
			el->u.c[0] = dx;
			el->u.c[1] = dy;
			el->v.c[0] = dy;
			el->v.c[1] = -dx;
			r2_dir (&(el->u), &(el->u));
			r2_dir (&(el->v), &(el->v));
			el->a = sqrt( (dx)*(dx) + (dy)*(dy));
			el->b = sqrt( (dx)*(dx) + (dy)*(dy));
			break;
	}
}


void HandleMouseMotion_Sphere(SDL_Event* event, workspace_t* wspc){
	int mx = event->button.x;
	int my = event->button.y;
	int x = mx - wspc->pan_x;
	int y = my - wspc->pan_y;
	if(wspc->status_action  == 0) return;
	if(wspc->ellipse_count[wspc->current_image] <= 0) return ;
	int current_el = wspc->ellipse_count[wspc->current_image]-1;
	ellipse_ouv_t* el = &(wspc->ellipse_list[wspc->current_image][current_el].el);

	switch(wspc->status_action ){
		case 2:
			{
				r2_t stretch;
				r2_t gab_center = wspc->ellipse_list[wspc->current_image][current_el].center;
				double dx = x - gab_center.c[0];
				double dy = y - gab_center.c[1];
				double gab_radius = hypot(dx,dy);
				float_image_t* current_image = wspc->float_image_list[wspc->current_image];
				r2_t cam_center =(r2_t) {{ current_image->sz[1]/2.0, current_image->sz[2]/2.0}};
				computeViewDirAndStretch(gab_center,gab_radius,cam_center,wspc->opts->cameraDist,NULL,&stretch);
				double st = r2_dir(&stretch,&el->u);
				if( st ==  0 ){
					el->u = (r2_t) {{ dx,dy}};
					r2_dir(&el->u,&el->u);
				}
				el->v = (r2_t){{el->u.c[1],-el->u.c[0]}};
				el->a = gab_radius + st;
				el->b = gab_radius;
			}
		break;
	}
}



void AdjustEllipseToFloatImage(workspace_t* wspc){
	int current_image = wspc->current_image;
	float_image_t* IMG = wspc->float_image_list[current_image];
	int index_el;
	int NC,NX,NY;
	NC = IMG->sz[0];  /* Image channels. */
    	NX = IMG->sz[1];  /* Image width (pixels). */
    	NY = IMG->sz[2];  /* Image height (pixels). */
	fprintf(stderr,"Fitting Ellipses\n");
	sprintf(debug," ");
	int gauge_sufix = 'A';
	for(index_el = 0; index_el < wspc->ellipse_count[wspc->current_image]; index_el++){
		ellipse_crs_t E;
		ellipse_t* el = &(wspc->ellipse_list[wspc->current_image][index_el]);
		
		char* gauge_prefix = NULL;
		char *gauge_prefix = jsprintf("%s-%c",wspc->opts->outPrefix,gauge_sufix);
		ellipse_ouv_to_crs(&(el->el), &(el->center), &E);

		fep_fit_ellipse(IMG, &E,wspc->opts->image_noise[wspc->current_image],el->type);
		
  		


		ellipse_crs_to_ouv(&E,&(el->el));
		//adjust center
		fprintf(stderr,"Parameters Adjusted\n");
		el->center = E.ctr;
		r2_t cam_center =  (r2_t){{0.5*NX,0.5*NY}};
		/*Compute view-dir and recompute strech if appropriate*/
		r3_t view_dir; 
		switch(el->type){
			case 1: 
				/*Ellipse*/
				{
					double d = r2_norm(&E.str);
					if(d == 0) { view_dir = (r3_t){{0,0,1}}; }
					else{
						//adjust signs of stretch so that it points towards optical center
						r2_t dgo;
						r2_sub(&cam_center,&E.ctr,&dgo);
						double dot = r2_dot(&E.str,&dgo);
						if(dot < 0 ) {
							r2_neg(&E.str,&E.str);
						}
						view_dir=(r3_t){{E.str.c[0],E.str.c[1],sqrt(2.0*E.rad*d)}};
						r3_dir(&view_dir,&view_dir);
						
					}
				}
				break;
			case 2:
				/*Circle*/
				assert( r2_norm(&E.str) == 0);
				view_dir = (r3_t){{0,0,1}};
				break;
			case 3:
				/*Sphere*/
				computeViewDirAndStretch(E.ctr,E.rad,cam_center, wspc->opts->cameraDist,&view_dir,&E.str);
				break;
		};
		if(!wspc->opts->yflip){
		  ellipse_crs_t Efl = E;
		  r3_t view_dir_fl = view_dir;
		  view_dir_fl.c[1] = -view_dir_fl.c[1];
		  Efl.str.c[1] = -Efl.str.c[1];
		  Efl.ctr.c[1]  = NX - Efl.ctr.c[1] - 1;
		  fep_write_ellipse_params(gauge_prefix, "full",Efl.ctr,Efl.rad,Efl.str,view_dir_fl);
		  
		}else{
		  fep_write_ellipse_params(gauge_prefix, "full",E.ctr,E.rad,E.str,view_dir);
		}

		/* Get a bounding box with some safety margin: */
    		double mrg = 2.0;
    		int xLo, xHi, yLo, yHi;
    		ellipse_crs_int_bbox(&E, mrg, &xLo, &xHi, &yLo, &yHi);

    		fprintf(stderr, "gauge cropping commands:\n"); 
    		fprintf(stderr, "  convert {IN} -crop '%dx%d+%d+%d' {OUT}\n", xHi-xLo, yHi-yLo, xLo, yLo);
    		fprintf(stderr, "  pnmcut %d %d %d %d < {IN} > {OUT}\n", xLo, yLo, xHi-xLo, yHi-yLo);
    		printCutParams(gauge_prefix, "cut",xHi-xLo, yHi-yLo, xLo, yLo); 
    		/* Write the clipped versions of the image and mask: */
		float_image_t *MSK = fep_compute_msk_image(&E, NX, NY);
       		fep_write_sub_image_as_pnm(gauge_prefix, "msk", MSK, 0, NX, 0, NY);
    		fep_write_sub_image_as_pnm(gauge_prefix, "clip-img", IMG, xLo, xHi, yLo, yHi);
		fep_write_sub_image_as_pnm(gauge_prefix, "clip-msk", MSK, xLo, xHi, yLo, yHi);
    		 //fep_write_sub_image_as_pnm(o->outPrefix, "clip-msk", MSK, xLo, xHi, yLo, yHi);
   
    		/* Compute the ellipse parameters {EC} relative to the clipped images: */
    		ellipse_crs_t EC = E; 
    		r2_t shift = (r2_t){{ xLo, yLo }};
    		r2_sub(&(E.ctr), &shift, &(EC.ctr));
    
    
    		/* Write the ellipse and camera parameters relative to the clipped image: */
		if(!wspc->opts->yflip){
		  ellipse_crs_t Efl = EC;
		  r3_t view_dir_fl = view_dir;
		  view_dir_fl.c[1] = -view_dir_fl.c[1];
		  Efl.str.c[1] = -Efl.str.c[1];
		  Efl.ctr.c[1]  = (yHi-yLo) - Efl.ctr.c[1] - 1;
		  fep_write_ellipse_params(gauge_prefix, "clip",Efl.ctr,Efl.rad,Efl.str,view_dir_fl);
		  
		}else{
		  fep_write_ellipse_params(gauge_prefix, "clip", EC.ctr,EC.rad,E.str,view_dir);
		}
		float_image_free(MSK);
		free(gauge_prefix);
		gauge_sufix++;
		
	}
	
	
}


void ClearScreen(SDL_Surface *screen)
{
	int i;
	/* Set the screen to black */
	if ( SDL_LockSurface(screen) == 0 ) {
		Uint32 black;
		Uint8 *pixels;
		black = SDL_MapRGB(screen->format, 0, 0, 0);
		pixels = (Uint8 *)screen->pixels;
		for ( i=0; i<screen->h; ++i ) {
			memset(pixels, black,
				screen->w*screen->format->BytesPerPixel);
			pixels += screen->pitch;
		}
		SDL_UnlockSurface(screen);
	}
}






void GetOptions(int argc, char** argv,SDL_options_t* sdl_opts, workspace_options_t* wspc_opts){
	
	argparser_t *pp = argparser_new(stderr, argc, argv);
	argparser_set_help(pp, PROG_HELP);
  	argparser_set_info(pp, PROG_INFO);
  	argparser_process_help_info_options(pp);

	/*set default values*/
	sdl_opts->w = 640;
	sdl_opts->h = 480;
	sdl_opts->desired_bpp = 0;
	sdl_opts->video_flags = 0;
	sdl_opts->screen_auto_size = 0;

	
	
	if (argparser_keyword_present(pp, "-width")) {
		sdl_opts->w = argparser_get_next_int(pp, 0, 2048);
	}
	if (argparser_keyword_present(pp, "-height")) {
		sdl_opts->h = argparser_get_next_int(pp, 0, 2048);
	}

	if (argparser_keyword_present(pp, "-bpp")) {
		sdl_opts->desired_bpp = argparser_get_next_int(pp, 8, 32);
	}
	if (argparser_keyword_present(pp, "-warp")) {
		sdl_opts->video_flags |= SDL_HWPALETTE ;
	}
	if (argparser_keyword_present(pp, "-hw")) {
		fprintf(stderr,"Trying alocate on HARDWARE\n");
		sdl_opts->video_flags |= SDL_HWSURFACE ;
	}
	if (argparser_keyword_present(pp, "-fullScreen")) {
		sdl_opts->video_flags |= SDL_FULLSCREEN ;
	}
	if (argparser_keyword_present(pp, "-autoSize")) {
		sdl_opts->screen_auto_size = 1;
	}
		

	sdl_opts->video_flags |= SDL_DOUBLEBUF;

	argparser_get_keyword(pp, "-cameraDist");
	wspc_opts->cameraDist = argparser_get_next_double(pp, 1.0, +INF);

	//now get workspace parameters
	argparser_get_keyword(pp, "-nImages");
  	wspc_opts->image_count = argparser_get_next_int(pp, 1, 1000);
	wspc_opts->image_noise = (double*)malloc(sizeof(double)*(wspc_opts->image_count));	
	int ind;
	float noise = 0.05;
	if (argparser_keyword_present(pp, "-noise")) {
		noise = argparser_get_next_double(pp, 0.0, 100.0);
	}
	for (ind = 0; ind < (wspc_opts->image_count); ind++){
		wspc_opts->image_noise[ind] = noise;
	}

	argparser_get_keyword(pp, "-prefix");
	wspc_opts->outPrefix = argparser_get_next(pp);
	
	argparser_get_keyword(pp, "-images");
  	fprintf(stderr, "Images:\n");
  	wspc_opts->float_image_filenames = (char**) malloc(sizeof(char*)*(wspc_opts->image_count));
  	
  	for (ind = 0; ind < (wspc_opts->image_count); ind++){
    		wspc_opts->float_image_filenames[ind] = argparser_get_next(pp);
    		fprintf(stderr, "    S[%02d] = %s\n", ind, wspc_opts->float_image_filenames[ind]);
  	}

	wspc_opts->gamma = 1.0;
	if (argparser_keyword_present(pp, "-gamma")) {
		wspc_opts->gamma = argparser_get_next_double(pp, 0.0, 100.0);
	}
	
	//Fliping goes here
	wspc_opts->yflip = argparser_keyword_present(pp, "-yflip");
	argparser_finish(pp);
	
}






/*ellipse fitting functions*/

void fep_fit_ellipse
  ( float_image_t *IMG, 
    ellipse_crs_t *E,
    double noise,
    int type
  )
  {
    fprintf(stderr, "finding the best-fit ellipse...\n");
    
    /* Get/check image dimensions: */
    int NX = IMG->sz[1];
    int NY = IMG->sz[2];

    /*Adjustement parameters*/
    double ctrAdj = 3.0; //center
    double radAdj = 3.0; //radius
    double strAdj = (type == 1 ? 3.0: 0.0); //stretch
    
    double minRadius = 5.0; //minumum radius for recursive call
    /* Print the user's guess: */
    fprintf(stderr,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(stderr, "initial guess for ellipse and adjustment amounts:\n");
    ellipse_crs_args_adjust_print(stderr, E, ctrAdj, radAdj, strAdj, "%8.2f");
    
    
    /* Check whether the initial guess is within the image bounds: */
    int xLo, xHi, yLo, yHi;
    double mrg = 0.0;
    ellipse_crs_int_bbox(E, mrg, &xLo, &xHi, &yLo, &yHi);
    fprintf(stderr, "initial bbox: [%d _ %d] �[%d _ %d]\n", xLo, xHi, yLo, yHi);
    demand((xLo >= 0) && (xHi <= NX), "guess extends outside valid X range");
    demand((yLo >= 0) && (yHi <= NY), "guess extends outside valid Y range");

    /* Do the fitting */
    ellipse_crs_t EFit = *E; /* Fitted ellipse. */
    double H = +INF;                /* Mismatch of {E}. */
    int maxIts = 15; /* Max iterations of the minimizer at coarsest scale. */
    H = pst_fit_ellipse_multiscale(
             IMG, noise, 
              &EFit, ctrAdj, radAdj, strAdj,
              minRadius, maxIts
            );
    /* Print the adjusted ellipse and camera: */
    fprintf(stderr, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(stderr, "adjusted ellipse:\n");
    ellipse_crs_args_adjust_print(stderr, &EFit, 0.0, 0.0, 0.0, "%8.2f");
    fprintf(stderr, "mismatch H = %12.6f\n", H);  
    fprintf(stderr, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    
    /* Return the fitted camera and ellipse: */
    (*E) = EFit;
}

void computeViewDirAndStretch(r2_t gab_center,double radius,r2_t cam_center, double cam_dist,r3_t* view_dirP, r2_t* stretchP){
	
	r2_t gcv; /*gabarito-camera vetor*/
	r2_sub(&gab_center,&cam_center,&gcv);
	/*Calcula view_dir e distancia h da camera ao gabarito*/
	r3_t view_dir = (r3_t){{ -gcv.c[0], -gcv.c[1], cam_dist}};
	double h = r3_dir(&view_dir,&view_dir);
	if(view_dirP != NULL) *view_dirP = view_dir;
	if(stretchP != NULL){
		/*calcula stretch*/
		double f = radius*((h/cam_dist) - 1.0 )/r2_norm(&gcv);
		r2_scale(f,&gcv,stretchP);
	}
}

void printGaugeParams(FILE* arq,r2_t gab_center,double gab_radius,  r2_t stretch, r3_t view_dir){
	
	fprintf(arq,"-gaugeCenter %lf %lf\n",gab_center.c[0],gab_center.c[1]);
	fprintf(arq,"-gaugeRadius %lf\n",gab_radius);
	fprintf(arq,"-gaugeStretch %lf %lf\n",stretch.c[0],stretch.c[1]);
	fprintf(arq,"-gaugeViewDir %lf %lf %lf\n",view_dir.c[0],view_dir.c[1],view_dir.c[2]);
}

void printCutParams(char *prefix, char* tag, int W, int H, int X, int Y){
	char *fname = NULL;
    	char *fname = jsprintf("%s-%s.cut", prefix, tag);
    	FILE *wr = open_write(fname, TRUE);
	fprintf(wr, "%dx%d+%d+%d", W, H, X, Y);
	fclose(wr);
    	free(fname);
}

void fep_write_ellipse_params(char *prefix, char *tag, r2_t gab_center,double gab_radius,  r2_t stretch, r3_t view_dir)
  { char *fname = NULL;
    char *fname = jsprintf("%s-%s.epar", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    printGaugeParams(wr,gab_center,gab_radius,stretch,view_dir);
    fclose(wr);
    free(fname);
  }

void fep_write_sub_image_as_pnm
  ( char *prefix, 
    char *tag, 
    float_image_t *A,
    int xLo,
    int xHi,
    int yLo, 
    int yHi
  )
  { /* Do we really need to crop? */
    int NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    bool_t crop = ((xLo != 0) || (xHi != NX) || (yLo != 0) || (yHi != NY));
    /* Get the (sub-)image {S} to write: */
    float_image_t *S = (crop ? float_image_crop(A, 0, NC, xLo, xHi, yLo, yHi, 0) : A);
    /* Write it: */
    char *name = jsprintf("%s-%s", prefix, tag);
    fep_write_image_as_pnm(name, S);
    free(name);
  }

void fep_write_image_as_pnm(char *name, float_image_t *A)
  { int NC = A->sz[0];
    char *ext = (NC == 1 ? "pgm" : "ppm");
    char *fname = jsprintf("%s.%s", name, ext);
    bool_t yup = FALSE; /* Use Y axis down, as in {pnmcut} and {convert}. */
    float_pnm_image_write(fname, A,FALSE, VIEW_GAMMA, VIEW_BIAS, yup, TRUE, FALSE);
    free(fname);
  }

float_image_t *fep_compute_msk_image(ellipse_crs_t *E, int NX, int NY)
  { 
    /* Make an image with the same size as the original: */
    float_image_t *MSK = float_image_new(1, NX, NY);
    float_image_fill(MSK, 0.0);
    
    /* Paint the precise filled ellipse into it: */
    float_image_paint_ellipse_crs(MSK, 0, E, 0.0, 1.0, NAN, 3 );
    
    return MSK;
  }


