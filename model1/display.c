#include <SDL/SDL.h>
#include "cairosdl.h"
#include <stdio.h>

cairo_t *cr; // Global cairo drawing context
int width = 600; //Width of drawing window
int height = 600; //Height of drawing window

void init_video(void) {
	if(SDL_Init(SDL_INIT_VIDEO) < 0) {
		fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	SDL_Surface *screen = SDL_SetVideoMode(width,height,32,SDL_HWSURFACE|SDL_ASYNCBLIT);
	SDL_Surface *canvas = SDL_CreateRGBSurface(SDL_SWSURFACE,width,height,32,CAIROSDL_RMASK,CAIROSDL_GMASK,CAIROSDL_BMASK,0);
	cr = cairosdl_create(canvas);

}

