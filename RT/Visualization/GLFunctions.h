#ifndef __GLFUNCTIONS_H__
#define __GLFUNCTIONS_H__

#include <algorithm>
#include <vector>

#include "load_nbin.h"
#include "vec3.h"
#include "RXTX.h"

struct inter_point {
    int type ;
    vec3 pos;
};

void draw_reference_frame();

void draw_plane(int Nx,int Nz);

void draw_segment(const vec3 A, const vec3 B);

void keyboardUp(unsigned char, int, int);
void keyboardDown(unsigned char, int, int);

void loop(int);

void Visualize(int *, char **);

//void * input_loop(void *);
//void input_loop();

void display(void);

void draw_model(void);
void draw_rays();
void draw_rx(const RX &);
void draw_tx(const TX &);

void reshape(int, int);

void draw_segment(const vec3, const vec3);

void draw_plane(int, int);

extern nvmesh mesh;
extern bool mesh_loaded;
extern std::vector<TX> txs;
extern std::vector<RX> rxs;

#endif
