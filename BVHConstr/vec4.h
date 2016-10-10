#ifndef __VEC4_H__
#define __VEC4_H__
#include <iostream>

struct vec4i{
	int a,b,c,d;
};


vec4i make_vec4i(int a, int b, int c, int d);

std::ostream& operator<< (std::ostream &os, const vec4i & v);


#endif
