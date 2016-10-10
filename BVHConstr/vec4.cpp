#include "vec4.h"

vec4i make_vec4i(int a, int b, int c, int d){
    vec4i v = {a,b,c,d};
    return v;
}


std::ostream & operator << (std::ostream &os, const vec4i & v){
	return os << v.a << " " << v.b << " " << v.c << " " << v.d; 
}
