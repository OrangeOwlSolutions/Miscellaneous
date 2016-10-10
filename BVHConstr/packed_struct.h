#ifndef __PACKED_STRUCT__H__
#define __PACKED_STRUCT__H__

/**************************/
/* VERTEX3F PACKED STRUCT */
/**************************/
#include "packed.h"
struct vertex3f {
    
	float x, y, z;
    
    void operator = (const vec3d &v) { x = v.x; y = v.y; z = v.z; }
} PACKED;
#include "endpacked.h"

/**************************/
/* VERTEX4F PACKED STRUCT */
/**************************/
#include "packed.h"
struct vertex4f {

	float x, y, z, w;

} PACKED;
#include "endpacked.h"

/************************/
/* FACE4I PACKED STRUCT */
/************************/
#include "packed.h"
struct face4i {

	int a, b, c, d;
    
    void operator = (const vec4i &f) { a = f.a; b = f.b; c = f.c; d = f.d; }
} PACKED;
#include "endpacked.h"

#endif
