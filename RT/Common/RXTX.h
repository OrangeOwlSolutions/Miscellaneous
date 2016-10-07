#ifndef __RXTX__H__
#define __RXTX__H__

#include "vec3.h"
#include "pattern.h"

/**********************/
/* TRANSMITTER STRUCT */
/**********************/
struct TX {

    vec3 	x_v;
    vec3 	y_v;
	vec3 	z_v;

    vec3 	pos;

    double 	azimuth;
    double 	zenith;
    double 	roll;

    pattern *patt;

	void update_unit_vectors();

};

/*******************/
/* RECEIVER STRUCT */
/*******************/
struct RX {

    vec3 	pos;

    double 	radius;

    vec3	x_v;
    vec3 	y_v;
    vec3 	z_v;

    double 	azimuth;
    double 	zenith;
    double 	roll;

    pattern *patt;

	void update_unit_vectors();

};

#endif
