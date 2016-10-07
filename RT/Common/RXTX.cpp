#include "RXTX.h"
#include "mat4.h"

void RX::update_unit_vectors() {

    x_v = normalize(make_vec3(cos(zenith) * cos(azimuth), cos(zenith) * sin(azimuth), -sin(zenith)));

    y_v = normalize(make_vec3(             -sin(azimuth),               cos(azimuth),            0));

    z_v = normalize(make_vec3(sin(zenith) * cos(azimuth),   sin(zenith)*sin(azimuth),  cos(zenith)));

    mat4 m = make_mat4_rotate_axis(z_v,roll);

    x_v = m * x_v;
    y_v = m * y_v;
    z_v = m * z_v;

}

void TX::update_unit_vectors(){

    x_v = normalize(make_vec3(cos(zenith) * cos(azimuth), cos(zenith) * sin(azimuth), -sin(zenith)));

    y_v = normalize(make_vec3(             -sin(azimuth),               cos(azimuth),            0));

    z_v = normalize(make_vec3(  sin(zenith)*cos(azimuth),   sin(zenith)*sin(azimuth),   cos(zenith)));

    mat4 m = make_mat4_rotate_axis(z_v,roll);

    x_v = m * x_v;
    y_v = m * y_v;
    z_v = m * z_v;

}
