#ifndef __VEC3_H__
#define __VEC3_H__

#include <iostream>
#include <algorithm>

struct vec3f{
	float x,y,z;	
};

struct vec3d{
	double x,y,z;	
};

struct vec3i{
	int a,b,c;	
};



std::ostream& operator<< (std::ostream &os, const vec3f & v);
std::ostream& operator<< (std::ostream &os, const vec3d & v);
std::ostream& operator<< (std::ostream &os, const vec3i & v);




//vec3f make_vec3f(float x,float y, float z);

//vec3d make_vec3d(double x,double y, double z);

//vec3i make_vec3i(int a, int b, int c);

inline
vec3f make_vec3f(float x,float y, float z){
    vec3f v = {x,y,z};
    return v;
}

inline
vec3d make_vec3d(double x,double y, double z){
    vec3d v = {x,y,z};
    return v;
}

inline
vec3i make_vec3i(int a, int b, int c){
    vec3i v = {a,b,c};
    return v;
}


vec3f operator+(const vec3f &a,const vec3f &b);

vec3f operator-(const vec3f &a,const vec3f &b);

vec3f operator+(const vec3f &a);

vec3f operator-(const vec3f &a);

vec3f operator*(const vec3f &a,const vec3f &b);

vec3f operator*(const float a,const vec3f &b);

vec3f operator*(const vec3f &a,const float b);

vec3f operator/(const vec3f &a,const vec3f &b);

float dot(const vec3f &a,const vec3f &b);

vec3f cross(const vec3f &a,const vec3f &b);

vec3f normalize(const vec3f &a);

float norm(const vec3f &a);

float dist(const vec3f &a,const vec3f &b);

float min(const vec3f &a);

float max(const vec3f &a);

/*vec3f cwisemin(const vec3f &a,const vec3f &b);

vec3f cwisemax(const vec3f &a,const vec3f &b);
*/

inline
vec3d cwisemin(const vec3d &a,const vec3d &b){
    return make_vec3d(std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z));
}


inline
vec3d cwisemax(const vec3d &a,const vec3d &b){
    return make_vec3d(std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z));
}




vec3d operator+(const vec3d &a,const vec3d &b);

vec3d operator-(const vec3d &a,const vec3d &b);

vec3d operator+(const vec3d &a);

vec3d operator-(const vec3d &a);

vec3d operator*(const vec3d &a,const vec3d &b);

vec3d operator*(const double a,const vec3d &b);

vec3d operator*(const vec3d &a,const double b);

vec3d operator/(const vec3d &a,const vec3d &b);

double dot(const vec3d &a,const vec3d &b);

vec3d cross(const vec3d &a,const vec3d &b);

vec3d normalize(const vec3d &a);

double norm(const vec3d &a);

double dist(const vec3d &a,const vec3d &b);

double min(const vec3d &a);

double max(const vec3d &a);

vec3d cwisemin(const vec3d &a,const vec3d &b);

vec3d cwisemax(const vec3d &a,const vec3d &b);


#endif
