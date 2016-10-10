#include "vec3.h"

#include <cmath>

#define maximum(a, b) (((a) > (b)) ? (a) : (b))
#define minimum(a, b) (((a) < (b)) ? (a) : (b))

std::ostream & operator << (std::ostream &os, const vec3f & v){
	return os << v.x << " " << v.y << " " << v.z ; 
}

std::ostream & operator << (std::ostream &os, const vec3d & v){
	return os << v.x << " " << v.y << " " << v.z ; 
}

std::ostream & operator << (std::ostream &os, const vec3i & v){
	return os << v.a << " " << v.b << " " << v.c ; 
}

//vec3f make_vec3f(float x,float y, float z){
    //vec3f v = {x,y,z};
    //return v;
//}

//vec3d make_vec3d(double x,double y, double z){
    //vec3d v = {x,y,z};
    //return v;
//}


//vec3i make_vec3i(int a, int b, int c){
    //vec3i v = {a,b,c};
    //return v;
//}


vec3f operator+(const vec3f &a,const vec3f &b){
    return make_vec3f(a.x+b.x,a.y+b.y,a.z+b.z);
}


vec3f operator-(const vec3f &a,const vec3f &b){
    return make_vec3f(a.x-b.x,a.y-b.y,a.z-b.z);
}


vec3f operator+(const vec3f &a){
    return a;
}


vec3f operator-(const vec3f &a){
    return make_vec3f(-a.x,-a.y,-a.z);
}



vec3f operator*(const vec3f &a,const vec3f &b){
    return make_vec3f(a.x*b.x,a.y*b.y,a.z*b.z);
}


vec3f operator*(const float a,const vec3f &b){
    return make_vec3f(a*b.x,a*b.y,a*b.z);
}


vec3f operator*(const vec3f &a,const float b){
    return make_vec3f(a.x*b,a.y*b,a.z*b);
}


vec3f operator/(const vec3f &a,const vec3f &b){
    return make_vec3f(a.x/b.x,a.y/b.y,a.z/b.z);
}


float dot(const vec3f &a,const vec3f &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}


vec3f cross(const vec3f &a,const vec3f &b){
    return make_vec3f(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}


vec3f normalize(const vec3f &a){
    return a*(1.0f/sqrt(dot(a,a)));
}


float norm(const vec3f &a){
    return sqrt(dot(a,a));
}


float dist(const vec3f &a,const vec3f &b){
    return norm(b-a);
}


float min(const vec3f &a){
#ifdef __linux__
    return fminf(fminf(a.x,a.y),a.z);
#elif _WIN32 || _WIN64
    return minimum(minimum(a.x,a.y),a.z);
#endif
}


float max(const vec3f &a){
#ifdef __linux__
	return fmaxf(fmaxf(a.x,a.y),a.z);
#elif _WIN32 || _WIN64
	return maximum(maximum(a.x,a.y),a.z);
#endif
}


vec3f cwisemin(const vec3f &a,const vec3f &b){
#ifdef __linux__
    return make_vec3f(fminf(a.x,b.x),fminf(a.y,b.y),fminf(a.z,b.z));
#elif _WIN32 || _WIN64
    return make_vec3f(minimum(a.x,b.x),minimum(a.y,b.y),minimum(a.z,b.z));
#endif
}


vec3f cwisemax(const vec3f &a,const vec3f &b){
#ifdef __linux__
    return make_vec3f(fmaxf(a.x,b.x),fmaxf(a.y,b.y),fmaxf(a.z,b.z));
#elif _WIN32 || _WIN64
    return make_vec3f(maximum(a.x,b.x),maximum(a.y,b.y),maximum(a.z,b.z));
#endif
}



//vec3d

vec3d operator+(const vec3d &a,const vec3d &b){
    return make_vec3d(a.x+b.x,a.y+b.y,a.z+b.z);
}


vec3d operator-(const vec3d &a,const vec3d &b){
    return make_vec3d(a.x-b.x,a.y-b.y,a.z-b.z);
}


vec3d operator+(const vec3d &a){
    return a;
}


vec3d operator-(const vec3d &a){
    return make_vec3d(-a.x,-a.y,-a.z);
}



vec3d operator*(const vec3d &a,const vec3d &b){
    return make_vec3d(a.x*b.x,a.y*b.y,a.z*b.z);
}


vec3d operator*(const double a,const vec3d &b){
    return make_vec3d(a*b.x,a*b.y,a*b.z);
}


vec3d operator*(const vec3d &a,const double b){
    return make_vec3d(a.x*b,a.y*b,a.z*b);
}


vec3d operator/(const vec3d &a,const vec3d &b){
    return make_vec3d(a.x/b.x,a.y/b.y,a.z/b.z);
}


double dot(const vec3d &a,const vec3d &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}


vec3d cross(const vec3d &a,const vec3d &b){
    return make_vec3d(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}


vec3d normalize(const vec3d &a){
    return a*(1.0/sqrt(dot(a,a)));
}


double norm(const vec3d &a){
    return sqrt(dot(a,a));
}


double dist(const vec3d &a,const vec3d &b){
    return norm(b-a);
}


double min(const vec3d &a){
    return std::min(std::min(a.x,a.y),a.z);
}


double max(const vec3d &a){
    return std::max(std::max(a.x,a.y),a.z);
}


//vec3d cwisemin(const vec3d &a,const vec3d &b){
    //return make_vec3d(std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z));
//}


//vec3d cwisemax(const vec3d &a,const vec3d &b){
    //return make_vec3d(std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z));
//}

