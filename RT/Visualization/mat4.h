#include "vec3.h"

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else
#include <cmath>
#define HOST_DEVICE_INLINE  inline
#endif

struct mat4 {
    float _00, _01, _02, _03;
    float _10, _11, _12, _13;
    float _20, _21, _22, _23;
    float _30, _31, _32, _33;
};

HOST_DEVICE_INLINE
mat4 make_mat4_const(float c){
    mat4 m;
    m._00 = c; m._01 = c;           m._02 = c;           m._03 = c;
    m._10 = c; m._11 = c;  m._12 = c; m._13 = c;
    m._20 = c; m._21 = c;  m._22 = c;  m._23 =c;
    m._30 = c; m._31 = c;           m._32 = c;           m._33 = c;
    return m;
}


HOST_DEVICE_INLINE
mat4 make_mat4_identity(){
    mat4 m  =make_mat4_const (0);
    m._00 = 1;
        m._11 = 1;
            m._22 = 1;
                m._33 = 1;
    return m;
}



HOST_DEVICE_INLINE
float& elem(mat4 &a,int row, int col){
    return ((float *)&a) [col + row*4] ;
}

HOST_DEVICE_INLINE
float elem(const mat4 &a,int row, int col){
    return ((float *)&a) [col + row*4] ;

}

HOST_DEVICE_INLINE
mat4 operator*(const mat4 &lhs,const mat4 &rhs){
    mat4 m = make_mat4_const(0);

    for(int row=0; row < 4; row++)
        for(int col=0; col < 4; col++)
            for(int c=0; c < 4; c++)
                elem(m,row,col) += elem(lhs,row,c) * elem(rhs,c,col);
    return m;

}

HOST_DEVICE_INLINE
mat4 operator*(const float c,const mat4 &b){

    mat4 m=b;
    m._00 *= c; m._01 *= c;           m._02 *= c;           m._03 *= c;
    m._10 *= c; m._11 *= c;  m._12 *= c; m._13 *= c;
    m._20 *= c; m._21 *= c;  m._22 *= c;  m._23 *=c;
    m._30 *= c; m._31 *= c;           m._32 *= c;           m._33 *= c;
    return m;

}

HOST_DEVICE_INLINE
mat4 operator*(const mat4 &a,const float b){
    return b*a;
}

/*************************************************************************/
/* OPERATOR OVERLOAD OF * FOR A MULTIPLICATION BETWEEN A mat4 AND A vec3 */
/*************************************************************************/
HOST_DEVICE_INLINE vec3 operator*(const mat4 &a, const vec3 &b) {
    
	return make_vec3(a._00 * b.x + a._01 * b.y + a._02 * b.z + a._03,
                     a._10 * b.x + a._11 * b.y + a._12 * b.z + a._13,
                     a._20 * b.x + a._21 * b.y + a._22 * b.z + a._23);
}

HOST_DEVICE_INLINE
mat4 make_mat4_rotate_x(float theta){
    mat4 m;
    m._00 = 1; m._01 = 0;           m._02 = 0;           m._03 = 0;
    m._10 = 0; m._11 = cos(theta);  m._12 = -sin(theta); m._13 = 0;
    m._20 = 0; m._21 = sin(theta);  m._22 = cos(theta);  m._23 = 0;
    m._30 = 0; m._31 = 0;           m._32 = 0;           m._33 = 1;
    return m;
}

HOST_DEVICE_INLINE
mat4 make_mat4_rotate_y(float theta){
    mat4 m;
    m._00 = cos(theta);     m._01 = 0;      m._02 = sin(theta);  m._03 = 0;
    m._10 = 0;              m._11 = 1;      m._12 = 0;           m._13 = 0;
    m._20 = -sin(theta);    m._21 = 0;      m._22 = cos(theta);  m._23 = 0;
    m._30 = 0;              m._31 = 0;      m._32 = 0;           m._33 = 1;
    return m;
}


HOST_DEVICE_INLINE
mat4 make_mat4_rotate_z(float theta){
    mat4 m;
    m._00 = cos(theta);     m._01 = -sin(theta);  m._02 = 0;  m._03 = 0;
    m._10 = sin(theta);     m._11 = cos(theta);   m._12 = 0;  m._13 = 0;
    m._20 = 0;              m._21 = 0;            m._22 = 1;  m._23 = 0;
    m._30 = 0;              m._31 = 0;            m._32 = 0;  m._33 = 1;
    return m;

}

/*************************/
/* ROTATE AROUND AN AXIX */
/*************************/
HOST_DEVICE_INLINE mat4 make_mat4_rotate_axis(const vec3 &u, float theta) {
    
	mat4 m;
    
	float c = cosf(theta);
    float s = sinf(theta);
    
	float p = 1 - c;

    m._00 = c + u.x * u.x * p;          m._01 = u.x * u.y * p - u.z * s;	m._02 = u.x * u.z * p + u.y * s;	m._03 = 0;
    m._10 = u.y * u.x * p + u.z * s;    m._11 = c + u.y * u.y * p;			m._12 = u.y * u.z * p - u.x * s;	m._13 = 0;
    m._20 = u.z * u.x * p - u.y * s;    m._21 = u.z * u.y * p + u.x * s;	m._22 = c + u.z * u.z * p;			m._23 = 0;
    m._30 = 0;							m._31 = 0;							m._32 = 0;							m._33 = 1;
    
	return m;

}


HOST_DEVICE_INLINE
mat4 make_mat4_translate(const vec3 &v){
    mat4 m;
    m._00 = 1;     m._01 = 0;  m._02 = 0;  m._03 = v.x;
    m._10 = 0;     m._11 = 1;  m._12 = 0;  m._13 = v.y;
    m._20 = 0;     m._21 = 0;  m._22 = 1;  m._23 = v.z;
    m._30 = 0;     m._31 = 0;  m._32 = 0;  m._33 = 1;
    return m;


}


HOST_DEVICE_INLINE
mat4 make_mat4_rotate_xyz(float theta_x, float theta_y, float theta_z){
    return make_mat4_rotate_z(theta_z)  *  //...
           make_mat4_rotate_y(theta_y)  *  //...
           make_mat4_rotate_x(theta_x);
}



