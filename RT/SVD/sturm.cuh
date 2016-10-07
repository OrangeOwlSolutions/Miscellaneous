#ifndef __STURM_H__
#define __STURM_H__

// -- RETURN THE FIRST INTEGER MULTIPLE -------------------------------------

inline int iDivUp(int a, int b); 

// ---------------------------------------------------------------------------

#include <float.h>
#include "utils.h"



// -----------------------------------------------------------------------------------------------------------------------

// -- Build tridiagonal matrix from bidiagonal --------------------------------------------------------------------------

template<class real_type, int N>
__host__ __device__ void build_tri(const real_type *d, const real_type *e, real_type *alpha, real_type *beta, uint batch_size){
    
    alpha[0] = d[0]*d[0];
    
    // -- Build alpha
    for(unsigned int i = 1; i < N; i++){
        alpha[i*batch_size] = (e[(i-1)*batch_size]*e[(i-1)*batch_size]) + d[i*batch_size]*d[i*batch_size];
    }
    
    // -- Build beta
    for(unsigned int i = 0; i < N - 1; i++)
        beta[i*batch_size] = d[i*batch_size]*e[i*batch_size];        

}

// -----------------------------------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------------------------------

// -- Evaluate Sturm Sequences
// -- STURM evaluates the Sturm succession in x ---------------------------------------------------------

template<class real_type>
__host__ __device__ void sturm(const real_type *d, const real_type *b,
                               real_type x, real_type *p, unsigned int length){

    unsigned int n = length;
    p[0] = 1;
    p[1] = d[0] - x;
    for(unsigned int i = 1; i < n; i++)
        p[i+1] = (d[i] - x)*p[i] - (b[i-1]*b[i-1])*p[i-1];
}
                               
// -------------------------------------------------------------------------------------------------------

// -- CHCKSIGN calcola il cambio di segni nella successione di Sturm -------------------------------------

template<class real_type>
__host__ __device__ void chcksign(const real_type *d, const real_type *b,
                                  real_type *p, real_type x, unsigned int &nch, unsigned int length){
    
    unsigned int n = length;
    unsigned int s = 0; // -- Count the number of zeros
    sturm(d, b, x, p, n);
    nch = 0;
    for(unsigned int i = 1; i <= n; i++){
        if(p[i]*p[i-1] <=0)
            nch += 1;
        else if(p[i] == 0)
            s += 1;    
    }
    
    nch = nch - s;
}
                                  
// --------------------------------------------------------------------------------------------------------        

// -- CHCKSIGN calcola il cambio di segni nella successione di Sturm -------------------------------------

template<class real_type>
__host__ __device__ void chcksign_stable(const real_type *d, const real_type *b,
                                         real_type x, uint &nch, uint length, uint batch_size){
    
    unsigned int n = length;
    nch = 0;
    
    real_type q = d[0] - x;
    if(q < 0.f){
        nch = nch + 1;
    }
    
    for(unsigned int i = 2; i <= n; i++){
    
        if( q != 0.f ){
            q = ( d[(i-1)*batch_size]  - ((b[(i-2)*batch_size]*b[(i-2)*batch_size])/q) ) - x; // -- This order of operation preserve monotonicity 
            if(q < 0.f){
                nch = nch + 1;
            }            
        } 
        
        else{
            q = ( d[(i-1)*batch_size] - (fabs(b[(i-2)*batch_size])/DBL_MIN) ) - x;
            if(q < 0.f){
                nch = nch + 1;
            }            
        }
    }
    
}
                                  
// -------------------------------------------------------------------------------------------------------- 


                                  
// -------------------------------------------------------------------------------------------------------- 

template<class real_type>
__host__ __device__ void chcksign_SignBit(const real_type *d, const real_type *b,
                                          real_type x, uint &nch, uint length, uint batch_size){
    
    unsigned int n = length;
    nch = 0;
    
    real_type q = d[0] - x;
    nch = nch + signbit(q);
    
    for(unsigned int i = 2; i <= n; i++){
        q = ( d[(i-1)*batch_size]  - ((b[(i-2)*batch_size]*b[(i-2)*batch_size])/q) ) - x; // -- This order of operation preserve monotonicity 
        nch = nch + signbit(q);
    }
    
}
                                  
// --------------------------------------------------------------------------------------------------------        

// -- Build WORKLIST ---------------------------------------------------------------------------------------------------

template<class real_type>
__host__ __device__ void build_worklist(uint subintervals, uint sub_index, const real_type *d, const real_type *b,
                                        uint length, uint batch_size,const real_type alpha,const  real_type beta,
                                        real_type *pivmin, uint *worklist){
    
    uint p = subintervals;
    real_type alpha_temp = 0;
    real_type beta_temp = 0;
    real_type width = (beta - alpha)/p;
    alpha_temp = alpha + sub_index * width;
    if(sub_index == p - 1){
        beta_temp = beta;
    }
    else{
        beta_temp = alpha + (sub_index+1) * width;
    }
    uint nch_temp1 = 0, nch_temp2 = 0;
    //chcksign_stable(d, b, alpha_temp, nch_temp1, length, batch_size);
    //chcksign_stable(d, b, beta_temp, nch_temp2, length, batch_size);
    chcksign_pivmin(d, b, pivmin[0], alpha_temp, nch_temp1, length, batch_size);
    chcksign_pivmin(d, b, pivmin[0], beta_temp, nch_temp2, length, batch_size);
    *worklist = nch_temp2 - nch_temp1;
}

// ---------------------------------------------------------------------------------------------------------------------                 
                                  
// -- CHCKSIGN calcola il cambio di segni nella successione di Sturm --------------------------------------------------------------

template<class real_type, int N>
__host__ __device__ void chcksign_optimized(const real_type *d, const real_type *b, real_type x, uint &nch, uint batch_size){
    
    real_type p_im1 = 1;
    real_type p_i = d[0] - x;
    real_type p_ip1 = 0;
    
    unsigned int s = 0; // -- Count the number of zeros

    nch = 0;
    
    #pragma unroll     
    for(unsigned int i = 2; i <= N; i++){
        p_ip1 = (d[(i-1)*batch_size] - x)*p_i - (b[(i-2)*batch_size]*b[(i-2)*batch_size])*p_im1;
        if(p_i*p_im1 <=0)
            nch += 1;
        else if(p_i == 0)
            s += 1;    
            
        p_im1 = p_i;
        p_i = p_ip1;    
    }
    
    if(p_i*p_im1 <=0)
        nch += 1;
    else if(p_i == 0)
        s += 1;      
    
    nch = nch - s;
}

// ---------------------------------------------------------------------------------------------------------------------------------

// -- CHCKSIGN calcola il cambio di segni nella successione di Sturm --------------------------------------------------------------

template<class real_type>
__host__ __device__ void chcksign_optimized_def(const real_type *d, const real_type *b, real_type x, uint &nch, uint length, uint batch_size){
    
    real_type p_im1 = 1;
    real_type p_i = d[0] - x;
    real_type p_ip1 = 0;
    
    unsigned int n = length;
    unsigned int s = 0; // -- Count the number of zeros

    nch = 0;
    
         
    for(unsigned int i = 2; i <= n; i++){
        p_ip1 = (d[(i-1)*batch_size] - x)*p_i - (b[(i-2)*batch_size]*b[(i-2)*batch_size])*p_im1;
        if(p_i*p_im1 <=0)
            nch += 1;
        else if(p_i == 0)
            s += 1;    
            
        p_im1 = p_i;
        p_i = p_ip1;    
    }
    
    if(p_i*p_im1 <=0)
        nch += 1;
    else if(p_i == 0)
        s += 1;      
    
    nch = nch - s;
}

// ---------------------------------------------------------------------------------------------------------------------------------
                                  

// ----------------------------------------------------------------------------------------------------------------------------------

// -- GIVSTURM Metodo di Givens basato sulla successione di Sturm -------------------------------------------------------------------

template<class real_type>
__host__ __device__ void givsturm(const real_type *d, const real_type *b, real_type *p, real_type *tol, real_type *eig, uint length){
    
    const unsigned int max_iter = 100; // -- Set the maximum number of iterations
    
    unsigned int n = length;
    unsigned int niter = 0; // -- Number of iteration of bisection
    unsigned int nch = 0; // -- Number of sign changes in the Sturm sequence
    double alpha = 0, beta = 0, alpha_temp = 0, beta_temp = 0, dist_temp = 0, s_temp = 0;
    double c = 0;
    bound(d, b, alpha, beta, length); // -- Evaluate minimum interval
    
    double dist = fabs(beta - alpha);
    double s = fabs(beta) + fabs(alpha);
    
    for(unsigned int index = 1; index <= n; index++){
    
        if(index == 1){
            alpha_temp = alpha;
            beta_temp = beta;
            dist_temp = dist;
            s_temp = s;
        }
        
        else{
            alpha_temp = alpha;
            beta_temp = c;
            dist_temp = fabs(beta_temp - alpha_temp);
            s_temp = fabs(beta_temp) + fabs(alpha_temp);        
        }

        niter = 0; // -- Reset niter at each index value
        while( (dist_temp > tol[0]*s_temp) && niter < max_iter ){
            niter += 1;
            c = (alpha_temp + beta_temp)/2;
            chcksign(d, b, p, c, nch, n);
            if(nch > (n - index)){
                beta_temp = c;
            }    
            else{
                alpha_temp = c;
            }    
            dist_temp = fabs(beta_temp - alpha_temp);
            s_temp = fabs(beta_temp) + fabs(alpha_temp);            
        }
        
        //eig[index - 1] = sqrt(ck[niter]);  
        eig[index - 1] = sqrt(c);  
    }
}

// -------------------------------------------------------------------------------------------------------------------------------

// -- GIVSTURM Metodo di Givens basato sulla successione di Sturm (Esecuzione parallela per blocchi) -----------------------------

template<class real_type>
__global__ void givsturm_parallel(real_type *d, real_type *b, real_type *dd, real_type *bb,
                                  real_type *eig, const uint cols, const uint batch_size){

    int tid_x = threadIdx.x + blockDim.x * blockIdx.x;
    int tid_y = threadIdx.y;
    
    //real_type tol = 0.00000001;
    //real_type tol = 0.00000000000001;
    real_type tol = 0.000000001;
    
    __shared__ real_type alpha_shared[NUM_THREADS_PER_BLOCK_STURM];
    __shared__ real_type beta_shared[NUM_THREADS_PER_BLOCK_STURM];
    
    const unsigned int n = cols;
    //const unsigned int max_iter = 200; // -- Set the maximum number of iterations
    
    if(tid_x < batch_size && tid_y == 0){
       build_tri(d + tid_x, b + tid_x, dd + tid_x, bb + tid_x, n, batch_size);
       bound(dd + tid_x, bb + tid_x, alpha_shared[threadIdx.x], beta_shared[threadIdx.x], cols, batch_size);

    }
    
    __syncthreads();
    
    if(tid_x < batch_size && tid_y < cols){
        

        uint eig_index = tid_y%n;
        //unsigned int eig_index = tid_y&(n - 1);
           
        real_type c = 0;
        unsigned int nch = 0;
                     
        real_type alpha , beta;
        
        alpha = alpha_shared[threadIdx.x];
        beta = beta_shared[threadIdx.x];
        
        real_type dist = fabs(beta - alpha); // -- Ampiezza iniziale dell'intervallo di ricerca diverso per ogni thread lungo X
        real_type s = fabs(beta) + fabs(alpha);
                
        unsigned int niter = 0; // -- Number of iteration of bisection

        //while( (dist > tol*s) && niter < max_iter ){
        while( dist > (tol + 2*DBL_MIN*s) ){
            niter += 1;
            c = (alpha + beta)/2;
            chcksign_stable<real_type>(dd + tid_x, bb + tid_x, c, nch, n, batch_size);  
            //chcksign_SignBit<real_type>(dd + tid_x, bb + tid_x, c, nch, n, batch_size);          
//            switch(n){
//                case 2:
//                    chcksign_optimized<real_type, 2>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 3:
//                    chcksign_optimized<real_type, 3>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 4:
//                    chcksign_optimized<real_type, 4>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 5:
//                    chcksign_optimized<real_type, 5>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 6:
//                    chcksign_optimized<real_type, 6>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 7:
//                    chcksign_optimized<real_type, 7>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 8:
//                    chcksign_optimized<real_type, 8>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 9:
//                    chcksign_optimized<real_type, 9>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;     
//                case 10:
//                    chcksign_optimized<real_type, 10>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 11:
//                    chcksign_optimized<real_type, 11>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                case 12:
//                    chcksign_optimized<real_type, 12>(dd + tid_x, bb + tid_x, c, nch, batch_size);
//                    break;
//                default:
//                    chcksign_optimized_def<real_type>(dd + tid_x, bb + tid_x, c, nch, n, batch_size);
//                    break;                                                                                                                                                                               
//            }
            if( nch > (n - ( (eig_index) + 1) ) ){
                beta = c;
            }    
            else{
                alpha = c;
            }    
            dist = fabs(beta - alpha);
            s = fabs(beta) + fabs(alpha);            
        }  
        
        eig[eig_index + tid_x*n] = sqrt(fabs(c));
        //eig[eig_index*batch_size + tid_x] = sqrt(c);    
       
    }
    
}
 
// -------------------------------------------------------------------------------------------------------------------------------

// -- GIVSTURM Metodo di Givens basato sulla successione di Sturm (Esecuzione parallela per blocchi) -----------------------------

template<class real_type>
__global__ void givsturm_parallel_worklist(real_type *d, real_type *b, real_type *dd, real_type *bb,
                                           real_type *eig, const uint cols, const uint batch_size){

    int tid_x = threadIdx.x + blockDim.x * blockIdx.x;
    int tid_y = threadIdx.y;
    const uint p = 3;
    
    //uint worklist[20];
    
    //real_type tol = 0.00000001;
    //real_type tol = 0.00000000000001;
    real_type tol = 0.000000001;
    
    __shared__ real_type alpha_shared[NUM_THREADS_PER_BLOCK_STURM];
    __shared__ real_type beta_shared[NUM_THREADS_PER_BLOCK_STURM];
    __shared__ real_type pivmin_shared[NUM_THREADS_PER_BLOCK_STURM];
    __shared__ uint worklist[p];
    
    const unsigned int n = cols;
    //const unsigned int max_iter = 200; // -- Set the maximum number of iterations
    
    if(tid_x < batch_size && tid_y == 0){
       build_tri(d + tid_x, b + tid_x, dd + tid_x, bb + tid_x, n, batch_size);
       build_pivmin(bb + tid_x, &pivmin_shared[threadIdx.x], cols, batch_size);
       bound(dd + tid_x, bb + tid_x, alpha_shared[threadIdx.x], beta_shared[threadIdx.x], cols, batch_size);

    }
    
    __syncthreads();
    
    if(tid_x < batch_size && tid_y < p){

        build_worklist(p, tid_y, dd + tid_x, bb + tid_x, cols, batch_size,
                       alpha_shared[threadIdx.x], beta_shared[threadIdx.x], &pivmin_shared[threadIdx.x], &worklist[threadIdx.y]);   
                       
    }                                         
    
    if(tid_y == 0 && tid_x == 0){
        for(uint i = 0; i < p; i++)
            printf("Work[%i] = %i \n",i, worklist[i]);
    }
    
    if(tid_y == 0 && tid_x == 3){
        printf("\n\n");
        for(uint i = 0; i < p; i++)
            printf("Work[%i] = %i \n",i, worklist[i]);
    }
        
//    if(tid_x < batch_size && tid_y < cols){
//        

//        uint eig_index = tid_y%n;
//        //unsigned int eig_index = tid_y&(n - 1);
//           
//        real_type c = 0;
//        unsigned int nch = 0;
//                     
//        real_type alpha , beta;
//        
//        alpha = alpha_shared[threadIdx.x];
//        beta = beta_shared[threadIdx.x];
//        
//        real_type dist = fabs(beta - alpha); // -- Ampiezza iniziale dell'intervallo di ricerca diverso per ogni thread lungo X
//        real_type s = fabs(beta) + fabs(alpha);
//                
//        unsigned int niter = 0; // -- Number of iteration of bisection

//        //while( (dist > tol*s) && niter < max_iter ){
//        while( dist > (tol + 2*DBL_MIN*s) ){
//            niter += 1;
//            c = (alpha + beta)/2;
//            chcksign_stable<real_type>(dd + tid_x, bb + tid_x, c, nch, n, batch_size);  
//            //chcksign_SignBit<real_type>(dd + tid_x, bb + tid_x, c, nch, n, batch_size);          
//            if( nch > (n - ( (eig_index) + 1) ) ){
//                beta = c;
//            }    
//            else{
//                alpha = c;
//            }    
//            dist = fabs(beta - alpha);
//            s = fabs(beta) + fabs(alpha);            
//        }  
//        
//        eig[eig_index + tid_x*n] = sqrt(fabs(c));
//        //eig[eig_index*batch_size + tid_x] = sqrt(c);    
//       
//    }
    
}
 
// -------------------------------------------------------------------------------------------------------------------------------

 
// -------------------------------------------------------------------------------------------------------------------------------


//// -- GIVSTURM Metodo di Givens basato sulla successione di Sturm (Esecuzione parallela per blocchi) -----------------------------

//template<class real_type>
//__global__ void givsturm_parallel(real_type *d, real_type *b, real_type *dd, real_type *bb,
//                                  real_type *eig, const unsigned int cols, const unsigned int batch_size){

//    int tid_x = threadIdx.x + blockDim.x * blockIdx.x;
//    int tid_y = threadIdx.y;
//    
//    double tol = 0.00000001;
//    
//    __shared__ real_type alpha_shared[NUM_THREADS_PER_BLOCK_STURM];
//    __shared__ real_type beta_shared[NUM_THREADS_PER_BLOCK_STURM];
//    
//    const unsigned int n = cols;
//    const unsigned int max_iter = 100; // -- Set the maximum number of iterations
//    
//    if(tid_x < batch_size && tid_y == 0){

//       build_tri(d + tid_x, b + tid_x, dd + tid_x, bb + tid_x, n, batch_size);
//       bound(dd + tid_x, bb + tid_x, alpha_shared[threadIdx.x], beta_shared[threadIdx.x], cols, batch_size);

//    }
//    
//    __syncthreads();
//    
//    if(tid_x < batch_size && tid_y < cols){
//        

//        unsigned int eig_index = tid_y%n;
//        //unsigned int eig_index = tid_y&(n - 1);
//           
//        real_type c = 0;
//        
//                            
//        real_type alpha , beta;
//        
//        alpha = alpha_shared[threadIdx.x];
//        beta = beta_shared[threadIdx.x];
//        
//        real_type dist = fabs(beta - alpha); // -- Ampiezza iniziale dell'intervallo di ricerca diverso per ogni thread lungo X
//        real_type s = fabs(beta) + fabs(alpha);
//                
//        uint niter = 0; // -- Number of iteration of bisection
//        uint nch = 0;
//       

//        
//        while( (dist > tol*s) && niter < max_iter ){
//            niter += 1;
//            c = (alpha + beta)/2;
//            
//            {
//                const real_type* DD = dd + tid_x;
//                const real_type* BB =  bb + tid_x; 
//                real_type p_im1 = 1;
//                real_type p_i = DD[0] - c;
//                real_type p_ip1 = 0;
//                uint s_zero = 0;
//                nch = 0;
//                uint i = 2;

//                
//                #define evalute_poly_x \
//                          p_ip1 = (DD[(i-1)*batch_size] - c)*p_i - (BB[(i-2)*batch_size]*BB[(i-2)*batch_size])*p_im1; \
//                          if(p_i*p_im1 <=0) \
//                             nch += 1; \
//                          else if(p_i == 0) \
//                             s_zero += 1; \
//                          p_im1 = p_i; \
//                          p_i = p_ip1; \
//                          i++;           
//                                 
//                switch(n){
//                    case 12:
//                    {
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        }
//                        break;
//                    case 11:
//                    {
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        }
//                        break;                        
//                    case 10:
//                    {
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        }
//                        break;
//                    case 9:
//                    {
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        }
//                        break;
//                    case 8:
//                    {
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                        evalute_poly_x
//                    }
//                        break;
//                    case 7:{
//                        evalute_poly_x
//                        }
//                        break;
//                    case 6:
//                    {
//                        evalute_poly_x
//                        }
//                        break;
//                    case 5:
//                    {
//                        evalute_poly_x
//                        }
//                        break;
//                    case 4:
//                    {
//                        evalute_poly_x
//                        
//                        }
//                        break;
//                    case 3:
//                    {
//                        evalute_poly_x
//                        }
//                        break;
//                        
//                    case 2:
//                    {
//                        evalute_poly_x
//                        }
//                        break;
//                    default:
//                        chcksign_optimized<real_type>(DD, BB, c, nch, n, batch_size);                                                                                                                                                                                                                                                                            
//                        break;
//                }

//                if(p_i*p_im1 <=0)
//                    nch += 1;
//                else if(p_i == 0)
//                    s_zero += 1;      

//                nch = nch - s_zero;                
//                   
//            }
//            
//            #undef evaluate_poly_x 
//          

//            if( nch > (n - ( (eig_index) + 1) ) ){
//                beta = c;
//            }    
//            else{
//                alpha = c;
//            }    
//            dist = fabs(beta - alpha);
//            s = fabs(beta) + fabs(alpha);            
//        }  
//        
//        eig[eig_index + tid_x*n] = sqrt(c);    
//       
//    }
//    
//}
// 
//// ------------------------------------------------------------------------------------------------------------------------------- 





// -------------------------------------------------------------------------------------------------------------------------------  
   

#endif
