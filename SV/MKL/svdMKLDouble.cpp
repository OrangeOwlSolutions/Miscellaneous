// --- This approach has been successfully checked only against the calculation of the singular values only. The check against the calculation of the
//     singular vectors has been unsuccessful for square matrices for reasons to be explored.
//     The approach has been tested only against square matrices.

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <mkl.h>
#include "mkl_lapacke.h"

#include "TimingCPU.h"
#include "InputOutput.h"

using namespace std;

//icpc -lpthread -lmkl_sequential -lmkl_core -lmkl_intel_lp64 -lrt -openmp svdMKLDouble.cpp TimingCPU.cpp InputOutput.cpp -O3 -o svdMKLDouble

/**********************************/
/* RANDOM MATRIX ENTRY GENERATION */
/**********************************/
double generateRandomMatrixEntry(int N) { return 2000. * ((double)rand() / (double)(RAND_MAX - 0.2) * (1. / (double)N)) + 100.*((double)rand() / (double)(RAND_MAX - 0.2) * (1. / (double)N)) + 24.; }

/********/
/* MAIN */
/********/
int main(int argc, char **argv) {

    	int Nrows, Ncols, K, numExecutions;
    
	char *jobu = "N", *jobvt = "N";		// --- "N" means not to calculate the corresponding singular vectors

	TimingCPU timer;
    
	if (argc != 5) { printf("usage: svdMKLDouble Nrows Ncols K numExecutions\n"); exit(0); }
    	{ std::stringstream ss(argv[1]); ss >> Nrows; }
    	{ std::stringstream ss(argv[2]); ss >> Ncols; }
    	{ std::stringstream ss(argv[3]); ss >> K; }    
    	{ std::stringstream ss(argv[4]); ss >> numExecutions; }

    	int lda = Nrows, ldu = Nrows, ldvt = Ncols;
    
    	double TotalTime = 0., AverageTime = 0.;

	double *A       = (double*)mkl_malloc(Nrows * Ncols * K * sizeof(double), 64); 	// --- Input matrices
	double *A_QUERY = (double*)mkl_malloc(Nrows * Ncols *     sizeof(double), 64); 	// --- Input matrix for MKL query
	double *U       = (double*)mkl_malloc(Nrows * Nrows * K * sizeof(double), 64); 	// --- Matrices of LEFT singular vectors
	double *S       = (double*)mkl_malloc(        Ncols * K * sizeof(double), 64); 	// --- Singular values
	double *V       = (double*)mkl_malloc(Ncols * Ncols * K * sizeof(double), 64); 	// --- Matrices of RIGHT singular vectors
	
	srand(time(NULL));
	//srand(0);
	int N = 5;

	/*************************************************/
	/* QUERY AND ALLOCATION OF THE OPTIMAL WORKSPACE */
	/*************************************************/
    	for (int h = 0; h < Nrows * Ncols; h++) A_QUERY[h] = generateRandomMatrixEntry(N);
        
	int info, lwork = -1;          
	double wkopt;
    	dgesvd(jobu, jobvt, &Nrows, &Ncols, A_QUERY, &lda, S, U, &ldu, V, &ldvt, &wkopt, &lwork, &info);
    	lwork = (int)wkopt;
    	double *work = (double *)malloc(lwork * K * sizeof(double));                
    
	/**************/
	/* DUMMY CALL */
	/**************/
	for (int i = 0; i < K; i++)
        	dgesvd(jobu, jobvt, &Nrows, &Ncols, A + (i * Nrows * Ncols), &lda, S + (i * Ncols), U + (i * Nrows * Nrows), &ldu, V + (i * Ncols * Ncols), &ldvt, work + i * lwork, &lwork, &info);
	
	/*********************/
	/* TIME MEASUREMENTS */
	/*********************/
	for (int k = 0; k < numExecutions; k++) {
	
	    	for (int h = 0; h < Nrows * Ncols * K; h++) { A[h] = generateRandomMatrixEntry(N); }
	
        	timer.StartCounter();
		#pragma omp parallel for
		for (int i = 0; i < K; i++) {
			dgesvd(jobu, jobvt, &Nrows, &Ncols, A + (i * Nrows * Ncols), &lda, S + (i * Ncols), U + (i * Nrows * Nrows), &ldu, V + (i * Ncols * Ncols), &ldvt, work + i * lwork, &lwork, &info);
			// --- Convergence check			
			if (info > 0) { printf("Convergence failure.\n"); exit(1); }
		}
		#pragma omp barrier   
        	TotalTime += timer.GetCounter(); // --- TotalTime in ms
		    
	}
	
	AverageTime = TotalTime / numExecutions;
	std::cout << std::scientific << Nrows << " " << Ncols << " " << K << " " << numExecutions << " " << AverageTime << " " << std::endl;
    
	/*****************************/
	/* SAVINGS FOR RESULT CHECKS */
	/*****************************/
	// --- Saving the last matrix
	saveCPUrealtxt<double>(A + ((K - 1) * Nrows * Ncols), "Source_Matrix.txt", Nrows * Ncols);	
	// --- Saving the SVD
	saveCPUrealtxt<double>(U + ((K - 1) * Nrows * Nrows), "U.txt", Nrows * Nrows);	
	saveCPUrealtxt<double>(V + ((K - 1) * Ncols * Ncols), "V.txt", Ncols * Ncols);	
	saveCPUrealtxt<double>(S + ((K - 1) * Ncols), 	      "S.txt", Ncols        );	

	return 0;
}
