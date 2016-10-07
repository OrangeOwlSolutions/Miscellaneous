#include<iostream>
#include<sstream>
#include<string>
#include<iomanip>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>

using namespace std;

//#define USE_DOUBLE
#define RELEASE
//#define uint unsigned int 

// -- RETURN THE FIRST INTEGER MULTIPLE -------------------------------------

//inline int iDivUp(int a, int b) { return ((a % b) != 0) ? (a / b + 1) : (a / b); }

// ---------------------------------------------------------------------------

// -- Stampa le matrici su FILE --------------------------------------------------------------------------------------------------

template<class real_type>
void print_matrix_on_file(real_type *in_mat, unsigned int rows, unsigned int cols, unsigned int batch_size){

	ofstream mat;
	mat.open("matrices.txt");
	mat << "ROWS = " << rows << endl;
	mat << "COLS = " << cols << endl;
	mat << "BATCH_SIZE = " << batch_size << endl;
	
	for(unsigned int num = 0; num < batch_size; num++){
		mat << "\n\n -------------------------- MATRIX "<<num + 1<<" ------------------------------------------------\n\n";
		for(unsigned int i = 0; i < rows; i++){
			for(unsigned int j = 0; j < cols; j++){
				mat << setprecision(3) << in_mat[i*batch_size*cols + j*batch_size + num] << "  ";
			}
			mat << "\n";
		}
		mat << "\n -------------------------------------------------------------------------------------------------\n\n";
	}
	mat.close();
}

// ---------------------------------------------------------------------------------------------------------------------------------------

// -- Stampa diagonale e prima diagonale superiore su FILE --------------------------------------------------------------------------------------------------

template<class real_type>
void print_bidiag_on_file(real_type *d, real_type *e, unsigned int cols, unsigned int batch_size){

	ofstream diag;
	ofstream supdiag;
	diag.open("d.txt");
	supdiag.open("e.txt");

	diag << "BATCH_SIZE = " << batch_size << endl;
	supdiag << "BATCH_SIZE = " << batch_size << endl;
		
	for(unsigned int num = 0; num < batch_size; num++){
		diag << "\n\n -------------------------- DIAG "<<num + 1<<" ------------------------------------------------\n\n";
		supdiag << "\n\n -------------------------- SUPDIAG "<<num + 1<<" ------------------------------------------------\n\n";

		for(unsigned int i = 0; i < cols; i++)
			diag << setprecision(3) << d[i + num*cols] << "  ";
		for(unsigned int j = 0; j < cols - 1; j++)
			supdiag << setprecision(3) << e[j + num*(cols-1)] << "  ";			
	}
	diag << "\n\n -------------------------------------------------------------------------------------------------\n\n";
	supdiag << "\n\n -------------------------------------------------------------------------------------------------\n\n";
	
		
	diag.close();
	supdiag.close();
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/***** TRIM WHITESPACES ****/
//std::string trim(const std::string& str, const std::string& whitespace = " \t"){
//    const size_t strBegin = str.find_first_not_of(whitespace);
//    if (strBegin == std::string::npos)
//        return ""; // no content

//    const size_t strEnd = str.find_last_not_of(whitespace);
//    const size_t strRange = strEnd - strBegin + 1;

//    return str.substr(strBegin, strRange);
//}


//void get_valid_line(std::istream &ifile,std::string &line){
//    do{
//        std::getline(ifile, line);
//        line = trim(line); // trim leading and trailing white spaces
//    }while(ifile.good() && (line[0]=='#'  || line[0]=='*' || line.size()==0));
//    
//}

// ***** COMPUTE ALLOCATED MEMORY ******************************************************* //
template<class real_type>
void malloc_size(uint rows, uint cols, uint batch_size){

    real_type singular_values; // -- CONTARE 1 VOLTA
    real_type input_matrices; // -- CONTARE 2 VOLTE
    real_type diagonal;       // -- CONTARE 2 VOLTE
    real_type off_diagonal;   // -- CONTARE 2 VOLTE
    
    real_type total_memory;
    
    singular_values = cols*batch_size*sizeof(real_type);
    input_matrices = 2*rows*cols*batch_size*sizeof(real_type);
    diagonal = 2*cols*batch_size*sizeof(real_type);
    off_diagonal = 2*(cols-1)*batch_size*sizeof(real_type);
    
    total_memory = singular_values + input_matrices + diagonal + off_diagonal;
    std::cout << "Memory allocated = " << total_memory/(1E6) << " MB" << std::endl;
}
// ************************************************************************************** //


/*****************************/
/* WRITE MATRIX TO FILE REAL */
/*****************************/
template<typename real_type>
void write_matrix_real(real_type* A, int dim, const char* filenameA) {

	ofstream wfile;
	wfile.precision(12);
	wfile.open(filenameA);

		for(int i=0;i<dim;i++){
			wfile << A[i] << " ";
		}


	wfile.close();
}
