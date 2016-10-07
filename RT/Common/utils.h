#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <string>
#include <vector>

// Common host and device function 
//inline 
//int iDivUp(const int a,const  int b){
//    return ((a % b) != 0) ? (a / b + 1) : (a / b);
//}

inline 
int iDivDown(const  int a,const  int b){
    return a / b;
}
 
inline  
int iAlignUp(const  int a,const  int b){
    return ((a % b) != 0) ?  (a - a % b + b) : a;
}

inline 
int iAlignDown(const  int a,const  int b){
    return a - a % b;
}




std::string trim(const std::string& str, const std::string& whitespace = " \t");

bool get_valid_line(std::istream &ifile,std::string &line);


std::string string_to_upper( std::string strToConvert);

//void linspace(
//    std::vector<double>::iterator begin,
//    std::vector<double>::iterator end,
//    double minv,
//    double maxv    
//);

//double deg2rad(double deg);


#endif
