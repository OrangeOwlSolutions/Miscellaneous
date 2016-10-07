#include <algorithm>

#include "math_and_physical_constants.h"
#include "utils.h"

/***** TRIM ****/
std::string trim(const std::string& str, const std::string& whitespace){
    const size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const size_t strEnd = str.find_last_not_of(whitespace);
    const size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

/***************************/
/* GET VALID LINE FUNCTION */
/***************************/
bool get_valid_line(std::istream &ifile, std::string &line) {
    
	do {
        std::getline(ifile, line);
		line = trim(line); // trim leading and trailing white spaces
    } while (ifile.good() && (line[0]=='#' || line.size()==0));
    
    return ifile.good();
}


std::string string_to_upper(std::string strToConvert){
    
    std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::toupper);

    return strToConvert;
}


//void linspace(
//    std::vector<double>::iterator begin,
//    std::vector<double>::iterator end,
//    double minv,
//    double maxv    
//){
//    
//    if(end-begin > 1){
//        double delta = (maxv-minv)/(end-begin-1);
//        for(std::vector<double>::iterator it = begin;
//            it != end;
//            ++it
//            ){
//            *it = minv+(it-begin)*delta;
//        }
//    }else{
//        *begin = minv;
//    }
//}


//double deg2rad(double deg){
//    return deg*PI_R/180;
//}
