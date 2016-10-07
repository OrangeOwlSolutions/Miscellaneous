//#include <sstream> 
//#include <fstream> 
// 
//#include "edge.h"
//#include "utils.h"
//
///*******************/
///* PRINT EDGE LIST */
///*******************/
//std::ostream &operator << (std::ostream & os, const edge& e){
//    return os << e.ID << " " <<  e.a << "    " << e.b;
//}
//
///***************************/
///* LOAD EDGE LIST FUNCTION */
///***************************/
//void load_edge_list(const std::string &filename, const unsigned int first_ID, std::vector<edge> &edge_list) {
//    
//	std::string		line;
//    std::ifstream	ifile(filename.c_str());
//    
//    edge_list.clear();
//    
//    // --- read the number of edges
//	size_t num_edges;
//    {
//        get_valid_line(ifile,line);
//        std::stringstream ss(line);        
//        ss >>num_edges;
//    }
//    
//    // --- Read the edges
//	{
//        for(size_t i =0; i < num_edges; i++){
//            get_valid_line(ifile,line);
//            std::stringstream ss(line);        
//            edge ed;                
//            ss >> ed.a.x >> ed.a.y >> ed.a.z >> ed.b.x >> ed.b.y >> ed.b.z;
//            
//            ed.ID = first_ID + i;
//            edge_list.push_back(ed);                              
//        }
//    }
//    
//
//}
//
//
//
////TODO use newton raphson 
//bool get_diffraction_point(const vec3 &tx_pos,const vec3 &rx_pos, const edge &ed, vec3 &diff_point){
//    
//    
//    vec3 bma = ed.b - ed.a;
//    
//    double uo =  0;
//    double uoo = 1;
//    
//    
//    //metodo delle secanti con 10 iterazioni fisse ma in genere finisce molto prima
//    for(int i =0; i < 10; i ++){
//        vec3 Po =  ed.a+uo*bma;  //p old
//        vec3 Poo =  ed.a+uoo*bma; //p old old
//        
//        //std::cerr << "Po " << Po << " Poo " << Poo << std::endl;
//    
//        double fo = dot(bma,(normalize(Po-tx_pos)+normalize(Po-rx_pos)));
//        double foo = dot(bma,(normalize(Poo-tx_pos)+normalize(Poo-rx_pos)));
//        //std::cerr << "fo " << fo << " foo " << foo << std::endl;
//        
//        if(fo*foo >= 0 && i == 0)// se non c'e il nullo non continuare
//            return false;
//        
//        if(fabs(fo-foo)/norm(bma) < 0.001){ //quando stai vicino alla soluzione  esci
//            break;
//        }
//        double un = uo - fo * (uo-uoo)/(fo - foo);
//        //std::cerr << "uo " << uo << " uoo " << uoo << " un " << un << std::endl;
//        
//        uoo = uo;
//        uo = un;        
//    }
//    
//    if(uo <0 || uo > 1){ //questo non dovrebbe verificarsi mai
//        return false ;
//    }
//    
//    diff_point = ed.a+uo*bma;
//    
//    return true;
//}
#include <sstream> 
#include <fstream> 
 
#include "edge.h"
#include "utils.h"

/*******************/
/* PRINT EDGE LIST */
/*******************/
std::ostream &operator << (std::ostream & os, const edge& e){
    return os << e.ID << " " <<  e.a << "    " << e.b;
}

/***************************/
/* LOAD EDGE LIST FUNCTION */
/***************************/
void load_edge_list(const std::string &filename, const unsigned int first_ID, std::vector<edge> &edge_list) {
    
	std::string		line;
    std::ifstream	ifile(filename.c_str());
    
    edge_list.clear();
    
    // --- read the number of edges
	size_t num_edges;
    {
        get_valid_line(ifile, line);
        std::stringstream ss(line);        
        ss >>num_edges;
    }
    
    // --- Read the edges
	{
        for(size_t i =0; i < num_edges; i++){
            get_valid_line(ifile, line);
            std::stringstream ss(line);        
            edge ed;                
            ss >> ed.a.x >> ed.a.y  >> ed.a.z  >> ed.b.x  >> ed.b.y  >> ed.b.z;
            ss >> ed.n   >> ed.no.x >> ed.no.y >> ed.no.z >> ed.nn.x >> ed.nn.y >> ed.nn.z;
            
            ed.ID = first_ID + i;
            edge_list.push_back(ed);                              
        }
    }
    

}

//TODO use newton raphson 
bool get_diffraction_point(const vec3 &tx_pos,const vec3 &rx_pos, const edge &ed, vec3 &diff_point){
    
    
    vec3 bma = ed.b - ed.a;
    
    double uo =  0;
    double uoo = 1;
    
    
    //metodo delle secanti con 10 iterazioni fisse ma in genere finisce molto prima
    for(int i =0; i < 10; i ++){
        vec3 Po =  ed.a+uo*bma;  //p old
        vec3 Poo =  ed.a+uoo*bma; //p old old
        
        //std::cerr << "Po " << Po << " Poo " << Poo << std::endl;
    
        double fo = dot(bma,(normalize(Po-tx_pos)+normalize(Po-rx_pos)));
        double foo = dot(bma,(normalize(Poo-tx_pos)+normalize(Poo-rx_pos)));
        //std::cerr << "fo " << fo << " foo " << foo << std::endl;
        
        if(fo*foo >= 0 && i == 0)// se non c'e il nullo non continuare
            return false;
        
        if(fabs(fo-foo)/norm(bma) < 0.001){ //quando stai vicino alla soluzione  esci
            break;
        }
        double un = uo - fo * (uo-uoo)/(fo - foo);
        //std::cerr << "uo " << uo << " uoo " << uoo << " un " << un << std::endl;
        
        uoo = uo;
        uo = un;        
    }
    
    if(uo <0 || uo > 1){ //questo non dovrebbe verificarsi mai
        return false ;
    }
    
    diff_point = ed.a+uo*bma;
    
    return true;
}
