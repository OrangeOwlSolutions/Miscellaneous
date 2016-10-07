#include <vector>

#include "setAntennas.h"
#include "Utilities.cuh"

extern "C" void setAntennas(segment *TXline, segment *RXline, vec3 origintx, vec3 originrx, std::vector<vec3> &azr_tx, std::vector<vec3> &azr_rx, 
	             TX *tx, RX *rx, std::vector<std::string> &pattern_tx, std::vector<std::string> &pattern_rx, vec3 dtx, vec3 drx, float radius,
				 std::string txPatternFileName, std::string rxPatternFileName, const int num_TX, const int num_RX) {

    //azr_tx.resize(num_TX);
    for (int j = 0; j < num_TX; j++) {
		azr_tx[j].x = deg2rad(0.);
		azr_tx[j].y = deg2rad(0.);
		azr_tx[j].z = deg2rad(0.);
    }

    //azr_rx.resize(num_RX);
    for (int j = 0; j < num_RX; j++) {
		azr_rx[j].x = deg2rad(0.);
		azr_rx[j].y = deg2rad(0.);
		azr_rx[j].z = deg2rad(0.);  
    } 

	float Ltx	=  1.;													// --- Transmitting antenna segment half-orientation
	float Lrx	=  1.;													// --- Receiving antenna segment half-orientation
	
	// --- Transmitting antenna segment
	TXline->orig.x	= origintx.x;
	TXline->orig.y	= origintx.y;
	TXline->orig.z	= origintx.z;
	TXline->dir_v.x = dtx.x;
	TXline->dir_v.y = dtx.y;
	TXline->dir_v.z = dtx.z;
	TXline->L = Ltx;

	// --- Receiving antenna segment
	RXline->orig.x	= originrx.x;
	RXline->orig.y	= originrx.y;
	RXline->orig.z	= originrx.z;
	RXline->dir_v.x = drx.x;
	RXline->dir_v.y = drx.y;
	RXline->dir_v.z = drx.z;
	RXline->L = Lrx;

	tx->patt = new pattern(txPatternFileName);

	tx->azimuth  = deg2rad(0.f);
	tx->zenith   = deg2rad(0.f);
	tx->roll     = deg2rad(0.f);

	rx->azimuth  = deg2rad(0.f);
	rx->zenith 	 = deg2rad(0.f);
	rx->roll     = deg2rad(0.f);

	rx->radius	= radius;													// --- Radius of the detector

	rx->patt = new pattern(rxPatternFileName);

    pattern_tx.resize(num_TX);
    //for (uint j = 0; j < num_TX; j++) pattern_tx[j].assign("./patterns/tx_1.txt");  
	pattern_tx[0].assign("antenna_pattern_3.txt");
	pattern_tx[1].assign("antenna_pattern_1.txt");

    pattern_rx.resize(num_RX); 
	//for (uint j = 0; j < num_RX; j++) pattern_rx[j].assign("./patterns/rx_1.txt");  
	pattern_rx[0].assign("antenna_pattern_4.txt");
	pattern_rx[1].assign("antenna_pattern_2.txt");

}

/*****************************************/
/* SET THE TX AND RX PATTERNS OF A GROUP */
/*****************************************/
void set_pattern(const std::vector<std::string> &pattern, std::string *pattern_x, unsigned int NX) {
    
    for(unsigned int j = 0; j < NX; j++) pattern_x[j] = pattern[j]; 
}

/********************************************/
/* SET AZIMUTH, ZERNITH AND ROLL OF A GROUP */
/********************************************/
extern "C" void set_azr(const std::vector<vec3> &azr, vec3 *azr_x, unsigned int NX) {

    for(unsigned int j = 0; j < NX; j++) azr_x[j] = azr[j];

}
