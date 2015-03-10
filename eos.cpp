#include "eos.h"
#include <iostream>
#include <cassert>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// Start of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////

double PolytropicGasEOS::getEnergy(double pressure, double density) {
	if(((m_fGamma - 1.) * density) != 0) 
		return ( pressure / ((m_fGamma - 1.) * density) );
	else { // divide by zero
		std::cout<<"Error (Divide by zero)! Computing energy by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", density = "<<density<<std::endl;
		assert(false);
	}
}

double PolytropicGasEOS::getSoundSpeed(double pressure, double density) {
	double cs;
	if(density != 0)
		cs = m_fGamma * pressure / density;
	else {
		std::cout<<"Error (Divide by zero)! Computing sound speed by EOS: "<<std::endl;
		std::cout<<"density = "<<density<<std::endl;
		assert(false);
	}
	if(cs > 0) 
		return sqrt(cs);
	else { // taking square root of a negative number
		std::cout<<"Error (Taking suqre root of 0 or a negative number)! Computing sound speed by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", pressure = "<<pressure<<", density = "<<density<<std::endl;
		assert(false);
	}
}

////////////////////////////////////////////////////////////////////////////////
// End of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////










////////////////////////////////////////////////////////////////////////////////
// Start : StiffPolytropicEOS
////////////////////////////////////////////////////////////////////////////////

double StiffPolytropicGasEOS::getEnergy(double pressure, double density) {
	if(((m_fGamma - 1.) * density) != 0) 
		return ( (pressure + m_fGamma * m_fPinf) / ((m_fGamma - 1.) * density) - m_fEinf);
	else {
		std::cout<<"Error (Divide by zero)! Computing energy by EOS: "<<std::endl;
		std::cout<<"gamma = "<<m_fGamma<<", density = "<<density<<std::endl;
		assert(false);
	}
  
}


double StiffPolytropicGasEOS::getSoundSpeed(double pressure, double density) {
	double cs;
	if(density != 0)
		cs = m_fGamma * (pressure + m_fPinf) / density;
	else {
		std::cout<<"Error (Divide by zero)! Computing sound speed by EOS: "<<std::endl;
		std::cout<<"density = "<<density<<std::endl;
		assert(false);
	}
	if(cs > 0)
		return sqrt(cs);
	else { 
		std::cout<<"Error (Taking suqre root of 0 or a negative number)! Computing sound speed by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", pinf = "<<m_fPinf<<", pressure = "<<pressure<<", density = "<<density<<std::endl;
		assert(false);
	}  
}

////////////////////////////////////////////////////////////////////////////////
// End of StiffPolytropicEOS
////////////////////////////////////////////////////////////////////////////////
