#include "state_1d.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of GaussianPressure1DState
////////////////////////////////////////////////////////////////////////////////////////

GaussianPressure1DState::GaussianPressure1DState(): m_fDen(0.01), m_fVelX(0), m_fPCenX(0), m_fPPeak(2.), m_fPCoeff(-100.) {} 


double GaussianPressure1DState::pressure(double x, double y, double z) {
	return 5. + m_fPPeak * exp(m_fPCoeff*( (x-m_fPCenX)*(x-m_fPCenX) ) );	
}

double GaussianPressure1DState::density(double x, double y, double z) {
	return m_fDen;
}

void GaussianPressure1DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = m_fVelX;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of GaussianPressure1DState
////////////////////////////////////////////////////////////////////////////////////////
