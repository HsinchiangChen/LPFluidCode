#include "state_collision.h"
#include <iostream>
#include <cmath>



////////////////////////////////////////////////////////////////////////////////////////
// Start of LeftUniformVelocityState
////////////////////////////////////////////////////////////////////////////////////////

LeftUniformVelocityState::LeftUniformVelocityState(): 
m_fDen(1.), m_fVelX(10.), m_fVelY(0), m_fVelZ(0), m_fPressure(0) {} 


double LeftUniformVelocityState::pressure(double x, double y, double z) {
	return m_fPressure;	
}

double LeftUniformVelocityState::density(double x, double y, double z) {
	return m_fDen;
}

void LeftUniformVelocityState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = m_fVelX;
	vY = m_fVelY;
	vZ = m_fVelZ;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of LeftUniformVelocityState
////////////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////////////
// Start of RightUniformVelocityState
////////////////////////////////////////////////////////////////////////////////////////

RightUniformVelocityState::RightUniformVelocityState(): 
m_fDen(1.), m_fVelX(-10.), m_fVelY(0), m_fVelZ(0), m_fPressure(0) {} 


double RightUniformVelocityState::pressure(double x, double y, double z) {
	return m_fPressure;	
}

double RightUniformVelocityState::density(double x, double y, double z) {
	return m_fDen;
}

void RightUniformVelocityState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = m_fVelX;
	vY = m_fVelY;
	vZ = m_fVelZ;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of RightUniformVelocityState
////////////////////////////////////////////////////////////////////////////////////////
