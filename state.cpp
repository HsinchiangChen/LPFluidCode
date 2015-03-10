#include "state.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of GaussianPressureState
////////////////////////////////////////////////////////////////////////////////////////

GaussianPressureState::GaussianPressureState(): m_fDen(1.), m_fVelX(0), m_fVelY(0), m_fVelZ(0), m_fPCenX(0), m_fPCenY(0), m_fPCenZ(0), m_fPPeak(7.), m_fPCoeff(-100.) {} 


double GaussianPressureState::pressure(double x, double y, double z) {
	return m_fPPeak * exp(m_fPCoeff*( (x-m_fPCenX)*(x-m_fPCenX)+(y-m_fPCenY)*(y-m_fPCenY)+(z-m_fPCenZ)*(z-m_fPCenZ) ) );	
}

double GaussianPressureState::density(double x, double y, double z) {
	return m_fDen;
}

void GaussianPressureState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = m_fVelX;
	vY = m_fVelY;
	vZ = m_fVelZ;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of GaussianPressureState
////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////
// Start of UniformVelocityState
////////////////////////////////////////////////////////////////////////////////////////

UniformVelocityState::UniformVelocityState(): m_fDen(1.66e-5), m_fPressure(1e-9), m_fVelocity(1e4), m_fCenX(0), m_fCenY(0), m_fCenZ(0) {} 


double UniformVelocityState::pressure(double x, double y, double z) {
	return m_fPressure;	
}

double UniformVelocityState::density(double x, double y, double z) {
	return m_fDen;
}

void UniformVelocityState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	double dist = sqrt((x-m_fCenX)*(x-m_fCenX)+(y-m_fCenY)*(y-m_fCenY)+(z-m_fCenZ)*(z-m_fCenZ));
	vX = -x/dist*m_fVelocity;
	vY = -y/dist*m_fVelocity;
	vZ = -z/dist*m_fVelocity;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of UniformVelocityState
////////////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////////////
// Start of StateFactory
////////////////////////////////////////////////////////////////////////////////////////

StateFactory& StateFactory::instance() { // singleton
	static StateFactory stateFactory;
	return stateFactory;
}

State* StateFactory::createState(std::string name) {
	const auto result = stateTable.find(name);
	if(result==stateTable.end()) {
		std::cout<<"This state class name is not registered!!!"<<std::endl;
		return nullptr;
	}
	return (result->second)();
}

void StateFactory::registerState(std::string name, StateCreateFunc func) {
	stateTable.insert({name,func});
}

////////////////////////////////////////////////////////////////////////////////////////
// End of StateFactory
////////////////////////////////////////////////////////////////////////////////////////


