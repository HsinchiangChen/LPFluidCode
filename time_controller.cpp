#include "time_controller.h"
#include "lp_solver.h"
#include "particle_viewer.h"
#include "initializer.h"
#include <cassert>
#include <iostream>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
// Start of TimeController
////////////////////////////////////////////////////////////////////////////////

bool TimeController::adjustDtByWriteTimeInterval() {
	if(m_fTime+m_fDt >= m_fNextWriteTime) {
		m_fDt = m_fNextWriteTime - m_fTime;
		m_fNextWriteTime += m_fWriteTimeInterval;
		if(m_fNextWriteTime > m_fEndTime) m_fNextWriteTime = m_fEndTime;
		assert(m_fDt >= 0);
		//cout<<"-------TimeController::adjustDtByWriteTimeInterval()-------"<<endl;
		//cout<<"m_fDt="<<m_fDt<<endl;
		//cout<<"-----------------------------------------------------------"<<endl;
		return true; // m_fDt adjusted
	}
	return false; // m_fDt did not get adjusted
}

////////////////////////////////////////////////////////////////////////////////
// End of TimeController
////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////
// Start of DefaultTimeController
////////////////////////////////////////////////////////////////////////////////

DefaultTimeController::DefaultTimeController(const Initializer& init, LPSolver* solver, const vector<ParticleViewer*>& viewers) {
	
	// The following are protected data members in the base class
	m_pSolver = solver;
	m_vViewers = viewers;

	m_fTime    = init.getStartTime();
	m_fEndTime = init.getEndTime();	
	m_fWriteTimeInterval = init.getWriteTimeInterval();
	m_fNextWriteTime = m_fTime + m_fWriteTimeInterval;
	if(m_fNextWriteTime > m_fEndTime) m_fNextWriteTime = m_fEndTime;
	m_fDt = 0;
	m_iWriteStep = 0;
	
	// private data member in this class
	m_fCFLCoeff = init.getCFLCoeff();
	
	m_iIfDebug = init.getIfDebug();
	debug.open(init.getDebugfileName(), std::ofstream::out | std::ofstream::app);
}


int DefaultTimeController::solve() {
	
	// visualization at zero time step
	for(auto pViewer:m_vViewers) {
		pViewer->writeResult(m_fTime, m_iWriteStep);
	}
	m_iWriteStep++;
	
	size_t iterationStep = 0; // counter of the number of iterations done
	// start iterations
	while(m_fTime < m_fEndTime) {
		
		// compute dt based on CFL condition
		computeDtByCFL();
		
		// adjust dt to fit into writeTimeInterval if needed
		bool isWriteTime = adjustDtByWriteTimeInterval();
		
		// print out the time and dt before solving
		printf("Time=%.16g, dt=%.16g, iterationStep=%ld\n",m_fTime, m_fDt, iterationStep);
		
		// call LPSolver to solve for this time step	
		int isIterSuccess = m_pSolver->solve(m_fDt); // 0 = success	
		
		if(isIterSuccess!=0) {
			cout<<"Time = "<<m_fTime<<": solver fails!!!"<<endl;
			return 1;
		}
		
		// increment time
		m_fTime += m_fDt;
		
		// increment the counter of number of iterations
		iterationStep++;
		
		// write results if necessary
		if(isWriteTime) { 
			printf("Time=%.16g, writeStep=%ld\n",m_fTime, m_iWriteStep); // check the write time is correct
			for(auto pViewer:m_vViewers) {
				pViewer->writeResult(m_fTime, m_iWriteStep);
			}
			m_iWriteStep++;	
		}

	}
	
	return 0;
}


void DefaultTimeController::computeDtByCFL() {
	double minParticleSpacing = m_pSolver->getMinParticleSpacing();
	double maxSoundSpeed = m_pSolver->getMaxSoundSpeed();
	double maxFluidVelocity = m_pSolver->getMaxFluidVelocity();
	double maxS = max(maxFluidVelocity,maxSoundSpeed);
	assert(maxS!=0);
	m_fDt = m_fCFLCoeff * minParticleSpacing / maxS;
	if(m_fDt > m_fWriteTimeInterval) { 
		cout<<"time = "<<m_fTime<<": m_fDt > m_fWriteTimeInterval!!!"<<endl; 
	}
	if(m_iIfDebug) {
		debug.precision(16);
		//debug<<"-------DefaultTimeController::computeDtByCFL()-------"<<endl;
		debug<<"m_fDt="<<m_fDt<<endl;
		//debug<<"-----------------------------------------------------"<<endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
// End of DefaultTimeController
////////////////////////////////////////////////////////////////////////////////




