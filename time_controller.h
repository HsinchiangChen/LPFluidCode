/**
 * \file   time_controller.h
 *
 * \brief  This header file contains classes for simulation time controls 
 *
 * The task of time controllers is to determine the length of time between any two iterations of 
 * the simulation and call the main solver. For example, a time controller which uses a 
 * constant time stepping independent of the results of the main solver, or a time stepping
 * calculated by the CFL condition, can be constructed by inheriting the TimeController class.
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/07/12
 *
 * Created on: 2014/07/02 
 *
 */


#ifndef __TIME_CONTROLLER_H__
#define __TIME_CONTROLLER_H__

#include <cstddef>
#include <vector>
#include <fstream>

class LPSolver; 
class ParticleViewer; 

/**
 * \class TimeController
 * 
 * \brief An abstract class for different types of simulation time controllers
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/07/12 
 *
 * Created on: 2014/07/02 
 *
 */
class TimeController {
public:
	/**
	 * \brief   Calls main solvers of the simulation and determines the time stepping between iterations 
	 *
	 * The simulation is performed by calling the main LP solver for each iteration 
	 * and the time length between iterations is determined by the type of time controller adopted 
	 *
	 * \param   None
	 * \return  0 if the entire simulation runs successfully to the end; 1 if main solver fails during any iteration
	 */
	virtual int solve() = 0;
protected:
	LPSolver* m_pSolver; ///< A pointer to the Lagrangian Particle solver
	std::vector<ParticleViewer*> m_vViewers; ///< A vector containing pointers to different types of particle viewers	
	double m_fTime; ///< Current physical time of simulation
	double m_fEndTime; ///< End physical time of simulation
	double m_fWriteTimeInterval; ///< The physical time interval between two writting of results
	double m_fNextWriteTime; ///< The next physical time point to write results
	double m_fDt; ///< The physical time interval between two iterations of the simulation
	std::size_t m_iWriteStep; ///< The number of times results are written to the particle veiwer
	bool m_iIfDebug;///< if true then print debug info
	std::ofstream debug;///< output information for debugging

	/**
	 * \brief   Adjusts the iteration time step when current time approaches the time for results writting 
	 *
	 * m_fDt will be adjusted(shrinked) when m_fTime + m_fDt > m_fNextWriteTime, 
	 * to match the exact time of results writting
	 *
	 * \param   None
	 * \return  \c true if m_fDt is adjusted; \c false otherwise
	 */
	bool adjustDtByWriteTimeInterval();
};









class Initializer;


/**
 * \class DefaultTimeController
 * 
 * \brief This time controller calculates the time stepping based on the CFL condition
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/07/12 
 *
 * Created on: 2014/07/02 
 *
 */
class DefaultTimeController : public TimeController {

public:		
	/**
	 * \brief              Constructor
	 *
	 * The constructor initializes parameters based on \init and also initializes the main solver and viewer
	 * based on \e solver and \viewers used to perform the main iterations and to write results for viewing
	 *
	 * \param [in] init    Const reference to an object of the Initializer class
	 * \param [in] solver  Pointer to an object of the LPSolver family
	 * \param [in] viewers Vector containing pointers to objects of the ParticleViewer family 
	 *
	 */
	DefaultTimeController(const Initializer& init, LPSolver* solver, const std::vector<ParticleViewer*>& viewers);
	
	/**
	 * \brief   Calls main solvers of the simulation and determines the time stepping between iterations 
	 *
	 * The simulation is performed by calling the main LP solver for each iteration 
	 * and the time length between iterations is determined by the CFL condition 
	 *
	 * \param   None
	 * \return  0 if the entire simulation runs successfully to the end; 1 if main solver fails during any iteration
	 */
	virtual int solve();

private:
	double m_fCFLCoeff; ///< a multiplier term which is between 0 and 1	used to shrink time stepping
	void computeDtByCFL(); ///< computes m_fDt based on the CFL condition

};


#endif // __TIME_CONTROLLER_H__
