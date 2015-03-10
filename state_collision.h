/**
 * \file   state_collision.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects in the 2D collision simulation
 *
 * This file serves as an example that different files (other than state.h and state.cpp) should be used to 
 * model the state of fluid objects in different simulations
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2014/09/08
 *
 * Created on: 2014/09/07 
 *
 */

#ifndef __STATE_COLLISION_H__
#define __STATE_COLLISION_H__

#include "state.h"




/**
 * \class LeftUniformVelocityState
 * 
 * \brief A class that implements the left uniform velocity state
 *
 * This state specifies a uniform velocity toward the \b right hand side along the x-coordinate,
 * with uniform density and uniform zero pressure. The keyword \b Left in the name indicates 
 * it is for the object on the left in the collision simulation
 *   
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/09/08 
 *
 * Created on: 2014/09/07
 * 
 *
 */
class LeftUniformVelocityState: public State {
public:
	/// constructor
	LeftUniformVelocityState();

	/// destructor
	virtual ~LeftUniformVelocityState() {};
	
	/**
	 * \brief         Specifies constant pressure as specified in construtor implementation 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        A constnat pressure value 
	 */
	virtual double pressure(double x, double y, double z);
	
	/**
	 * \brief         Specifies a constant value as specified in constructor implementation
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        A constant density value
	 */
	virtual double density(double x, double y, double z);
	
	/**
	 * \brief           Specifies constant velocity to the right-hand-side along the x-coordinate 
	 *                  with magnitude specified in constructor implementation
	 * \param [in]  x   The x-coordinate
	 * \param [in]  y   The y-coordinate
	 * \param [in]  z   The z-coordinate
	 * \param [out] vX  A constant value (>0)   
	 * \param [out] vY  zero
	 * \param [out] vZ  zero
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
	double m_fDen;
	double m_fVelX, m_fVelY, m_fVelZ;
	double m_fPressure;	
};










/**
 * \class RightUniformVelocityState
 * 
 * \brief A class that implements the right uniform velocity state
 *
 * This state specifies a uniform velocity toward the \b left hand side along the x-coordinate,
 * with uniform density and uniform zero pressure. The keyword \b Right in the name indicates 
 * it is for the object on the right in the collision simulation
 *   
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/09/08 
 *
 * Created on: 2014/09/07
 * 
 *
 */
class RightUniformVelocityState: public State {
public:
	/// construtor
	RightUniformVelocityState();
	
	/// destructor
	virtual ~RightUniformVelocityState() {};
	
	/**
	 * \brief         Specifies constant pressure as specified in construtor implementation 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        A constnat pressure value 
	 */
	virtual double pressure(double x, double y, double z);
	
	/**
	 * \brief         Specifies a constant value as specified in constructor implementation
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        A constant density value
	 */
	virtual double density(double x, double y, double z);
	
	/**
	 * \brief           Specifies constant velocity to the right-hand-side along the x-coordinate 
	 *                  with magnitude specified in constructor implementation
	 * \param [in]  x   The x-coordinate
	 * \param [in]  y   The y-coordinate
	 * \param [in]  z   The z-coordinate
	 * \param [out] vX  A constant value (>0)   
	 * \param [out] vY  zero
	 * \param [out] vZ  zero
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);	
private:
	double m_fDen;
	double m_fVelX, m_fVelY, m_fVelZ;
	double m_fPressure;	
};

#endif //__STATE_COLLISION_H__
