/**
 * \file   state_1d.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of 1D fluid objects  
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2015/03/13
 *
 * Created on: 2015/03/13 
 *
 */

#ifndef __STATE_1D_H__
#define __STATE_1D_H__

#include "state.h"




/**
 * \class GaussianPressure1DState
 * 
 * \brief A class that implements the Gaussian-pressure state on a line
 *
 * This state specifies a Gaussian-profile pressure state with uniform density and uniform zero velocity 
 *   
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/03/13
 *
 * Created on: 2015/03/13
 * 
 *
 */
class GaussianPressure1DState: public State {
public:
	/** 
	 * \brief Constructor
	 *
	 * a 1D Gaussian-profile pressure state with uniform density and uniform zero velocity
	 *
	 *
	 *
	 */
	GaussianPressure1DState();

	/// destructor
	virtual ~GaussianPressure1DState() {};
	
	/**
	 * \brief         Calculates pressure based on the Cartesian coordinate x of a particle and Gaussian distribution 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        The calculated pressure value 
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
	 * \brief           Specifies uniform zero velocity as specified in constructor implementation
	 * \param [in]  x   The x-coordinate
	 * \param [in]  y   The y-coordinate
	 * \param [in]  z   The z-coordinate
	 * \param [out] vX  A constant velocity value 
	 * \param [out] vY  A constant velocity value
	 * \param [out] vZ  A constant velocity value
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
	double m_fDen; ///< density
	double m_fVelX; ///< velocity
	double m_fPCenX; ///< center of the Gaussian pressure profile
	double m_fPPeak; ///< peak of pressure value
	double m_fPCoeff; ///< Gaussian coefficient	
};




#endif //__STATE_1D_H__
