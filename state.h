/**
 * \file   state.h
 *
 * \brief  This header file contains classes for the initialization of the state of fluid objects
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2014/08/08
 *
 * Created on: 2014/06/07 
 *
 */

#ifndef __STATE_H__
#define __STATE_H__

#include <unordered_map>


/**
 * \class State
 * 
 * \brief An abstract class for the initialization of the state of fluid objects
 *
 * The specified states include pressure, density, and velocity. 
 * These states are assigned based on the initial Cartesian coordinate of particles (\f$x\f$, \f$y\f$, \f$z\f$).
 *  
 * \note 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/08/08 
 *
 * Created on: 2014/06/07 
 *
 */
class State {
public:
	/**
	 * \brief      Calculates pressure based on the Cartesian coordinate (x,y,z) of a particle 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     The calculated pressure value
	 */
	virtual double pressure(double x, double y, double z)=0;
	/**
	 * \brief      Calculates density based on the Cartesian coordinate (x,y,z) of a particle 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     The calculated density value
	 */
	virtual double density(double x, double y, double z)=0;
	/**
	 * \brief       Calculates velocity based on the Cartesian coordinate (x,y,z) of a particle 
	 * \param [in]  x  The x-coordinate
	 * \param [in]  y  The y-coordinate
	 * \param [in]  z  The z-coordinate
	 * \param [out] vX  The calculated velocity value in x-coordinate
	 * \param [out] vY  The calculated velocity value in y-coordinate
	 * \param [out] vZ  The calculated velocity value in z-coordinate
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ)=0;
	/**
	 * \brief Destructor 
	 */
	virtual ~State() {};
};









/**
 * \class GaussianPressureState
 * 
 * \brief A class that implements the Gaussian pressure state
 *
 * The Gaussian state specifies the pressure distribution as Gaussian,
 * with uniform density and uniform zero velocities in all x, y, and z-coordinates
 *   
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/08/08 
 *
 * Created on: 2014/06/07
 * 
 *
 */
class GaussianPressureState: public State {
public:
	/// Constructor 
	GaussianPressureState();
	
	/// Destructor 
	virtual ~GaussianPressureState() {};
	
	/**
	 * \brief         Calculates pressure based on the Cartesian coordinate (x,y,z) of a particle and Gaussian distribution 
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
	 * \brief           Specifies uniform velocity as specified in constructor implementation
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
	double m_fVelX, m_fVelY, m_fVelZ; ///< velocities
	double m_fPCenX, m_fPCenY, m_fPCenZ; ///< center of the Gaussian pressure profile
	double m_fPPeak; ///< peak of pressure value
	double m_fPCoeff; ///< Gaussian coefficient	
};







/**
 * \class UniformVelocityState
 * 
 * \brief A class that implements the uniform velocity state
 *
 * The uniform velocity state specifies the velocity as inward uniform radial velocity,
 * with uniform density and uniform pressure everywhere
 *  
 * \note 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/08/08 
 *
 * Created on: 2014/06/07
 * 
 *
 */
class UniformVelocityState: public State {
public:
	/// Constructor	
	UniformVelocityState();
	/// Destructor
	virtual ~UniformVelocityState() {};
	
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
	 * \brief           Specifies uniform radial inward velocity with magnitude specified in constructor implementation
	 * \param [in]  x   The x-coordinate
	 * \param [in]  y   The y-coordinate
	 * \param [in]  z   The z-coordinate
	 * \param [out] vX  A uniform radial inward velocity value  
	 * \param [out] vY  A uniform radial inward velocity value
	 * \param [out] vZ  A uniform radial inward velocity value
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);	
private:
	double m_fDen; ///< density
	double m_fPressure; ///< pressure
	double m_fVelocity; ///< velocity
	double m_fCenX, m_fCenY, m_fCenZ; ///< center of disk/sphere
};








/**
 * \class StateFactory
 * 
 * \brief The abstract factory class for creating objects in the State family 
 *   
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/08/08 
 *
 * Created on: 2014/06/07
 * 
 *
 */
class StateFactory {
public:
	/**
     * \brief Defines a function pointer pointing to a function which creates objects in the State family 
	 */
	typedef State* (*StateCreateFunc)();
	
	/**
     * \brief   Returns reference to a Singleton object of class StateFactory
	 * 
	 * \param   None 
     *
	 * \return  Reference to a static object of class StateFactory 
	 * 
	 *
	 * Example usage: 
	 * \code
	 *          StateFactory& factory = StateFactory::instance();
	 * \endcode
	 *
	 * \note    This function is implemented based on the lazy Singleton design pattern; 
	 *          therefore, only one StateFactory instance is allowed in each program
	 */
	static StateFactory& instance(); 
	
	/**
     * \brief      Registers (links) the state name \e name 
	 *		       with the function \e func for creating objects in the State family
	 *			   
	 *             After registration, \e name can be used as an argument in the createState member function
	 *             for creating objects of the linked type in the State family
	 *	           
	 *  
	 * \param [in] name the state name 
	 * \param [in] func the function pointer pointing to the function that creates objects of a specific type
	 *             in the State family
	 * 
	 * \return     None  
	 *
	 * \note       Instead of using this function directly, consider using the StateRegistrar class for
	 *		       the purpose of linking a state name and a specific class in the State family. 
	 *             The function is kept public in case one wants to use it directly
	 *		  
	 *
	 */
	void registerState(std::string name, StateCreateFunc func);

	/**
     * \brief      This function creates an object of the class linked to the \name  
	 *
	 * \param [in] name the name linked to a specific class in the State family via the 
	 *				    registrerState member function
	 *  
	 * \return     A State * pointer pointing to an object of a specific class in the State family 
	 * 
	 * Example usage: 
	 * \code
	 *            StateFactory& factory = StateFactory::instance();
	 *            State* newState = factory.createState(name);
	 * \endcode
	 *
	 */
	State* createState(std::string name);	

private:	
	std::unordered_map<std::string, StateCreateFunc> stateTable; ///< hash table for the (name,creatFunction) pair	
	StateFactory() {}; ///< for singleton design pattern
	StateFactory(const StateFactory& other); ///< for singleton design pattern (Don't implement)
	StateFactory& operator=(const StateFactory& other); ///< for singleton design pattern (Don't implement)

};

#endif // __STATE_H__
