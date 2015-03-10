/**
 * \file   registrar.h
 *
 * \brief  This header file contains template registrar classes for the Geometry and State family
 *
 * 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/08/10 
 *
 * Created on: 2014/06/12 
 *
 */

#ifndef __REGISTRAR_H__
#define __REGISTRAR_H__

#include "geometry.h"
#include "state.h"



/**
 * \class GeometryRegistrar
 * 
 * \brief A template class for the registration of classes in the Geometry family
 *
 * The purpose of this class is to have a unified place to register objects in the Geometry family.
 * Registration means linking a geometry name (given in the argument list of the constructor) 
 * and a class in the Geometry family.
 * It is implemented as a template class to accomondate various child classes in the Geometry family. 
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
template<typename Derived>
class GeometryRegistrar {
public:	
	/**
     * \brief  constructor  
	 *
	 * Registers the \e name given in the argument list. 
	 * Registration means linking \e name to a class in the Geometry family
	 *
	 * \param  name a geometry name           
	 *
	 * Example usage: Registers \e a_name with class \e SomeGeometry
	 * \code
	 *         GeometryRegistrar<SomeGeometry> g("a_name");
	 * \endcode
	 */
	GeometryRegistrar(std::string name);
	
	/**
     * \brief   This function creates an object of the \e Derived class and returns a pointer to the Geometry class  
	 *
	 * \param   None   
	 *  
	 * \return  A Geoemtry * pointer that points to an object of \e Derived class  
	 *
	 */
	static Geometry* createFunc();
};



////////////////////////////////////////////////////////////////////////////////////////
// Start of GeometryRegistrar
////////////////////////////////////////////////////////////////////////////////////////

template<typename Derived>
Geometry* GeometryRegistrar<Derived>::createFunc() {
	return new Derived();
}


template<typename Derived>
GeometryRegistrar<Derived>::GeometryRegistrar(std::string name) {
	GeometryFactory& factory = GeometryFactory::instance(); // & cannot be omitted
	factory.registerGeometry(name, GeometryRegistrar<Derived>::createFunc); // second arg is function pointer
}

////////////////////////////////////////////////////////////////////////////////////////
// Start of GeometryRegistrar
////////////////////////////////////////////////////////////////////////////////////////








/**
 * \class StateRegistrar
 * 
 * \brief A template class for the registration of classes in the State family
 *
 * The purpose of this class is to have a unified place to register objects in the State family.
 * Registration means linking a state name (given in the argument list of the constructor) 
 * and a class in the State family.
 * It is implemented as a template class to accomondate various child classes in the State family. 
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
template<typename Derived>
class StateRegistrar {
public:
	/**
     * \brief  constructor 
	 *
	 * Registers the \e name given in the argument list.
	 * Registration means linking \e name to a class in the Geometry family.
	 *
	 * \param  name a state name           
	 *
	 * Example usage: Registers \e a_name with class \e SomeState
	 * \code
	 *         StateRegistrar<SomeState> s("a_name");
	 * \endcode
	 */
	StateRegistrar(std::string name);
	/**
     * \brief   This function creates an object of the \e Derived class and returns a pointer to the State class 
	 *
	 * \param   None   
	 *  
	 * \return  A State * pointer that points to an object of \e Derived class  
	 *
	 */
	static State* createFunc();	
};


////////////////////////////////////////////////////////////////////////////////////////
// Start of StateRegistrar
////////////////////////////////////////////////////////////////////////////////////////

template<typename Derived>
State* StateRegistrar<Derived>::createFunc() {
	return new Derived();
}

template<typename Derived>
StateRegistrar<Derived>::StateRegistrar(std::string name) {
	StateFactory& factory = StateFactory::instance(); // note here is the difference from the GeometryRegistrar class
	factory.registerState(name, StateRegistrar<Derived>::createFunc);
}

////////////////////////////////////////////////////////////////////////////////////////
// End of StateRegistrar
////////////////////////////////////////////////////////////////////////////////////////

#endif // __REGISTRAR_H__

