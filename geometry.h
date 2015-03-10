/**
 * \file   geometry.h
 *
 * \brief  This header file contains classes for the initialization of the geometry of fluid objects
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

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <unordered_map>

/**
 * \class Geometry
 * 
 * \brief An abstract class for the initialization of the geometry of fluid objects
 *
 * This class can be used as a function object: the operator() is overloaded
 * as the level set function of particle geometry 
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
class Geometry {
public:
	/// Destructor
	virtual ~Geometry() {}
	
	/**
	 * \brief      Level set function of particle geometry 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     \c true if (x,y,z) is inside the level set function; \c false otherwise
	 */
	virtual bool operator()(double x, double y, double z) const=0; 
	
	/**
	 * \brief       Calculates the bounding box of a geometric shape
	 *	
     *              For example, for a ball in 3D with center (xcen, ycen, zcen) and radius r
	 *			    the bounding box is calculated as 
	 *			    ([xcen-r,xcen+r],[ycen-r,ycen+r],[zcen-r,zcen+r]) 
	 *
	 * \param [out] xmin  The minimum in x-coordinate
	 * \param [out] xmax  The maximum in x-coordinate
	 * \param [out] ymin  The minimum in y-coordinate
	 * \param [out] ymax  The maximum in y-coordinate
	 * \param [out] zmin  The minimum in z-coordinate
	 * \param [out] zmax  The maximum in z-coordinate
	 * \return      None
	 */
	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax)=0;
	
};









/**
 * \class Ball
 * 
 * \brief Supply functions for generating a 3D ball geometry
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
class Ball: public Geometry {
public:
	/// Constructor
	Ball();
	
	/// Destructor
	virtual ~Ball() {}
	
	/**
	 * \brief         Level set function of a 3D ball 
	 *
	 * The level set function is 
	 * \f$
	 *   (x-\mbox{xcen})^2+(y-\mbox{ycen})^2+(z-\mbox{zcen})^2 \leq r^2   
	 * \f$
	 * , where (xcen,ycen,zcen) and \e r refer to the center and the radius of the ball
	 *
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     \c true if (x,y,z) is inside the level set function; \c false otherwise
	 */
	virtual bool operator()(double x, double y, double z) const;
	
	/**
	 * \brief       Calculates the bounding box of 3D ball 
	 *
	 * \param [out] xmin  The minimum in x-coordinate
	 * \param [out] xmax  The maximum in x-coordinate
	 * \param [out] ymin  The minimum in y-coordinate
	 * \param [out] ymax  The maximum in y-coordinate
	 * \param [out] zmin  The minimum in z-coordinate
	 * \param [out] zmax  The maximum in z-coordinate
	 * \return      None
	 */
	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);	
private:
	double radius; ///< ball radius
	double xCen; ///< ball center in x-coordinate
	double yCen; ///< ball center in y-coordinate
	double zCen; ///< ball center in z-coordinate
};









/**
 * \class Disk
 * 
 * \brief Supply functions for generating a 2D disk geometry
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
class Disk: public Geometry {
public:
	/// Constructor	
	Disk();
	
	/// Desturctor
	virtual ~Disk() {}
	
	/**
	 * \brief         Level set function of a 2D disk 
	 *
	 * The level set function is 
	 * \f$
	 *   (x-\mbox{xcen})^2+(y-\mbox{ycen})^2 \leq r^2   
	 * \f$
	 * , where (xcen,ycen) and \e r refer to the center and the radius of the disk 
	 *
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     \c true if (x,y,z) is inside the level set function; \c false otherwise
	 */
	virtual bool operator()(double x, double y, double z) const;	
	
	/**
	 * \brief       Calculates the bounding box of 2D disk 
	 *
	 * \param [out] xmin  The minimum in x-coordinate
	 * \param [out] xmax  The maximum in x-coordinate
	 * \param [out] ymin  The minimum in y-coordinate
	 * \param [out] ymax  The maximum in y-coordinate
	 * \param [out] zmin  The minimum in z-coordinate
	 * \param [out] zmax  The maximum in z-coordinate
	 * \return      None
	 */	
	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);	
private:
	double radius; ///< disk radius
	double xCen; ///< disk center in x-coordinate
	double yCen; ///< disk center in y-coordinate
};











/**
 * \class GeometryFactory
 * 
 * \brief The abstract factory class for creating objects in the Geometry family 
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
class GeometryFactory {
public:
	/**
     * \brief Defines a function pointer pointing to a function which creates objects in the Geometry family 
	 */
	typedef Geometry* (*GeoCreateFunc)();
	
	/**
     * \brief   Returns reference to a Singleton object of class GeometryFactory
	 * 
	 * \param   None 
     *
	 * \return  Reference to a static object of class GeometryFactory 
	 * 
	 *
	 * Example usage: 
	 * \code
	 *          GeometryFactory& factory = GeometryFactory::instance();
	 * \endcode
	 *
	 * \note    This function is implemented based on the lazy Singleton design pattern; 
	 *          therefore, only one GeometryFactory instance is allowed in each program
	 */
	static GeometryFactory& instance(); 
	
	/**
     * \brief      Registers (links) the geometry name \e name 
	 *		       with the function \e func for creating objects in the Geometry family
	 *			   
	 *             After registration, \e name can be used as an argument in the createGeometry member function
	 *             for creating objects of the linked type in the Geometry family
	 *	           
	 *  
	 * \param [in] name the geometry name 
	 * \param [in] func the function pointer pointing to the function that creates objects of a specific type
	 *             in the Geometry family
	 * 
	 * \return     None  
	 *
	 * \note       Instead of using this function directly, consider using the GeometryRegistrar class for
	 *		       the purpose of linking a geometry name and a specific class in the Geometry family. 
	 *             The function is kept public in case one wants to use it directly
	 *		  
	 *
	 */
	void registerGeometry(std::string name, GeoCreateFunc func);
	
	/**
     * \brief      This function creates an object of the class linked to the \name  
	 *
	 * \param [in] name the name linked to a specific class in the Geometry family via the 
	 *				    registrerGeometry member function
	 *  
	 * \return     A Geometry * pointer pointing to an object of a specific class in the Geometry family 
	 * 
	 * Example usage: 
	 * \code
	 *            GeometryFactory& factory = GeometryFactory::instance();
	 *            Geometry* newGeometry = factory.createGeometry(name);
	 * \endcode
	 *
	 */
	Geometry* createGeometry(std::string name);
private:
	std::unordered_map<std::string,GeoCreateFunc> geoTable; ///< hash table for the (name,creatFunction) pair
	GeometryFactory() {} ///< for singleton design pattern
	GeometryFactory(const GeometryFactory& other); ///< for singleton design pattern (Don't implement)
	GeometryFactory& operator=(const GeometryFactory& other); ///< for singleton design pattern (Don't implement)
};


#endif // __GEOMETRY_H__
