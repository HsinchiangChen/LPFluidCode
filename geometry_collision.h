/**
 * \file   geometry_collision.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the geometry of fluid objects in the collision simulation
 *
 * This file serves as an example that different files (other than geometry.h and geometry.cpp) should be used to 
 * model the geometry of fluid objects in different simulations
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

#ifndef __GEOMETRY_COLLISION_H__
#define __GEOMETRY_COLLISION_H__

#include "geometry.h"



/**
 * \class DiskLeft
 * 
 * \brief Supply functions for generating a 2D disk geometry for the 2D collision simulation
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
class DiskLeft: public Geometry {
public:
	/// constructor
	DiskLeft();

	/// destructor
	virtual ~DiskLeft() {}
	
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
	double radius;
	double xCen;
	double yCen;
};










/**
 * \class DiskRight
 * 
 * \brief Supply functions for generating a 2D disk geometry for the 2D collision simulation
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
class DiskRight: public Geometry {
public:
	/// constructor
	DiskRight();
	
	/// destructor
	virtual ~DiskRight() {}
	
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
	double radius;
	double xCen;
	double yCen;
};


#endif // __GEOMETRY_COLLISION_H__
