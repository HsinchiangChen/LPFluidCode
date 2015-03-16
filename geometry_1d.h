/**
 * \file   geometry_1d.h
 *
 * \brief  This header file contains classes for the initialization of 1d fluid objects
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

#ifndef __GEOMETRY_1D_H__
#define __GEOMETRY_1D_H__

#include "geometry.h"

/**
 * \class Line
 * 
 * \brief Supply functions for generating a 1D line geometry
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
class Line: public Geometry {
public:
	/** 	
	 * \brief Constructor
	 *
	 * This line object is centered at zero with length on the left being \e leftLen and 
	 * that on the right being \e rightLen
	 */
	Line();
	
	/// Desturctor
	virtual ~Line() {}
	
	/**
	 * \brief Level set function of a 1D line with center at zero
	 *
	 * The level set function is 
	 * \f$
	 *    -\mbox{rightLen} \leq x \leq \mbox{leftLen}    
	 * \f$
	 *
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return     \c true if (x,y,z) is inside the level set function; \c false otherwise
	 *  
	 */
	virtual bool operator()(double x, double y, double z) const;	
	
	/**
	 * \brief       Calculates the bounding box of 1D line centered at zero 
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
	double leftLen; ///< length on the left of zero
	double rightLen; ///< length on the right of zero
};



#endif //__GEOMETRY_1D_H__ 
