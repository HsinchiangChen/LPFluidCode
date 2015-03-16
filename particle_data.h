/**
 * \file   particle_data.h
 *
 * \brief  This header file contains classes for storing information and data of particles, 
 *         such as x, y, and z coordinates, neighbour lists, and bounding boxes   
 *      
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on interface and data structure design 
 *                           
 *
 * \version 1.0 
 *
 * \date 2014/7/12
 *
 * Created on: 2014/5/23 
 *
 */

#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include <cstddef>
#include <vector>
#include <string>

class LPSolver;
class Initializer;
class BoundingBox;



/**
 * \class ParticleData
 * 
 * \brief A class that stores all information and data of particles, 
 *        such as x, y, and z coordinates, neighbour lists, and bounding boxes
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface and data structure design
 *
 * \version 1.0 
 *
 * \date 2014/07/12 
 *
 * Created on: 2014/05/23 
 *
 */
class ParticleData {
public:	
	/**
	 * \brief             Constructor
	 * 
	 * Use \e init to retrieve information from the Initializer class, including getting scalar-valued parameters 
	 * and ownership of data arrays. This method also allocates space for all other necessary particle data arrays 
	 *
	 * \param [in] init   To retrieve information from \e init   
	 * \warning           Obtains ownership of some data arrays (see implementation for details) 
	 *                    from the Initializer class. These data arrays are deleted in the Destructor 
	 *                    of this class and should not be deleted anywhere else
	 */
	ParticleData(const Initializer& init);
	
	/**
	 * \brief Destructor 
	 *
	 * Deletes memory of data arrays
	 * 
	 */
	virtual ~ParticleData();

	/**
	 * \brief   Getter function of the specified array length for particle data arrays 
	 *          such as x, y, and z-coordinates, pressure, volume, and sound speeds    
	 * \param   None
	 * \return  The specified array length for particle data arrays 
	 *          such as x, y, and z-coordinates, pressure, volume, and sound speeds 
	 * \note    Capacity is a number larger than the total number of \e all types of particles
	 */
	std::size_t getCapacity() {return m_iCapacity;}

	/**
	 * \brief   Getter function of the number of \e all types of particles  
	 * \param   None
	 * \return  The number of \e all types of particles 
	 */
	std::size_t getTotalNum() {return m_iTotalNum;}
	
	/**
	 * \brief   Getter function of the number of \e fluid particles  
	 * \param   None
	 * \return  The number of \e fluid particles 
	 */
	std::size_t getFluidNum() {return m_iFluidNum;}
	
	/**
	 * \brief   Getter function of the number of \e boundary particles  
	 * \param   None
	 * \return  The number of \e boundary particles 
	 */
	std::size_t getBoundaryNum() {return m_iBoundaryNum;}
	
	/**
	 * \brief   Getter function of the number of \e ghost particles  
	 * \param   None
	 * \return  The number of \e ghost particles 
	 */
	std::size_t getGhostNum() {return m_iGhostNum;}
	
	/**
	 * \brief   Getter function of the specified maximum number of neighbours of a particle  
	 * \param   None
	 * \return  The specified maximum number of neighbours of a particle 
	 */
	std::size_t getMaxNeighbourNum() {return m_iMaxNeighbourNum;}
	
	/**
	 * \brief   Getter function of the specified maximum number of neighbours in one direction 
	 *          (eg. right-hand-side in the x-coordinate) of a particle  
	 * \param   None
	 * \return  The specified maximum number of neighbours in one direction 
	 *          (eg. right-hand-side in the x-coordinate) of a particle
	 */
	std::size_t getMaxNeighbourNumInOneDir() {return m_iMaxNeighbourNumInOneDir;}

	/**
	 * \brief   Getter function of the start index of fluid particles in the particle data arrays  
	 * \param   None
	 * \return  The start index of fluid particles in the particle data arrays
	 */
	std::size_t getFluidStartIndex() {return m_iFluidStartIndex;}
	
	/**
	 * \brief   Getter function of the start index of boundary particles in the particle data arrays  
	 * \param   None
	 * \return  The start index of boundary particles in the particle data arrays
	 */
	std::size_t getBoundaryStartIndex() {return m_iBoundaryStartIndex;}
	
	/**
	 * \brief   Getter function of the start index of ghost particles in the particle data arrays  
	 * \param   None
	 * \return  The start index of ghost particles in the particle data arrays
	 */
	std::size_t getGhostStartIndex() {return m_iGhostStartIndex;} 
	
	/**
	 * \brief   Getter function of the specified dimension of simulation 
	 * \param   None
	 * \return  The specified dimension of simulation
	 */
	int getDimension() {return m_iDimension;}
	
	
	/**
	 * \brief   Getter function of the array of x-coordinates of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of x-coordinates of the initialized particles
	 * 
	 */ 
	const double* const getPositionX() const {return m_vPositionX;}
	
	/**
	 * \brief   Getter function of the array of y-coordinates of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of y-coordinates of the initialized particles
	 * 
	 */
	const double* const getPositionY() const {return m_vPositionY;}
	
	/**
	 * \brief   Getter function of the array of z-coordinates of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of z-coordinates of the initialized particles
	 * 
	 */
	const double* const getPositionZ() const {return m_vPositionZ;}
	
	/**
	 * \brief   Getter function of the array of velocity in x-coordinate of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of velocity in x-coordinate of the initialized particles
	 *
	 */
	const double* const getVelocityU() const {return m_vVelocityU;}
	
	/**
	 * \brief   Getter function of the array of velocity in y-coordinate of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of velocity in y-coordinate of the initialized particles
	 *
	 */
	const double* const getVelocityV() const {return m_vVelocityV;}
	
	/**
	 * \brief   Getter function of the array of velocity in z-coordinate of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of velocity in z-coordinate of the initialized particles
	 *
	 */
	const double* const getVelocityW() const {return m_vVelocityW;}
	
	/**
	 * \brief   Getter function of the array of volume of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of volume of the initialized particles
	 *
	 */
	const double* const getVolume() const {return m_vVolume;}
	
	/**
	 * \brief   Getter function of the array of pressure of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of pressure of the initialized particles
	 *
	 */
	const double* const getPressure() const {return m_vPressure;}	
	
	/**
	 * \brief   Getter function of the array of sound speed of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of sound speed of the initialized particles
	 *
	 */
	const double* const getSoundSpeed() const {return m_vSoundSpeed;}	
	
	
	//const double* const getEnergy() const {return m_vEnergy;}

	
	/**
	 * \brief   Getter function of the array of "object tags" of the initialized particles  
	 * \param   None
	 * \return  Pointer pointing to the array of "object tags" of the initialized particles
	 * \note    The tag of a particle in a fluid object is the index of the fluid object 
	 *          (the index of fluid object is 1,2,...). 
	 *          The tag of a particle \b not in fluid object is 0
	 */
	const int* const getObjectTag() const {return m_vObjectTag;}
	

	/**
	 * \brief   Getter function of the array of divided difference of order 1   
	 * \param   None
	 * \return  Pointer pointing to the array of divided difference of order 1
	 *
	 */
	const double* const getDD1() const {return m_vDD1;}
	
	/**
	 * \brief   Getter function of the array of divided difference of order 2   
	 * \param   None
	 * \return  Pointer pointing to the array of divided difference of order 2
	 *
	 */
	const double* const getDD2Left() const {return m_vDD2Left;}
	
	/**
	 * \brief   Getter function of the array of divided difference of order 2   
	 * \param   None
	 * \return  Pointer pointing to the array of divided difference of order 2
	 *
	 */
	const double* const getDD2Right() const {return m_vDD2Right;}
	
	/**
	 * \brief   Getter function of the array of divided difference of order 3   
	 * \param   None
	 * \return  Pointer pointing to the array of divided difference of order 3
	 *
	 */
	const double* const getDD3Left() const {return m_vDD3Left;}
	
	/**
	 * \brief   Getter function of the array of divided difference of order 3   
	 * \param   None
	 * \return  Pointer pointing to the array of divided difference of order 3
	 *
	 */
	const double* const getDD3Right() const {return m_vDD3Right;}


private:
	
	friend class HyperbolicLPSolver; ///< To facilitate access of data in this class for HyperbolicLPSolver
	friend class HyperbolicLPSolver1D; ///< To facilitate access of data in this class for HyperbolicLPSolver1D

	int m_iDimension; ///< dimension
	std::size_t m_iCapacity;///< Maximum length of particle arrays (> the total number of all particles)
	std::size_t m_iTotalNum;///< Number of all types of particles
	std::size_t m_iFluidNum;///< Number of fluid particles
	std::size_t m_iBoundaryNum;///< Number of boundary particles
	std::size_t m_iGhostNum;///< Number of ghost particles
	std::size_t m_iFluidStartIndex;///< Start index of fluid particles in the particle array
	std::size_t m_iBoundaryStartIndex;///< Start index of boundary particles in the particle array
	std::size_t m_iGhostStartIndex;///< Start index of ghost particles in the particle array 
	std::size_t m_iMaxNeighbourNum;///< maximum number of neighbours of a particle
	std::size_t m_iMaxNeighbourNumInOneDir;///< maximum number of neighbours of a particle in one direction

	double* m_vPositionX;///< x
	double* m_vPositionY;///< y
	double* m_vPositionZ;///< z
	double* m_vVelocityU;///< velocity in x-coordinate
	double* m_vVelocityV;///< velocity in y-coordinate
	double* m_vVelocityW;///< velocity in z-coordinate
	double* m_vVolume;///< volume
	double* m_vPressure;///< pressure	
	double* m_vSoundSpeed;///< sound speed
	//double* m_vEnergy;///< energy
	
	double* m_vTemp1VelocityU;///< Temp array 1 of velocity in x-coordinate (for Stang Splitting)
	double* m_vTemp1VelocityV;///< Temp array 1 of velocity in y-coordinate (for Stang Splitting)
	double* m_vTemp1VelocityW;///< Temp array 1 of velocity in z-coordinate (for Stang Splitting)
	double* m_vTemp1Volume;///< Temp array 1 of volume (for Stang Splitting)
	double* m_vTemp1Pressure;///< Temp array 1 of pressure (for Stang Splitting)
	double* m_vTemp1SoundSpeed;///< Temp array 1 of sound speed (for Stang Splitting)

	double* m_vTemp2VelocityU;///< Temp array 2 of velocity in x-coordinate (for Stang Splitting)
	double* m_vTemp2VelocityV;///< Temp array 2 of velocity in y-coordinate (for Stang Splitting)
	double* m_vTemp2VelocityW;///< Temp array 2 of velocity in z-coordinate (for Stang Splitting)
	double* m_vTemp2Volume;///< Temp array 2 of volume (for Stang Splitting)
	double* m_vTemp2Pressure;///< Temp array 2 of pressure (for Stang Splitting)
	double* m_vTemp2SoundSpeed;///< Temp array 2 of sound speed (for Stang Splitting)

	int* m_vLPFOrderRight;///< The order of local polynomial fitting in \e right direction (in x-coordinate)
	int* m_vLPFOrderLeft;///< The order of local polynomial fitting in \e left direction (in x-coordinate) 
	int* m_vLPFOrderNorth;///< The order of local polynomial fitting in \e north direction (in y-coordinate) 
	int* m_vLPFOrderSouth;///< The order of local polynomial fitting in \e south direction (in y-coordinate) 
	int* m_vLPFOrderUp;///< The order of local polynomial fitting in \e up direction (in z-coordinate) 
	int* m_vLPFOrderDown;///< The order of local polynomial fitting in \e down direction (in z-coordinate) 
	
	int* m_vNeighbourList;///< The entire neighbour list
	int* m_vNeighbourListRight;///< The one-sided (right) neighbour list (\f& x > x_0 \f$)
	int* m_vNeighbourListLeft;///< The one-sided (left) neighbour list (\f& x < x_0 \f$)  
	int* m_vNeighbourListNorth;///< The one-sided (north) neighbour list (\f& y > y_0 \f$) 
	int* m_vNeighbourListSouth;///< The one-sided (south) neighbour list (\f& y < y_0 \f$) 
	int* m_vNeighbourListUp;///< The one-sided (up) neighbour list (\f& z > z_0 \f$)    
	int* m_vNeighbourListDown;///< The one-sided (down) neighbour list (\f& z < z_0 \f$)  
	
	int* m_vNeighbourListSize;///< The size of the entire neighbour list
	int* m_vNeighbourListRightSize;///< The size of the one-sided (right) neighbour list (\f& x > x_0 \f$) 
	int* m_vNeighbourListLeftSize;///< The size of the one-sided (left) neighbour list (\f& x < x_0 \f$)  
	int* m_vNeighbourListNorthSize;///< The size of the one-sided (north) neighbour list (\f& y > y_0 \f$) 
	int* m_vNeighbourListSouthSize;///< The size of the one-sided (south) neighbour list (\f& y < y_0 \f$) 
	int* m_vNeighbourListUpSize;///< The size of the one-sided (up) neighbour list (\f& z > z_0 \f$)    
	int* m_vNeighbourListDownSize;///< The size of the one-sided (down) neighbour list (\f& z < z_0 \f$)  

	 
	int* m_vObjectTag;///< tag=1,2,3,...:fluid objects; otherwise: boundary or ghost particles

	std::vector<BoundingBox*> m_vFluidBoundingBox;///< The vector containing all bounding box of fluid objects	
	std::vector<std::string> m_vBoundaryObjTypes; ///< The vector containing all boundary particle types

	double* m_vDD1; ///< divided difference order 1 (for 1D limiter)
	double* m_vDD2Left;///< divided difference order 2 (for 1D limiter)
	double* m_vDD2Right;///< divided difference order 2 (for 1D limiter)
	double* m_vDD3Left;///< divided difference order 3 (for 1D limiter)
	double* m_vDD3Right;///< divided difference order 3 (for 1D limiter)
	double* m_vCumP;///< cumulative pressure (for 1D limiter)
	double* m_vPositionXm;/// middle of positionX (for 1D limiter)

};


#endif
