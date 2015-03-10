/**
 * \file   lp_solver.h
 *
 * \brief  This header file contains classes of the main Lagrangian Particle solvers such as 
 *         the hyperbolic solvers and the elliptic solvers
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design 
 *            and the design of data pointer swaping algorithms in the Strang splitting method               
 *
 *
 * \version 1.0 
 *
 * \date 2014/10/09
 *
 * Created on: 2014/9/20 
 *
 */


#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__

#include <cstddef>
#include <vector>
#include <string>


class Initializer;
class ParticleData;
class NeighbourSearcher;
class EOS;


/**
 * \class LPSolver
 * 
 * \brief An abstract class for the family of Lagrangian Particle solvers
 *
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/10/09 
 *
 * Created on: 2014/09/20 
 *
 */
class LPSolver {

public:
	/// Destructor
	virtual ~LPSolver() {}

	/**
	 * \brief         The black box main Lagrangian particle solver for one iteration step
	 * 
	 * The method should be called in the context of TimeController
	 *
	 * \param [in] dt The length of physical time for this iteration
	 * \return        0 if the iteration is success 
	 * \warning       The function should always return 0 because all exceptions should be handled inside this class
	 */ 
	virtual int solve(double dt) = 0;
	
	/**
	 * \brief   Getter function of the minimum inter-particle distance among all fluid particles 
	 * \param   None
	 * \return  The minimum inter-particle distance among all fluid particles
	 */		
	virtual double getMinParticleSpacing() const {return m_fMinParticleSpacing;}
	
	/**
	 * \brief   Getter function of the maximum sound speed among all fluid particles 
	 * \param   None
	 * \return  The maximum sound speed among all fluid particles
	 */	
	virtual double getMaxSoundSpeed() const {return m_fMaxSoundSpeed;}

	/**
	 * \brief   Getter function of the maximum absolute value velocity among all fluid particles 
	 * \param   None
	 * \return  The maximum absolute value velocity among all fluid particles
	 */	
	virtual double getMaxFluidVelocity() const {return m_fMaxFluidVelocity;}

protected:

	double m_fMinParticleSpacing; ///< Minimum inter-particle spacing among fluid particles at a time step		
	double m_fMaxSoundSpeed; ///< Maximum sound speed of fluid particles at a time step	
	double m_fMaxFluidVelocity; ///< Maximum absolute value velocity of fluid particles at a time step
};


/**
 * \class HyperbolicLPSolver
 * 
 * \brief The default Lagrangian Particle solver for the compressible Euler's equation
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *            and the design of data pointer swaping algorithms in the Strang splitting method
 *
 * \version 1.0 
 *
 * \date 2014/10/09 
 *
 * Created on: 2014/09/20 
 *
 */
class HyperbolicLPSolver : public LPSolver {

public:
	/**
	 * \brief       Constructor
	 * 
	 * Get and set up parameters and obtain access to objects needed for the main solver 
	 *
	 * \param [in] init   To retrieve information from \e init   
	 * \param [in] pData  To obtain access to an object of the PaticleData clas
	 * \param [in] ns     To obtain access to an object in the NeighbourSearcher class
	 */
	HyperbolicLPSolver(const Initializer& init, ParticleData* pData, NeighbourSearcher* ns);
	
	/**
	 * \brief         The Lagrangian particle solver for the compressible Euler's equations for one iteration step
	 * 
	 * The method should be called in the context of TimeController
	 *
	 * \param [in] dt The length of physical time for this iteration
	 * \return        0 if the iteration is success 
	 * \warning       The function should always return 0 because all exceptions should be handled inside this class
	 */
	virtual int solve(double dt);	

private:

	//-----------------------------------------Data----------------------------------------

	//--------------------------Info got from input argument list---------------------------

	ParticleData* m_pParticleData; ///< Pointer to the object containing major particle data arrays 	
	NeighbourSearcher* m_pNeighbourSearcher; ///< Pointer to the object for the neighbour search task
	
	//--------------------------------------------------------------------------------------	

	//--------------------------Info get from Initializer class------------------------------
	
	EOS* m_pEOS; ///< Pointer to the object for computing eos
	int m_iNumThreads; ///< Number of threads	
	bool m_iIfMultiThreads;	///< true if use multithreads
	int m_iDimension; ///< dimension
	bool m_iRandomDirSplitOrder; ///< if true then the order of directional splitting is randomly set
	int m_iLPFOrder; ///< the order of Local Polynomial Fitting (LPF)
	std::size_t m_iNumRow2ndOrder; ///< the smallest number of rows of A to solve 2nd order LPF
	std::size_t m_iNumRow1stOrder; ///< the smallest number of rows of A to solve 1st order LPF
	std::size_t m_iNumCol2ndOrder; ///< the number of columns of A when solving 2nd order LPF
	std::size_t m_iNumCol1stOrder; ///< the number of columns of A when solving 1st order LPF	
	bool m_iMovingBoxForGhostParticle; ///< if true then the fluid box will be updated
	double m_fInitParticleSpacing; ///< the initial particle spacing
	double m_fGravity; ///< Gravity 
	double m_fInvalidPressure; ///< if p < invalid pressure => invalid state
	double m_fInvalidVolume; ///< volume cannot be negative: if volume < invalid volume => invalid state
	double m_fNeiSearchRadius; ///< the radius for neighbour search
	size_t m_iNumParticleWithinSearchRadius; ///< the number of particles within the search radius at initialization time
	double m_fContactLength; ///< defined length such that for two fluid particles from different fluid object, if the distance from each other is shorter than the length the two fluid particles start to interact with each other
	
	//--------------------------------------------------------------------------------------
	
	//---------------------------------Other parameters-------------------------------------

	int m_iNumPhase; ///< number of phases in directional splitting
	/**
	 *\brief 2D: A 2X3 table which maps direction split order and split phase to 0(x) or 1(y)\n
		     3D: A 6X5 table which maps direction split order and split phase to 0(x), 1(y), or 2(z)
	*/
	std::vector<std::vector<int> > m_vDirSplitTable; 

	int m_iDirSplitOrder;///< In 3D: 0=xyzyx, 1=xzyzx, 2=yxzxy, 3=yzxzy, 4=zxyxz, 5=zyx. In 2D: 0=xyx, 1=yxy	
	
	double m_fDt; ///< the time length of this iteration 
	
	int m_iContactAlert; ///< If true two or more fluid objects may soon in contact with each other; the value is determined by whether the bounding boxes of these fluid objects are overlap or not

	//-------------------------------------------------------------------------------------




	//-------------------------------------Methods-----------------------------------------
	
	/**  
	 * \brief A composite function that calls a bunch of other methods in order to set up
	 *        the environment for the next iteration based on the update of this iteration
	 *
	 * This function calls the following methods (names only for clarity)\n
	 * 
	 * updateFluidBoundingBox();\n
	 * checkForContactAlert();\n
	 * generateGhostParticle();\n 
	 * searchNeighbourForAllParticle();\n 	
	 * setUpwindNeighbourList();\n
	 * resetLPFOrder();\n
	 * computeMinParticleSpacing();\n
	 * computeMaxSoundSpeed();\n
	 * computeMaxFluidVelocity();\n
	 * setBoundaryPressureAndVelocity(int phase);
	 * setGhostPressureAndVelocity(int phase);
	 */	
	void computeSetupsForNextIteration(); 
	
		
	/**  
	 * \brief Update the bounding boxes of fluid objects based on the updated particle location  
	 *
	 *
	 *
	 */
	void updateFluidBoundingBox();
	
	/**  
	 * \brief Checks if two fluid objects in the simulation is in contact or not 
	 *
	 * The check is based on if the bounding box of the two fluid objects overlap or not 
	 *
	 *
	 */
	void checkForContactAlert();
	
	/**  
	 * \brief Generates ghost particles  
	 *  
	 *
	 */
	void generateGhostParticle();
	
	/**  
	 * \brief Determines if a ghost candidate particle is added as a ghost particle based on some criteria
	 *
	 *  
	 *
	 */
	bool isValidGhostParticle(double x, double y, double z, int* neiList, size_t numNei, int objectTag); 
	
	/**  
	 * \brief Searches neighbours for all particles including fluid, boundary, and ghost particles 
	 *        based on octree neighobur search 
	 *
	 *  
	 *
	 */
	void searchNeighbourForAllParticle(); 	
	
	/**  
	 * \brief Set up one-sided neighbour list based on the entire neighbour list (for fluid particles only) 
	 *
	 *  
	 *
	 */
	void setUpwindNeighbourList();
	
	/**  
	* \brief Reset the order of local polynomial fitting to the pre-set value \e m_iLPFOrder  
	*
	* The order of local polynomial fitting are multiple arrays, each represent a direction 
	* (eg. two arrays for right and left in the x-coordinate) 
	*
	*/	
	void resetLPFOrder();
	
	/**  
	 * \brief Compute the minimum inter-particle spacing among fluid particles
	 *
	 * For the computation of next dt  
	 *
	 */
	void computeMinParticleSpacing();
	
	/**  
	 * \brief Computes the maximum sound speed of fluid particles
	 *
	 * For the computation of next dt 
	 *
	 */
	void computeMaxSoundSpeed();
	
	/**  
	 * \brief Computes the maximum absolute value velocity of fluid particles
	 *
	 * For the computation of next dt 
	 *
	 */
	void computeMaxFluidVelocity();
	
	
	/**  
	 * \brief Changes m_iObjectTag of a fluid object if it is in contact with particles in other fluid object 
	 *
	 * If a fluid particle have a fluid neighbour from other fluid object within m_fContactLength
     * then set the m_vObjectTag of this particle to be the negative of original
	 *
	 */
	void changeFluidObjectTagIfInContact(int index, size_t numNeiFound, const double* neiListDist);
	
	/**  
	 * \brief Change the neighbour list of a fluid particle 
	 * 
	 * Change the neighbour list to only include
	 * 1. non-fluid particles and\n 
	 * 2. fluid particles in the same fluid object  
	 *
	 * \note This is done only when this fluid particle is \e not in contact with other fluid particles
	 *
	 */
	void changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(int index, size_t numNeiFound);
	
	/**  
	 * \brief Changes the neighbour list to include only fluid particles 
	 *
	 * This function is now only used by boundary and ghost particles
	 *
	 */
	void changeNeighbourhoodToIncludeOnlyFluidNei(int index, size_t numNeiFound);

	
	/**  
	 * \brief Calculates the spatial derivatives and performs time integration by the Strang splitting method
	 *
	 * Calls the following functions (names only for clarity):\n
	 * setNeighbourListPointers();\n
     * setInAndOutDataPointers();\n
	 * setLPFOrderPointers();\n
     * computeSpatialDer();\n 
	 * timeIntegration();\n 
	 * printInvalidState();\n
	 * lowerLPFOrder();
	 *
	 */
	bool directionalSplitting(int phase);
	
	/**  
	 * \brief Alias two one-sided neighbour lists based on if this round is on x, y, or z 
	 *
	 *
	 */	
	void setNeighbourListPointers(int dir, // input
		const int **neighbourList0, const int **neighbourList1, // output
		const int **neighbourListSize0, const int **neighbourListSize1);
	
	/**  
	 * \brief Alias input and output data pointers for this round 
	 *
	 * \note There are 3 rounds in total for 2D and 5 rounds in toal for 3D
	 */
	void setInAndOutDataPointers(int phase, int dir,
		const double** inVelocity, const double** inPressure, const double** inVolume, const double** inSoundSpeed, 
		double** outVelocity, double** outPressure, double** outVolume, double** outSoundSpeed);
		
	/**  
	 * \brief Alias the pointers to the arrays of the order of local polynomial fitting
	 *
	 *
	 */
	void setLPFOrderPointers(int dir, // input
		int** LPFOrder0, int** LPFOrder1, std::vector<int*>& LPFOrderOther); // output

	
	/**  
	 * \brief Compute the spatial derivatives by solving least squares problem for this round
	 *
	 *
	 */
	void computeSpatialDer(int dir, size_t index, // input
		int offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*),
		const double* inPressure, const double* inVelocity,
		const int *neighbourList, const int *neighbourListSize, 
		int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd); // output
	
	
	/**  
	 * \brief Performs time integration for this round
	 *
	 * \note There are 3 rounds in total for 2D and 5 rounds in toal for 3D
	 */
	void timeIntegration(double real_dt, double multiplier1st, double multiplier2nd, 
		double gravity, double inVolume, double inVelocity, double inPressure,
		double inSoundSpeed, 
		double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
		double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
		double* outVolume, double* outVelocity, double* outPressure); // output
	
	
	/**  
	 * \brief Print out info about the particle which evolve into invalid states 
	 *
	 *
	 */
	void printInvalidState(int phase, int dir,  int index, double positionX, double positionY, double positionZ,
		double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0, 
		double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1);
	
	/*!
	\brief lower the order the LPF in a direction if get invalid state
		   then see if all LPFOrder are zero, if it is the case, return false and go back to phase 0
		   otherwise, return true and go to next phase
	*/
	
	/**  
	 * \brief Lower the order of local polynomial fitting for a particle in a specific direction (x,y, or z)
	 *
	 * \return \c true if at least one of the LPFOrder array for this particle in a direction is not zero;
	 *         \c false otherwise
	 *
	 * \note   If return \c then we will redo this particle with the lowered order of local polynomial fitting;
	 *         if return \c false we will go back to phase for all particles, with this particle basically not 
	 *         updated for the following one entirew iteration
	 *
	 */
	bool lowerLPFOrder(int index, const std::vector<int*>& LPFOrderOther, // input
		int* LPFOrder0, int* LPFOrder1); // output	
	
	/**  
	 * \brief Computes the number of rows and columns used for a particle at this round 
	 *
	 * \note A helper function of computeSpatialDer()
	 *
	 */
	void computeNumRowAndNumColAndLPFOrder(size_t index, // input
		const int *neighbourList, const int *neighbourListSize, size_t numRow2nd, size_t numRow1st,
		int* LPFOrder, size_t *numRow, size_t *numCol); // output


	/**  
	 * \brief Computes the matrix A in the least squares problem Ax~b in the 2D context 
	 *
	 * \note A helper function of computeSpatialDer()
	 */
	void computeA2D(size_t index, const int *neighbourList, const int* LPFOrder, size_t numRow, size_t numCol, // input
					double *A); // output 
	
	/**  
	 * \brief Computes the matrix A in the least squares problem Ax~b in the 3D context
	 *
	 * \note A helper function of computeSpatialDer()
	 */
	void computeA3D(size_t index, const int *neighbourList, const int* LPFOrder, size_t numRow, size_t numCol, // input
					double *A); // output
	
	/**  
	 * \brief Computes the vector b in the least squares problem Ax~b
	 *
	 * \note A helper function of computeSpatialDer()
	 */
	void computeB(size_t index, const int *neighbourList, size_t numRow, const double* inData, 
				  double *b); // output
	
	
	/**  
	 * \brief Compute the pressure and velocity by 0-th order local polynomial fitting for boundary particles
	 *
	 * \note calls the function computeOthOrderWeightedLPF()
	 */
	void setBoundaryPressureAndVelocity(int phase);
	
	/**  
	 * \brief Compute the pressure and velocity by 0-th order local polynomial fitting for ghost particles
	 *
	 * \note calls the function computeOthOrderWeightedLPF()
	 */
	void setGhostPressureAndVelocity(int phase); 
	
	/**  
	 * \brief Computes the 0-th order local polynomial fitting
	 *
	 *  
	 */
	void computeOthOrderWeightedLPF(std::vector<const double*>& position, 
									std::size_t startIndex, std::size_t numParticle,
									std::vector<double*>& data);
	
	/**  
	 * \brief Updates the states of fluid particles at the end of one iteration by swapping pointers 
	 *
	 *
	 */
	void updateFluidState();
	
	/**  
	 * \brief Updates the location of fluid particles based on velocities at the end of one iteration
	 *
	 * Based on a combination of forward and backward Euler's method
	 */
	void moveFluidParticle();
	
	/**  
	 * \brief Update the velocities of fluid particles at the end of one iteration by swapping pointers
	 *
	 * 
	 */
	void updateFluidVelocity();


	/**  
	 * \brief Tests on the neighbour search methods 
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void testNeighbourSearch();
	
	/**  
	 * \brief Use brute force neighbour search to search neighbours for fluid particles
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void searchNeighbourBruteForce(int, int* neighbourList, int* neighbourListSize); 
	
	/**  
	 * \brief Use brute force neighbour search to search neighbours for ghost particles
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void searchNeighbourBruteForce(double x, double y, double z, int* neighbourList, size_t& numNeiFound); 
	
	/**  
	 * \brief Generate ghost particles by the brute force neighbour search
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void generateGhostParticleByBruteForceNeighbourSearch();
	
	/**  
	 * \brief A wrapper function which uses brute force neighbour search to search neighbours for all particles 
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void searchNeighbourForAllParticleByBruteForceNeighbourSearch();
	
	/**  
	 * \brief Check the validity of all one-sided neighbour lists
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	bool checkUpwindNeighbourList();
	
	/**  
	 * \brief helper function of writeResult()
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	std::string rightFlush(size_t writeStep, size_t numDigits);	
	
	/**  
	 * \brief Writes results for debugging purposes
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	int writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle);
};



#endif // __LP_SOLVER_H__
