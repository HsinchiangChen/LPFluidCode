/**
 * \file   neighbour_searcher.h
 *
 * \brief  This header file contains classes for searching nearest neighbours of particles 
 *         using different algorithms such as octrees   
 *      
 *
 * \author  Yu, Kwangmin (yukwangmin@gmail.com)
 *
 *  
 *                           
 *
 * \version 1.0 
 *
 * \date 2014/10/15
 *
 * Created on: 2014/8/22 
 *
 */


#ifndef __NEIGHBOUR_SEARCHER_H__
#define __NEIGHBOUR_SEARCHER_H__



#ifdef _MSC_VER  
typedef __int32 int32_t; 
typedef unsigned __int32 uint32_t; 
typedef __int64 int64_t; 
typedef unsigned __int64 uint64_t;  
#else 
#include <stdint.h> 
#endif 


#include <cstdlib>
#include <vector>
#include <deque>
#include "octree.h"


class Initializer;


/**
 * \class NeighbourSearcher
 * 
 * \brief An abstract class for the family of classes that performs the task of
 *        nearest neighbour search for particles
 *
 * \author  Yu, Kwangmin (yukwangmin@gmail.com)
 *
 *  
 *                           
 *
 * \version 1.0 
 *
 * \date 2014/10/14
 *
 * Created on: 2014/8/22 
 *
 */
class NeighbourSearcher {

public:
  
	/// Destructor
	virtual ~NeighbourSearcher() {} 

    /**
	 * \brief Build the data structure for the nearest neighbour search, based on 
	 *        the input particle location ayyays in the x, y, and z-coordinates and the number of particles 
	 * 
	 * The function uses the information of the input arrays in the range [0,\e numParticles)
	 *
	 * \param [in] x            The particle location array in the x-coordinate
	 * \param [in] y            The particle location array in the y-coordinate 
	 * \param [in] z            The particle location array in the z-coordinate
	 * \param [in] numParticles The number of particles used in the particle location arrays 
	 * \return                  0, if success        
	 */
	virtual int buildSearchStructure(const double* x, const double* y, const double* z, size_t numParticles) = 0;
	
	/**
	 * \brief Build the data structure for the nearest neighbour search, based on 
	 *        the input particle location arrays in the x, y, and z-coordinates, 
	 *        the starting point in the arrays and the number of particles 
	 * 
	 * The function uses the information of the input arrays in the range [\e begin,\e numParticles)
	 *
	 * \param [in] x            The particle location array in the x-coordinate
	 * \param [in] y            The particle location array in the y-coordinate 
	 * \param [in] z            The particle location array in the z-coordinate
	 * \param [in] begin        The starting point in the particle array to use to build the search structure
	 * \param [in] numParticles The number of particles used in the particle location arrays 
	 * \return                  0, if success        
	 */
	virtual int buildSearchStructure(const double* x, const double* y, const double* z, size_t begin, size_t numParticles) = 0;
 
#ifdef _OPENMP
	/**
	 * \brief Searches the nearest neighbours for a particle on the built search data structure  
	 *
	 * \param [in]  x             The location in the x-coordinate of the target particle
	 * \param [in]  y             The location in the y-coordinate of the target particle 
	 * \param [in]  z             The location in the z-coordinate of the target particle
	 * \param [in]  radius        The neighbour search radius 
	 * \param [out] result        An array saving the results of neighbour search (the neighbour list), 
	 *                            which is sorted based on the 
	 *                            the \e distance (in the \b ascending order of \e distance) 
	 * \param [out] distance      An array saving the distance between the target particle and its neighbours
	 * \param [out] result_length The number of neighbours obtained for the target particle
	 * \param [in]  tid           The thread id
	 * \param [in]  index         The index of target particle; this index is used to remove the target particle
	 *                            itself from its own neighbour list. The default argument is -1,
	 *                            which should only be used when the target particle itself is \b not
	 *                            in the search structure
	 * \return                    0, if success        
	 *
	 * \note This is the \b multithread version of this function  
	 *       
	 */	
	virtual int searchNeighbour(const double x, const double y, const double z, const double radius, int* result, double* distance, size_t& result_length, int tid, int index = -1) = 0;
#else
	/**
	 * \brief Searches the nearest neighbours for a particle on the built search data structure  
	 *
	 * \param [in]  x             The location in the x-coordinate of the target particle
	 * \param [in]  y             The location in the y-coordinate of the target particle 
	 * \param [in]  z             The location in the z-coordinate of the target particle
	 * \param [in]  radius        The neighbour search radius 
	 * \param [out] result        An array saving the results of neighbour search (the neighbour list), 
	 *                            which is sorted based on the 
	 *                            the \e distance (in the \b ascending order of \e distance) 
	 * \param [out] distance      An array saving the distance between the target particle and its neighbours
	 * \param [out] result_length The number of neighbours obtained for the target particle
	 * \param [in]  index         The index of target particle; this index is used to remove the target particle
	 *                            itself from its own neighbour list. The default argument is -1,
	 *                            which should only be used when the target particle itself is \b not
	 *                            in the search structure
	 * \return                    0, if success        
	 *
	 * \note This is the \b single-thread version of this function  
	 *       
	 */
	virtual int searchNeighbour(const double x, const double y, const double z, const double radius, int* result, double* distance, size_t& result_length, int index = -1) = 0;
#endif

//virtual int searchNeighbour(const double x, const double y, const double z, const double radius, int* result, size_t& result_length, int index = -1) = 0;

  //virtual int searchNeighbour(const double x, const double y, const double z, const double radius, int* result, size_t& result_length) = 0;


};

/**
 * \class OctreeSearcher
 * 
 * \brief A class that searches nearest neighbours of particles based on the Octree data structure and algorithm
 *
 * \author  Yu, Kwangmin (yukwangmin@gmail.com)
 *
 *  
 *                           
 *
 * \version 1.0 
 *
 * \date 2014/10/14
 *
 * Created on: 2014/8/22 
 *
 */
class OctreeSearcher : public NeighbourSearcher {
public:
	/**
	 * \brief Constructor
	 * 
	 * Obtain all necessary information for the Initializer object and allocate memory for 
	 * the Octree object and array of SearchResult
	 *
	 * \param [in] init An object of the Initializer class
	 *
	 *
	 */
	OctreeSearcher(const Initializer& init);
	

	/**
	 * \brief Constructor
	 * 
	 * Obtain scalar parameters and allocate memory for the Octree object and array of SearchResult
	 *
	 * \param [in] maxParticleNum    The capacity of the particle array (> the total number of particles)
	 * \param [in] maxNeighbourNum   The maximum number of neighbours of a particle
	 * \param [in] treedepth         The depth of octrees
	 *
	 */		
	OctreeSearcher(size_t maxParticleNum, size_t maxNeighbourNum, int treedepth);
	

	/**
	 * \brief Destructor
	 * 
	 * Delete the Octree object and the array of SearchResult
	 *
	 *
	 *
	 */
	virtual ~OctreeSearcher() {
		delete m_pOctree;
#ifdef _OPENMP
		for(int i=0 ; i<m_iTheadNum ; i++) {
			delete[] m_pSearchResult[i];
		}
		delete[] m_pSearchResult;
#else
		delete[] m_pSearchResult;
#endif
	}

public:
	
	/**
	 * \brief Build an Octree object for the nearest neighbour search, based on 
	 *        the input particle location ayyays in the x, y, and z-coordinates and the number of particles 
	 * 
	 * The function uses the information of the input arrays in the range [0,\e numParticles)
	 *
	 * \param [in] x            The particle location array in the x-coordinate
	 * \param [in] y            The particle location array in the y-coordinate 
	 * \param [in] z            The particle location array in the z-coordinate
	 * \param [in] numParticles The number of particles used in the particle location arrays 
	 * \return                  0, if success        
	 */
	virtual int buildSearchStructure(const double* x, const double* y, const double* z, size_t numParticles) {
		return buildSearchStructure(x, y, z, 0, numParticles);
	}
	

	/**
	 * \brief Build an Octree object for the nearest neighbour search, based on 
	 *        the input particle location arrays in the x, y, and z-coordinates, 
	 *        the starting point in the arrays and the number of particles 
	 * 
	 * The function uses the information of the input arrays in the range [\e begin,\e numParticles)
	 *
	 * \param [in] x            The particle location array in the x-coordinate
	 * \param [in] y            The particle location array in the y-coordinate 
	 * \param [in] z            The particle location array in the z-coordinate
	 * \param [in] begin        The starting point in the particle array to use to build the search structure
	 * \param [in] numParticles The number of particles used in the particle location arrays 
	 * \return                  0, if success        
	 */	
	virtual int buildSearchStructure(const double* x, const double* y, const double* z, size_t begin, size_t numParticles);
   
   
#ifdef _OPENMP
	/**
	 * \brief Searches the nearest neighbours for a particle on the built Octree object  
	 *
	 * \param [in]  x             The location in the x-coordinate of the target particle
	 * \param [in]  y             The location in the y-coordinate of the target particle 
	 * \param [in]  z             The location in the z-coordinate of the target particle
	 * \param [in]  radius        The neighbour search radius 
	 * \param [out] result        An array saving the results of neighbour search (the neighbour list), 
	 *                            which is sorted based on the 
	 *                            the \e distance (in the \b ascending order of \e distance) 
	 * \param [out] distance      An array saving the distance between the target particle and its neighbours
	 * \param [out] result_length The number of neighbours obtained for the target particle
	 * \param [in]  tid           The thread id
	 * \param [in]  index         The index of target particle; this index is used to remove the target particle
	 *                            itself from its own neighbour list. The default argument is -1,
	 *                            which should only be used when the target particle itself is \b not
	 *                            in the search structure
	 * \return                    0, if success        
	 *
	 * \note This is the \b multithread version of this function  
	 *       
	 */
	virtual int searchNeighbour(const double x, const double y, const double z, const double radius, int* result, double* distance, size_t& result_length, int tid, int index = -1);
#else
	/**
	 * \brief Searches the nearest neighbours for a particle on the built Octree object  
	 *
	 * \param [in]  x             The location in the x-coordinate of the target particle
	 * \param [in]  y             The location in the y-coordinate of the target particle 
	 * \param [in]  z             The location in the z-coordinate of the target particle
	 * \param [in]  radius        The neighbour search radius 
	 * \param [out] result        An array saving the results of neighbour search (the neighbour list), 
	 *                            which is sorted based on the 
	 *                            the \e distance (in the \b ascending order of \e distance) 
	 * \param [out] distance      An array saving the distance between the target particle and its neighbours
	 * \param [out] result_length The number of neighbours obtained for the target particle
	 * \param [in]  index         The index of target particle; this index is used to remove the target particle
	 *                            itself from its own neighbour list. The default argument is -1,
	 *                            which should only be used when the target particle itself is \b not
	 *                            in the search structure
	 * \return                    0, if success        
	 *
	 * \note This is the \b single-thread version of this function  
	 *       
	 */
	virtual int searchNeighbour(const double x, const double y, const double z, const double radius, int* result, double* distance, size_t& result_length, int index = -1);
#endif



private:
	Octree* m_pOctree;///< The octree
	size_t m_iMaxParticleNum;///< The capacity of the particle array (> the total number of particles)
	size_t m_iMaxNeighborNum;///< The maximum number of neighbours of a particle

#ifdef _OPENMP
	int m_iTheadNum;///< The number of threads 
	SearchResult** m_pSearchResult;///< The search result for each thread id  
	//vector<SearchResult*> m_pSearchResult;
#else
	SearchResult* m_pSearchResult;/// The search result (of length \e m_iMaxNeighborNum)
#endif

};

#endif
