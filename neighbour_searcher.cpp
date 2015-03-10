/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Oct. 15 2014
 * 
 * \brief
 *      
 */


#include <algorithm>
#include <iostream>
#include "neighbour_searcher.h"
#include "initializer.h"
#include <omp.h>


////////////////////////////////////////////////////////////////////////////////////////
// Start of OctreeSearcher
////////////////////////////////////////////////////////////////////////////////////////

OctreeSearcher::OctreeSearcher(const Initializer& init) {
	m_iMaxParticleNum = init.getCapacity();
	m_iMaxNeighborNum = init.getMaxNeighbourNum();
	int treeDepth = init.getTreeDepth(); 

    //m_iMaxParticleNum = 1000000;
	//m_iMaxNeighborNum = 1000;
	//int treeDepth = 6;

	m_pOctree = new Octree(treeDepth, m_iMaxParticleNum);
	//m_pSearchResult = new SearchResult[m_iMaxNeighborNum];
	
#ifdef _OPENMP
	m_iTheadNum = omp_get_max_threads();
	
	m_pSearchResult = new SearchResult*[m_iTheadNum];
	for(int i=0 ; i<m_iTheadNum ; i++) {
		m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
	}
#else
	m_pSearchResult = new SearchResult[m_iMaxNeighborNum];
#endif


	// print debug info
	std::cout<<"-------OctreeSearcher::OctreeSearcher()-------"<<std::endl;
	std::cout<<"m_iMaxParticleNum = "<<m_iMaxParticleNum<<std::endl;
	std::cout<<"m_iMaxNeighborNum = "<<m_iMaxNeighborNum<<std::endl;
	std::cout<<"treeDepth = "<<treeDepth<<std::endl;
	std::cout<<"----------------------------------------------"<<std::endl;
}



OctreeSearcher::OctreeSearcher(size_t maxParticleNum, size_t maxNeighbourNum, int treedepth)
: m_iMaxParticleNum(maxParticleNum), m_iMaxNeighborNum(maxNeighbourNum)
{
  m_pOctree = new Octree(treedepth, m_iMaxParticleNum);

#ifdef _OPENMP
  m_iTheadNum = omp_get_max_threads();
  /*
  m_pSearchResult.resize(m_iTheadNum);
  for(int i=0 ; i<m_iTheadNum ; i++) {
    //m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
    m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
  }
  */
  m_pSearchResult = new SearchResult*[m_iTheadNum];
  for(int i=0 ; i<m_iTheadNum ; i++) {
    m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
  }
#else
  m_pSearchResult = new SearchResult[m_iMaxNeighborNum];
#endif
}



int OctreeSearcher::buildSearchStructure(const double* x, const double* y, const double* z, size_t begin, size_t numParticles) {
  return m_pOctree->buildOctree(x+begin, y+begin, z+begin, numParticles);
}

/*
int OctreeSearcher::searchNeighbour(const double x, const double y, const double z, const double radius, int* result, size_t& result_length) {
  return m_pOctree->searchNeighbor(x, y, z, radius, result, result_length);
}
*/


#ifdef _OPENMP
int OctreeSearcher::searchNeighbour(const double x, const double y, const double z, const double radius, int* result, double* distance, size_t& result_length, int tid, int index) {

  m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult[tid], result_length);

  sort(m_pSearchResult[tid], m_pSearchResult[tid]+result_length);
  
  result_length = min(result_length, m_iMaxNeighborNum);

  if(index > -1) {
    size_t i=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[tid][i].index) {
        i++;
        break;
      } else {
        result[i] = m_pSearchResult[tid][i].index;
        distance[i] = m_pSearchResult[tid][i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[tid][i].index;
      distance[i-1] = m_pSearchResult[tid][i].distance;
    }

    result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[tid][i].index;
      distance[i] = m_pSearchResult[tid][i].distance;
    }
  }



  return 0;
}

#else
int OctreeSearcher::searchNeighbour(const double x, const double y, const double z, const double radius, int* result, double* distance, size_t& result_length, int index) {

  m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult, result_length);


  sort(m_pSearchResult, m_pSearchResult+result_length);
  
  result_length = min(result_length, m_iMaxNeighborNum);

  if(index > -1) {
    size_t i=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[i].index) {
        i++;
        break;
      } else {
        result[i] = m_pSearchResult[i].index;
        distance[i] = m_pSearchResult[i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[i].index;
      distance[i-1] = m_pSearchResult[i].distance;
    }

    result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[i].index;
      distance[i] = m_pSearchResult[i].distance;
    }
  }



  return 0;
}

#endif

////////////////////////////////////////////////////////////////////////////////////////
// End of OctreeSearcher
////////////////////////////////////////////////////////////////////////////////////////

  



  
