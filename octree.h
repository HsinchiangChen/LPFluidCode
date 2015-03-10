/**
 * \file   octree.h
 *
 * \brief  This header file contains classes for the Octree data structure and nearest neighbour search algorithms   
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

#ifndef __OCTREE_H__
#define __OCTREE_H__


#ifdef _MSC_VER  
typedef __int32 int32_t; 
typedef unsigned __int32 uint32_t; 
typedef __int64 int64_t; 
typedef unsigned __int64 uint64_t;  
#else 
#include <stdint.h> 
#endif 

#include <deque>
#include <vector>

#ifdef _OPENMP
//#include <thrust/sort.h>
//#include <thrust/host_vector.h>
//#include <thrust/system/omp/vector.h>
#endif


//#include "particle_bunch.h"


using namespace std;

/**
 * \brief A simple struct for the neighbour search result
 *
 */
struct SearchResult {
  double distance; ///< The distance between the target particle and its neighbour
  int index; ///< The index of the neighbour particle
  /// Overloaded operator <
  bool operator < (const SearchResult& sr) const {
    return (distance < sr.distance);
  }
};





/*!
\brief This struct contains the morton key and the index of a particle. 
       a vector of such structs is then sorted by the "key" field. 
 */
struct KeyIndex
{
  int key;//Should be a 32/64 bit integer. With a 32 bit integer we can use Octrees of depth 10. With a 64 bit integer we can use octrees of depth 20.
  int index;
  bool operator < (const KeyIndex& a) const
  {
    return (key < a.key);
  }
};



/**
 * \brief This class contains the Octree data structure and algorithm for nearest neighbour search of particles
 *      
 * This class is very strongly related with its implementation and the data structure of DefaultParticles class. 
 * So do not use this class with the general purpose. This class does not have octree search algorithm. 
 * The search algorithm is implemented in DefaultParticles and its subclasses.
 * The member fields starting from "m_v" must have same length with DefaultParticles::m_iNumOfParticles 
 * and m_iTotalNumberOfParticles. 
 * That is, the length of these member field means the total number of particles.
 * Each value of the array represent octree nodes.
 * The location information can be fetched by decoding m_vKey value.
 * \see DefaultParticles
 */
class Octree {

private:

  size_t m_iMaxParticleNum;

  //vector<ParticleBunch*> *m_vParticleBunchList;

  /*!
  \brief array of x-coordinates of particles.
  */
  const double *m_vCoordX;
  /*!
  \brief array of y-coordinates of particles.
  */
  const double *m_vCoordY;
  /*!
  \brief array of z-coordinates of particles.
  */
  const double *m_vCoordZ;

  /*!
  \brief Total number of particles. The value of this field must be same with DefaultParticles::m_iNumOfParticles.
  */
  int m_iTotalNumberOfParticles;
 
  /*!
  \brief Array of size m_iTotalNumberOfParticles which contains the Morton key and Index of each particle. 
  */
#ifdef _OPENMP
  //thrust::omp::vector<KeyIndex> m_vParticleKeyIndex;
  //thrust::host_vector<KeyIndex> m_vParticleKeyIndex;
  //std::vector<KeyIndex> m_vParticleKeyIndex;
  KeyIndex *m_vParticleKeyIndex;
#else
  KeyIndex *m_vParticleKeyIndex;
#endif
                        
  /*!
  \brief The maximum depth of this octree. If this value is 1, then the octree hase one root node and 8 children nodes.
  \see m_vDepth
  */
  int m_iMaxDepth;

  /*!
  \brief This field is for search algorithm. Do not use this field except the implementor of octree search algorithm.
  */
  int m_iLinearSearchDepth;
  //int maxParticleNumInCell;

  /*!
  \brief Minimum X-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_min_x ;

  /*!
  \brief Maximum X-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_max_x ;
 
  /*!
  \brief Minimum Y-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_min_y ;

  /*!
  \brief Maximum Y-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_max_y ;

  /*!
  \brief Minimum Z-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_min_z ;

  /*!
  \brief Maximum Z-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_max_z ;


  /*!
  \brief This is the length of the octree. Specifically it is the length of the arrays listed below
  */
  int m_iTreeLength;
  


  /*!
  \brief Morton Key. 
         This key has octal base number representation.
         So the total length having meaning is 3 * (m_iMaxDepth - 1) bits.
         The first 3 bits represent the first region divided by 2.
         For example,
         000 means x <= 1/2, y <= 1/2, z <= 1/2
         001 means x <= 1/2, y <= 1/2, z > 1/2
         010 means x <= 1/2, y > 1/2, z <= 1/2
         ...
         111 means x > 1/2, y > 1/2, z > 1/2
         The second successive three bits have the same meaning above in the region the first three bits represent.
         So the root node of octree must have 0 key value.
         In the second depth of octree, all bits execpt the first three bits have 0 values.
  */
  int* m_vNodeKey;
  
  /*!
  \brief Represent the depth of the node which is designated by the array index.
         The value of the root node is 0.
  */
  int* m_vDepth;

  /*!
  \brief The first index of particles of the nodes. In combination with m_vNumberOfParticle, we can obtain all indexes of particles in the node which is designated by the array index.
             In order to travel all particles in the node which is designated by the index, reference from m_vFirstParticleIndex to m_vFirstParticleIndex + m_vNumberOfParticle - 1.
  \see m_vNumberOfParticle
  */
  int* m_vFirstParticleIndex;
  
  /*!
  \brief The total number of particles including in the node which is designated by the array index.
  \see m_vFirstParticleIndex
  */
  int* m_vNumberOfContainedParticles;

  /*!
  \brief The first index of child node. In combination with m_vNumberOfChild, we can obtain all indexes of children nodes in the node which is designated by the array index.
         In order to travel all nodes in the node which is designated by the index, reference from m_vFirstChildIndex to m_vFirstChildIndex + m_vNumberOfChild - 1.
         If the node designated tye the array index is a leaf node, this field is set by -1.
  \see m_vNumberOfChild
  */
  int* m_vFirstChildIndex;
  
  /*!
  \brief The total number of children nodes including in the node which is designated by the array index.
         If the node designated tye the array index is a leaf node, this field is set by 0.
  \see m_vFirstChildIndex
  */
  int* m_vNumberOfChildren;

  /*!
  \brief The lower limit of x coordinates of the node which is designated by the array index.
  */
  double* m_vLowerLimitOfX;
  
  /*!
  \brief The lower limit of y coordinates of the node which is designated by the array index.
  */
  double* m_vLowerLimitOfY;
  
  /*!
  \brief The lower limit of z coordinates of the node which is designated by the array index.
  */
  double* m_vLowerLimitOfZ;

  /*!
  \brief The upper limit of x coordinates of the node which is designated by the array index.
  */
  double* m_vUpperLimitOfX;
  
  /*!
  \brief The upper limit of y coordinates of the node which is designated by the array index.
  */
  double* m_vUpperLimitOfY;
  
  /*!
  \brief The upper limit of z coordinates of the node which is designated by the array index.
  */
  double* m_vUpperLimitOfZ;


  //vector<int>* m_vParticleNumList;

public:
  //Octree(int treedepth, int numOfParticles);
  Octree(int treedepth, size_t maxParticleNum);
  virtual ~Octree();


  //int initialize();
  /*!
  brief Building octree data structure.
  */
  int buildOctree(const double* x, const double* y, const double* z, size_t numParticles);
  //int buildOctree(const double *xp, const double *yp, const double *zp);

  int searchNeighbor(const double search_x, const double search_y, const double search_z, const double radius, int* result, size_t& result_length);

  int searchNeighbor(const double search_x, const double search_y, const double search_z, const double radius, SearchResult* result, size_t& result_length);

private:
  /*!
  \brief This method computes the Morton Key for a particle whose coordinates are x,y and z.
   */
  inline uint32_t computeKey(const double& x, const double& y, const double& z);

  //void oneDimIndexTo2DimPair(vector<int>* numParticles, int index, int& clusterIndex, int& particleIndex);
};


#endif // __OCTREE_H__


