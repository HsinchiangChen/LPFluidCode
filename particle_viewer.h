/**
 * \file   particle_viewer.h
 *
 * \brief  This header file contains classes for output simulation results in various formats such as .vtk and .txt  
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

#ifndef __PARTICLE_VIEWER_H__
#define __PARTICLE_VIEWER_H__

#include <cstddef>
#include <string>


class ParticleData;


/**
 * \class ParticleViewer
 * 
 * \brief An abstract class for classes that output simulations results 
 *        
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
class ParticleViewer {

public:
	/**
	 * \brief                Write simulation results to the output file
	 * \param [in] time      The physical output time
	 * \param [in] writeStep The number of times of output
	 * \return               0 if output success; 1 otherwise
	 */
	virtual int writeResult(double time, std::size_t writeStep) = 0;

protected:		
	ParticleData* m_pParticleData;///< A pointer to the object which holds particle information and data
	std::string m_sOutputfileName;///< The name of the output file 
	int m_iNumDigits;///< The number of digits for the indexing of output file name 
	std::string m_sParticleType;///< The type of particle data to output (all, fluid, boundary, etc)

	/**
	 * \brief                A small helper function that organizes the format of output filename
	 * \param [in] writeStep The number of times of output
	 * \param [in] numDigits The number of digits for the indexing of output file name
	 *
	 */
	std::string rightFlush(std::size_t writeStep, std::size_t numDigits);

};



/**
 * \class VTKParticleViewer
 * 
 * \brief A class that output simulation results in the .vtk format 
 *        
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
class VTKParticleViewer : public ParticleViewer {

public:
	/**
	 * \brief                      Constructor 
	 *
	 * \param [in] data            A pointer to the object which holds particle information and data
	 * \param [in] particleType    The type of particle to output (all, fluid, boundary, ghost)
	 * \param [in] outputfileName  The name of the output file
	 * \param [in] numDigits       The number of digits for the indexing of output file name
	 *                
	 */
	VTKParticleViewer(ParticleData* data, const std::string& particleType, const std::string& outputfileName="", int numDigits=7);
	
	/**
	 * \brief                Write simulation results to the output file in the .vtk format
	 * \param [in] time      The physical output time
	 * \param [in] writeStep The number of times of output
	 * \return               0 if output success; 1 otherwise
	 */
	virtual int writeResult(double time, std::size_t writeStep);
	
};







/**
 * \class TXTParticleViewer1D
 * 
 * \brief A class that output 1D simulation results in the .txt format 
 *        
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 *
 * \version 1.0 
 *
 * \date 2015/03/14 
 *
 * Created on: 2015/03/14 
 *
 */
class TXTParticleViewer1D : public ParticleViewer {

public:
	/**
	 * \brief                      Constructor 
	 *
	 * \param [in] data            A pointer to the object which holds particle information and data
	 * \param [in] particleType    The type of particle to output (all, fluid, boundary, ghost)
	 * \param [in] outputfileName  The name of the output file
	 * \param [in] numDigits       The number of digits for the indexing of output file name
	 *                
	 */
	TXTParticleViewer1D(ParticleData* data, const std::string& particleType, const std::string& outputfileName="", int numDigits=7);
	
	/**
	 * \brief                Write 1D simulation results to the output file in the .txt format
	 * \param [in] time      The physical output time
	 * \param [in] writeStep The number of times of output
	 * \return               0 if output success; 1 otherwise
	 */
	virtual int writeResult(double time, std::size_t writeStep);
	
};

#endif // __PARTICLE_VIEWER_H__
