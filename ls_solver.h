/**
 * \file   ls_solver.h
 *
 * \brief  This header file contains classes for solvers of the least squares problem Ax = b 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design                
 *
 *
 * \version 1.0 
 *
 * \date 2014/10/09
 *
 * Created on: 2014/9/20 
 *
 */


#ifndef __LS_SOLVER_H__
#define __LS_SOLVER_H__

#include <cstddef>
#include <vector>



/**
 * \class LSSolver
 * 
 * \brief An abstract class for the family of solvers for the least squares problem
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
class LSSolver {
public:
	/**
	 * \brief Destructor
	 *
	 * \note The data members \e m_vA and \e m_vb will point to some memory outside this class
	 *       , as a result no memory should be deleted in this destructor
	 *
	 */
	virtual ~LSSolver() {}
	
	/**
	 * \brief Solves the least squares problem Ax = b 
	 *
	 * A is a mxn matrix, b is a mx1 vector, and x is the solution which is a nx1 vector. 
	 *
	 * \param [out] result The \e x in Ax = b, which is a nx1 vector
	 * \param [in]  b      The right-hand-size, which is a mx1 vector 
	 * \return             0 if success; otherwise if failure
	 * \warning The input array \e b \b may be modified in this method           
	 */
	virtual int solve(double* result, double* b) = 0; 

protected:	
	std::size_t m_iNumRow;///< The number of rows in matrix A	
	std::size_t m_iNumCol;///< The number of columns in matrix A	
	double *m_vA;///< The matrix A
	double *m_vb;///< The vector b (the right-hand-side)		
};






/**
 * \class QRSolver
 * 
 * \brief A class which solves the least squares problem by the QR decomposition method
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
class QRSolver : public LSSolver {
public:	
	/**
	 * \brief Constructor
	 *
	 * Set up the matrix A and other scalar parameters of this class
	 *
	 * \param [in] numRow The number of rows in matrix A 
	 * \param [in] numCol The number of columns in matrix A
	 * \param [in] A      The matrix A
	 * \param [in] limitR A specified scalar multiplier (<1) for trimming diagonal entry of R in the QR decomposition
	 *
	 */
	QRSolver(std::size_t numRow, std::size_t numCol, double *A, double limitR=1e-3) {  
		m_iNumRow = numRow;
		m_iNumCol = numCol;
		m_vA = A;
		isDecomposed = false; // The matrix A is not decomposed into QR at initialization
		m_fLimitR = limitR;
	}

	/*!
	\brief returns 0 if success and result will be non-empty
	       returns k<0 if the kth argument is illegal, result is empty
		   returns k>0, k=rank in this case, result is empty
	*/
	/**
	 * \brief Solves the least squares problem Ax = b by the QR decomposition method 
	 *
	 * A is a mxn matrix, b is a mx1 vector, and x is the solution which is a nx1 vector
	 * This method solves for \f$ x \f$ by the following steps:\n
	 * 1. Call the function \e LAPACKE_dgeqp3() in LAPACK 
	 * to perform QR decomposition with column pivoting to obtain \f[ QRPx = b \f]
	 * where \f$ P \f$ is a permutation matrix\n 
	 * 2. Multiply both sides on the left by \f$ Q^T  \f$ to get \f[ RPx = Q^T b \f]\n
	 * 3. Use back substitution to get \f$ y = Px \f$\n
	 * 4. Solves \f$ Px = y \f$ for \f$ x \f$\n
	 *
	 *
	 * \param [out] result The \e x in Ax = b, which is a nx1 vector
	 * \param [in]  b      The right-hand-size, which is a mx1 vector 
	 * \return             0 if success; otherwise if failure
	 * \warning The input array \e b \b is modified in this method           
	 */
	virtual int solve(double* result, double* b);

private:
	bool isDecomposed;///< if true then the matrix A is decomposed into QR already; if false otherwise
	double m_fLimitR;///< A specified scalar multiplier (<1) for trimming diagonal entry of R in the QR decomposition
	std::vector<int> m_vJPVT;///< The pivoting matrix 	
	std::vector<double> m_vTAU;///< The scalar factors of elementary reflectors
};


#endif // __LS_SOLVER_H__
