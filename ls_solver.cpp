#include "ls_solver.h"
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
using namespace std;

#if defined(cray) || defined(_AIX)
#	define 	FORTRAN_NAME(a)		a
#else
#	define	FORTRAN_NAME(a)		a ## _
#endif

#if defined(__cplusplus)
extern "C" {
#endif

    //int FORTRAN_NAME(LAPACKE_dgeqp3)(int, int, int, double*, int, int*, double*);
	int LAPACKE_dgeqp3(int, int, int, double*, int, int*, double*);
	void FORTRAN_NAME(dgels)( char* trans, int* m, int* n, int* nrhs,
                              double* a, int* lda, double* b, int* ldb,
                              double* work, int* lwork, int *info );
	
#if defined(__cplusplus)
}
#endif


////////////////////////////////////////////////////////////////////////////////
// Start : QRSolver
////////////////////////////////////////////////////////////////////////////////

int QRSolver::solve(double* result, double* b) {
	//cout<<"-------------QRSolver::solve()--------------"<<endl;
	
	m_vb = b;
	
	//cout<<"Before solving. A="<<endl;
	//for(size_t i=0; i<m_iNumRow; i++) {
	//	for(size_t j=0; j<m_iNumCol; j++) 
	//		cout<<m_vA[i+j*m_iNumRow]<<"	";
	//	cout<<endl;
	//}
	//cout<<"b="<<endl;
	//for(size_t i=0; i<m_iNumRow; i++) cout<<m_vb[i]<<endl;


	if(!isDecomposed) {

		//---------------------------Start QR decomposition----------------------------------
		
		int JPVT_[m_iNumCol]; // Pivoting matrix
		double TAU_[m_iNumCol]; // Scalar factors of elementary reflectors	
		
		//FOR STORING RETURNING INFO BY DGEQP3 SUBROUTINE
		int info = 0;	
		
		// Perform QR decomposition
		//               CALLING LAPACK QR \W PIVOTING SUBROUTINE
		//-----------1. LAPACK_ROW_MAJOR (101) INDICATES MATRIX STORED IN ROW MAJOR
		//-----------   LAPACK_COL_MAJOR (102) INDICATES MATRIX STORED IN COL MAJOR
		//-----------2. m > n
		info = LAPACKE_dgeqp3(102, m_iNumRow, m_iNumCol, m_vA, m_iNumRow, JPVT_, TAU_);	
		
		//cout<<"After solving. A="<<endl;
		//for(size_t i=0; i<m_iNumRow; i++) {
		//	for(size_t j=0; j<m_iNumCol; j++) 
		//		cout<<m_vA[i+j*m_iNumRow]<<"	";
		//	cout<<endl;
		//}		
		//cout<<"info="<<info<<endl;	
		
		// return value will be negative if the -info th argument is illegal
		if(info != 0) { 
			//std::cout << "The " << -info << "th argument had an illegal value." <<std::endl;
			return info;
		}

		size_t rank = m_iNumCol;
		for(size_t i=1; i<m_iNumCol; i++) {
			if(fabs(m_vA[i*m_iNumRow+i]) <= m_fLimitR*fabs(m_vA[0])) {
				rank = i;
				break;
			}
		}
		//cout<<"rank="<<rank<<endl;
		

		// return value will be >0 (==rank) 
		if(rank<m_iNumCol) return rank;
		
		// so that next time will not do decomposition again
		isDecomposed = true;

		m_vTAU.assign(TAU_, TAU_+m_iNumCol);
		m_vJPVT.assign(JPVT_, JPVT_+m_iNumCol);
			
		//------------------------------End QR decomposition---------------------------------
	
	}
	
//	cout<<"m_vTAU:"<<endl;
//	for(size_t k=0; k<m_iNumCol; k++) 
//		cout<<m_vTAU[k]<<endl;
//	cout<<"m_vJPVT:"<<endl;
//	for(size_t k=0; k<m_iNumCol; k++) 
//		cout<<m_vJPVT[k]<<endl;

	//------------------------------Start Multiply QTb-----------------------------------
	
	// b will become QTb
	double v[m_iNumRow];	
	for(size_t i=0; i<m_iNumCol; i++) {
		
		// assign v
		for(size_t j=0; j<m_iNumRow; j++) {
			if(j < i) v[j] = 0.;
			else if(j == i) v[j] = 1.;
			else v[j] = m_vA[j+i*m_iNumRow];
		}

		double v_times_b = 0.;
		for(size_t j=0; j<m_iNumRow; j++)
			v_times_b += v[j]*m_vb[j];

		v_times_b *= m_vTAU[i];

		for(size_t j=0; j<m_iNumRow; j++)
			m_vb[j] -= v_times_b *v[j];      
	}
	//cout<<"QTb computed"<<endl;

	//--------------------------------End Multiply QTb-----------------------------------
	
	//------------------------------Start Backsubstitution-------------------------------		

	// compute PTx
	double PTx[m_iNumCol];
	fill_n(PTx,m_iNumCol,0);
	// original code starts from i=rank
	for(int i=m_iNumCol-1; i>=0; i--) { // Note that do not use size_t here because here is i-- 
		PTx[i] = m_vb[i]/m_vA[i*m_iNumRow+i];
		for(int j=0; j<i; j++) 
			m_vb[j] -= m_vA[j+i*m_iNumRow]*PTx[i];	  
	} 
	//cout<<"PTx computed"<<endl;

	// consider permutation to get the right result
	for(size_t i=1; i<=m_iNumCol; i++) {
		bool found = false;
		for(size_t j=0; j<m_iNumCol; j++) {
			if(m_vJPVT[j] == (int)i) {
				result[i-1] = PTx[j];
				found = true;
				break;
			}
		}
		if(found == false) assert(false);
		//cout<<"result["<<i-1<<"]="<<result[i-1]<<endl;
	}
	//cout<<"result obtained"<<endl;

	//--------------------------------End Backsubstitution-------------------------------
	
	//cout<<"--------------------------------------------"<<endl;
	return 0;
}


////////////////////////////////////////////////////////////////////////////////
// End : QRSolver
////////////////////////////////////////////////////////////////////////////////
