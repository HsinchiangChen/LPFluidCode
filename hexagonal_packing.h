/**
 * \file   hexagonal_packing.h
 *
 * \brief  This header file contains classes for creating particles based on the 2D and 3D hexagonal close packing  
 *
 * Refer to http://en.wikipedia.org/wiki/Close-packing_of_equal_spheres for the hexagonal close packing 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 *
 * \version 1.0 
 *
 * \date 2014/05/09
 *
 * Created on: 2014/05/01 
 *
 */


#ifndef __HEXAGONALPACKING_H__
#define __HEXAGONALPACKING_H__

#include <cstddef>
#include <cmath>


/**
 * \class HexagonalPacking2D
 * 
 * \brief Computes the 2D Cartesian coordinates of a particle based on the hexagonal close packing 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/05/09 
 *
 * Created on: 2014/05/01 
 *
 */
class HexagonalPacking2D {
public:
	/**
	 * \brief             Constructor
	 *
	 * The constructor initializes parameters needed to compute the 2D Cartesian coordinate of particles
	 * based on the hexagonal close packing using the input spatial domain and the input radius of the 
	 * constructing equal sphere. All particles will be located inside the input spatial domain. 
	 *
	 * \param [in] xmin_  The minimum value of the input spatial domain in x-coordinate
	 * \param [in] xmax_  The maximum value of the input spatial domain in x-coordinate
	 * \param [in] ymin_  The minimum value of the input spatial domain in y-coordinate
	 * \param [in] ymax_  The maximum value of the input spatial domain in y-coordinate
	 * \param [in] h_r_   The radius of the constructing sphere
	 *
	 * \note The parameter h_r_ should be set to be one half of the initial inter-particle distance
	 *
	 */
	HexagonalPacking2D(double xmin_, double xmax_, double ymin_, double ymax_, double h_r_);
	
	/**
	 * \brief        Getter function to retrieve parameters of the boundaries of this hexagonal close packing  
	 *  
	 * \param [out]  m0_       The minimum index of rows
	 * \param [out]  m1_       The maximum index of rows
	 * \param [out]  n0_odd_   The minimum index of columns in odd-numbered rows
	 * \param [out]  n1_odd_   The maximum index of columns in odd-numbered rows
	 * \param [out]  n0_even_  The minimum index of columns in even-numbered rows
	 * \param [out]  n1_even_  The maximum index of columns in even-numbered rows
	 *
	 * \note         The information is retrieved to compute the Cartesian coordinate of particles one-by-one outside
	 *               this class. This enables the caller to use different criterions to screen each particle
	 *               to decide if this particle should be included into the computational context or not.
	 *               This will save a lot of computational time compared to if all particles within the 
	 *               input spatial domain is included in the computational context first, and then
	 *               screen and delete unnecessary particles.
	 */
	void getParameters(size_t& m0_, size_t& m1_, 
				       size_t& n0_odd_, size_t& n1_odd_, 
				       size_t& n0_even_, size_t& n1_even_) {
		m0_ = m0; 
		m1_ = m1;
		n0_odd_ = n0_odd; 
		n1_odd_ = n1_odd;
		n0_even_ = n0_even;
		n1_even_ = n1_even;
	}
	
	/**
	 * \brief           Computes the particle location in the x-coordinate        
     * \param [in] tag  Value 0 indicates odd-numbered rows; otherwise indicates even-numbered rows  
	 * \param [in] k    The k-th column 
	 * \return          Particle location in the x-coordinate
	 *
	 */
	double computeX(int tag, size_t k) {
		if(tag==0) return xmin+2.0*h_r+k*2.0*h_r; //odd-numbered rows
		else return  xmin+h_r+k*2.0*h_r; //even-numbered rows
	}
	
	/**
	 * \brief           Computes the particle location in the y-coordinate          
	 * \param [in] j    The j-th row 
	 * \return          Particle location in the y-coordinate
	 *
	 */
	double computeY(size_t j) {
		return ymin+h_r+j*sqrt(3.0)*h_r;
	}

private:
	double xmin, xmax, ymin, ymax, h_r;
	size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;		
	//HexagonalPacking2D() {}
	//HexagonalPacking2D(const HexagonalPacking2D& other) {}
	//HexagonalPacking2D& operator=(const HexagonalPacking2D& other) {return *this;}
};










/**
 * \class HexagonalPacking3D
 * 
 * \brief Computes the 3D Cartesian coordinates of a particle based on the hexagonal close packing 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2014/05/09 
 *
 * Created on: 2014/05/01 
 *
 */
class HexagonalPacking3D {
public:
	/**
	 * \brief             Constructor
	 *
	 * The constructor initializes parameters needed to compute the 3D Cartesian coordinate of particles
	 * based on the hexagonal close packing using the input spatial domain and the input radius of the 
	 * constructing equal sphere. All particles will be located inside the input spatial domain. 
	 *
	 * \param [in] xmin_  The minimum value of the input spatial domain in x-coordinate
	 * \param [in] xmax_  The maximum value of the input spatial domain in x-coordinate
	 * \param [in] ymin_  The minimum value of the input spatial domain in y-coordinate
	 * \param [in] ymax_  The maximum value of the input spatial domain in y-coordinate
	 * \param [in] zmin_  The minimum value of the input spatial domain in z-coordinate
	 * \param [in] zmax_  The maximum value of the input spatial domain in z-coordinate
	 * \param [in] h_r_   The radius of the constructing sphere
	 *
	 * \note The parameter h_r_ should be set to be one half of the initial inter-particle distance
	 *
	 */
	HexagonalPacking3D(double xmin_, double xmax_, double ymin_, double ymax_, double zmin_, double zmax_, double h_r_);
	
	/**
	 * \brief        Getter function to retrieve parameters of the boundaries of this hexagonal close packing  
	 *  
	 * \param [out]  l0_        The minimum index of layers
	 * \param [out]  l1_        The maximum index of layers
	 * \param [out]  m0_odd_    The minimum index of rows in odd-numbered layers
	 * \param [out]  m1_odd_    The maximum index of rows in odd-numbered layers
	 * \param [out]  m0_even_   The minimum index of rows in even-numbered layers
	 * \param [out]  m1_even_   The maximum index of rows in even-numbered layers
	 * \param [out]  n0_odd_    The minimum index of columns in odd rows & odd layers 
	 * \param [out]  n1_odd_    The maximum index of columns in odd rows & odd layers
	 * \param [out]  n0_even_   The minimum index of columns in even rows & odd layers
	 * \param [out]  n1_even_   The maximum index of columns in even rows & odd layers
	 * \param [out]  nn0_odd_   The minimum index of columns in odd rows & even layers
	 * \param [out]  nn1_odd_   The maximum index of columns in odd rows & even layers
	 * \param [out]  nn0_even_  The minimum index of columns in even rows & even layers
	 * \param [out]  nn1_even_  The maximum index of columns in even rows & even layers
	 *
	 * \note         The information is retrieved to compute the Cartesian coordinate of particles one-by-one outside
	 *               this class. This enables the caller to use different criterions to screen each particle
	 *               to decide if this particle should be included into the computational context or not.
	 *               This will save a lot of computational time compared to if all particles within the 
	 *               input spatial domain is included in the computational context first, and then
	 *               screen and delete unnecessary particles.
	 */
	void getParameters(size_t& l0_, size_t& l1_, 
					   size_t& m0_odd_, size_t& m1_odd_, 
					   size_t& m0_even_, size_t& m1_even_,
					   size_t& n0_odd_, size_t& n1_odd_,
					   size_t& n0_even_, size_t& n1_even_,
					   size_t& nn0_odd_, size_t& nn1_odd_,
					   size_t& nn0_even_, size_t& nn1_even_) {
		l0_ = l0; 
		l1_ = l1; 
	    m0_odd_ = m0_odd; 
		m1_odd_ = m1_odd; 
	    m0_even_ = m0_even; 
		m1_even_ = m1_even;
	    n0_odd_ = n0_odd; 
		n1_odd_ = n1_odd;
	    n0_even_ = n0_even; 
		n1_even_ = n1_even;
	    nn0_odd_ = nn0_odd; 
		nn1_odd_ = nn1_odd;
	    nn0_even_ = nn0_even; 
		nn1_even_ = nn1_even;
	}
	

	/**
	 * \brief           Computes the particle location in the x-coordinate        
     * \param [in] tag  Value 0 indicates: 1. odd layers and odd rows, or 2. even layers and even rows 
	 *                  ; otherwise indicates: 1. odd layers and even rows, or 2. even layers and odd rows  
	 * \param [in] k    The k-th column 
	 * \return          Particle location in the x-coordinate
	 *
	 */
	double computeX(int tag, size_t k) {
		if(tag==0) return xmin+2.0*h_r+k*2.0*h_r; //(odd layers && odd rows) || (even layers && even rows)
		else return  xmin+h_r+k*2.0*h_r; // (odd layers && even rows) || (even layers && odd rows)
	}

	/**
	 * \brief           Computes the particle location in the y-coordinate          
	 * \param [in] tag  Value 0 indicates: odd-numbered layers; otherwise indicates even-numbered layers
	 * \param [in] j    The j-th row 
	 * \return          Particle location in the y-coordinate
	 *
	 */
	double computeY(int tag, size_t j) {
		if(tag==0) return ymin+h_r+j*sqrt(3.0)*h_r;	//odd-numbered layers
		else return ymin+sqrt(3.0)/3.0*h_r+j*sqrt(3.0)*h_r; //even-numbered layers
	}


	/**
	 * \brief           Computes the particle location in the z-coordinate           
	 * \param [in] j    The i-th layer 
	 * \return          Particle location in the z-coordinate
	 *
	 */
	double computeZ(size_t i) {
		return zmin+h_r+i*2.0*sqrt(6.0)/3.0*h_r;
	}
	
private:
	
	double xmin, xmax, ymin, ymax, zmin, zmax, h_r;
	size_t l0,l1;
	size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
	size_t nn0_odd, nn1_odd, nn0_even, nn1_even;		
	//HexagonalPacking3D() {}
	//HexagonalPacking3D(const HexagonalPacking3D& other) {}
	//HexagonalPacking3D& operator=(const HexagonalPacking3D& other) {return *this;}
};


#endif // __HEXAGONALPACKING_H__



