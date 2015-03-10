#include "hexagonal_packing.h"
#include <cmath>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
// Start of HexagonalPacking2D
////////////////////////////////////////////////////////////////////////////////

HexagonalPacking2D::HexagonalPacking2D(double xmin_, double xmax_, double ymin_, double ymax_, double h_r_) 
:xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), h_r(h_r_) {
	
	//---------------compute parameters of hexagonal packing--------------
	m0 = n0_odd = n0_even = 0;
	//rows (for odd-numbered layers)
	m1 = (size_t) ((ymax-ymin-h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1; 
	//columns (for odd-numbered rows)
	n1_odd = (size_t) ((xmax-xmin-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
	//columns (for even-numbered rows)       
	n1_even = (size_t) ((xmax-xmin-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;  
	//------------Adjust the number of layers, rows, and columns-----------
	//------------to make consistent periodic buffers----------------------
	if(m1%2 == 0) m1--; //if there are odd number of rows, make it even
	if(n1_odd>n1_even) n1_even = n1_odd;
	else n1_odd = n1_even;	
	//---------------------------------------------------------------------

}

////////////////////////////////////////////////////////////////////////////////
// End of HexagonalPacking2D
////////////////////////////////////////////////////////////////////////////////








////////////////////////////////////////////////////////////////////////////////
// Start of HexagonalPacking3D
////////////////////////////////////////////////////////////////////////////////

HexagonalPacking3D::HexagonalPacking3D(double xmin_, double xmax_, double ymin_, double ymax_, 
									   double zmin_, double zmax_, double h_r_) 
: xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), zmin(zmin_), zmax(zmax_), h_r(h_r_) {


	//---------------compute parameters of hexagonal packing--------------		
	l0 = m0_odd = m0_even = n0_odd = n0_even = nn0_odd = nn0_even = 0;
	//layers
	l1 = (int) ((zmax-zmin-h_r-0.4999999*(2.0*sqrt(6.0)/3.0*h_r))/(2.0*sqrt(6.0)/3.0*h_r))+1; 
	//rows (for odd-numbered layers)
	m1_odd = (int) ((ymax-ymin-h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1; 
	//columns (for odd-numbered layers & odd-numbered rows)
	n1_odd = (int) ((xmax-xmin-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
	//columns (for odd-numbered layers & even-numbered rows)       
	n1_even = (int) ((xmax-xmin-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
	//rows (for even-numbered layers)
	m1_even = (int) ((ymax-ymin-sqrt(3.0)/3.0*h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1; 
	//columns (for even-numbered layers & odd-numbered rows)       
	nn1_odd = (int) ((xmax-xmin-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
	//columns (for even-numbered layers & even-numbered rows)
	nn1_even = (int) ((xmax-xmin-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
	//------------Adjust the number of layers, rows, and columns-----------
	//------------to make consistent periodic buffers----------------------
	if(l1%2 == 0) l1--; //if there are odd number of layers, make it even
	if(m1_odd>m1_even) m1_even = m1_odd;
	else m1_odd = m1_even;
	if(m1_odd%2 == 0) m1_odd--; //if there are odd number of rows, make it even
	if(m1_even%2 == 0) m1_even--; //if there are odd number of rows, make it even
	//odd layers
	if(n1_odd>n1_even) n1_even = n1_odd;
	else n1_odd = n1_even;
	//even layers
	if(nn1_odd>nn1_even) nn1_even = nn1_odd;
	else nn1_odd = nn1_even;
	if(n1_odd>nn1_odd) nn1_odd = n1_odd;
	else n1_odd = nn1_odd;
	n1_even = n1_odd;
	nn1_even = nn1_odd;	
	//-----------------------------------------------------------------


}

////////////////////////////////////////////////////////////////////////////////
// End of HexagonalPacking3D
////////////////////////////////////////////////////////////////////////////////
