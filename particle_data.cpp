#include "particle_data.h"
#include "initializer.h"
#include <algorithm>
//#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Start : ParticleData
////////////////////////////////////////////////////////////////////////////////


ParticleData::ParticleData(const Initializer& init) {
	
	// get the scalar parameters
	m_iDimension = init.getDimension(); 
	m_iCapacity = init.getCapacity();
	m_iMaxNeighbourNum = init.getMaxNeighbourNum(); 
	m_iMaxNeighbourNumInOneDir = init.getMaxNeighbourNumInOneDir();	
	
	m_iFluidNum = init.getFluidNum();
	m_iBoundaryNum = init.getBoundaryNum();
	m_iGhostNum = 0;
	m_iTotalNum = m_iFluidNum + m_iBoundaryNum + m_iGhostNum;
	
	m_iFluidStartIndex = init.getFluidStartIndex();
	m_iBoundaryStartIndex = init.getBoundaryStartIndex();
	m_iGhostStartIndex = init.getGhostStartIndex();

	// get the resources allocated on heap by the Initializer class
	m_vPositionX = init.getPositionX();
	m_vPositionY = init.getPositionY();	
	m_vPositionZ = init.getPositionZ();
	m_vVelocityU = init.getVelocityU();
	m_vVelocityV = init.getVelocityV();
	m_vVelocityW = init.getVelocityW();
	m_vVolume = init.getVolume();
	m_vPressure = init.getPressure();
	m_vSoundSpeed = init.getSoundSpeed();
	
	m_vObjectTag = init.getObjectTag();

	const vector<BoundingBox*>& tmp = init.getFluidBoundingBox();	
	for(size_t i=0; i<tmp.size(); i++) 
		m_vFluidBoundingBox.push_back(new BoundingBox(*tmp[i]));

	// allocate memory on the heap
	m_vTemp1VelocityU = new double[m_iCapacity];
	m_vTemp1VelocityV = new double[m_iCapacity];
	if(m_iDimension==3) m_vTemp1VelocityW = new double[m_iCapacity];
	m_vTemp1Volume = new double[m_iCapacity];
	m_vTemp1Pressure = new double[m_iCapacity];
	m_vTemp1SoundSpeed = new double[m_iCapacity];
	
	fill_n(m_vTemp1VelocityU,m_iCapacity,0);
	fill_n(m_vTemp1VelocityV,m_iCapacity,0);	
	if(m_iDimension==3) fill_n(m_vTemp1VelocityW,m_iCapacity,0);
	fill_n(m_vTemp1Volume,m_iCapacity,0);
	fill_n(m_vTemp1Pressure,m_iCapacity,0);
	fill_n(m_vTemp1SoundSpeed,m_iCapacity,0);

	m_vTemp2VelocityU = new double[m_iCapacity];
	m_vTemp2VelocityV = new double[m_iCapacity];
	if(m_iDimension==3) m_vTemp2VelocityW = new double[m_iCapacity];
	m_vTemp2Volume = new double[m_iCapacity];
	m_vTemp2Pressure = new double[m_iCapacity];
	m_vTemp2SoundSpeed = new double[m_iCapacity];

	fill_n(m_vTemp2VelocityU,m_iCapacity,0);
	fill_n(m_vTemp2VelocityV,m_iCapacity,0);	
	if(m_iDimension==3) fill_n(m_vTemp2VelocityW,m_iCapacity,0);
	fill_n(m_vTemp2Volume,m_iCapacity,0);
	fill_n(m_vTemp2Pressure,m_iCapacity,0);
	fill_n(m_vTemp2SoundSpeed,m_iCapacity,0);

	m_vLPFOrderRight = new int[m_iCapacity];
	m_vLPFOrderLeft = new int[m_iCapacity];
	m_vLPFOrderNorth = new int[m_iCapacity];
	m_vLPFOrderSouth = new int[m_iCapacity]; 
	if(m_iDimension==3) m_vLPFOrderUp = new int[m_iCapacity];
	if(m_iDimension==3) m_vLPFOrderDown = new int[m_iCapacity];	
	
	fill_n(m_vLPFOrderRight,m_iCapacity,0);
	fill_n(m_vLPFOrderLeft,m_iCapacity,0);	
	fill_n(m_vLPFOrderNorth,m_iCapacity,0);
	fill_n(m_vLPFOrderSouth,m_iCapacity,0);
	if(m_iDimension==3) fill_n(m_vLPFOrderUp,m_iCapacity,0);
	if(m_iDimension==3) fill_n(m_vLPFOrderDown,m_iCapacity,0);
	
	m_vNeighbourList = new int[m_iCapacity*m_iMaxNeighbourNum]; // whole
	m_vNeighbourListRight = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // right
	m_vNeighbourListLeft = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // left
	m_vNeighbourListNorth = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // north
	m_vNeighbourListSouth = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // south
	if(m_iDimension==3) m_vNeighbourListUp = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // up
	if(m_iDimension==3)	m_vNeighbourListDown = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // down		
	
	fill_n(m_vNeighbourList,m_iCapacity,0);
	fill_n(m_vNeighbourListRight,m_iCapacity,0);
	fill_n(m_vNeighbourListLeft,m_iCapacity,0);
	fill_n(m_vNeighbourListNorth,m_iCapacity,0);
	fill_n(m_vNeighbourListSouth,m_iCapacity,0);
	if(m_iDimension==3) fill_n(m_vNeighbourListUp,m_iCapacity,0);
	if(m_iDimension==3) fill_n(m_vNeighbourListDown,m_iCapacity,0);	
	
	m_vNeighbourListSize = new int [m_iCapacity];
	m_vNeighbourListRightSize = new int [m_iCapacity];
	m_vNeighbourListLeftSize = new int [m_iCapacity];
	m_vNeighbourListNorthSize = new int [m_iCapacity];
	m_vNeighbourListSouthSize = new int [m_iCapacity];
	if(m_iDimension==3) m_vNeighbourListUpSize = new int [m_iCapacity];
	if(m_iDimension==3) m_vNeighbourListDownSize = new int [m_iCapacity];
	
	fill_n(m_vNeighbourListSize,m_iCapacity,0);
	fill_n(m_vNeighbourListRightSize,m_iCapacity,0);
	fill_n(m_vNeighbourListLeftSize,m_iCapacity,0);
	fill_n(m_vNeighbourListNorthSize,m_iCapacity,0);
	fill_n(m_vNeighbourListSouthSize,m_iCapacity,0);
	if(m_iDimension==3) fill_n(m_vNeighbourListUpSize,m_iCapacity,0);
	if(m_iDimension==3) fill_n(m_vNeighbourListDownSize,m_iCapacity,0);
	
}



ParticleData::~ParticleData() {
  
	delete[] m_vPositionX;
	delete[] m_vPositionY;
	delete[] m_vPositionZ;
	//if(m_iDimension==3) delete[] m_vPositionZ;

	delete[] m_vVelocityU;
	delete[] m_vVelocityV;
	if(m_iDimension==3) delete[] m_vVelocityW;

	delete[] m_vVolume;
	delete[] m_vPressure;
	//delete[] m_vEnergy;
	delete[] m_vSoundSpeed;

	delete[] m_vTemp1VelocityU;
	delete[] m_vTemp1VelocityV;
	if(m_iDimension==3) delete[] m_vTemp1VelocityW;
	delete[] m_vTemp1Volume;
	delete[] m_vTemp1Pressure;
	delete[] m_vTemp1SoundSpeed;

	delete[] m_vTemp2VelocityU;
	delete[] m_vTemp2VelocityV;
	if(m_iDimension==3) delete[] m_vTemp2VelocityW;
	delete[] m_vTemp2Volume;
	delete[] m_vTemp2Pressure;
	delete[] m_vTemp2SoundSpeed;

	delete[] m_vLPFOrderRight;
	delete[] m_vLPFOrderLeft;
	delete[] m_vLPFOrderNorth;
	delete[] m_vLPFOrderSouth;
	if(m_iDimension==3) delete[] m_vLPFOrderUp;
	if(m_iDimension==3) delete[] m_vLPFOrderDown;
	
	
	delete[] m_vNeighbourList;
	delete[] m_vNeighbourListRight;
	delete[] m_vNeighbourListLeft;
	delete[] m_vNeighbourListNorth;
	delete[] m_vNeighbourListSouth;	
	if(m_iDimension==3) delete[] m_vNeighbourListUp;	
	if(m_iDimension==3)	delete[] m_vNeighbourListDown;
	

	delete[] m_vNeighbourListSize;
	delete[] m_vNeighbourListRightSize;
	delete[] m_vNeighbourListLeftSize;
	delete[] m_vNeighbourListNorthSize;
	delete[] m_vNeighbourListSouthSize;
	delete[] m_vNeighbourListUpSize;
	delete[] m_vNeighbourListDownSize;
	
	delete[] m_vObjectTag;
}


////////////////////////////////////////////////////////////////////////////////
// End : ParticleData
////////////////////////////////////////////////////////////////////////////////
