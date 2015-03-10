#include "initializer.h"
#include "geometry.h"
#include "state.h"
#include "eos.h"
#include "hexagonal_packing.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of Initializer
////////////////////////////////////////////////////////////////////////////////////////

Initializer::Initializer(const string& inputFileName) {
	
	ifstream ifs(inputFileName);
	this->ofs.open("inputdata_verification");

	vector<string> lines;
	string s;
	while(getline(ifs, s)) {
		cout<<s<<endl;
		lines.push_back(s);
	}

	istringstream iss;
	size_t i = 0; // input file line number
	
	iss.str(lines[i++]);
	iss>>m_iNumThreads;
	ofs<<"m_iNumThreads = "<<m_iNumThreads<<endl;

	iss.str(lines[i++]);
	iss>>m_fStartTime;
	ofs<<"m_fStartTime = "<<m_fStartTime<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fEndTime;
	ofs<<"m_fEndTime = "<<m_fEndTime<<endl;

	iss.str(lines[i++]);
	iss>>m_fWriteTimeInterval;
	ofs<<"m_fWriteTimeInterval = "<<m_fWriteTimeInterval<<endl;

	iss.str(lines[i++]);
	iss>>m_fCFLCoeff;
	ofs<<"m_fCFLCoeff = "<<m_fCFLCoeff<<endl;

	iss.str(lines[i++]);
	iss>>m_iDimension;
	ofs<<"m_iDimension = "<<m_iDimension<<endl;

	iss.str(lines[i++]);
	iss>>m_iFluidObjNum;
	ofs<<"m_iFluidObjNum= "<<m_iFluidObjNum<<endl;

	vector<string> fluidObjNames; // temp vector
	for(int j=0; j<m_iFluidObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		fluidObjNames.push_back(tmpS);
		ofs<<"fluidObjectName"<<j<<"="<<fluidObjNames[j]<<endl;
	}
	
	vector<string> fluidObjStateNames; // temp vector
	for(int j=0; j<m_iFluidObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		fluidObjStateNames.push_back(tmpS);
		ofs<<"fluidObjectStateName"<<j<<"="<<fluidObjStateNames[j]<<endl;
	}

	iss.str(lines[i++]);
	iss>>m_iBoundaryObjNum;
	ofs<<"m_iBoundaryObjNum= "<<m_iBoundaryObjNum<<endl;

	vector<string> boundaryObjNames; // temp vector
	for(int j=0; j<m_iBoundaryObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		boundaryObjNames.push_back(tmpS);
		ofs<<"boundaryObjectName"<<j<<"="<<boundaryObjNames[j]<<endl;
	}

	iss.str(lines[i++]);
	iss>>m_iRandomDirSplitOrder;
	ofs<<"m_iRandomDirSplitOrder = "<<m_iRandomDirSplitOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iLPFOrder;
	ofs<<"m_iLPFOrder = "<<m_iLPFOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iEOSChoice;
	ofs<<"m_iEOSChoice = "<<m_iEOSChoice<<endl;

	iss.str(lines[i++]);
	iss>>m_fGamma;
	ofs<<"m_fGamma = "<<m_fGamma<<endl;

	iss.str(lines[i++]);
	iss>>m_fPinf;
	ofs<<"m_fPinf = "<<m_fPinf<<endl;

	iss.str(lines[i++]);
	iss>>m_fEinf;
	ofs<<"m_fEinf = "<<m_fEinf<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fInitParticleSpacing;
	ofs<<"m_fInitParticleSpacing = "<<m_fInitParticleSpacing<<endl;	

	iss.str(lines[i++]);
	iss>>m_fGravity;
	ofs<<"m_fGravity = "<<m_fGravity<<endl;

	iss.str(lines[i++]);
	iss>>m_iMovingBoxForGhostParticle;
	ofs<<"m_iMovingBoxForGhostParticle = "<<m_iMovingBoxForGhostParticle<<endl;	
	
	ofs<<"-------------------Input Data as verified above-----------------------"<<endl;


	// create EOS object and get a pointer to it
	if(m_iEOSChoice==1) m_pEOS = new PolytropicGasEOS(m_fGamma);
	else if(m_iEOSChoice==2) m_pEOS = new StiffPolytropicGasEOS(m_fGamma,m_fPinf,m_fEinf);
	else assert(false);

	// create fluid objects and get pointers to them
	for(int j=0; j<m_iFluidObjNum; j++) 
		m_vFluidObj.push_back(GeometryFactory::instance().createGeometry(fluidObjNames[j]));
	
	// create fluid state objects and get pointers to them
	for(int j=0; j<m_iFluidObjNum; j++) 
		m_vFluidObjState.push_back(StateFactory::instance().createState(fluidObjStateNames[j]));

	// create boundary objects and get pointers to them
	for(int j=0; j<m_iBoundaryObjNum; j++) 
		m_vBoundaryObj.push_back(GeometryFactory::instance().createGeometry(boundaryObjNames[j]));
	
	// bounding box of fluid and boundary objects
	computeInitBoundingBox();	
	
	initGeometryAndState();	
	
	setBoundingBoxStartIndex();

	setObjectTag();

	setParams();	

}

Initializer::~Initializer() {
	for(auto obj:m_vFluidObj) delete obj;
	for(auto box:m_vFluidBoundingBox) delete box;
	for(auto obj:m_vBoundaryObj) delete obj;
	for(auto box:m_vBoundaryBoundingBox) delete box;
	for(auto state:m_vFluidObjState) delete state;
}

void Initializer::computeInitBoundingBox() {
	
		
	for(size_t i=0; i<m_vFluidObj.size(); i++) {
		double xmin, xmax, ymin, ymax, zmin, zmax;
		m_vFluidObj[i]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		xmin -= 3.*m_fInitParticleSpacing;
		xmax += 3.*m_fInitParticleSpacing;
		ymin -= 3.*m_fInitParticleSpacing;
		ymax += 3.*m_fInitParticleSpacing;
		if(m_iDimension==3) {
			zmin -= 3.*m_fInitParticleSpacing;
			zmax += 3.*m_fInitParticleSpacing;
		}	

		BoundingBox* box = new BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		m_vFluidBoundingBox.push_back(box);

		ofs<<"m_vFluidBoundingBox("<<i<<"):"<<endl;
		ofs<<"m_fXmin="<<m_vFluidBoundingBox[i]->getXmin()<<endl;
		ofs<<"m_fXmax="<<m_vFluidBoundingBox[i]->getXmax()<<endl;
		ofs<<"m_fYmin="<<m_vFluidBoundingBox[i]->getYmin()<<endl;
		ofs<<"m_fYmax="<<m_vFluidBoundingBox[i]->getYmax()<<endl;
		ofs<<"m_fZmin="<<m_vFluidBoundingBox[i]->getZmin()<<endl;
		ofs<<"m_fZmax="<<m_vFluidBoundingBox[i]->getZmax()<<endl;
		ofs<<"-----------------"<<endl;
	}

	for(size_t i=0; i<m_vBoundaryObj.size(); i++) {
		double xmin, xmax, ymin, ymax, zmin, zmax;
		m_vBoundaryObj[i]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);

		xmin -= 3.*m_fInitParticleSpacing;
		xmax += 3.*m_fInitParticleSpacing;
		ymin -= 3.*m_fInitParticleSpacing;
		ymax += 3.*m_fInitParticleSpacing;
		if(m_iDimension==3) {
			zmin -= 3.*m_fInitParticleSpacing;
			zmax += 3.*m_fInitParticleSpacing;
		}	

		BoundingBox* box = new BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		m_vBoundaryBoundingBox.push_back(box);

		ofs<<"m_vBoundaryBoundingBox("<<i<<"):"<<endl;
		ofs<<"m_fXmin="<<m_vBoundaryBoundingBox[i]->getXmin()<<endl;
		ofs<<"m_fXmax="<<m_vBoundaryBoundingBox[i]->getXmax()<<endl;
		ofs<<"m_fYmin="<<m_vBoundaryBoundingBox[i]->getYmin()<<endl;
		ofs<<"m_fYmax="<<m_vBoundaryBoundingBox[i]->getYmax()<<endl;
		ofs<<"m_fZmin="<<m_vBoundaryBoundingBox[i]->getZmin()<<endl;
		ofs<<"m_fZmax="<<m_vBoundaryBoundingBox[i]->getZmax()<<endl;
		ofs<<"-----------------"<<endl;
	}
	
}

void Initializer::setParams() {
		
	if(m_iLPFOrder==2)	
		m_fNeiSearchRadius = 2.2*m_fInitParticleSpacing;
	else if(m_iLPFOrder==1)
		m_fNeiSearchRadius = 2.*m_fInitParticleSpacing;
	
	m_fContactLength = 1.1*m_fInitParticleSpacing;

	// compute the number of particles within m_fNeiSearchRadius based on the initial packing of particles
	computeNumParticleWithinSearchRadius();
	
	if(m_iDimension==3) {
		m_iMaxNeighbourNum = 2*m_iNumParticleWithinSearchRadius;
		m_iMaxNeighbourNumInOneDir = m_iMaxNeighbourNum/3;
				
		m_iNumRow2ndOrder = 20; //TODO
		m_iNumRow1stOrder = 5; //TODO
		m_iNumCol2ndOrder = 9; 
		m_iNumCol1stOrder = 3; 
	}
	else if(m_iDimension==2) {
		m_iMaxNeighbourNum = 7*m_iNumParticleWithinSearchRadius; 
		m_iMaxNeighbourNumInOneDir = m_iMaxNeighbourNum/2; 
				
		m_iNumRow2ndOrder = 7; 
		m_iNumRow1stOrder = 3; 
		m_iNumCol2ndOrder = 5; 
		m_iNumCol1stOrder = 2;	
	}

	// (eos chocie=1:poly 2:spoly)	
	if(m_iEOSChoice==1) 
		m_fInvalidPressure = 0;
	else if(m_iEOSChoice==2) 	
		m_fInvalidPressure = -m_fPinf;
	m_fInvalidVolume = 0;
	
	m_iTreeDepth = 5;
	
	ofs<<"m_fNeiSearchRadius = "<<m_fNeiSearchRadius<<endl;
	ofs<<"m_iNumParticleWithinSearchRadius = "<<m_iNumParticleWithinSearchRadius<<endl;
	ofs<<"m_iMaxNeighbourNum = "<<m_iMaxNeighbourNum<<endl;
	ofs<<"m_iMaxNeighbourNumInOneDir = "<<m_iMaxNeighbourNumInOneDir<<endl;
	ofs<<"m_iNumRow2ndOrder = "<<m_iNumRow2ndOrder<<endl;
	ofs<<"m_iNumRow1stOrder = "<<m_iNumRow1stOrder<<endl;
	ofs<<"m_iNumCol2ndOrder = "<<m_iNumCol2ndOrder<<endl;
	ofs<<"m_iNumCol1stOrder = "<<m_iNumCol1stOrder<<endl;	
	ofs<<"m_fInvalidPressure = "<<m_fInvalidPressure<<endl;
	ofs<<"m_fInvalidVolume = "<<m_fInvalidVolume<<endl;
	ofs<<"m_iTreeDepth = "<<m_iTreeDepth<<endl;	
	ofs<<"-------------------Prespecified params as above-----------------------"<<endl;
	
}


void Initializer::initGeometryAndState() {
	
	bool saveData = false;
	m_iFluidNum = initGeometryAndStateOnHexPacking(saveData, "fluid");
	m_iBoundaryNum = initGeometryAndStateOnHexPacking(saveData, "boundary");

	m_iCapacity = (size_t)(1.5*(m_iFluidNum+m_iBoundaryNum));
	
	// use m_iCapacity to create memory
	initParticleDataMemory();
	
	saveData = true;
	initGeometryAndStateOnHexPacking(saveData, "fluid");
	initGeometryAndStateOnHexPacking(saveData, "boundary");

	// assign the start index of fluid, boundary, and ghost particles
	m_iFluidStartIndex = 0;
	m_iBoundaryStartIndex = m_iFluidStartIndex + m_iFluidNum;
	m_iGhostStartIndex = m_iFluidStartIndex + m_iFluidNum + m_iBoundaryNum;
	
	/* // TO DELETE
	size_t tmp=0;
	for(size_t index=0; index<m_iCapacity; index++) m_vObjectTag[index] = -1; // init object tag
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) {
		// set the start index of fluid bound box
		m_vFluidBoundingBox[p]->setStartIndex(m_iFluidStartIndex+tmp); 
		size_t num = m_vFluidBoundingBox[p]->getNumber();
		// set the fluid object tag
		for(size_t index=m_iFluidStartIndex+tmp; index<m_iFluidStartIndex+tmp+num; index++) {
			m_vObjectTag[index] = p;
		}
		tmp += num;
		ofs<<"m_vFluidBoundingBox["<<p<<"]->m_iStartIndex="<<m_vFluidBoundingBox[p]->getStartIndex()<<endl;
		ofs<<"m_vFluidBoundingBox["<<p<<"]->m_iNumber="<<num<<endl;
	}		
	*/

	ofs<<"m_iFluidNum = "<<m_iFluidNum<<endl;
	ofs<<"m_iBoundaryNum = "<<m_iBoundaryNum<<endl;
	ofs<<"m_iCapacity = "<<m_iCapacity<<endl;
	ofs<<"m_iFluidStartIndex = "<<m_iFluidStartIndex<<endl;
	ofs<<"m_iBoundaryStartIndex = "<<m_iBoundaryStartIndex<<endl;
	ofs<<"m_iGhostStartIndex = "<<m_iGhostStartIndex<<endl;
	
	ofs<<"-----------------Initialized geometry and state info as above-------------------"<<endl;
	
}


void Initializer::setBoundingBoxStartIndex() {
	size_t tmp=0;	
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) {
		// set the start index of fluid bound box
		m_vFluidBoundingBox[p]->setStartIndex(m_iFluidStartIndex+tmp); 
		size_t num = m_vFluidBoundingBox[p]->getNumber();	
		tmp += num;
		ofs<<"m_vFluidBoundingBox["<<p<<"]->m_iStartIndex="<<m_vFluidBoundingBox[p]->getStartIndex()<<endl;
		ofs<<"m_vFluidBoundingBox["<<p<<"]->m_iNumber="<<num<<endl;
	}
	tmp=0;	
	for(size_t p=0; p<m_vBoundaryBoundingBox.size(); p++) {
		// set the start index of boundary bound box
		m_vBoundaryBoundingBox[p]->setStartIndex(m_iBoundaryStartIndex+tmp); 
		size_t num = m_vBoundaryBoundingBox[p]->getNumber();	
		tmp += num;
		ofs<<"m_vBoundaryBoundingBox["<<p<<"]->m_iStartIndex="<<m_vBoundaryBoundingBox[p]->getStartIndex()<<endl;
		ofs<<"m_vBoundaryBoundingBox["<<p<<"]->m_iNumber="<<num<<endl;
	}
	ofs<<"-----------------Initialized Bounding Box info as above-------------------"<<endl;
}

//intial tags: fluid tag = the object num: 1, 2, 3,...; non-fluid: 0
void Initializer::setObjectTag() {
	size_t tmp=0;
	//for(size_t index=0; index<m_iCapacity; index++) m_vObjectTag[index] = 0; // init object tag
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) { 
		size_t num = m_vFluidBoundingBox[p]->getNumber();
		// set the fluid object tag
		for(size_t index=m_iFluidStartIndex+tmp; index<m_iFluidStartIndex+tmp+num; index++) {
			m_vObjectTag[index] = p+1;
		}
		tmp += num;
	}
	//cout<<"-------Initializer::setObjectTag()-------"<<endl;
	//for(size_t index=0; index<m_iCapacity; index++) cout<<"m_vObjectTag["<<index<<"]="<<m_vObjectTag[index]<<endl; 
	//cout<<"-----------------------------------------"<<endl;
}


void Initializer::initParticleDataMemory() {
	
	m_vPositionX = new double[m_iCapacity];
	m_vPositionY = new double[m_iCapacity];
	m_vPositionZ = new double[m_iCapacity];
	fill_n(m_vPositionX,m_iCapacity,0);
	fill_n(m_vPositionY,m_iCapacity,0);
	if(m_iDimension==2) fill_n(m_vPositionZ,m_iCapacity,0);	

	m_vVelocityU = new double[m_iCapacity];
	m_vVelocityV = new double[m_iCapacity];
	fill_n(m_vVelocityU,m_iCapacity,0);
	fill_n(m_vVelocityV,m_iCapacity,0);
	if(m_iDimension==3) {
		m_vVelocityW = new double[m_iCapacity];
		fill_n(m_vVelocityW,m_iCapacity,0);
	}
	else	m_vVelocityW = nullptr; 

	m_vVolume = new double[m_iCapacity];
	fill_n(m_vVolume,m_iCapacity,0);
	
	m_vPressure = new double[m_iCapacity];
	fill_n(m_vPressure,m_iCapacity,0);	
	
	m_vSoundSpeed = new double[m_iCapacity];
	fill_n(m_vSoundSpeed,m_iCapacity,0);
	
	m_vObjectTag = new int[m_iCapacity];
	fill_n(m_vObjectTag,m_iCapacity,0);
	
}


size_t Initializer::initGeometryAndStateOnHexPacking(bool saveData, const string& kind) {

	const double h_r = 0.5*m_fInitParticleSpacing;
	
	vector<BoundingBox*>& boxes = kind=="fluid"? m_vFluidBoundingBox:m_vBoundaryBoundingBox;
	const vector<Geometry*>& objs = kind=="fluid"? m_vFluidObj:m_vBoundaryObj;
	const vector<State*>& states = m_vFluidObjState;

	size_t numParticle = 0;

	if(m_iDimension==2) {

		for(size_t p=0; p<objs.size(); p++) {
		
			size_t numParticleI = 0;

			double xmin = boxes[p]->getXmin();	
			double xmax = boxes[p]->getXmax();
			double ymin = boxes[p]->getYmin();	
			double ymax = boxes[p]->getYmax();
				
			HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	

			// compute the location of particles
			// and compute the number of fluid particles
			for(size_t j=m0; j<=m1; j++) { 
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2D.computeX(0,k);
						double y = hex2D.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; // level function of fluid object	
						if(saveData) {
							m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;
							if(kind=="fluid") {
								m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
								m_vPressure[numParticle] = states[p]->pressure(x,y,0);
								double tmpZ;
								states[p]->velocity(x,y,0,
									m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
								m_vSoundSpeed[numParticle] = 
									m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							}
						}	
						numParticle++;
						numParticleI++;
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2D.computeX(1,k);
						double y = hex2D.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; 	
						if(saveData) {
							m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;
							if(kind=="fluid") {
								m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
								m_vPressure[numParticle] = states[p]->pressure(x,y,0);
								double tmpZ;
								states[p]->velocity(x,y,0,
									m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
								m_vSoundSpeed[numParticle] = 
									m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							}
						}	
						numParticle++;
						numParticleI++;
					}
				}
			}
			boxes[p]->setNumber(numParticleI);
		}	
	}
	else if(m_iDimension==3) {
	
		for(size_t p=0; p<objs.size(); p++) {
			size_t numParticleI = 0;

			double xmin = boxes[p]->getXmin();	
			double xmax = boxes[p]->getXmax();
			double ymin = boxes[p]->getYmin();	
			double ymax = boxes[p]->getYmax();
			double zmin = boxes[p]->getZmin();	
			double zmax = boxes[p]->getZmax();

			HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
			//get parameters of hexagonal packing
			size_t l0,l1;
			size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
			size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
			hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
								n0_odd, n1_odd, n0_even, n1_even, 
								nn0_odd, nn1_odd, nn0_even, nn1_even);	
			
			// compute the location of particles
			// and compute the number of fluid particles
			for(size_t i=l0; i<=l1; i++) { 
				if((i+1)%2 != 0) { //odd-numbered layers
					for(size_t j=m0_odd; j<=m1_odd; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows 
							for(size_t k=n0_odd; k<=n1_odd; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue;	
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									if(kind=="fluid") {
										m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
										m_vPressure[numParticle] = states[p]->pressure(x,y,z);
										states[p]->velocity(x,y,z,
											m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
										m_vSoundSpeed[numParticle] = 
											m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									}
								}	
								numParticle++;
								numParticleI++;
								
							}
						} 
						else{ //even-numbered rows
							for(size_t k=n0_even; k<=n1_even; k++) {
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue;	
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									if(kind=="fluid") {
										m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
										m_vPressure[numParticle] = states[p]->pressure(x,y,z);
										states[p]->velocity(x,y,z,
											m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
										m_vSoundSpeed[numParticle] = 
											m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									}
								}
								numParticle++;
								numParticleI++;
								
							}
						}
					}
						
				} 
				else { //even-numbered layers
					for(size_t j=m0_even; j<=m1_even; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows
							for(size_t k=nn0_odd; k<=nn1_odd; k++) { 
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue;
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									if(kind=="fluid") {
										m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
										m_vPressure[numParticle] = states[p]->pressure(x,y,z);
										states[p]->velocity(x,y,z,
											m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
										m_vSoundSpeed[numParticle] = 
											m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									}
								}	
								numParticle++;
								numParticleI++;	
							}
						} 
						else { //even-numbered rows
							for(size_t k=nn0_even; k<=nn1_even; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue; 	
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									if(kind=="fluid") {
										m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
										m_vPressure[numParticle] = states[p]->pressure(x,y,z);
										states[p]->velocity(x,y,z,
											m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
										m_vSoundSpeed[numParticle] = 
											m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									}
								}
								numParticle++;
								numParticleI++;	
							}
						}
					}	
				}    
			}
			boxes[p]->setNumber(numParticleI);
		}	
	
	}

	return numParticle;
}


void Initializer::computeNumParticleWithinSearchRadius() {
	
	size_t result = 0;

	double h_r = 0.5*m_fInitParticleSpacing;

	if(m_iDimension==2) {
		
		// a small rectangular space 
		double xmin = -1.5*m_fNeiSearchRadius;
		double xmax =  1.5*m_fNeiSearchRadius;
		double ymin = -1.5*m_fNeiSearchRadius;
		double ymax =  1.5*m_fNeiSearchRadius;		
		
		HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
		// get parameters of hexagonal packing
		size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
		hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	
			
		//compute the location of particles 
		for(size_t j=m0; j<=m1; j++) { 
			if((j+1)%2 != 0) { //odd-numbered rows 
				for(size_t k=n0_odd; k<=n1_odd; k++) { 
					double x = hex2D.computeX(0,k);
					double y = hex2D.computeY(j);	
					if(sqrt(x*x+y*y)<=m_fNeiSearchRadius) result++;
				}
			} 
			else{ //even-numbered rows
				for(size_t k=n0_even; k<=n1_even; k++) {
					double x = hex2D.computeX(1,k);
					double y = hex2D.computeY(j);	
					if(sqrt(x*x+y*y)<=m_fNeiSearchRadius) result++;
				}
			}
		}
	}
	else if(m_iDimension==3) {
		
		// a small cube 
		double xmin = -1.5*m_fNeiSearchRadius;
		double xmax =  1.5*m_fNeiSearchRadius;
		double ymin = -1.5*m_fNeiSearchRadius;
		double ymax =  1.5*m_fNeiSearchRadius;
		double zmin = -1.5*m_fNeiSearchRadius;
		double zmax =  1.5*m_fNeiSearchRadius;
		
		HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
		//get parameters of hexagonal packing
		size_t l0,l1;
		size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
		size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
		hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
						    n0_odd, n1_odd, n0_even, n1_even, 
							nn0_odd, nn1_odd, nn0_even, nn1_even);	
		
		// compute the location of particles
		// and compute the number of fluid particles
		for(size_t i=l0; i<=l1; i++) { 
			if((i+1)%2 != 0) { //odd-numbered layers
				for(size_t j=m0_odd; j<=m1_odd; j++) { 
					if((j+1)%2 != 0) { //odd-numbered rows 
						for(size_t k=n0_odd; k<=n1_odd; k++) {
							double x = hex3D.computeX(0,k);
							double y = hex3D.computeY(0,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					} 
					else{ //even-numbered rows
						for(size_t k=n0_even; k<=n1_even; k++) {
							double x = hex3D.computeX(1,k);
							double y = hex3D.computeY(0,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					}
				}	
			} 
			else { //even-numbered layers
				for(size_t j=m0_even; j<=m1_even; j++) { 
					if((j+1)%2 != 0) { //odd-numbered rows
						for(size_t k=nn0_odd; k<=nn1_odd; k++) { 
							double x = hex3D.computeX(1,k);
							double y = hex3D.computeY(1,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					} 
					else { //even-numbered rows
						for(size_t k=nn0_even; k<=nn1_even; k++) {
							double x = hex3D.computeX(0,k);
							double y = hex3D.computeY(1,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					}
				}	
			}    
		}	
	}

	m_iNumParticleWithinSearchRadius = result;

	cout<<"-------Initializer::computeNumParticleWithinSearchRadius()-------"<<endl;
	cout<<"m_iNumParticleWithinSearchRadius="<<m_iNumParticleWithinSearchRadius<<endl;
	cout<<"-----------------------------------------------------------------"<<endl<<endl;
}


////////////////////////////////////////////////////////////////////////////////////////
// End of Initializer
////////////////////////////////////////////////////////////////////////////////////////
