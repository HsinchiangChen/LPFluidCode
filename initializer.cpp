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

Initializer::Initializer(const string& inputfileName, const string& debugfileName, bool ifDebug)
: m_sDebugfileName(debugfileName), m_iIfDebug(ifDebug) {
	
	ifstream ifs(inputfileName);
	this->debug.open(m_sDebugfileName);

	vector<string> lines;
	string s;
	debug<<"--------------The original input file-------------------"<<endl;
	while(getline(ifs, s)) {
		debug<<s<<endl;
		lines.push_back(s);
	}
	debug<<"--------------End of original input file----------------"<<endl<<endl;

	istringstream iss;
	size_t i = 0; // input file line number
	
	debug<<"--------------Variables assigned by input file-------------------"<<endl;
	iss.str(lines[i++]);
	iss>>m_iNumThreads;
	debug<<"m_iNumThreads = "<<m_iNumThreads<<endl;

	iss.str(lines[i++]);
	iss>>m_fStartTime;
	debug<<"m_fStartTime = "<<m_fStartTime<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fEndTime;
	debug<<"m_fEndTime = "<<m_fEndTime<<endl;

	iss.str(lines[i++]);
	iss>>m_fWriteTimeInterval;
	debug<<"m_fWriteTimeInterval = "<<m_fWriteTimeInterval<<endl;

	iss.str(lines[i++]);
	iss>>m_fCFLCoeff;
	debug<<"m_fCFLCoeff = "<<m_fCFLCoeff<<endl;

	iss.str(lines[i++]);
	iss>>m_iDimension;
	debug<<"m_iDimension = "<<m_iDimension<<endl;

	iss.str(lines[i++]);
	iss>>m_iFluidObjNum;
	debug<<"m_iFluidObjNum= "<<m_iFluidObjNum<<endl;
 
	for(int j=0; j<m_iFluidObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		m_vFluidObjNames.push_back(tmpS);
		debug<<"m_vFluidObjectNames"<<j<<"="<<m_vFluidObjNames[j]<<endl;
	}
	
	for(int j=0; j<m_iFluidObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		m_vFluidObjStateNames.push_back(tmpS);
		debug<<"m_vFluidObjectStateNames"<<j<<"="<<m_vFluidObjStateNames[j]<<endl;
	}

	iss.str(lines[i++]);
	iss>>m_iBoundaryObjNum;
	debug<<"m_iBoundaryObjNum= "<<m_iBoundaryObjNum<<endl;
	
	if(m_iDimension == 1) { // In the 1D case only boundary type is relevant (only one type is allowed)
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		m_vBoundaryObjTypes.push_back(tmpS);
		debug<<"m_vBoundaryObjectTypes[0]="<<m_vBoundaryObjTypes[0]<<endl;
	}
	else { // In 2D & 3D both boundary geometry name and boundary type should be specified
		for(int j=0; j<m_iBoundaryObjNum; j++) {
			string tmpS;
		
			iss.str(lines[i++]);
			iss>>tmpS;
			m_vBoundaryObjNames.push_back(tmpS);
			debug<<"m_vBoundaryObjectNames"<<j<<"="<<m_vBoundaryObjNames[j]<<endl;
			
			iss.str(lines[i++]);
			iss>>tmpS;
			m_vBoundaryObjTypes.push_back(tmpS);
			debug<<"m_vBoundaryObjectTypes"<<j<<"="<<m_vBoundaryObjTypes[j]<<endl;
		}
	}

	iss.str(lines[i++]);
	iss>>m_iRandomDirSplitOrder;
	debug<<"m_iRandomDirSplitOrder = "<<m_iRandomDirSplitOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iLPFOrder;
	debug<<"m_iLPFOrder = "<<m_iLPFOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iEOSChoice;
	debug<<"m_iEOSChoice = "<<m_iEOSChoice<<endl;

	iss.str(lines[i++]);
	iss>>m_fGamma;
	debug<<"m_fGamma = "<<m_fGamma<<endl;

	iss.str(lines[i++]);
	iss>>m_fPinf;
	debug<<"m_fPinf = "<<m_fPinf<<endl;

	iss.str(lines[i++]);
	iss>>m_fEinf;
	debug<<"m_fEinf = "<<m_fEinf<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fInitParticleSpacing;
	debug<<"m_fInitParticleSpacing = "<<m_fInitParticleSpacing<<endl;	

	iss.str(lines[i++]);
	iss>>m_fGravity;
	debug<<"m_fGravity = "<<m_fGravity<<endl;

	iss.str(lines[i++]);
	iss>>m_iMovingBoxForGhostParticle;
	debug<<"m_iMovingBoxForGhostParticle = "<<m_iMovingBoxForGhostParticle<<endl;	
	
	iss.str(lines[i++]);
	iss>>m_iUseLimiter;
	debug<<"m_iUseLimiter = "<<m_iUseLimiter<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fThresholdP;
	debug<<"m_fThresholdP = "<<m_fThresholdP<<endl;

	debug<<"--------------End Variables assigned by input file----------------"<<endl<<endl;


	// create EOS object and get a pointer to it
	if(m_iEOSChoice==1) m_pEOS = new PolytropicGasEOS(m_fGamma);
	else if(m_iEOSChoice==2) m_pEOS = new StiffPolytropicGasEOS(m_fGamma,m_fPinf,m_fEinf);
	else {
		debug<<"The EOS does not exist!!!"<<endl;
		assert(false);
	}

	// create fluid objects and get pointers to them
	for(int j=0; j<m_iFluidObjNum; j++) 
		m_vFluidObj.push_back(GeometryFactory::instance().createGeometry(m_vFluidObjNames[j]));
	
	// create fluid state objects and get pointers to them
	for(int j=0; j<m_iFluidObjNum; j++) 
		m_vFluidObjState.push_back(StateFactory::instance().createState(m_vFluidObjStateNames[j]));

	
	if(m_iDimension != 1) { // Only do these when in 2D & 3D	
		// create boundary objects and get pointers to them
		for(int j=0; j<m_iBoundaryObjNum; j++) 
			m_vBoundaryObj.push_back(GeometryFactory::instance().createGeometry(m_vBoundaryObjNames[j]));
		
		// bounding box of fluid and boundary objects
		computeInitBoundingBox();
	}	
	
	initGeometryAndState();	
	
	if(m_iDimension != 1) { // Only do these when in 2D & 3D
		setBoundingBoxStartIndex();
		setObjectTag();
	}

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

		debug<<"m_vFluidBoundingBox("<<i<<"):"<<endl;
		debug<<"m_fXmin="<<m_vFluidBoundingBox[i]->getXmin()<<endl;
		debug<<"m_fXmax="<<m_vFluidBoundingBox[i]->getXmax()<<endl;
		debug<<"m_fYmin="<<m_vFluidBoundingBox[i]->getYmin()<<endl;
		debug<<"m_fYmax="<<m_vFluidBoundingBox[i]->getYmax()<<endl;
		debug<<"m_fZmin="<<m_vFluidBoundingBox[i]->getZmin()<<endl;
		debug<<"m_fZmax="<<m_vFluidBoundingBox[i]->getZmax()<<endl;
		debug<<"-----------------"<<endl;
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

		debug<<"m_vBoundaryBoundingBox("<<i<<"):"<<endl;
		debug<<"m_fXmin="<<m_vBoundaryBoundingBox[i]->getXmin()<<endl;
		debug<<"m_fXmax="<<m_vBoundaryBoundingBox[i]->getXmax()<<endl;
		debug<<"m_fYmin="<<m_vBoundaryBoundingBox[i]->getYmin()<<endl;
		debug<<"m_fYmax="<<m_vBoundaryBoundingBox[i]->getYmax()<<endl;
		debug<<"m_fZmin="<<m_vBoundaryBoundingBox[i]->getZmin()<<endl;
		debug<<"m_fZmax="<<m_vBoundaryBoundingBox[i]->getZmax()<<endl;
		debug<<"-----------------"<<endl;
	}
	
}

void Initializer::setParams() {
	debug<<"-------------------Initializer::setParams()-----------------------"<<endl;	
	
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
				
		m_iNumRow2ndOrder = 20; 
		m_iNumRow1stOrder = 5; 
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
	else { // 1D
		m_iMaxNeighbourNum = m_iLPFOrder==1? 2:4; 
		m_iMaxNeighbourNumInOneDir = m_iLPFOrder==1? 1:2;
				
		m_iNumRow2ndOrder = 2; 
		m_iNumRow1stOrder = 1; 
		m_iNumCol2ndOrder = 2; 
		m_iNumCol1stOrder = 1;	
	
	}

	// (eos chocie=1:poly 2:spoly)	
	if(m_iEOSChoice==1) 
		m_fInvalidPressure = 0;
	else if(m_iEOSChoice==2) 	
		m_fInvalidPressure = -m_fPinf;
	m_fInvalidVolume = 0;
	
	m_iTreeDepth = 5;
	
	debug<<"m_fNeiSearchRadius = "<<m_fNeiSearchRadius<<endl;
	debug<<"m_iNumParticleWithinSearchRadius = "<<m_iNumParticleWithinSearchRadius<<endl;
	debug<<"m_iMaxNeighbourNum = "<<m_iMaxNeighbourNum<<endl;
	debug<<"m_iMaxNeighbourNumInOneDir = "<<m_iMaxNeighbourNumInOneDir<<endl;
	debug<<"m_iNumRow2ndOrder = "<<m_iNumRow2ndOrder<<endl;
	debug<<"m_iNumRow1stOrder = "<<m_iNumRow1stOrder<<endl;
	debug<<"m_iNumCol2ndOrder = "<<m_iNumCol2ndOrder<<endl;
	debug<<"m_iNumCol1stOrder = "<<m_iNumCol1stOrder<<endl;	
	debug<<"m_fInvalidPressure = "<<m_fInvalidPressure<<endl;
	debug<<"m_fInvalidVolume = "<<m_fInvalidVolume<<endl;
	debug<<"m_iTreeDepth = "<<m_iTreeDepth<<endl;		
	debug<<"-------------------End Initializer::setParams()-------------------"<<endl<<endl;
	
}


void Initializer::initGeometryAndState() {
	
	debug<<"-----------------Initializer::initGeometryAndState()-------------------"<<endl;

	if(m_iDimension == 1) {
		
		bool saveData = false;
		// Compute the number of particles to initialize the particle array	
		m_iBoundaryNum = 0; 
		m_iFluidNum = initGeometryAndState1D(saveData);
		
		m_iCapacity = m_iFluidNum;
		
		// use m_iCapacity to create memory
		initParticleDataMemory();
		
		saveData = true;
		// Really initialize particle location and state this time	
		initGeometryAndState1D(saveData);	
		
		// NOT relevant in 1D; In 1D only m_iFluidNum is used
		m_iFluidStartIndex = 0;
		m_iBoundaryStartIndex = 0;
		m_iGhostStartIndex = 0;

		debug<<"m_iFluidNum = "<<m_iFluidNum<<endl;
		//debug<<"m_iBoundaryNum = "<<m_iBoundaryNum<<endl;
		debug<<"m_iCapacity = "<<m_iCapacity<<endl;
		//debug<<"m_iFluidStartIndex = "<<m_iFluidStartIndex<<endl;
		//debug<<"m_iBoundaryStartIndex = "<<m_iBoundaryStartIndex<<endl;
		//debug<<"m_iGhostStartIndex = "<<m_iGhostStartIndex<<endl;

	}
	else {
		bool saveData = false;
		// Compute the number of particles to initialize the particle array
		m_iFluidNum = initGeometryAndStateOnHexPacking(saveData, "fluid");
		m_iBoundaryNum = initGeometryAndStateOnHexPacking(saveData, "boundary");
		
		m_iCapacity = (size_t)(1.5*(m_iFluidNum+m_iBoundaryNum));
		
		// use m_iCapacity to create memory
		initParticleDataMemory();
		
		saveData = true;
		// Really initialize particle location and state this time
		initGeometryAndStateOnHexPacking(saveData, "fluid");
		initGeometryAndStateOnHexPacking(saveData, "boundary");

		// assign the start index of fluid, boundary, and ghost particles
		// this variables are only relevant in 2D & 3D; In 1D only m_iFluidNum is used
		m_iFluidStartIndex = 0;
		m_iBoundaryStartIndex = m_iFluidStartIndex + m_iFluidNum;
		m_iGhostStartIndex = m_iFluidStartIndex + m_iFluidNum + m_iBoundaryNum;
			
		debug<<"m_iFluidNum = "<<m_iFluidNum<<endl;
		debug<<"m_iBoundaryNum = "<<m_iBoundaryNum<<endl;
		debug<<"m_iCapacity = "<<m_iCapacity<<endl;
		debug<<"m_iFluidStartIndex = "<<m_iFluidStartIndex<<endl;
		debug<<"m_iBoundaryStartIndex = "<<m_iBoundaryStartIndex<<endl;
		debug<<"m_iGhostStartIndex = "<<m_iGhostStartIndex<<endl;	
	}	
	
	debug<<"-----------------End Initializer::initGeometryAndState()---------------"<<endl<<endl;
	
}


void Initializer::setBoundingBoxStartIndex() {
	size_t tmp=0;	
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) {
		// set the start index of fluid bound box
		m_vFluidBoundingBox[p]->setStartIndex(m_iFluidStartIndex+tmp); 
		size_t num = m_vFluidBoundingBox[p]->getNumber();	
		tmp += num;
		debug<<"m_vFluidBoundingBox["<<p<<"]->m_iStartIndex="<<m_vFluidBoundingBox[p]->getStartIndex()<<endl;
		debug<<"m_vFluidBoundingBox["<<p<<"]->m_iNumber="<<num<<endl;
	}
	tmp=0;	
	for(size_t p=0; p<m_vBoundaryBoundingBox.size(); p++) {
		// set the start index of boundary bound box
		m_vBoundaryBoundingBox[p]->setStartIndex(m_iBoundaryStartIndex+tmp); 
		size_t num = m_vBoundaryBoundingBox[p]->getNumber();	
		tmp += num;
		debug<<"m_vBoundaryBoundingBox["<<p<<"]->m_iStartIndex="<<m_vBoundaryBoundingBox[p]->getStartIndex()<<endl;
		debug<<"m_vBoundaryBoundingBox["<<p<<"]->m_iNumber="<<num<<endl;
	}
	debug<<"-----------------Initialized Bounding Box info as above-------------------"<<endl;
}


//intial tags: fluid tag = the object num: 1, 2, 3,...; non-fluid: 0
void Initializer::setObjectTag() {
	// set the fluid object tag
	size_t tmp=0;
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) { 
		size_t num = m_vFluidBoundingBox[p]->getNumber();	
		for(size_t index=m_iFluidStartIndex+tmp; index<m_iFluidStartIndex+tmp+num; index++) {
			m_vObjectTag[index] = p+1;
		}
		tmp += num;
	}
	 
	//cout<<"-------Initializer::setObjectTag()-------"<<endl;
	//for(size_t index=0; index<m_iCapacity; index++) cout<<"m_vObjectTag["<<index<<"]="<<m_vObjectTag[index]<<endl; 
	//cout<<"-----------------------------------------"<<endl;

	//TODO Set the tags for boundary particles
}


void Initializer::initParticleDataMemory() {
	
	// location
	m_vPositionX = new double[m_iCapacity];
	m_vPositionY = new double[m_iCapacity];
	m_vPositionZ = new double[m_iCapacity];
	fill_n(m_vPositionX,m_iCapacity,0);
	fill_n(m_vPositionY,m_iCapacity,0);
	fill_n(m_vPositionZ,m_iCapacity,0);
	//if(m_iDimension==2) fill_n(m_vPositionZ,m_iCapacity,0);	
	
	// velocity
	m_vVelocityU = new double[m_iCapacity];
	fill_n(m_vVelocityU,m_iCapacity,0);
	
	if(m_iDimension==2 || m_iDimension==3) {
		m_vVelocityV = new double[m_iCapacity];
		fill_n(m_vVelocityV,m_iCapacity,0);
	}
	else m_vVelocityV = nullptr;
		
	if(m_iDimension==3) {
		m_vVelocityW = new double[m_iCapacity];
		fill_n(m_vVelocityW,m_iCapacity,0);
	}
	else m_vVelocityW = nullptr; 
	
	// states
	m_vVolume = new double[m_iCapacity];
	fill_n(m_vVolume,m_iCapacity,0);
	
	m_vPressure = new double[m_iCapacity];
	fill_n(m_vPressure,m_iCapacity,0);	
	
	m_vSoundSpeed = new double[m_iCapacity];
	fill_n(m_vSoundSpeed,m_iCapacity,0);
	
	// object tags
	m_vObjectTag = new int[m_iCapacity];
	fill_n(m_vObjectTag,m_iCapacity,0);
	
}


size_t Initializer::initGeometryAndState1D(bool saveData) {
	
	// alias
	const vector<Geometry*>& objs = m_vFluidObj;
	const vector<State*>& states = m_vFluidObjState;	

	double xmin, xmax, ymin, ymax, zmin, zmax;
	objs[0]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
	size_t numParticle = (size_t)((xmax - xmin)/m_fInitParticleSpacing) + 1;
		
	if(saveData) { 	
		for(size_t i=0; i<numParticle; i++) { 
			double x = xmin + i*m_fInitParticleSpacing;
			m_vPositionX[i] = x; 
			m_vVolume[i]    = 1./states[0]->density(x,0,0);
			m_vPressure[i]  = states[0]->pressure(x,0,0);
			double tmpY, tmpZ;
			states[0]->velocity(x,0,0,m_vVelocityU[i],tmpY,tmpZ);
			m_vSoundSpeed[i] = m_pEOS->getSoundSpeed(m_vPressure[i],1./m_vVolume[i]);
		}
	}

	//cout<<"xmin="<<xmin<<"	xmax="<<xmax<<"	numParticle="<<numParticle<<endl;
	return numParticle;
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

	debug<<"-------Initializer::computeNumParticleWithinSearchRadius()-------"<<endl;
	debug<<"m_iNumParticleWithinSearchRadius="<<m_iNumParticleWithinSearchRadius<<endl;
	debug<<"-----------------------------------------------------------------"<<endl<<endl;
}


////////////////////////////////////////////////////////////////////////////////////////
// End of Initializer
////////////////////////////////////////////////////////////////////////////////////////
