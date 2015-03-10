#include "lp_solver.h"
#include "neighbour_searcher.h"
#include "eos.h"
#include "particle_data.h"
#include "initializer.h"
#include "ls_solver.h"
#include "hexagonal_packing.h"
#include "omp.h"
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of HyperbolicLPSolver
////////////////////////////////////////////////////////////////////////////////////////


HyperbolicLPSolver::HyperbolicLPSolver(const Initializer& init, ParticleData* pData, NeighbourSearcher* ns) {
	
	srand(time(0));

	m_pParticleData = pData; 
	m_pNeighbourSearcher = ns;
	m_pEOS = init.getEOS();
	
	// get parameters from init
	m_iNumThreads = init.getNumThreads();
	m_iDimension = init.getDimension();
	m_iNumPhase = m_iDimension==3? 5:3;
	m_iRandomDirSplitOrder = init.getRandomDirSplitOrder();
	m_iLPFOrder = init.getLPFOrder(); 
	m_iNumRow1stOrder = init.getNumRow1stOrder();
	m_iNumRow2ndOrder = init.getNumRow2ndOrder();
	m_iNumCol1stOrder = init.getNumCol1stOrder(); 	
	m_iNumCol2ndOrder = init.getNumCol2ndOrder();
	m_iMovingBoxForGhostParticle = init.getMovingBoxForGhostParticle(); 
	m_fInitParticleSpacing = init.getInitParticleSpacing();
	m_fGravity = init.getGravity(); 
	m_fInvalidPressure = init.getInvalidPressure(); 
	m_fInvalidVolume = init.getInvalidVolume();
	m_fNeiSearchRadius = init.getNeiSearchRadius(); 
	m_iNumParticleWithinSearchRadius = init.getNumParticleWithinSearchRadius(); 
	m_fContactLength = init.getContactLength();
	
	// all fluid objects should be distant at initialization
	// this variable should always be false in the single fluid object case
	m_iContactAlert = false; 

	// set OpenMP environment
	m_iIfMultiThreads = false;
	if(m_iNumThreads > 1) {
		// set the number of threads
		omp_set_num_threads(min(omp_get_max_threads(), m_iNumThreads));	
		m_iIfMultiThreads = true;	
		
		cout<<"-------HyperbolicLPSolver::HyperbolicLPSolver()-------"<<endl;
		cout<<"m_iNumThreads = "<<m_iNumThreads<<endl;
		cout<<"omp_get_num_procs() = "<<omp_get_num_procs()<<endl;
		cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<" (after omp_set_num_threads())"<<endl;
		cout<<"------------------------------------------------------"<<endl;
	}

	if(m_iDimension==2)
		m_vDirSplitTable = vector<vector<int> >({{0,1,0},{1,0,1}});
	else if(m_iDimension==3)
		m_vDirSplitTable = vector<vector<int> >
		({{0,1,2,1,0},
		  {0,2,1,2,0},
		  {1,0,2,0,1},
		  {1,2,0,2,1},
		  {2,0,1,0,2},
		  {2,1,0,1,2}});
	
	// for completeness initialize to zero
	m_iDirSplitOrder = 0;
	m_fDt = 0;
	
	computeSetupsForNextIteration();
		

}

void HyperbolicLPSolver::computeSetupsForNextIteration() {
		
//testNeighbourSearch();
//assert(false);
	//static int counter=0;

	// update the fluid bounding box based on new location of fluid particles
	if(m_iMovingBoxForGhostParticle) 
		updateFluidBoundingBox();

	// will check for contact alert when having multiple fluid objects and when it is false
	// when m_iContactAlert becomes true it remains true for good
	if(m_pParticleData->m_vFluidBoundingBox.size() > 1 && !m_iContactAlert)
		checkForContactAlert();

	// generate ghost particles (based on fluid particles)
	generateGhostParticleByBruteForceNeighbourSearch();
	//generateGhostParticle();
	//writeResult(0,0,0,m_pParticleData->getTotalNum());


	// search neighbours and put them into the neighbour list
	searchNeighbourForAllParticle();
	//searchNeighbourForAllParticleByBruteForceNeighbourSearch();


	// set neihbour lists in directions right, left, north, and south
	setUpwindNeighbourList();


	// initialize the LPF order (1 or 2) in directions right, left, north, and south
	//writeResult(0,counter++,0,m_pParticleData->getTotalNum());
	resetLPFOrder();


	// set boundary and ghost states 
	setBoundaryPressureAndVelocity(-1);
	setGhostPressureAndVelocity(-1);
	//static int sstep = 0;
	//writeResult(sstep, sstep, 0, m_pParticleData->getTotalNum());
	//sstep++;

	// to determine the dt for next step
	computeMinParticleSpacing();
	computeMaxSoundSpeed();
	computeMaxFluidVelocity();

}

int HyperbolicLPSolver::solve(double dt) {	
	//cout<<"--------------HyperbolicLPSolver::solve()--------------"<<endl;
	
	// dt for this time step 
	m_fDt = dt;

	if(m_iRandomDirSplitOrder) {// get a random order of split
		if(m_iDimension==2) m_iDirSplitOrder = rand()%2; // 0, 1
		else if(m_iDimension==3) m_iDirSplitOrder = rand()%6; // 0, 1, 2, 3, 4, 5
	}
	//cout<<"m_iDirSplitOrder="<<m_iDirSplitOrder<<endl;  

	
	for(int phase=0; phase<m_iNumPhase; ) {
		//cout<<"phase="<<phase<<endl;
		bool phase_success = directionalSplitting(phase);

		// in directionalSplitting, only LPFOrder will be changed
		// so if phase does not succeed, can go back to phase 0 and do again
		// note that phase fails if a particle has LPFOrder in all 4 directions
		if(!phase_success) {phase=0; continue;} 
		
		setBoundaryPressureAndVelocity(phase);
		setGhostPressureAndVelocity(phase);
		
		phase++;
	}
	
	updateFluidState();
	moveFluidParticle();	
	updateFluidVelocity();

	computeSetupsForNextIteration();
	
	//cout<<"-------------------------------------------------------"<<endl;

	return 0;
}

void HyperbolicLPSolver::updateFluidBoundingBox() {
	// PRINT DEBUG INFO
	//cout<<"-------HyperbolicLPSolver::updateFluidBoundingBox()-------"<<endl;

	double *positionX = m_pParticleData->m_vPositionX;
	double *positionY = m_pParticleData->m_vPositionY;
	double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;
	for(size_t p=0; p<fluidBoxes.size(); p++) {
		size_t startIndex = fluidBoxes[p]->getStartIndex();
		size_t number = fluidBoxes[p]->getNumber();
			
		// PRINT DEBUG INFO
		//cout<<"----------------------------------------------------------"<<endl;
		//cout<<"Before update:"<<endl;
		//cout<<"fluidBoxes[p]->getStartIndex()="<<fluidBoxes[p]->getStartIndex()<<endl;
		//cout<<"fluidBoxes[p]->getNumber()="<<fluidBoxes[p]->getNumber()<<endl;
		//cout<<"fluidBoxes[p]->getXmin()="<<fluidBoxes[p]->getXmin()<<endl;
		//cout<<"fluidBoxes[p]->getXmax()="<<fluidBoxes[p]->getXmax()<<endl;
		//cout<<"fluidBoxes[p]->getYmin()="<<fluidBoxes[p]->getYmin()<<endl;
		//cout<<"fluidBoxes[p]->getYmax()="<<fluidBoxes[p]->getYmax()<<endl;
		//cout<<"fluidBoxes[p]->getZmin()="<<fluidBoxes[p]->getZmin()<<endl;
		//cout<<"fluidBoxes[p]->getZmax()="<<fluidBoxes[p]->getZmax()<<endl;
		//cout<<"----------------------------------------------------------"<<endl;	


		// TODO THIS IS A PLACE TO USE OPENMP LOOP LEVEL PARALLELIZATION!!!
		fluidBoxes[p]->setXmin(
		*min_element(positionX+startIndex, positionX+startIndex+number) - 3.*m_fInitParticleSpacing);
		
		fluidBoxes[p]->setXmax(
		*max_element(positionX+startIndex, positionX+startIndex+number) + 3.*m_fInitParticleSpacing);
		
		fluidBoxes[p]->setYmin(
		*min_element(positionY+startIndex, positionY+startIndex+number) - 3.*m_fInitParticleSpacing);
		
		fluidBoxes[p]->setYmax(
		*max_element(positionY+startIndex, positionY+startIndex+number) + 3.*m_fInitParticleSpacing);
		
		if(m_iDimension==3) {
			fluidBoxes[p]->setZmin(
			*min_element(positionZ+startIndex, positionZ+startIndex+number) - 3.*m_fInitParticleSpacing);
			
			fluidBoxes[p]->setZmax(
			*max_element(positionZ+startIndex, positionZ+startIndex+number) + 3.*m_fInitParticleSpacing);
		}

		// PRINT DEBUG INFO
		//cout<<"----------------------------------------------------------"<<endl;
		//cout<<"After update:"<<endl;
		//cout<<"fluidBoxes[p]->getStartIndex()="<<fluidBoxes[p]->getStartIndex()<<endl;
		//cout<<"fluidBoxes[p]->getNumber()="<<fluidBoxes[p]->getNumber()<<endl;
		//cout<<"fluidBoxes[p]->getXmin()="<<fluidBoxes[p]->getXmin()<<endl;
		//cout<<"fluidBoxes[p]->getXmax()="<<fluidBoxes[p]->getXmax()<<endl;
		//cout<<"fluidBoxes[p]->getYmin()="<<fluidBoxes[p]->getYmin()<<endl;
		//cout<<"fluidBoxes[p]->getYmax()="<<fluidBoxes[p]->getYmax()<<endl;
		//cout<<"fluidBoxes[p]->getZmin()="<<fluidBoxes[p]->getZmin()<<endl;
		//cout<<"fluidBoxes[p]->getZmax()="<<fluidBoxes[p]->getZmax()<<endl;
		//cout<<"----------------------------------------------------------"<<endl;	

	}	

}


void HyperbolicLPSolver::checkForContactAlert() {
	cout<<"-------HyperbolicLPSolver::checkForContactAlert()-------"<<endl;	
	
	vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;

	if(fluidBoxes.size() <= 1) return;

	for(size_t p=0; p<fluidBoxes.size(); p++) {
		double xmin1 = fluidBoxes[p]->getXmin();
		double xmax1 = fluidBoxes[p]->getXmax();
		double ymin1 = fluidBoxes[p]->getYmin();
		double ymax1 = fluidBoxes[p]->getYmax();
		double zmin1 = fluidBoxes[p]->getZmin();
		double zmax1 = fluidBoxes[p]->getZmax();
		for(size_t q=p+1; q<fluidBoxes.size(); q++) {
			double xmin2 = fluidBoxes[q]->getXmin();
			double xmax2 = fluidBoxes[q]->getXmax();
			double ymin2 = fluidBoxes[q]->getYmin();
			double ymax2 = fluidBoxes[q]->getYmax();
			double zmin2 = fluidBoxes[q]->getZmin();
			double zmax2 = fluidBoxes[q]->getZmax();
			
			bool x_intersect = (xmax1>=xmax2 && xmin1<=xmax2) || (xmax2>=xmax1 && xmin2<=xmax1);
			bool y_intersect = (ymax1>=ymax2 && ymin1<=ymax2) || (ymax2>=ymax1 && ymin2<=ymax1);
			bool z_intersect = (zmax1>=zmax2 && zmin1<=zmax2) || (zmax2>=zmax1 && zmin2<=zmax1);

			if(x_intersect && y_intersect && z_intersect) { 
				cout<<"fluid object "<<p+1<<" is approaching fluid object "<<q+1<<endl;
				m_iContactAlert = true;
				return;
			}
		}	
	}
	cout<<"All fluid objects are distant so far"<<endl;
	
}





void HyperbolicLPSolver::generateGhostParticle() {
	cout<<"-------HyperbolicLPSolver::generateGhostParticle()-------"<<endl;

	double *positionX = m_pParticleData->m_vPositionX;
	double *positionY = m_pParticleData->m_vPositionY;
	double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	const size_t fluidNum = m_pParticleData->getFluidNum();	
	
	//cout<<"fluidNum = "<<fluidNum<<endl;
	//double startTime = omp_get_wtime();

	// build the search structure on "fluid particles ONLY"
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, fluidStartIndex, fluidNum);
	
	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Building search structure takes %.16g seconds\n", elipsedTime);
	
	// use each fluid bounding box to generate ghost particles
	const vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;
	
	const double h_r = 0.5*m_fInitParticleSpacing;
	const size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;
	const size_t capacity = m_pParticleData->m_iCapacity;
	
	int neiListTemp[maxNeiNum]; // the index of the neighbours of the ghost particle
	double neiListDistTemp[maxNeiNum]; // the distance between the ghost particle and its neighbours //TODO not used
	
	size_t ghostIndex = m_pParticleData->getGhostStartIndex();
	if(m_iDimension==2) {
		//int ttmp = m_pParticleData->getGhostStartIndex();
		for(size_t p=0; p<fluidBoxes.size(); p++) {

			double xmin = fluidBoxes[p]->getXmin();	
			double xmax = fluidBoxes[p]->getXmax();
			double ymin = fluidBoxes[p]->getYmin();	
			double ymax = fluidBoxes[p]->getYmax();
				
			HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);		

			// compute the location of particles	
			#ifdef _OPENMP
			
			#pragma omp parallel  
			{
			
			int tid = omp_get_thread_num();	
			omp_lock_t lockForAddGhostPar; // lock when adding ghost particle into particle array
			omp_init_lock(&lockForAddGhostPar);
			#pragma omp for 
			
			#endif
			for(size_t j=m0; j<=m1; j++) { 
				
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2D.computeX(0,k);
						double y = hex2D.computeY(j);
					
						size_t numNeiFound;
						#ifdef _OPENMP				
						m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
															  neiListTemp, neiListDistTemp, numNeiFound, tid); // output	
						#else
						m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
															  neiListTemp, neiListDistTemp, numNeiFound); // output
						
						#endif
						
						//cout<<numNeiFound<<endl;
						if(!isValidGhostParticle(x,y,0,neiListTemp,numNeiFound,p+1)) continue;		
						if(ghostIndex>=capacity) assert(false); // exceed array size capacity
						
						// This is a valid ghost particle so add this to the particle array and increment the index
						#ifdef _OPENMP
						omp_set_lock(&lockForAddGhostPar);		
						#endif

						positionX[ghostIndex] = x;
						positionY[ghostIndex] = y;	
						ghostIndex++;
						
						#ifdef _OPENMP
						omp_unset_lock(&lockForAddGhostPar);
						#endif
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2D.computeX(1,k);
						double y = hex2D.computeY(j);
						
						size_t numNeiFound;
						#ifdef _OPENMP				
						m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
															  neiListTemp, neiListDistTemp, numNeiFound, tid); // output	
						#else
						m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
															  neiListTemp, neiListDistTemp, numNeiFound); // output
						
						#endif	
						//cout<<numNeiFound<<endl;
						if(!isValidGhostParticle(x,y,0,neiListTemp,numNeiFound,p+1)) continue;		
						if(ghostIndex>=capacity) assert(false); // exceed array size capacity
						
						// This is a valid ghost particle so add this to the particle array and increment the index
						#ifdef _OPENMP
						omp_set_lock(&lockForAddGhostPar);
						#endif

						positionX[ghostIndex] = x;
						positionY[ghostIndex] = y;	
						ghostIndex++;
						
						#ifdef _OPENMP
						omp_unset_lock(&lockForAddGhostPar);
						#endif
					}
				}
			} 
			//cout<<"p="<<p<<" ghostNum="<<ghostIndex-ttmp<<endl;
			//ttmp = ghostIndex;
			
			#ifdef _OPENMP
			omp_destroy_lock(&lockForAddGhostPar);
			} // omp parallel section
			#endif
		}	
	}
	else if(m_iDimension==3) {
	
		for(size_t p=0; p<fluidBoxes.size(); p++) {
			
			double xmin = fluidBoxes[p]->getXmin();	
			double xmax = fluidBoxes[p]->getXmax();
			double ymin = fluidBoxes[p]->getYmin();	
			double ymax = fluidBoxes[p]->getYmax();
			double zmin = fluidBoxes[p]->getZmin();	
			double zmax = fluidBoxes[p]->getZmax();

			HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
			//get parameters of hexagonal packing
			size_t l0,l1;
			size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
			size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
			hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
								n0_odd, n1_odd, n0_even, n1_even, 
								nn0_odd, nn1_odd, nn0_even, nn1_even);	
			
			#ifdef _OPENMP

			#pragma omp parallel 
			{
			
			int tid = omp_get_thread_num();
			cout<<"tid="<<tid<<endl;
			omp_lock_t lockForAddGhostPar; // lock when adding ghost particle into particle array
			omp_init_lock(&lockForAddGhostPar);			
			#pragma omp for 

			#endif
			for(size_t i=l0; i<=l1; i++) { // compute the location of particles
				if((i+1)%2 != 0) { //odd-numbered layers
					for(size_t j=m0_odd; j<=m1_odd; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows 
							for(size_t k=n0_odd; k<=n1_odd; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);	
								
								size_t numNeiFound;
								#ifdef _OPENMP
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound,tid); // output	
								#else
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound); // output
								#endif
								
								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostIndex>=capacity) assert(false); // exceed array size capacity
								
								// This is a valid ghost particle so add this to the particle array and increment the index
								#ifdef _OPENMP
								omp_set_lock(&lockForAddGhostPar);		
								#endif
								
								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;
								
								#ifdef _OPENMP
								omp_unset_lock(&lockForAddGhostPar);
								#endif

							}
						} 
						else{ //even-numbered rows
							for(size_t k=n0_even; k<=n1_even; k++) {
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								#ifdef _OPENMP
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound,tid); // output	
								#else
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound); // output
								#endif	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostIndex>=capacity) assert(false); // exceed array size capacity
								
								// This is a valid ghost particle so add this to the particle array and increment the index
								#ifdef _OPENMP
								omp_set_lock(&lockForAddGhostPar);		
								#endif

								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;

								#ifdef _OPENMP
								omp_unset_lock(&lockForAddGhostPar);
								#endif

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
								
								size_t numNeiFound;
								#ifdef _OPENMP
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound,tid); // output	
								#else
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound); // output
								#endif	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostIndex>=capacity) assert(false); // exceed array size capacity
								
								
								// This is a valid ghost particle so add this to the particle array and increment the index
								#ifdef _OPENMP
								omp_set_lock(&lockForAddGhostPar);		
								#endif
								
								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;

								#ifdef _OPENMP
								omp_unset_lock(&lockForAddGhostPar);
								#endif

							}
						} 
						else { //even-numbered rows
							for(size_t k=nn0_even; k<=nn1_even; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								#ifdef _OPENMP
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound,tid); // output	
								#else
								m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																	  neiListTemp,neiListDistTemp,numNeiFound); // output
								#endif	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostIndex>=capacity) assert(false); // exceed array size capacity
								
								// This is a valid ghost particle so add this to the particle array and increment the index
								#ifdef _OPENMP
								omp_set_lock(&lockForAddGhostPar);		
								#endif

								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;	

								#ifdef _OPENMP
								omp_unset_lock(&lockForAddGhostPar);
								#endif

							}
						}
					}	
				}    
			}
			
			#ifdef _OPENMP
			omp_destroy_lock(&lockForAddGhostPar);
			} //omp parallel section
			#endif

		}	
	
	}
	
	m_pParticleData->m_iGhostNum = ghostIndex - m_pParticleData->m_iGhostStartIndex;
	m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
								   m_pParticleData->m_iBoundaryNum + 
								   m_pParticleData->m_iGhostNum;
	
	
	cout<<"m_pParticleData->m_iGhostStartIndex="<<m_pParticleData->m_iGhostStartIndex<<endl;
	cout<<"m_pParticleData->m_iGhostNum="<<m_pParticleData->m_iGhostNum<<endl;
	cout<<"m_pParticleData->m_iTotalNum="<<m_pParticleData->m_iTotalNum<<endl;
	cout<<"---------------------------------------------------------"<<endl;
}


//TODO
bool HyperbolicLPSolver::isValidGhostParticle(double x, double y, double z, int* neiList, size_t numNei, int objectTag) {

	// 1. have at leat one fluid neighbour
	if(numNei < 1) return false; 
	
	// 2. If neighbour comes from other fluid object than this fluid bounding box
	if(m_pParticleData->m_vFluidBoundingBox.size() > 1) {
	//if(m_iContactAlert) {
		for(size_t l=0; l<numNei; l++) {
			if(m_pParticleData->m_vObjectTag[neiList[l]] != objectTag &&
			   m_pParticleData->m_vObjectTag[neiList[l]] != -objectTag) 
				return false;
		}
	}

	// 3. the percentage of fluid is not too large
	if((double)numNei/(double)m_iNumParticleWithinSearchRadius > 0.7) return false;
			
	// 4. the distance from the nearest neighbour is not too small
	double xDiff = m_pParticleData->m_vPositionX[neiList[0]] - x;
	double yDiff = m_pParticleData->m_vPositionY[neiList[0]] - y;
	double zDiff = m_pParticleData->m_vPositionZ[neiList[0]] - z;
	double dist = sqrt(xDiff*xDiff+yDiff*yDiff+zDiff*zDiff);
	if(dist < 0.7*m_fInitParticleSpacing) return false;
	
	//TODO the ObjectTag may be replaced by info in bounding box
	// 4. cannot have neighbours from different fluid objects
	// only do when m_iContactAlert is true
	//if(m_iContactAlert) { 
	//	int firstTag = abs(m_pParticleData->m_vObjectTag[neiList[0]]); // first tag, +-num are from same object
	//	for(size_t l=1; l<numNei; l++) 
	//		if(abs(m_pParticleData->m_vObjectTag[neiList[l]]) != firstTag) return false; 	
	//}

	// passed all 4 criteria
	return true;
}


// TODO
void HyperbolicLPSolver::searchNeighbourForAllParticle() {
	//cout<<"-------HyperbolicLPSolver::searchNeighbourForAllParticle()-------"<<endl;

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t totalNumParticle = m_pParticleData->getTotalNum();
	
	//cout<<"totalNumParticle = "<<totalNumParticle<<endl;
	//double startTime = omp_get_wtime();
	
	// build the search structure
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, totalNumParticle);	
	
	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Building search structure takes %.16g seconds\n", elipsedTime);

	// the entire neighbour list (to be updated in the following loops)
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	size_t boundaryStartIndex = m_pParticleData->getBoundaryStartIndex();
	size_t boundaryEndIndex = boundaryStartIndex + m_pParticleData->getBoundaryNum();
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	size_t ghostEndIndex = ghostStartIndex + m_pParticleData->getGhostNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
	
	//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;	
	//double neiListDistTemp[omp_get_max_threads()][maxNeiNum];	
	
	#ifdef _OPENMP
	#pragma omp parallel
	{
	
	int tid = omp_get_thread_num();
	
	#endif	
	
	double neiListDist[maxNeiNum]; // a temp array for dist between a particle and its neighbours
	size_t numNeiFound;	

	// fluid
	//if(m_iContactAlert) { // multiple fluid objects (need to handle collisions)
	if(m_pParticleData->m_vFluidBoundingBox.size() > 1) {
		
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif

			neighbourListSize[index] = numNeiFound;		    
			changeFluidObjectTagIfInContact(index,numNeiFound,neiListDist);
			changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(index,numNeiFound);
		}

	}
	else { // single fluid object
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif	
			
			neighbourListSize[index] = numNeiFound;		    	
		}
	}
	

	// boundary
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;		
		
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		
		neighbourListSize[index] = numNeiFound;
		changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	}

	//for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;	
	//	size_t numNeiFound;
	//	m_pNeighbourSearcher->searchNeighbour(positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius,
	//	neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
	//	
	//	//neighbourListSize[index] = numNeiFound;	
	//	
	//	// ghost particle take only fluid particles as neighbours
	//	int incr = 0;
	//	for(size_t k=0; k<numNeiFound; k++) {
	//		size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
	//		if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
	//			neighbourList[neiListStartIndex+incr] = neiI;
	//			incr++;
	//		} 
	//	}
	//	neighbourListSize[index] = incr;		    
	//}

	// ghost
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;	
		
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		
		neighbourListSize[index] = numNeiFound;
	//	cout<<"index="<<index<<" numNeiFound="<<numNeiFound<<endl;
	//	for(size_t k=0; k<numNeiFound; k++) {
	//		int neiIndex = neighbourList[neiListStartIndex+k];
	//		cout<<"neiIndex="<<neiIndex<<" m_vObjectTag="<<m_pParticleData->m_vObjectTag[neiIndex]<<endl;
	//		cout<<"dist="<<neiListDist[k]<<endl;
	//	}	
		changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	}
	
	#ifdef _OPENMP
	}
	#endif

//	cout<<"totalNumParticle="<<totalNumParticle<<endl;
//	cout<<"maxNeiNum="<<maxNeiNum<<endl;
//	cout<<"Fluid particles:"<<endl;
//	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
//		
//		size_t neiListStartIndex = index*maxNeiNum;			
//		cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//		cout<<"-----------------------------"<<endl;
//		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
//			cout<<"neiIndex="<<neighbourList[k]<<endl;
//		cout<<"-----------------------------"<<endl;
//	}
//	cout<<"Boundary particles:"<<endl;
//	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
//		
//		size_t neiListStartIndex = index*maxNeiNum;			
//		cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//		cout<<"-----------------------------"<<endl;
//		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
//			cout<<"neiIndex="<<neighbourList[k]<<endl;
//		cout<<"-----------------------------"<<endl;
//	}
//	cout<<"Ghost particles:"<<endl;
//	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
//		
//		size_t neiListStartIndex = index*maxNeiNum;			
//		cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//		cout<<"-----------------------------"<<endl;
//		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
//			cout<<"neiIndex="<<neighbourList[k]<<endl;
//		cout<<"-----------------------------"<<endl;
//	}
	//cout<<"-----------------------------------------------------------------"<<endl;
	

}

//
void HyperbolicLPSolver::changeFluidObjectTagIfInContact(int index, size_t numNeiFound, const double* neiListDist) {
	int* objectTag = m_pParticleData->m_vObjectTag;
	
	// tag < 0  means have already contacted with other fluid object	
	if(objectTag[index] < 0) return;

	// tag == 0 means this is not a fluid particle (This function should not be called by non-fluid particle!)
	if(objectTag[index] == 0) assert(false);

	int *neighbourList = m_pParticleData->m_vNeighbourList;	
	 	
	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;	
	for(size_t k=0; k<numNeiFound; k++) {
		int neiIndex = neighbourList[neiListStartIndex+k];
		double neiDist = neiListDist[k];
		if(objectTag[neiIndex] == 0) continue; // 1. neighbour is not a fluid particle
		if(objectTag[neiIndex] == objectTag[index] || 
		   objectTag[neiIndex] == -objectTag[index]) // 2. neighbour is the same fluid object
			continue; 
		if(neiDist > m_fContactLength) continue; // 3. neighbour is in other fluid object but not within contact length
		
		objectTag[index] = -objectTag[index];
		return; // changed tag so returns - no need to check the remaining neighbours	
	}	
	// the end of function: if is not in contact with other fluid objects then its object tag is not modified
}

 
void HyperbolicLPSolver::changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(int index, size_t numNeiFound) {
	int* objectTag = m_pParticleData->m_vObjectTag;
	
	// only do this for fluid particles NOT in contact with fluid particles from other fluid object
	if(objectTag[index] < 0) return;
	
	// tag == 0 means this is not a fluid particle (This function should not be called by non-fluid particles)
	if(objectTag[index] == 0) assert(false);

	int *neighbourList = m_pParticleData->m_vNeighbourList;	
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;	
	int count = 0;
	for(size_t k=0; k<numNeiFound; k++) {
		int neiIndex = neighbourList[neiListStartIndex+k];
		// only admit neighbours that are 1: non-fluid 2: fluid from the same fluid object
		if(objectTag[neiIndex]==0 || 
		   objectTag[neiIndex] == objectTag[index] ||
		   objectTag[neiIndex] == -objectTag[index]) {  
			neighbourList[neiListStartIndex+count] = neiIndex;
			count++;
		}	
	}	
	// modify neighbourListSize
	neighbourListSize[index] = count;
}

void HyperbolicLPSolver::changeNeighbourhoodToIncludeOnlyFluidNei(int index, size_t numNeiFound) {
	int* objectTag = m_pParticleData->m_vObjectTag;

	if(objectTag[index] != 0) return; // This function should be called by non-fluid particles
	
	int *neighbourList = m_pParticleData->m_vNeighbourList;	
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;	
	int count = 0;
	for(size_t k=0; k<numNeiFound; k++) {
		int neiIndex = neighbourList[neiListStartIndex+k];
		// only admit neighbours that are fluid particles
		if(objectTag[neiIndex] != 0) {  
			neighbourList[neiListStartIndex+count] = neiIndex;
			count++;
		}	
	}
//	cout<<"-------HyperbolicLPSolver::changeNeighbourhoodToIncludeOnlyFluidNei()-------"<<endl;
//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//	cout<<"fluid neighbour count ="<<count<<endl;
//	for(int k=0; k<count; k++) {
//		int neiIndex = neighbourList[neiListStartIndex+k];
//		cout<<"neiIndex="<<neiIndex<<endl;	
//	}
//	cout<<"----------------------------------------------------------------------------"<<endl;
	
	// modify neighbourListSize
	neighbourListSize[index] = count;


}


// a1>a0  1 2 
// a1<a0  3 4
// b1>=b0 1 3 
// b1<b0  2 4 
void setListInOneDir2D(size_t neiIndex, double a0, double a1, double b0, double b1, 
					   int* arr1, size_t& n1, int* arr2, size_t& n2,
					   int* arr3, size_t& n3, int* arr4, size_t& n4) {
	
	if(a1 > a0) {
		//printf("a1 > a0\n");
		//printf("a1=%.16g\n",a1);
		//printf("a0=%.16g\n",a0);
		if(b1 >= b0) arr1[n1++]=neiIndex;
		else arr2[n2++]=neiIndex;
	}
	else if(a1 < a0) { 
		//printf("a1 < a0\n");
		//printf("a1=%.16g\n",a1);
		//printf("a0=%.16g\n",a0);
		if(b1 >= b0) arr3[n3++]=neiIndex;
		else arr4[n4++]=neiIndex;
	}

}


// a1>a0  1 2 3 4
// a1<a0  5 6 7 8
// b1>=b0 1 2 5 6
// b1<b0  3 4 7 8
// c1>=c0 1 3 5 7
// c1<c0  2 4 6 8
void setListInOneDir3D(size_t neiIndex, double a0, double a1, double b0, double b1, double c0, double c1,
					   int* arr1, size_t& n1, int* arr2, size_t& n2,
					   int* arr3, size_t& n3, int* arr4, size_t& n4,
					   int* arr5, size_t& n5, int* arr6, size_t& n6,
					   int* arr7, size_t& n7, int* arr8, size_t& n8) {
	
	if(a1 > a0) { 
		if(b1 >= b0) {
			if(c1 >= c0) arr1[n1++]=neiIndex;
			else arr2[n2++]=neiIndex;
		}
		else {
			if(c1 >= c0) arr3[n3++]=neiIndex;
			else arr4[n4++]=neiIndex;
		}
	
	}
	else if(a1 < a0) { 
		if(b1 >= b0) {
			if(c1 >= c0) arr5[n5++]=neiIndex;
			else arr6[n6++]=neiIndex;
		}
		else {
			if(c1 >= c0) arr7[n7++]=neiIndex;
			else arr8[n8++]=neiIndex;
		}	
	}

}


// helper function of HyperbolicLPSolver::setUpwindNeighbourList()
void setListInOneDir2D(size_t index, size_t maxNeiNumInOneDir, 
					   const int* list1, size_t sz1, const int* list2, size_t sz2,
					   int* upwindNeighbourList, int* upwindNeighbourListSize) { // output
	
	size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
	size_t numInOneDir = 0, n1 = 0, n2 = 0;
	while(numInOneDir < maxNeiNumInOneDir) {
		if(n1 < sz1) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list1[n1++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n2 < sz2) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list2[n2++];
			numInOneDir++;
		}
		if(n1==sz1 && n2==sz2) break;
	}
	upwindNeighbourListSize[index] = numInOneDir;
}


// helper function of HyperbolicLPSolver::setUpwindNeighbourList()
void setListInOneDir3D(size_t index, size_t maxNeiNumInOneDir, 
					   const int* list1, size_t sz1, const int* list2, size_t sz2,
					   const int* list3, size_t sz3, const int* list4, size_t sz4,
					   int* upwindNeighbourList, int* upwindNeighbourListSize) { // output
	
	size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
	size_t numInOneDir = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	while(numInOneDir < maxNeiNumInOneDir) {
		if(n1 < sz1) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list1[n1++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n2 < sz2) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list2[n2++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n3 < sz3) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list3[n3++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n4 < sz4) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list4[n4++];
			numInOneDir++;
		}	
		
		if(n1==sz1 && n2==sz2 && n3==sz3 && n4==sz4) break;
	}
	upwindNeighbourListSize[index] = numInOneDir;
}



// only fluid particles have upwind neighbour list
void HyperbolicLPSolver::setUpwindNeighbourList() {
	
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;		
	const double *positionZ = m_pParticleData->m_vPositionZ;
	if(m_iDimension==2) positionZ = nullptr;

	// get the upwind neighbour lists (to be updated in the following)
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListRight = m_pParticleData->m_vNeighbourListRight;
	int *neighbourListLeft = m_pParticleData->m_vNeighbourListLeft;
	int *neighbourListNorth = m_pParticleData->m_vNeighbourListNorth;
	int *neighbourListSouth = m_pParticleData->m_vNeighbourListSouth;
	int *neighbourListUp, *neighbourListDown;
	if(m_iDimension==3) {
		neighbourListUp = m_pParticleData->m_vNeighbourListUp;
		neighbourListDown = m_pParticleData->m_vNeighbourListDown;
	}	
	// get the size of neighbour lists (to be updated in the following)
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;
	int *neighbourListRightSize = m_pParticleData->m_vNeighbourListRightSize;
	int *neighbourListLeftSize = m_pParticleData->m_vNeighbourListLeftSize;
	int *neighbourListNorthSize = m_pParticleData->m_vNeighbourListNorthSize;
	int *neighbourListSouthSize = m_pParticleData->m_vNeighbourListSouthSize;
	int *neighbourListUpSize, *neighbourListDownSize;
	if(m_iDimension==3) {
		neighbourListUpSize = m_pParticleData->m_vNeighbourListUpSize;
		neighbourListDownSize = m_pParticleData->m_vNeighbourListDownSize;
	}
		
	// get the maximum size of neighbour lists	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
	size_t maxNeiNumInOneDir = m_pParticleData->m_iMaxNeighbourNumInOneDir;

	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	if(m_iDimension==2) {

		int rn[maxNeiNum]; int rs[maxNeiNum]; 
		int ln[maxNeiNum]; int ls[maxNeiNum]; 
		int nr[maxNeiNum]; int nl[maxNeiNum]; 
		int sr[maxNeiNum]; int sl[maxNeiNum]; 
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			size_t rnSz = 0, rsSz = 0, lnSz = 0, lsSz = 0;
			size_t nrSz = 0, nlSz = 0, srSz = 0, slSz = 0;
			
			size_t neiListStartIndex = index*maxNeiNum;
			size_t neiListEndIndex = neiListStartIndex + neighbourListSize[index];
			
			double x0 = positionX[index], y0 = positionY[index];
			for(size_t i=neiListStartIndex; i<neiListEndIndex; i++) {	
				size_t  neiIndex = neighbourList[i];
				double x1 = positionX[neiIndex], y1 = positionY[neiIndex];
				// right or left
				setListInOneDir2D(neiIndex, x0, x1, y0, y1,
								  rn, rnSz, rs, rsSz, ln, lnSz, ls, lsSz);
				
				// north or south
				setListInOneDir2D(neiIndex, y0, y1, x0, x1,
								  nr, nrSz, nl, nlSz, sr, srSz, sl, slSz);
	
			}
				
			//DEBUG INFO
			//cout<<"-------HyperbolicLPSolver::setUpwindNeighbourList()-------"<<endl;
			//cout<<"rnSz="<<rnSz<<endl;
			//for(size_t k=0; k<rnSz; k++) {
			//	cout<<"rn["<<k<<"]="<<rn[k]<<endl;
			//	printf("positionX[neiIndex]=%.16g\n",positionX[rn[k]]);
			//	printf("positionX[index]=%.16g\n",positionX[index]);
			//}
			//cout<<"rsSz="<<rsSz<<endl;
			//for(size_t k=0; k<rsSz; k++) {
			//	cout<<"rs["<<k<<"]="<<rs[k]<<endl;
			//	printf("positionX[neiIndex]=%.16g\n",positionX[rs[k]]);
			//	printf("positionX[index]=%.16g\n",positionX[index]);
			//}	
			//cout<<"lnSz="<<lnSz<<endl;
			//for(size_t k=0; k<lnSz; k++) cout<<"ln["<<k<<"]="<<ln[k]<<endl;
			//cout<<"lsSz="<<lsSz<<endl;
			//for(size_t k=0; k<lsSz; k++) cout<<"ls["<<k<<"]="<<ls[k]<<endl;
			//cout<<"nrSz="<<nrSz<<endl;
			//for(size_t k=0; k<nrSz; k++) cout<<"nr["<<k<<"]="<<nr[k]<<endl;
			//cout<<"nlSz="<<nlSz<<endl;
			//for(size_t k=0; k<nlSz; k++) cout<<"nl["<<k<<"]="<<nl[k]<<endl;
			//cout<<"srSz="<<srSz<<endl;
			//for(size_t k=0; k<srSz; k++) cout<<"sr["<<k<<"]="<<sr[k]<<endl;
			//cout<<"slSz="<<slSz<<endl;
			//for(size_t k=0; k<slSz; k++) cout<<"sl["<<k<<"]="<<sl[k]<<endl;
			//cout<<"----------------------------------------------------------"<<endl;


			// call helper function to set the lists
			//right
			setListInOneDir2D(index, maxNeiNumInOneDir, rn, rnSz, rs, rsSz, 
							  neighbourListRight, neighbourListRightSize); // output
			// left
			setListInOneDir2D(index, maxNeiNumInOneDir, ln, lnSz, ls, lsSz, 
							  neighbourListLeft, neighbourListLeftSize); // output
			//north
			setListInOneDir2D(index, maxNeiNumInOneDir, nr, nrSz, nl, nlSz, 
							  neighbourListNorth, neighbourListNorthSize); // output
			//south
			setListInOneDir2D(index, maxNeiNumInOneDir, sr, srSz, sl, slSz, 
							  neighbourListSouth, neighbourListSouthSize); // output	

		}
	
	}
	else if(m_iDimension==3) {

		int rnu[maxNeiNum]; int rsu[maxNeiNum]; int rnd[maxNeiNum]; int rsd[maxNeiNum]; // right
		int lnu[maxNeiNum]; int lsu[maxNeiNum]; int lnd[maxNeiNum]; int lsd[maxNeiNum]; // left
		int nru[maxNeiNum]; int nlu[maxNeiNum]; int nrd[maxNeiNum]; int nld[maxNeiNum]; // north
		int sru[maxNeiNum]; int slu[maxNeiNum]; int srd[maxNeiNum]; int sld[maxNeiNum]; // south
		int urn[maxNeiNum]; int uln[maxNeiNum]; int urs[maxNeiNum]; int uls[maxNeiNum]; // up
		int drn[maxNeiNum]; int dln[maxNeiNum]; int drs[maxNeiNum]; int dls[maxNeiNum]; // down 
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			size_t rnuSz = 0, rsuSz = 0, rndSz = 0, rsdSz = 0;
			size_t lnuSz = 0, lsuSz = 0, lndSz = 0, lsdSz = 0;
			size_t nruSz = 0, nluSz = 0, nrdSz = 0, nldSz = 0;
			size_t sruSz = 0, sluSz = 0, srdSz = 0, sldSz = 0;
			size_t urnSz = 0, ulnSz = 0, ursSz = 0, ulsSz = 0; 
			size_t drnSz = 0, dlnSz = 0, drsSz = 0, dlsSz = 0; 
			
			size_t neiListStartIndex = index*maxNeiNum;
			size_t neiListEndIndex = neiListStartIndex + neighbourListSize[index];
			
			double x0 = positionX[index], y0 = positionY[index], z0 = positionZ[index];
			for(size_t i=neiListStartIndex; i<neiListEndIndex; i++) {	
				size_t  neiIndex = neighbourList[i];
				double x1 = positionX[neiIndex], y1 = positionY[neiIndex], z1 = positionZ[neiIndex];
				// right or left
				setListInOneDir3D(neiIndex, x0, x1, y0, y1, z0, z1,
								  rnu, rnuSz, rnd, rndSz, rsu, rsuSz, rsd, rsdSz,
								  lnu, lnuSz, lnd, lndSz, lsu, lsuSz, lsd, lsdSz);
				// north or south
				setListInOneDir3D(neiIndex, y0, y1, x0, x1, z0, z1,
								  nru, nruSz, nrd, nrdSz, nlu, nluSz, nld, nldSz,
								  sru, sruSz, srd, srdSz, slu, sluSz, sld, sldSz);
				// up or down
				setListInOneDir3D(neiIndex, z0, z1, x0, x1, y0, y1,
								  urn, urnSz, urs, ursSz, uln, ulnSz, uls, ulsSz,	
								  drn, drnSz, drs, drsSz, dln, dlnSz, dls, dlsSz);	
			}
		

			//DEBUG INFO
			//if(index%100==0 && index<=1000) {
			//	cout<<"-------HyperbolicLPSolver::setUpwindNeighbourList()-------"<<endl;
			//	cout<<"index="<<index<<endl;
			//	cout<<"rnuSz="<<rnuSz<<endl;
			//	cout<<"rsuSz="<<rsuSz<<endl;
			//	cout<<"rndSz="<<rndSz<<endl;
			//	cout<<"rsdSz="<<rsdSz<<endl;	
			//	cout<<"right size = "<<rnuSz+rsuSz+rndSz+rsdSz<<endl;
			//	cout<<"lnuSz="<<lnuSz<<endl;
			//	cout<<"lsuSz="<<lsuSz<<endl;
			//	cout<<"lndSz="<<lndSz<<endl;
			//	cout<<"lsdSz="<<lsdSz<<endl;
			//	cout<<"left size = "<<lnuSz+lsuSz+lndSz+lsdSz<<endl;
			//	cout<<"nruSz="<<nruSz<<endl;
			//	cout<<"nluSz="<<nluSz<<endl;
			//	cout<<"nrdSz="<<nrdSz<<endl;
			//	cout<<"nldSz="<<nldSz<<endl;
			//	cout<<"north size = "<<nruSz+nluSz+nrdSz+nldSz<<endl;
			//	cout<<"sruSz="<<sruSz<<endl;
			//	cout<<"sluSz="<<sluSz<<endl;
			//	cout<<"srdSz="<<srdSz<<endl;
			//	cout<<"sldSz="<<sldSz<<endl;
			//	cout<<"south size = "<<sruSz+sluSz+srdSz+sldSz<<endl;
			//	cout<<"urnSz="<<urnSz<<endl;
			//	cout<<"ulnSz="<<ulnSz<<endl;
			//	cout<<"ursSz="<<ursSz<<endl;
			//	cout<<"ulsSz="<<ulsSz<<endl;
			//	cout<<"up size = "<<urnSz+ulnSz+ursSz+ulsSz<<endl;
			//	cout<<"drnSz="<<drnSz<<endl;
			//	cout<<"dlnSz="<<dlnSz<<endl;
			//	cout<<"drsSz="<<drsSz<<endl;
			//	cout<<"dlsSz="<<dlsSz<<endl;
			//	cout<<"down size = "<<drnSz+dlnSz+drsSz+dlsSz<<endl;
			//	cout<<"----------------------------------------------------------"<<endl;
			//}

			// call helper function to set the lists
			//right
			setListInOneDir3D(index, maxNeiNumInOneDir, rnu, rnuSz, rsu, rsuSz, rnd, rndSz, rsd, rsdSz,
							  neighbourListRight, neighbourListRightSize); // output
			// left
			setListInOneDir3D(index, maxNeiNumInOneDir, lnu, lnuSz, lsu, lsuSz, lnd, lndSz, lsd, lsdSz,
							  neighbourListLeft, neighbourListLeftSize); // output
			//north
			setListInOneDir3D(index, maxNeiNumInOneDir, nru, nruSz, nlu, nluSz, nrd, nrdSz, nld, nldSz,
							  neighbourListNorth, neighbourListNorthSize); // output
			//south
			setListInOneDir3D(index, maxNeiNumInOneDir, sru, sruSz, slu, sluSz, srd, srdSz, sld, sldSz,
							  neighbourListSouth, neighbourListSouthSize); // output
			//up
			setListInOneDir3D(index, maxNeiNumInOneDir, urn, urnSz, uln, ulnSz, urs, ursSz, uls, ulsSz,
							  neighbourListUp, neighbourListUpSize); // output
			//down
			setListInOneDir3D(index, maxNeiNumInOneDir, drn, drnSz, dln, dlnSz, drs, drsSz, dls, dlsSz,
							  neighbourListDown, neighbourListDownSize); // output
		}
	
	}

	// DEBUG (check if the upwind neighbour list is correct)
	checkUpwindNeighbourList();
}



void HyperbolicLPSolver::resetLPFOrder() {
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidNum = m_pParticleData->getFluidNum();
	fill_n(m_pParticleData->m_vLPFOrderRight+fluidStartIndex, fluidNum, m_iLPFOrder);
	fill_n(m_pParticleData->m_vLPFOrderLeft+fluidStartIndex,  fluidNum, m_iLPFOrder);	
	fill_n(m_pParticleData->m_vLPFOrderNorth+fluidStartIndex, fluidNum, m_iLPFOrder);
	fill_n(m_pParticleData->m_vLPFOrderSouth+fluidStartIndex, fluidNum, m_iLPFOrder);
	if(m_iDimension==3) {
		fill_n(m_pParticleData->m_vLPFOrderUp+fluidStartIndex, fluidNum, m_iLPFOrder);
		fill_n(m_pParticleData->m_vLPFOrderDown+fluidStartIndex, fluidNum, m_iLPFOrder);
	}
}



void HyperbolicLPSolver::computeMinParticleSpacing() {
	
	// set particle position pointers (to save typing)
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ;

	// whole neighbour list
	const int *neighbourList = m_pParticleData->m_vNeighbourList;
	const int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;
	
	// initial value
	m_fMinParticleSpacing=numeric_limits<double>::max();	

	size_t startIndex = m_pParticleData->getFluidStartIndex();
	size_t endIndex = startIndex + m_pParticleData->getFluidNum();
	if(m_iDimension==3) {
		for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle

			size_t totalNumNei = neighbourListSize[index];	
			for(size_t i=0; i<totalNumNei; i++) { // for all of its neighbours
				size_t neiIndex = neighbourList[index*maxNeiNum+i];
				if(neiIndex>=startIndex && neiIndex<endIndex) { // find one fluid neighbour, compare, and break
					double dx = positionX[neiIndex] - positionX[index];
					double dy = positionY[neiIndex] - positionY[index];
					double dz = positionZ[neiIndex] - positionZ[index];
					double dist  = sqrt(dx*dx + dy*dy + dz*dz);
					m_fMinParticleSpacing = min(m_fMinParticleSpacing,dist);
					break;
				}
			}
		}
	}
	else if(m_iDimension==2) {
		for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle

			size_t totalNumNei = neighbourListSize[index];	
			for(size_t i=0; i<totalNumNei; i++) { // for all of its neighbours
				size_t neiIndex = neighbourList[index*maxNeiNum+i];
				if(neiIndex>=startIndex && neiIndex<endIndex) { // find one fluid neighbour, compare, and break
					double dx = positionX[neiIndex] - positionX[index];
					double dy = positionY[neiIndex] - positionY[index];
					double dist  = sqrt(dx*dx + dy*dy);
					m_fMinParticleSpacing = min(m_fMinParticleSpacing,dist);
					break;
				}
			}
		}
	}

	// the corner case when no fluid particle has a fluid neighbour
	assert(m_fMinParticleSpacing != numeric_limits<double>::max());
	
	cout<<"-------HyperbolicLPSolver::computeMinParticleSpacing()-------"<<endl;
	cout<<"m_fMinParticleSpacing="<<m_fMinParticleSpacing<<endl;
	cout<<"-------------------------------------------------------------"<<endl;
}


void HyperbolicLPSolver::computeMaxSoundSpeed() {
	
	const double* soundSpeed = m_pParticleData->m_vSoundSpeed;

	//initial value
	//m_fMaxSoundSpeed = numeric_limits<double>::min();
	m_fMaxSoundSpeed = -1;

	size_t startIndex = m_pParticleData->getFluidStartIndex();
	size_t endIndex = startIndex + m_pParticleData->getFluidNum();	
	for(size_t index=startIndex; index<endIndex; index++) // for each fluid particle
		m_fMaxSoundSpeed = max(m_fMaxSoundSpeed,soundSpeed[index]);
	
	//assert(m_fMaxSoundSpeed != numeric_limits<double>::min());	
	assert(m_fMaxSoundSpeed != -1);

	cout<<"-------HyperbolicLPSolver::computeMaxSoundSpeed()-------"<<endl;
	cout<<"m_fMaxSoundSpeed="<<m_fMaxSoundSpeed<<endl;
	cout<<"--------------------------------------------------------"<<endl;
}


void HyperbolicLPSolver::computeMaxFluidVelocity() {
	
	const double* vU = m_pParticleData->m_vVelocityU;
	const double* vV = m_pParticleData->m_vVelocityV;
	const double* vW = (m_iDimension==3)? m_pParticleData->m_vVelocityW:nullptr;
	
	// initial value
	m_fMaxFluidVelocity = -1;

	size_t startIndex = m_pParticleData->getFluidStartIndex();
	size_t endIndex = startIndex + m_pParticleData->getFluidNum();	
	if(m_iDimension==3) {
		for(size_t index=startIndex; index<endIndex; index++) // for each fluid particle
			m_fMaxFluidVelocity = max(m_fMaxFluidVelocity,
				vU[index]*vU[index]+vV[index]*vV[index]+vW[index]*vW[index]);			
	}
	else if(m_iDimension==2) {
		for(size_t index=startIndex; index<endIndex; index++) // for each fluid particle
			m_fMaxFluidVelocity = max(m_fMaxFluidVelocity,vU[index]*vU[index]+vV[index]*vV[index]);		
	}

	assert(m_fMaxFluidVelocity != -1);

	m_fMaxFluidVelocity = sqrt(m_fMaxFluidVelocity);

	cout<<"-------HyperbolicLPSolver::computeMaxFluidVelocity()-------"<<endl;
	cout<<"m_fMaxFluidVelocity="<<m_fMaxFluidVelocity<<endl;
	cout<<"--------------------------------------------------------"<<endl;

}

bool HyperbolicLPSolver::directionalSplitting(int phase) {	
	//cout<<"--------------HyperbolicLPSolver::directionalSplitting()--------------"<<endl;
	
	// determine dir: x(0), y(1), or z(2)
	const int dir = m_vDirSplitTable[m_iDirSplitOrder][phase];
//cout<<"dir="<<dir<<endl;	
	// set particle position pointers (to save typing)
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ;
	
	// set neighbour list pointers by dir (dir=0->right/left, dir=1->north/south, dir=2->up/down)
	const int *neiList0=nullptr, *neiList1=nullptr;
	const int *neiListSize0=nullptr, *neiListSize1=nullptr;
	setNeighbourListPointers(dir, &neiList0, &neiList1, &neiListSize0, &neiListSize1);
	
	// input data pointers
	const double *inVelocity=nullptr, *inPressure=nullptr, *inVolume=nullptr, *inSoundSpeed=nullptr;
	// output data pointers
	double *outVelocity=nullptr, *outPressure=nullptr, *outVolume=nullptr, *outSoundSpeed=nullptr;
	// set data pointers by phase (for temp1 or temp2) and dir(for velocity U or V) 
	setInAndOutDataPointers(phase,dir,&inVelocity,&inPressure,&inVolume,&inSoundSpeed,
							&outVelocity,&outPressure,&outVolume,&outSoundSpeed);
	
	// set local polynomail order pointers
	// dir==0->right(0) & left(1), dir==1->north(0) & south(1), dir==2->up(0) & down(1)
	int *LPFOrder0=nullptr, *LPFOrder1=nullptr;
	vector<int*> LPFOrderOther;	
	setLPFOrderPointers(dir,&LPFOrder0,&LPFOrder1,LPFOrderOther);

	// phase_success will be false if one particle LPF order is zero in all directions (x,y, and z)
	bool phaseSuccess = true;
	
	// gravity is only on y(1) direction 
	double gravity;
	if(m_iDimension==2) gravity = dir==1? m_fGravity:0; // only in y direction
	else if(m_iDimension==3) gravity = dir==2? m_fGravity:0; // only in z direction

	// set real dt:  
	// 2D: phase:	0		1     2 
	//	 		   dt/4   dt/2   dt/4
	// 3D: phase:		0		1		2		3		4
	//				  dt/6    dt/6     dt/3    dt/6    dt/6
	double realDt;
	if(m_iDimension==2) realDt = phase==1? m_fDt/2.:m_fDt/4.;
	else if(m_iDimension==3) realDt = phase==2? m_fDt/3.:m_fDt/6.;
	
	// set the function pointer to compute the A matrix for QR solver
	void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*);
	int offset; // offset is used to get results of computed spatial derivatives from QR solver
	if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; offset=2;}
	else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; offset=3;}
	
	// the coeff before first and second order term during time integration
	double multiplier1st, multiplier2nd;
	if(m_iDimension==2) {multiplier1st=2; multiplier2nd=m_fDt/2.;}
	else if(m_iDimension==3) {multiplier1st=3; multiplier2nd=m_fDt*3./4.;}
	/*	
	if(m_iDimension==2) {
		multiplier1st=2;
		if(phase==1) multiplier2nd = m_fDt/4.;
		else multiplier2nd = m_fDt/8.;
		//if(phase==1) multiplier2nd = m_fDt/2.;
		//else multiplier2nd = m_fDt/4.;
	}
	else if(m_iDimension==3) {
		multiplier1st=3;
		if(phase==2) multiplier2nd = m_fDt/4.;
		else multiplier2nd = m_fDt/8.;
		//if(phase==2) multiplier2nd = m_fDt/2.;
		//else multiplier2nd = m_fDt/4.;
	}
	*/

	// iteration start index
	size_t startIndex = m_pParticleData->getFluidStartIndex();
	// iteration end index (1 after the last element)
	size_t endIndex = startIndex + m_pParticleData->getFluidNum(); 
	
	//cout<<"m_pParticleData->getFluidNum() = "<<m_pParticleData->getFluidNum()<<endl;
	//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;
	//double startTime = omp_get_wtime();	

	// iterate through fluid particles
	#ifdef _OPENMP
	#pragma omp parallel for 
	#endif
	for(size_t index=startIndex; index<endIndex; index++) {
		
		// spatial derivatives (0:right/north/up; 1:left/south/down)
		double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output	
		computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiList0, neiListSize0, 
						  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
		computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiList1, neiListSize1, 
						  LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output	
		
		// update outVolume, outVelocity, outPressure
		timeIntegration(realDt, multiplier1st, multiplier2nd, 
						gravity, inVolume[index], inVelocity[index], inPressure[index],
						inSoundSpeed[index], 
						vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1,
						&outVolume[index], &outVelocity[index], &outPressure[index]); // output	
		
		bool isInvalid = (outPressure[index]<m_fInvalidPressure || outVolume[index]<m_fInvalidVolume);
	
		if(isInvalid) {
			printInvalidState(phase,dir,index,positionX[index],positionY[index],positionZ[index], 
							  vel_d_0,vel_dd_0,p_d_0,p_dd_0,vel_d_1,vel_dd_1,p_d_1,p_dd_1); // input	
			bool atLeastOneNonzeroOrder = lowerLPFOrder(index,LPFOrderOther,LPFOrder0,LPFOrder1);	
			if(atLeastOneNonzeroOrder) { // then use the lowered order to recompute this partcle again
				index--;
				continue;
			}
			else {
				phaseSuccess = false; // mark this phase as false
				outVolume[index]   = inVolume[index];
				outPressure[index] = inPressure[index];
				outVelocity[index] = inVelocity[index];
			}
		}
		
		//cout<<"LPFOrder0["<<index<<"]="<<LPFOrder0[index]<<endl;
		//cout<<"LPFOrder1["<<index<<"]="<<LPFOrder1[index]<<endl;

		// update outSoundSpeed
		outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
	}	

	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Directional splitting takes %.16g seconds\n", elipsedTime);	

	//cout<<"----------------------------------------------------------------------"<<endl;

	// if(phaseSuccess==true) ==> continue to next phase
	// else                   ==> go back to phase 0 (all particles)
	return phaseSuccess;
}


void HyperbolicLPSolver::setNeighbourListPointers(int dir, // input
	const int **neiList0, const int **neiList1, // output
	const int **neiListSize0, const int **neiListSize1) { 
	
	if(dir==0) { // x
		*neiList0 = m_pParticleData->m_vNeighbourListRight; 
		*neiList1 = m_pParticleData->m_vNeighbourListLeft;
		*neiListSize0 = m_pParticleData->m_vNeighbourListRightSize; 
		*neiListSize1 = m_pParticleData->m_vNeighbourListLeftSize;
	}
	else if(dir==1) { // y
		*neiList0 = m_pParticleData->m_vNeighbourListNorth; 
		*neiList1 = m_pParticleData->m_vNeighbourListSouth;
		*neiListSize0 = m_pParticleData->m_vNeighbourListNorthSize; 
		*neiListSize1 = m_pParticleData->m_vNeighbourListSouthSize;
	}
	else if(dir==2) { // z (if m_iDimension==2, dir != 2 for sure)
		*neiList0 = m_pParticleData->m_vNeighbourListUp; 
		*neiList1 = m_pParticleData->m_vNeighbourListDown;
		*neiListSize0 = m_pParticleData->m_vNeighbourListUpSize; 
		*neiListSize1 = m_pParticleData->m_vNeighbourListDownSize;	
	}
	else 
		assert(false);

}


void HyperbolicLPSolver::setInAndOutDataPointers(int phase, int dir,
	const double** inVelocity, const double** inPressure, const double** inVolume, const double** inSoundSpeed, 
	double** outVelocity, double** outPressure, double** outVolume, double** outSoundSpeed) {
	
	// assign pressure, volume, and sound_speed pointers
	if(phase==0) { // input: original output:temp1
		*inPressure   = m_pParticleData->m_vPressure;
		*inVolume     = m_pParticleData->m_vVolume;
		*inSoundSpeed = m_pParticleData->m_vSoundSpeed;
		
		*outPressure   = m_pParticleData->m_vTemp1Pressure;
		*outVolume     = m_pParticleData->m_vTemp1Volume;
		*outSoundSpeed = m_pParticleData->m_vTemp1SoundSpeed;	
	}
	else if(phase==1 || phase==3) { // input:temp1 output:temp2
		*inPressure   = m_pParticleData->m_vTemp1Pressure;
		*inVolume     = m_pParticleData->m_vTemp1Volume;
		*inSoundSpeed = m_pParticleData->m_vTemp1SoundSpeed;
		
		*outPressure   = m_pParticleData->m_vTemp2Pressure;
		*outVolume     = m_pParticleData->m_vTemp2Volume;
		*outSoundSpeed = m_pParticleData->m_vTemp2SoundSpeed;	
	}
	else if(phase==2 || phase==4){ // input:temp2 output: temp1
		*inPressure   = m_pParticleData->m_vTemp2Pressure;
		*inVolume     = m_pParticleData->m_vTemp2Volume;
		*inSoundSpeed = m_pParticleData->m_vTemp2SoundSpeed;
		
		*outPressure   = m_pParticleData->m_vTemp1Pressure;
		*outVolume	  = m_pParticleData->m_vTemp1Volume;
		*outSoundSpeed = m_pParticleData->m_vTemp1SoundSpeed;	
	}	
	else assert(false);
	
	// assign velocity pointers
	if(m_iDimension==2) {
		if(phase==0 || phase==1) { // input: original output:temp2
			if(dir==0) {
				*inVelocity = m_pParticleData->m_vVelocityU;
				*outVelocity = m_pParticleData->m_vTemp2VelocityU;
			}
			else if(dir==1) {
				*inVelocity = m_pParticleData->m_vVelocityV;
				*outVelocity = m_pParticleData->m_vTemp2VelocityV;	
			}	
		}	
		else if(phase==2){ // input:temp2 output: temp1
			if(dir==0) { // u v u
				*inVelocity = m_pParticleData->m_vTemp2VelocityU;
				*outVelocity = m_pParticleData->m_vTemp1VelocityU;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);	
			}
			else if(dir==1) { // v u v
				*inVelocity = m_pParticleData->m_vTemp2VelocityV;
				*outVelocity = m_pParticleData->m_vTemp1VelocityV;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
			}		
		}	
		else assert(false);	
	}
	else if(m_iDimension==3) {
		if(phase==0 || phase==1 || phase==2) { // input: original output:temp1
			if(dir==0) {
				*inVelocity = m_pParticleData->m_vVelocityU;
				*outVelocity = m_pParticleData->m_vTemp1VelocityU;
			}
			else if(dir==1) {
				*inVelocity = m_pParticleData->m_vVelocityV;
				*outVelocity = m_pParticleData->m_vTemp1VelocityV;	
			}
			else if(dir==2) {
				*inVelocity = m_pParticleData->m_vVelocityW;
				*outVelocity = m_pParticleData->m_vTemp1VelocityW;	
			}
		}	
		else if(phase==3){ // input:temp1 output: temp2
			if(dir==0) { 
				*inVelocity = m_pParticleData->m_vTemp1VelocityU;
				*outVelocity = m_pParticleData->m_vTemp2VelocityU;
				// swap pointers so that temp2 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==1) { 
				*inVelocity = m_pParticleData->m_vTemp1VelocityV;
				*outVelocity = m_pParticleData->m_vTemp2VelocityV;
				// swap pointers so that temp2 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==2) { 
				*inVelocity = m_pParticleData->m_vTemp1VelocityW;
				*outVelocity = m_pParticleData->m_vTemp2VelocityW;
				// swap pointers so that temp2 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
			}
		}
		else if(phase==4){ // input:temp2 output: temp1
			if(dir==0) { 
				*inVelocity = m_pParticleData->m_vTemp2VelocityU;
				*outVelocity = m_pParticleData->m_vTemp1VelocityU;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==1) { 
				*inVelocity = m_pParticleData->m_vTemp2VelocityV;
				*outVelocity = m_pParticleData->m_vTemp1VelocityV;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==2) { 
				*inVelocity = m_pParticleData->m_vTemp2VelocityW;
				*outVelocity = m_pParticleData->m_vTemp1VelocityW;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
			}
		}
		else assert(false);	
	}

}


void HyperbolicLPSolver::setLPFOrderPointers(int dir, // input
	int** LPFOrder0, int** LPFOrder1, vector<int*>& LPFOrderOther) { // output
	
	if(dir==0) { // x
		*LPFOrder0 = m_pParticleData->m_vLPFOrderRight; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderLeft;
		
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderNorth); // other directions
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderSouth);
		if(m_iDimension==3) {
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderUp); 
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderDown);
		}
	}
	else if(dir==1) { // y
		*LPFOrder0 = m_pParticleData->m_vLPFOrderNorth; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderSouth;
		
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderRight); // other directions
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderLeft);
		if(m_iDimension==3) {
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderUp); 
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderDown);
		}
	}
	else if(dir==2) { // z
		*LPFOrder0 = m_pParticleData->m_vLPFOrderUp; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderDown;
		
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderRight); // other directions
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderLeft);
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderNorth); 
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderSouth);
	}
	else
		assert(false);

}



void HyperbolicLPSolver::computeSpatialDer(int dir, size_t index, // input 
	int offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*),
	const double* inPressure, const double* inVelocity,
	const int *neighbourList, const int *neighbourListSize,
	int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd) { // output
	
	bool sufficientRank = false;

	// the initial value of the try-error process of finding a good A matrix with sufficient rank
	// these two variables will be increased if necessary in the while loop
	size_t numRow2nd = m_iNumRow2ndOrder;
	size_t numRow1st = m_iNumRow1stOrder;
	
	//cout<<"-------HyperbolicLPSolver::computeSpatialDer()-------"<<endl;
	//cout<<"numRow2nd="<<numRow2nd<<endl;
	//cout<<"numRow1st="<<numRow1st<<endl;

	while(!sufficientRank) {
		
		// decide row and column number based on LPFOrder (and numRow2nd, numRow1st) and the total number of neighbours
		size_t numRow, numCol;
		computeNumRowAndNumColAndLPFOrder(index, neighbourList, neighbourListSize, numRow2nd, numRow1st, // input
										  LPFOrder, &numRow, &numCol); // output	

		if(LPFOrder[index] == 0) { 
			*vel_d  = 0; *vel_dd = 0; *p_d  = 0; *p_dd = 0;
			return;
		}
		
		// compute A
		double A[numRow*numCol];
		(this->*computeA)(index, neighbourList, LPFOrder, numRow, numCol, // input
						  A); // output
		
		double b[numRow]; 
		computeB(index, neighbourList, numRow, inPressure, // input: pressure first
				 b); // output

		QRSolver qrSolver(numRow,numCol,A);
		
		double result[numCol];
		int info = qrSolver.solve(result,b);
		
		if(info!=0) { // then need to recompute A
			//cout<<"rank="<<info<<endl;
			if(LPFOrder[index]==2) numRow2nd++;
			else if(LPFOrder[index]==1) numRow1st++;	
		}	
		else {
			sufficientRank = true;
			
			*p_d = result[dir]; //dir=0 (x), dir=1(y), dir=2(z)
			*p_dd = LPFOrder[index]==2? result[dir+offset]:0;

			computeB(index, neighbourList, numRow, inVelocity, // input: velocity comes second 
					 b); // output (rewrite array b)	
			
			qrSolver.solve(result,b);

			*vel_d = result[dir]; //dir=0 (x), dir=1(y), dir=2(z)
			*vel_dd = LPFOrder[index]==2? result[dir+offset]:0;
		}

	}
	
	//cout<<"p_d="<<*p_d<<endl;
	//cout<<"p_dd="<<*p_dd<<endl;
	//cout<<"vel_d="<<*vel_d<<endl;
	//cout<<"vel_dd="<<*vel_dd<<endl;
	//cout<<"-----------------------------------------------------"<<endl;

}



void HyperbolicLPSolver::computeNumRowAndNumColAndLPFOrder(size_t index, // input
    const int *neighbourList, const int *neighbourListSize, size_t numRow2nd, size_t numRow1st,
	int* LPFOrder, size_t *numRow, size_t *numCol) { // output
	
	// compute the numRow and numCol for matrix A
	size_t totalNeiNum = neighbourListSize[index];

	if(LPFOrder[index] == 2) {
		if(totalNeiNum >= numRow2nd) { // if the total number of neighbour >= current numRow2nd
			(*numRow) = numRow2nd;
			(*numCol) = m_iNumCol2ndOrder; 
		}
		else LPFOrder[index] = 1; // if no enough neighbour -> reduce LPFOrder to 1
	}

	if(LPFOrder[index] == 1) { 
		if(totalNeiNum >= numRow1st) { // if the total number of neighbour >= current numRow1st
			(*numRow) = numRow1st;
			(*numCol) = m_iNumCol1stOrder; 
		}
		else LPFOrder[index] = 0;
	}

	if(LPFOrder[index] == 0) {
		(*numRow) = 0;
		(*numCol) = 0;
	}

	//cout<<"-------HyperbolicLPSolver::computeNumRowAndNumColAndLPFOrder()-------"<<endl;
	//cout<<"index="<<index<<endl;
	//cout<<"totalNeiNum="<<totalNeiNum<<endl;
	//cout<<"LPFOrder="<<LPFOrder[index]<<endl;
	//cout<<"numRow="<<*numRow<<endl;
	//cout<<"numCol="<<*numCol<<endl;
	//cout<<"---------------------------------------------------------------------"<<endl;

}


void HyperbolicLPSolver::computeA2D(size_t index, const int *neighbourList, 
							  	    const int* LPFOrder, size_t numRow, size_t numCol,
								    double *A) { // output 	
	//cout<<"--------------HyperbolicLPSolver::computeA2D()--------------"<<endl;
	//cout<<"index="<<index<<endl;
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
	
	if(LPFOrder[index] == 1) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself 
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index];
			double k = m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index];

			A[i]            = h;
			A[i + 1*numRow] = k;	

			//cout<<A[i]<<"	"<<A[i + 1*numRow]<<endl;
		}
	}
	else if(LPFOrder[index] == 2) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index];
			double k = m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index];
			
			A[i]            = h;
			A[i + 1*numRow] = k;
			A[i + 2*numRow]	= 0.5*h*h;
			A[i + 3*numRow]	= 0.5*k*k;
			A[i + 4*numRow]	= h*k;
			

			//cout<<A[i]<<"	"<<A[i + 1*numRow]<<"	"<<A[i + 2*numRow]<<"	"
			//	<<A[i + 3*numRow]<<"	"<<A[i + 4*numRow]<<endl;
		}
	}	
	
	//cout<<"------------------------------------------------------------"<<endl;
}


void HyperbolicLPSolver::computeA3D(size_t index, const int *neighbourList, 
								    const int* LPFOrder, size_t numRow, size_t numCol,
								    double *A) { // output 	
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
	
	if(LPFOrder[index] == 1) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index];
			double k = m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index];
			double l = m_pParticleData->m_vPositionZ[neiIndex] - m_pParticleData->m_vPositionZ[index];

			A[i]            = h;
			A[i + 1*numRow] = k;
			A[i + 2*numRow] = l;
		}
	}
	else if(LPFOrder[index] == 2) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index];
			double k = m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index];
			double l = m_pParticleData->m_vPositionZ[neiIndex] - m_pParticleData->m_vPositionZ[index];

			A[i]            = h;
			A[i + 1*numRow] = k;
			A[i + 2*numRow] = l;
			A[i + 3*numRow] = 0.5*h*h;
			A[i + 4*numRow] = 0.5*k*k;
			A[i + 5*numRow] = 0.5*l*l;
			A[i + 6*numRow] = h*k;
			A[i + 7*numRow] = h*l;
			A[i + 8*numRow] = k*l;
		}
	}
}


void HyperbolicLPSolver::computeB(size_t index, const int *neighbourList, size_t numRow, const double* inData, 
								  double *b) { // output 	
	
	//cout<<"--------------HyperbolicLPSolver::computeB()--------------"<<endl;
	//cout<<"index="<<index<<endl;

	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
	for(size_t i=0; i<numRow; i++) { 
		int neiIndex = neighbourList[index*maxNeiNum+i];	
		b[i] = inData[neiIndex] - inData[index];

		//cout<<b[i]<<endl;
	}	

	//cout<<"----------------------------------------------------------"<<endl;

}



void HyperbolicLPSolver::timeIntegration(double realDt, double multiplier1st, double multiplier2nd, 
	double gravity, double inVolume, double inVelocity, double inPressure, 
	double inSoundSpeed, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
	double* outVolume, double* outVelocity, double* outPressure) { // output
	
	// TODO
	// Note that this coeff K only works for Poly and Spoly EOS!!!!!!!
	double K = inSoundSpeed*inSoundSpeed/inVolume/inVolume; 

	// Pt
	double Pt1st = -0.5*inVolume*K*(vel_d_0+vel_d_1) + 0.5*inVolume*sqrt(K)*(p_d_0-p_d_1);
	double Pt2nd = -inVolume*inVolume*pow(K,1.5)*(vel_dd_0-vel_dd_1) + inVolume*inVolume*K*(p_dd_0+p_dd_1);
	double Pt = multiplier1st*Pt1st + multiplier2nd*Pt2nd;
	
	//double Pt = multiplier1st*Pt1st + realDt*Pt2nd; ///<wrong coeff
	//double Pt = multiplier1st*Pt1st + 0.5*multiplier2nd*Pt2nd; ///<wrong coeff

	// Vt
	double Vt = -Pt/K;

	// VELt
	double VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0-vel_d_1) - 0.5*inVolume*(p_d_0+p_d_1);
	double VELt2nd = inVolume*inVolume*K*(vel_dd_0+vel_dd_1) - inVolume*inVolume*sqrt(K)*(p_dd_0-p_dd_1);
	double VELt = multiplier1st*VELt1st + multiplier2nd*VELt2nd;

	//double VELt = multiplier1st*VELt1st + realDt*VELt2nd; ///<wrong coeff
	//double VELt = multiplier1st*VELt1st + 0.5*multiplier2nd*VELt2nd; ///<wrong coeff

	// Note that the data pointers of in and out are different!!!!!!!
	(*outVolume)   = inVolume   + realDt*Vt;
	(*outPressure) = inPressure + realDt*Pt;
	(*outVelocity) = inVelocity + realDt*(VELt+gravity);
	
	/*
	double tmp = realDt; // coeff before second order term	
	if(m_iDimension==3) {
		tmp*=1.5;
	}	

	double tmp1 = -0.5*multiplier*inVolume*K;	
	double tmp2 = 0.5*multiplier*inVolume*sqrt(K);
	double tmp3 = -tmp*inVolume*inVolume*pow(K,1.5);
	double tmp4 = tmp*inVolume*inVolume*K;
	double Pt = tmp1*(vel_d_0+vel_d_1) + tmp2*(p_d_0-p_d_1) + tmp3*(vel_dd_0-vel_dd_1) + tmp4*(p_dd_0+p_dd_1);
	
	double Vt = -Pt/K;
	
	tmp1 = 0.5*multiplier*inVolume*sqrt(K);
	tmp2 = -0.5*multiplier*inVolume;
	tmp3 = tmp*inVolume*inVolume*K;
	tmp4 = -tmp*inVolume*inVolume*sqrt(K);
	double VELt = tmp1*(vel_d_0-vel_d_1) + tmp2*(p_d_0+p_d_1) + tmp3*(vel_dd_0+vel_dd_1) + tmp4*(p_dd_0-p_dd_1);
	*/

	//cout<<"-------HyperbolicLPSolver::timeIntegration()-------"<<endl;
	//cout<<"realDt="<<realDt<<endl;
	//cout<<"multiplier="<<multiplier<<endl;
	//cout<<"gravity="<<gravity<<endl;	
	//cout<<"vel_d_0="<<vel_d_0<<endl;
	//cout<<"vel_dd_0="<<vel_dd_0<<endl;
	//cout<<"p_d_0="<<p_d_0<<endl;
	//cout<<"p_dd_0="<<p_dd_0<<endl;
	//cout<<"vel_d_1="<<vel_d_1<<endl;
	//cout<<"vel_dd_1="<<vel_dd_1<<endl;
	//cout<<"p_d_1="<<p_d_1<<endl;
	//cout<<"p_dd_1="<<p_dd_1<<endl;
	//cout<<"Vt="<<Vt<<endl;
	//cout<<"Pt="<<Pt<<endl;
	//cout<<"VELt="<<VELt<<endl;
	//cout<<"inSoundSpeed="<<inSoundSpeed<<endl;
	//cout<<"inVolume="<<inVolume<<endl;
	//cout<<"inPressure="<<inPressure<<endl;
	//cout<<"inVelocity="<<inVelocity<<endl;
	//cout<<"outVolume="<<*outVolume<<endl;
	//cout<<"outPressure="<<*outPressure<<endl;
	//cout<<"outVelocity="<<*outVelocity<<endl;	
	//cout<<"---------------------------------------------------"<<endl;

}


void HyperbolicLPSolver::printInvalidState(int phase, int dir, int index, 
	double positionX, double positionY, double positionZ, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0, 
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1) {
	
	cout<<"---------------------Invalid state-----------------------"<<endl;
	cout<<"phase="<<phase<<", direction="<<dir<<endl;
	cout<<"x["<<index<<"]="<<positionX<<", y["<<index<<"]="<<positionY
		<<", z["<<index<<"]="<<positionZ<<endl;	
	cout<<"---------------First Order Derivatives-------------------"<<endl;
	cout<<"vel_d_0="<<vel_d_0<<", vel_d_1="<<vel_d_1<<endl;
	cout<<"p_d_0="<<p_d_0<<", p_d_1="<<p_d_1<<endl;
	cout<<"---------------Second Order Derivatives------------------"<<endl;
	cout<<"vel_dd_0="<<vel_dd_0<<", vel_dd_1="<<vel_dd_1<<endl;
	cout<<"p_dd_0="<<p_dd_0<<", p_dd_1="<<p_dd_1<<endl;
	cout<<"---------------------------------------------------------"<<endl;
}


bool HyperbolicLPSolver::lowerLPFOrder(int index, const vector<int*>& LPFOrderOther, // input
	int* LPFOrder0, int* LPFOrder1) { // output
	
	//TODO: Note that here I used different algorithms for the return value	
	// If the order in x are both zero, then check the order on the y direction
	// if both north and south derivatives are zero, then return false (then will go back to phase 0)
	// otherwise, return true (will then use the lowered order to redo this particle)
	
	// lower the order in this direction (both sides)
	if(LPFOrder0[index]>0) LPFOrder0[index]--;
	if(LPFOrder1[index]>0) LPFOrder1[index]--;
	
	for(size_t i=0; i<LPFOrderOther.size(); i++) {
		if(LPFOrderOther[i][index]>0) return true;
	}

	if(LPFOrder0[index]==0 && LPFOrder1[index]==0) return false;
	else return true;

}



void HyperbolicLPSolver::setBoundaryPressureAndVelocity(int phase) {
	
	// no need to update for the last phase
	if((m_iDimension==2 && phase==2) || (m_iDimension==3 && phase==4)) return;

	// iteration start and end index
	size_t startIndex = m_pParticleData->getBoundaryStartIndex();
	size_t numParticle = m_pParticleData->getBoundaryNum();
	if(numParticle==0) return;
	
	//------- specifications of pointers and parameters-------
	
	// determine dir: x(0), y(1), or z(2)
	const int dir = (phase==-1)? -1:m_vDirSplitTable[m_iDirSplitOrder][phase];
	
	vector<const double*> position;
	position.push_back(m_pParticleData->m_vPositionX);
	position.push_back(m_pParticleData->m_vPositionY);
	if(m_iDimension==3) position.push_back(m_pParticleData->m_vPositionZ);		

	// input and also output data
	vector<double*> data;

	//-------------------------------------------------------
	
	if(phase==-1) {		
		// set input/output data
		data.push_back(m_pParticleData->m_vPressure);
		data.push_back(m_pParticleData->m_vVelocityU);
		data.push_back(m_pParticleData->m_vVelocityV);
		if(m_iDimension==3) data.push_back(m_pParticleData->m_vVelocityW);	
	}
	else {
		if(m_iDimension==2) {
			if(phase==0) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);
				if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityV); // next step is V
				else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityU); // next step is U
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
			}
			else assert(false);
		}
		else if(m_iDimension==3) {
			if(phase==0) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);	
			}
			else if(phase==2) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);
				if(dir==0) { 
					if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp1VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp1VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp1VelocityW);
			}
			else if(phase==3) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);	
				if(dir==0) { 
					if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp2VelocityW);
			}
			else assert(false);
		}	
	}
	
	computeOthOrderWeightedLPF(position, startIndex, numParticle, data);	
	
	for(size_t k=1; k<data.size(); k++) // k starts from 1: change velocity to the opposite of fluid particles
		for(size_t index=startIndex; index<startIndex+numParticle; index++)
			data[k][index] *= -1;
}


void HyperbolicLPSolver::setGhostPressureAndVelocity(int phase) {
	
	// no need to update for the last phase
	if((m_iDimension==2 && phase==2) || (m_iDimension==3 && phase==4)) return;

	// iteration start and end index
	size_t startIndex = m_pParticleData->getGhostStartIndex();
	size_t numParticle = m_pParticleData->getGhostNum();
	if(numParticle==0) return;
//	cout<<"-------HyperbolicLPSolver::setGhostPressureAndVelocity()-------"<<endl;
//	cout<<"startIndex = "<<startIndex<<endl;
//	cout<<"numParticle = "<<numParticle<<endl;
//	cout<<"---------------------------------------------------------------"<<endl;
	//------- specifications of pointers and parameters-------
	
	// determine dir: x(0), y(1), or z(2)
	const int dir = (phase==-1)? -1:m_vDirSplitTable[m_iDirSplitOrder][phase];
	
	vector<const double*> position;
	position.push_back(m_pParticleData->m_vPositionX);
	position.push_back(m_pParticleData->m_vPositionY);
	if(m_iDimension==3) position.push_back(m_pParticleData->m_vPositionZ);		

	// input and also output data
	vector<double*> data;

	//-------------------------------------------------------
	
	if(phase==-1) {		
		// set input/output data
		data.push_back(m_pParticleData->m_vPressure);
		data.push_back(m_pParticleData->m_vVelocityU);
		data.push_back(m_pParticleData->m_vVelocityV);
		if(m_iDimension==3) data.push_back(m_pParticleData->m_vVelocityW);	
	}
	else {
		if(m_iDimension==2) {
			if(phase==0) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);
				if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityV); // next step is V
				else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityU); // next step is U	
			}	
			else assert(false);
		}
		else if(m_iDimension==3) {
			if(phase==0) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);	
			}
			else if(phase==2) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);
				if(dir==0) { 
					if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp1VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp1VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp1VelocityW);
			}
			else if(phase==3) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);	
				if(dir==0) { 
					if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp2VelocityW);
			}
			else assert(false);
		}	
	}
	
	computeOthOrderWeightedLPF(position, startIndex, numParticle, data);	
		
}


void HyperbolicLPSolver::computeOthOrderWeightedLPF(vector<const double*>& position, 
													size_t startIndex, size_t numParticle,
													vector<double*>& data) {
	//cout<<"-------HyperbolicLPSolver::computeOthOrderWeightedLPF()-------"<<endl;
	// whole neighbour list
	const int *neighbourList = m_pParticleData->m_vNeighbourList;
	const int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;

	// iteration end index (1 after the last element)
	size_t endIndex = startIndex + numParticle;	

	//cout<<"numParticle = "<<numParticle<<endl;
	//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;
	//double startTime = omp_get_wtime();
	
	//cout<<"startIndex = "<<startIndex<<endl;
	//cout<<"numParticle = "<<numParticle<<endl;
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(size_t index=startIndex; index<endIndex; index++) {	

		size_t totalNumNei = neighbourListSize[index];
		if(totalNumNei == 0) { // does not have neighbour
			for(size_t k=0; k<data.size(); k++) {data[k][index] = 0;}
			continue; 
		}

		// do 0th order weighted LPF
		double data_w2_total[data.size()];
		fill_n(data_w2_total,data.size(),0);
		double w2_total = 0;
		for(size_t i=0; i<totalNumNei; i++) { 
			int neiIndex = neighbourList[index*maxNeiNum+i];
			//assert(nei_index >= nb && nei_index < nf); // fluid neighbour only
			
			double dis2 = 0;
			for(size_t j=0; j<position.size(); j++) {
				double h = ( position[j][neiIndex]-position[j][index] );
				dis2 += (h*h);
			}	

			double weight = exp(-dis2);
			w2_total += (weight*weight);
			
			for(size_t k=0; k<data.size(); k++) 
				data_w2_total[k] += (data[k][neiIndex]*weight*weight);	
			
		}
		
		if(w2_total != 0) {
			for(size_t k=0; k<data.size(); k++)
				data[k][index] = data_w2_total[k]/w2_total;
		}
		else assert(false);
		
		//if(data.size()==3 && position[0][index]<-24 && data[1][index]<0) { 
		//	cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!index="<<index<<endl;
		//	cout<<"####### u="<<data[1][index]
		//		<<", x="<<position[0][index]<<", y="<<position[1][index]<<endl;
		//	cout<<"neighbours:"<<endl;
		//	for(size_t i=0; i<totalNumNei; i++) { 
		//		int neiIndex = neighbourList[index*maxNeiNum+i];
		//		double dis2 = 0;
		//		for(size_t j=0; j<position.size(); j++) {
		//			double h = ( position[j][neiIndex]-position[j][index] );
		//			dis2 += (h*h);
		//		}	
		//		cout<<"neiIndex="<<neiIndex<<", u="<<data[1][neiIndex]
		//			<<", x="<<position[0][neiIndex]<<", y="<<position[1][neiIndex]
		//			<<", dist="<<sqrt(dis2)<<endl;
		//	}
		//}
	}

	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Compute zeroth order LPF takes %.16g seconds\n", elipsedTime);	
	
	//cout<<"--------------------------------------------------------------"<<endl;

}



void HyperbolicLPSolver::updateFluidState() {
		
	swap(m_pParticleData->m_vTemp1Volume, m_pParticleData->m_vVolume);
	swap(m_pParticleData->m_vTemp1Pressure, m_pParticleData->m_vPressure);
	swap(m_pParticleData->m_vTemp1SoundSpeed, m_pParticleData->m_vSoundSpeed);	
	
}


void HyperbolicLPSolver::moveFluidParticle() {
		
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	if(m_iDimension==2) {
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			m_pParticleData->m_vPositionX[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp1VelocityU[index]); // 0.5 (old + new)	
			
			m_pParticleData->m_vPositionY[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityV[index] + m_pParticleData->m_vTemp1VelocityV[index]);
		}	
	}
	else if(m_iDimension==3) {
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			m_pParticleData->m_vPositionX[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp1VelocityU[index]); // 0.5 (old + new)	
			
			m_pParticleData->m_vPositionY[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityV[index] + m_pParticleData->m_vTemp1VelocityV[index]);
			
			m_pParticleData->m_vPositionZ[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityW[index] + m_pParticleData->m_vTemp1VelocityW[index]);
		}
	}

}


void HyperbolicLPSolver::updateFluidVelocity() {
	
	swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vVelocityU);
	swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vVelocityV);	
	if(m_iDimension==3)	swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vVelocityW);

}

/*****************************************************************************
*
* The following functions are for testing purposes only
*
*****************************************************************************/

bool compByDist(const pair<double,int>& l, const pair<double,int>& r) {
	return l.first < r.first;
}

void HyperbolicLPSolver::testNeighbourSearch() {
	
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t totalNumParticle = m_pParticleData->getTotalNum();
	
	// build the search structure
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, totalNumParticle);	
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;

	// the entire neighbour list (to be updated in the following loops)
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	
	//double neiListDistTemp[maxNeiNum]; // a temp array for dist between a particle and its neighbours
//cout<<"1"<<endl;
	double* neighbourListDist = new double[m_pParticleData->m_iCapacity*maxNeiNum];	
	//double neighbourListDist[m_pParticleData->m_iCapacity*maxNeiNum]; 
//cout<<"2"<<endl;
	//
	//int* bruteForceNeighbourList = new int[m_pParticleData->m_iCapacity*maxNeiNum];
	//int* bruteForceNeighbourListSize = new int[m_pParticleData->m_iCapacity];
	int bruteForceNeighbourList[m_pParticleData->m_iCapacity*maxNeiNum];
//cout<<"3"<<endl;
	int bruteForceNeighbourListSize[m_pParticleData->m_iCapacity];
//cout<<"4"<<endl;
    //double bruteForceNeighbourListDist[m_pParticleData->m_iCapacity*maxNeiNum];
	double* bruteForceNeighbourListDist = new double[m_pParticleData->m_iCapacity*maxNeiNum];
//cout<<"5"<<endl;
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();		
	
	// fluid
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;
		
		size_t count = 0;
		vector<pair<double,int>> neiDistIndex;
		
		//-------brute-force neighbour search-------
		for(size_t j=fluidStartIndex; j<fluidStartIndex+totalNumParticle; j++) {
			double xd = positionX[index]-positionX[j];
			double yd = positionY[index]-positionY[j];
			double zd = positionZ[index]-positionZ[j];
			double dist = sqrt(xd*xd+yd*yd+zd*zd);
			if(dist <= m_fNeiSearchRadius && j!=index) {	
				neiDistIndex.push_back({dist,j});
				count++; // found a neighbour
			}
			if(count == maxNeiNum) break;
		}
		sort(neiDistIndex.begin(),neiDistIndex.end(),compByDist);
		for(size_t k=0; k<count; k++) {
			bruteForceNeighbourList[neiListStartIndex+k] = neiDistIndex[k].second;
			bruteForceNeighbourListDist[neiListStartIndex+k] = neiDistIndex[k].first;
		}
		bruteForceNeighbourListSize[index] = count;
		//------------------------------------------	
	
		
		//-------begin octree search-------
		size_t numNeiFound;
		m_pNeighbourSearcher->searchNeighbour(positionX[index], positionY[index], positionZ[index], m_fNeiSearchRadius, 
											  neighbourList+neiListStartIndex,neighbourListDist+neiListStartIndex,
											  numNeiFound,index); // output
		
		neighbourListSize[index] = numNeiFound;
		//---------------------------------
	}	
	
	// print debug info
	ofstream ofs1("bruteForceSearch");
	ofstream ofs2("octreeSearch");
	//cout<<"-------HyperbolicLPSolver::testNeighbourSearch()-------"<<endl;
	//cout<<"totalNumParticle="<<totalNumParticle<<endl;
	//cout<<"maxNeiNum="<<maxNeiNum<<endl;
	//cout<<"The neighbour list for Fluid particles:"<<endl;
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 	
		size_t neiListStartIndex = index*maxNeiNum;
		
		//cout<<"-----Brute Force SEARCH------"<<endl;
		ofs1<<"neighbourListSize["<<index<<"]="<<bruteForceNeighbourListSize[index]<<endl;
		ofs1<<"-----------------------------"<<endl;
		for(size_t k=neiListStartIndex; k<neiListStartIndex+bruteForceNeighbourListSize[index]; k++) 
			ofs1<<"neiIndex="<<bruteForceNeighbourList[k]<<"   "<<bruteForceNeighbourListDist[k]<<endl;
		ofs1<<"-----------------------------"<<endl;

		//cout<<"---------OCTREE SEARCH-------"<<endl;
		ofs2<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
		ofs2<<"-----------------------------"<<endl;
		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++) 
			ofs2<<"neiIndex="<<neighbourList[k]<<"   "<<neighbourListDist[k]<<endl;
		ofs2<<"-----------------------------"<<endl;

	}	
	//cout<<"-----------------------------------------------------------------"<<endl;
	
	delete[] neighbourListDist;
	delete[] bruteForceNeighbourListDist;
}


void HyperbolicLPSolver::searchNeighbourBruteForce(int index, int* neighbourList, int* neighbourListSize) {

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t totalNumParticle = m_pParticleData->getTotalNum();			
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
		
	size_t count = 0;
	vector<pair<double,int>> neiDistIndex(maxNeiNum);
	
	// WILL SEARCH FOR ALL PARTICLES AS NEIGHBOURS!!!
	//-------brute-force neighbour search-------
	for(size_t j=0; j<totalNumParticle; j++) {
		if(j == (size_t)index) continue; // DO NOT INCLUDE THE PARTICLE ITSELF AS A NEIGHBOUR!!!

		double xd = positionX[index]-positionX[j];
		double yd = positionY[index]-positionY[j];
		double zd = positionZ[index]-positionZ[j];
		double dist = sqrt(xd*xd+yd*yd+zd*zd);
		if(dist <= m_fNeiSearchRadius) {	
			neiDistIndex[count] = {dist,j};	
			count++;
		}
		if(count == maxNeiNum) { 
			cout<<"-------HyperbolicSolver::searchNeighbourBruteForce()-------"<<endl;
			cout<<"particle "<<index<<" has more than "<<maxNeiNum<<" neighbours!!!"<<endl;
			cout<<"-----------------------------------------------------------"<<endl;
			break;
		}
	}
	sort(neiDistIndex.begin(),neiDistIndex.begin()+count,compByDist);
	for(size_t k=0; k<count; k++) {
		neighbourList[k] = neiDistIndex[k].second;
		//cout<<"neiDistIndex["<<k<<"].first="<<neiDistIndex[k].first<<endl;
		//cout<<"neiDistIndex["<<k<<"].second="<<neiDistIndex[k].second<<endl;
		//cout<<"(size_t)neighbourList+k="<<(size_t)neighbourList+k<<endl;
		//cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	}
	neighbourListSize[0] = count;
	//------------------------------------------			
	
	//cout<<"-------HyperbolicLPSolver::searchNeighbourBruteForce()-------"<<endl;
	//cout<<"neighbourList="<<neighbourList<<endl;
	//cout<<"neighbourListSize[0]="<<neighbourListSize[0]<<endl;
	//cout<<"-------------------------------------------------------------"<<endl<<endl;

}


void HyperbolicLPSolver::searchNeighbourBruteForce(double x0, double y0, double z0, 
	int* neighbourList, size_t& numNeiFound) {

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case			
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
		
	size_t count = 0;
	vector<pair<double,int>> neiDistIndex(maxNeiNum);
	
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	// WILL ONLY SEARCH FLUID PARTICLES AS NEIGHOBURS!!!
	//-------brute-force neighbour search-------
	for(size_t j=fluidStartIndex; j<fluidEndIndex; j++) {
		double xd = x0-positionX[j];
		double yd = y0-positionY[j];
		double zd = z0-positionZ[j];
		double dist = sqrt(xd*xd+yd*yd+zd*zd);
		if(dist <= m_fNeiSearchRadius) {	
			neiDistIndex[count] = {dist,j};
			count++;
		}
		if(count == maxNeiNum) { 
			cout<<"-------HyperbolicSolver::searchNeighbourBruteForce()-------"<<endl;
			cout<<"A ghost particle has more than "<<maxNeiNum<<" neighbours!!!"<<endl;
			cout<<"-----------------------------------------------------------"<<endl;
			break;
		}
	}
	sort(neiDistIndex.begin(),neiDistIndex.begin()+count,compByDist);
	for(size_t k=0; k<count; k++) 
		neighbourList[k] = neiDistIndex[k].second;
	numNeiFound = count;
	//------------------------------------------			

}


void HyperbolicLPSolver::generateGhostParticleByBruteForceNeighbourSearch() {		
	
	// use each fluid bounding box to generate ghost particles
	const vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;
	
	const double h_r = 0.5*m_fInitParticleSpacing;
	const size_t capacity = m_pParticleData->m_iCapacity;	
	
	// temp space for valid ghost particles
	const size_t n = capacity - (m_pParticleData->m_iFluidNum + m_pParticleData->m_iBoundaryNum);
	vector<double> xGhost(n,0), yGhost(n,0), zGhost(n,0);
		
	int neiListTemp[m_pParticleData->m_iMaxNeighbourNum]; // the index of the neighbours of the ghost particle	
	
	size_t ghostCount = 0;
	if(m_iDimension==2) {

		for(size_t p=0; p<fluidBoxes.size(); p++) {

			double xmin = fluidBoxes[p]->getXmin();	
			double xmax = fluidBoxes[p]->getXmax();
			double ymin = fluidBoxes[p]->getYmin();	
			double ymax = fluidBoxes[p]->getYmax();
				
			HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	

			// compute the location of particles	
			for(size_t j=m0; j<=m1; j++) { 
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2D.computeX(0,k);
						double y = hex2D.computeY(j);
						
						size_t numNeiFound;
						searchNeighbourBruteForce(x,y,0,neiListTemp,numNeiFound);
						
						//cout<<numNeiFound<<endl;
						if(!isValidGhostParticle(x,y,0,neiListTemp,numNeiFound,p+1)) continue;		
						if(ghostCount >= n) assert(false); // exceed array size capacity
						xGhost[ghostCount] = x;
						yGhost[ghostCount] = y;
						ghostCount++;
					
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2D.computeX(1,k);
						double y = hex2D.computeY(j);
						
						size_t numNeiFound;
						searchNeighbourBruteForce(x,y,0,neiListTemp,numNeiFound);

						//cout<<numNeiFound<<endl;
						if(!isValidGhostParticle(x,y,0,neiListTemp,numNeiFound,p+1)) continue;		
						if(ghostCount >= n) assert(false); // exceed array size capacity	
						xGhost[ghostCount] = x;
						yGhost[ghostCount] = y;
						ghostCount++;
					}
				}
			}	
		}	
	}
	else if(m_iDimension==3) {
	
		for(size_t p=0; p<fluidBoxes.size(); p++) {
			
			double xmin = fluidBoxes[p]->getXmin();	
			double xmax = fluidBoxes[p]->getXmax();
			double ymin = fluidBoxes[p]->getYmin();	
			double ymax = fluidBoxes[p]->getYmax();
			double zmin = fluidBoxes[p]->getZmin();	
			double zmax = fluidBoxes[p]->getZmax();

			HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
			//get parameters of hexagonal packing
			size_t l0,l1;
			size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
			size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
			hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
								n0_odd, n1_odd, n0_even, n1_even, 
								nn0_odd, nn1_odd, nn0_even, nn1_even);	
			
			// compute the location of particles
			for(size_t i=l0; i<=l1; i++) { 
				if((i+1)%2 != 0) { //odd-numbered layers
					for(size_t j=m0_odd; j<=m1_odd; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows 
							for(size_t k=n0_odd; k<=n1_odd; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);	
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;		

							}
						} 
						else{ //even-numbered rows
							for(size_t k=n0_even; k<=n1_even; k++) {
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;	
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
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;	
							}
						} 
						else { //even-numbered rows
							for(size_t k=nn0_even; k<=nn1_even; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,p+1)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;	
							}
						}
					}	
				}    
			}	
		}	
	
	}
	
	size_t ghostIndex = m_pParticleData->getGhostStartIndex();
	for(size_t count=0; count<ghostCount; count++) {
		m_pParticleData->m_vPositionX[ghostIndex] = xGhost[count];
		m_pParticleData->m_vPositionY[ghostIndex] = yGhost[count];
		m_pParticleData->m_vPositionZ[ghostIndex] = zGhost[count];
		ghostIndex++;
	}

	m_pParticleData->m_iGhostNum = ghostIndex - m_pParticleData->m_iGhostStartIndex;
	m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
								   m_pParticleData->m_iBoundaryNum + 
								   m_pParticleData->m_iGhostNum;
	
	cout<<"-------HyperbolicLPSolver::generateGhostParticleByBruteForceNeighbourSearch()-------"<<endl;
	cout<<"m_pParticleData->m_iGhostStartIndex="<<m_pParticleData->m_iGhostStartIndex<<endl;
	cout<<"m_pParticleData->m_iGhostNum="<<m_pParticleData->m_iGhostNum<<endl;
	cout<<"m_pParticleData->m_iTotalNum="<<m_pParticleData->m_iTotalNum<<endl;
	cout<<"------------------------------------------------------------------------------------"<<endl;
}


void HyperbolicLPSolver::searchNeighbourForAllParticleByBruteForceNeighbourSearch() {	
	
	// the entire neighbour list (to be updated in the following loops)
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	size_t boundaryStartIndex = m_pParticleData->getBoundaryStartIndex();
	size_t boundaryEndIndex = boundaryStartIndex + m_pParticleData->getBoundaryNum();
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	size_t ghostEndIndex = ghostStartIndex + m_pParticleData->getGhostNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
	
	//cout<<"---HyperbolicLPSolver::searchNeighbourForAllParticleByBruteForceNeighbourSearch()---"<<endl;

	// fluid
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;
		searchNeighbourBruteForce(index,neighbourList+neiListStartIndex,neighbourListSize+index);
			

		// PRINT DEBUG INFO
		//cout<<"index="<<index<<endl;
		//cout<<"neiListStartIndex="<<neiListStartIndex<<endl;
		//cout<<"neighbourList+neiListStartIndex="<<neighbourList+neiListStartIndex<<endl;
		//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
		//for(int k=0; k<neighbourListSize[index]; k++) 
		//	cout<<"neighbourList["<<neiListStartIndex+k<<"]="<<neighbourList[neiListStartIndex+k]<<endl;
		
	}
	
	// boundary
	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;			
		searchNeighbourBruteForce(index,neighbourList+neiListStartIndex,neighbourListSize+index);

		// boundary particle take only fluid particles as neighbours
		int incr = 0;
		for(int k=0; k<neighbourListSize[index]; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;

		// PRINT DEBUG INFO
		//cout<<"index="<<index<<endl;
		//cout<<"neiListStartIndex="<<neiListStartIndex<<endl;
		//cout<<"neighbourList+neiListStartIndex="<<neighbourList+neiListStartIndex<<endl;
		//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
		//for(int k=0; k<neighbourListSize[index]; k++) 
		//	cout<<"neighbourList["<<neiListStartIndex+k<<"]="<<neighbourList[neiListStartIndex+k]<<endl;
	}


	// ghost
	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;
		searchNeighbourBruteForce(index,neighbourList+neiListStartIndex,neighbourListSize+index);

		// boundary particle take only fluid particles as neighbours
		int incr = 0;
		for(int k=0; k<neighbourListSize[index]; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;

		// PRINT DEBUG INFO
		//cout<<"index="<<index<<endl;
		//cout<<"neiListStartIndex="<<neiListStartIndex<<endl;
		//cout<<"neighbourList+neiListStartIndex="<<neighbourList+neiListStartIndex<<endl;
		//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
		//for(int k=0; k<neighbourListSize[index]; k++) 
		//	cout<<"neighbourList["<<neiListStartIndex+k<<"]="<<neighbourList[neiListStartIndex+k]<<endl;		    
	}


	// PRINT DEBUG INFO
	//cout<<"-------HyperbolicLPSolver::searchNeighbourForAllParticle()-------"<<endl;
	//cout<<"m_pParticleData->getTotalNum()="<<m_pParticleData->getTotalNum()<<endl;
	//cout<<"maxNeiNum="<<maxNeiNum<<endl;
	//cout<<"Fluid particles:"<<endl;
	//for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;			
	//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//	for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
	//		cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//}
	//cout<<"Boundary particles:"<<endl;
	//for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;			
	//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//	for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
	//		cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//}
	//cout<<"Ghost particles:"<<endl;
	//for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;			
	//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//	for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
	//		cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//}
	cout<<"-----------------------------------------------------------------"<<endl;
	
}



bool HyperbolicLPSolver::checkUpwindNeighbourList() {
	
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;		
	const double *positionZ = m_pParticleData->m_vPositionZ;
	if(m_iDimension==2) positionZ = nullptr;

	// get the upwind neighbour lists (to be updated in the following)
	int *neighbourListRight = m_pParticleData->m_vNeighbourListRight;
	int *neighbourListLeft = m_pParticleData->m_vNeighbourListLeft;
	int *neighbourListNorth = m_pParticleData->m_vNeighbourListNorth;
	int *neighbourListSouth = m_pParticleData->m_vNeighbourListSouth;
	int *neighbourListUp, *neighbourListDown;
	if(m_iDimension==3) {
		neighbourListUp = m_pParticleData->m_vNeighbourListUp;
		neighbourListDown = m_pParticleData->m_vNeighbourListDown;
	}

	// get the size of neighbour lists (to be updated in the following)
	int *neighbourListRightSize = m_pParticleData->m_vNeighbourListRightSize;
	int *neighbourListLeftSize = m_pParticleData->m_vNeighbourListLeftSize;
	int *neighbourListNorthSize = m_pParticleData->m_vNeighbourListNorthSize;
	int *neighbourListSouthSize = m_pParticleData->m_vNeighbourListSouthSize;
	int *neighbourListUpSize, *neighbourListDownSize;
	if(m_iDimension==3) {
		neighbourListUpSize = m_pParticleData->m_vNeighbourListUpSize;
		neighbourListDownSize = m_pParticleData->m_vNeighbourListDownSize;
	}
	
	size_t maxNeiNumInOneDir = m_pParticleData->m_iMaxNeighbourNumInOneDir;

	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		
		size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
			
		// print the upwind neighbour lists
		//if(index%100==0 && index<=1000) cout<<"neighbourListRightSize["<<index<<"]="<<neighbourListRightSize[index]<<endl;
		for(int k=0; k<neighbourListRightSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListRight[j];
			assert(positionX[neiIndex] > positionX[index]);
			//cout<<"neighbourListRight["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionX[neiIndex]=%.16g\n",positionX[neiIndex]);
			//printf("positionX[index]=%.16g\n",positionX[index]);
		}
		//if(index%100==0 && index<=1000)cout<<"neighbourListLeftSize["<<index<<"]="<<neighbourListLeftSize[index]<<endl;
		for(int k=0; k<neighbourListLeftSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListLeft[j];
			assert(positionX[neiIndex] < positionX[index]);
			//cout<<"neighbourListLeft["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionX[neiIndex]=%.16g\n",positionX[neiIndex]);
			//printf("positionX[index]=%.16g\n",positionX[index]);
		}
		//if(index%100==0 && index<=1000) cout<<"neighbourListNorthSize["<<index<<"]="<<neighbourListNorthSize[index]<<endl;
		for(int k=0; k<neighbourListNorthSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListNorth[j];
			assert(positionY[neiIndex] > positionY[index]);
			//cout<<"neighbourListNorth["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionY[neiIndex]=%.16g\n",positionY[neiIndex]);
			//printf("positionY[index]=%.16g\n",positionY[index]);
		}
		//if(index%100==0 && index<=1000) cout<<"neighbourListSouthSize["<<index<<"]="<<neighbourListSouthSize[index]<<endl;
		for(int k=0; k<neighbourListSouthSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListSouth[j];
			assert(positionY[neiIndex] < positionY[index]);
			//cout<<"neighbourListSOuth["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionY[neiIndex]=%.16g\n",positionY[neiIndex]);
			//printf("positionY[index]=%.16g\n",positionY[index]);
		}
	}	
	if(m_iDimension==3) {
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		
			size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
			
			// print the upwind neighbour lists
			//if(index%100==0 && index<=1000) cout<<"neighbourListUpSize["<<index<<"]="<<neighbourListUpSize[index]<<endl;
			for(int k=0; k<neighbourListUpSize[index]; k++) {
				int j = neiListInOneDirStartIndex+k;
				int neiIndex = neighbourListUp[j];
				assert(positionZ[neiIndex] > positionZ[index]);
				//cout<<"neighbourListUp["<<j<<"]="<<neiIndex<<endl;	
				//printf("positionZ[neiIndex]=%.16g\n",positionZ[neiIndex]);
				//printf("positionZ[index]=%.16g\n",positionZ[index]);
			}
			//if(index%100==0 && index<=1000) cout<<"neighbourListDownSize["<<index<<"]="<<neighbourListDownSize[index]<<endl;
			for(int k=0; k<neighbourListDownSize[index]; k++) {
				int j = neiListInOneDirStartIndex+k;
				int neiIndex = neighbourListDown[j];
				assert(positionZ[neiIndex] < positionZ[index]);
				//cout<<"neighbourListDown["<<j<<"]="<<neiIndex<<endl;	
				//printf("positionZ[neiIndex]=%.16g\n",positionZ[neiIndex]);
				//printf("positionZ[index]=%.16g\n",positionZ[index]);
			}	
		}	
	}	

	return true;
}



string HyperbolicLPSolver::rightFlush(size_t writeStep, size_t numDigits) {
	
	assert(pow(10,numDigits) >= writeStep);

	string result;

	if(writeStep == 0) numDigits--;
	for(size_t i=writeStep; i!=0; i /= 10) numDigits--;

	for( ; numDigits>0; numDigits--) result.push_back('0');
	
	result += to_string(writeStep); 
	
	return result;

}

int HyperbolicLPSolver::writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle) {
	
	// Create an output file the name "filename"
	string filename = "vtkDebug" + rightFlush(writeStep, 7) + ".vtk";
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Error opening file: %s\n",filename.c_str()); 
		return 1;
	}
	size_t endIndex = startIndex + numParticle;


	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",numParticle);
	if(m_pParticleData->getDimension()==2) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vPositionX[i], m_pParticleData->m_vPositionY[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vPositionX[i], m_pParticleData->m_vPositionY[i], m_pParticleData->m_vPositionZ[i]);
	}
	fprintf(outfile,"POINT_DATA %ld\n",numParticle);
	
	fprintf(outfile,"VECTORS Velocity double\n");
	if(m_pParticleData->getDimension()==2) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vVelocityU[i], m_pParticleData->m_vVelocityV[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vVelocityU[i], m_pParticleData->m_vVelocityV[i], m_pParticleData->m_vVelocityW[i]);
	}
	
	fprintf(outfile,"SCALARS LPFOrder_right int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderRight[i]);
	
	fprintf(outfile,"SCALARS LPFOrder_left int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderLeft[i]);

	fprintf(outfile,"SCALARS LPFOrder_north int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderNorth[i]);

	fprintf(outfile,"SCALARS LPFOrder_south int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderSouth[i]);
		

	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS LPFOrder_up int\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderUp[i]);

		fprintf(outfile,"SCALARS LPFOrder_down int\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderDown[i]);	
	}


	fprintf(outfile,"SCALARS velocity_u double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vVelocityU[i]);
	
	fprintf(outfile,"SCALARS velocity_v double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vVelocityV[i]);
	
	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS velocity_w double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g\n",m_pParticleData->m_vVelocityW[i]);
	}
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vSoundSpeed[i]);
	
	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vPressure[i]);
		
	fprintf(outfile,"SCALARS density double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",1./m_pParticleData->m_vVolume[i]);
		

	fclose(outfile);
	
	return 0;
	
	
}

////////////////////////////////////////////////////////////////////////////////////////
// End of HyperbolicLPSolver
////////////////////////////////////////////////////////////////////////////////////////





