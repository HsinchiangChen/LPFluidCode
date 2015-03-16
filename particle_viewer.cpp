#include "particle_viewer.h"
#include "particle_data.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip> //setw

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Start of ParticleViewer
////////////////////////////////////////////////////////////////////////////////////////


string ParticleViewer::rightFlush(size_t writeStep, size_t numDigits) {
	
	assert(pow(10,numDigits) > writeStep);

	string result;

	if(writeStep == 0) numDigits--;
	for(size_t i=writeStep; i!=0; i /= 10) numDigits--;

	for( ; numDigits>0; numDigits--) result.push_back('0');
	
	result += to_string(writeStep); 
	
	return result;

}


////////////////////////////////////////////////////////////////////////////////////////
// End of ParticleViewer
////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////
// Start of VTKParticleViewer
////////////////////////////////////////////////////////////////////////////////////////

VTKParticleViewer::VTKParticleViewer(ParticleData* data, const std::string& particleType, 
const string& outputfileName, int numDigits) {
	m_pParticleData = data;
	m_sOutputfileName = outputfileName;
	m_iNumDigits = numDigits; 
	m_sParticleType = particleType;
}


int VTKParticleViewer::writeResult(double time, size_t writeStep) {
	 
	// alias pointers
	const double* positionX = m_pParticleData->getPositionX();
	const double* positionY = m_pParticleData->getPositionY();
	const double* positionZ = m_pParticleData->getPositionZ();

	const double* velocityU = m_pParticleData->getVelocityU();
	const double* velocityV = m_pParticleData->getVelocityV();
	const double* velocityW = m_pParticleData->getVelocityW();

	const double* volume = m_pParticleData->getVolume();
	const double* pressure = m_pParticleData->getPressure();
	const double* soundSpeed = m_pParticleData->getSoundSpeed();

	const int* objectTag = m_pParticleData->getObjectTag();

	// Create an output file the name "filename"
	string filename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".vtk";
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Unable to open file: %s\n",filename.c_str()); 
		return 1;
	}
	size_t startIndex, numParticle;
	if(m_sParticleType=="all") {
		startIndex = 0;
		numParticle = m_pParticleData->getTotalNum();
	}
	else if(m_sParticleType=="fluid") {
		startIndex = m_pParticleData->getFluidStartIndex();
		numParticle = m_pParticleData->getFluidNum();
	}

	size_t endIndex = startIndex + numParticle;


	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",numParticle);
	if(m_pParticleData->getDimension()==2) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], positionZ[i]);
	}
	fprintf(outfile,"POINT_DATA %ld\n",numParticle);
	
	fprintf(outfile,"VECTORS Velocity double\n");
	if(m_pParticleData->getDimension()==2) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], velocityW[i]);
	}

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",pressure[i]);
		
	fprintf(outfile,"SCALARS volume double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",volume[i]);
	
	fprintf(outfile,"SCALARS velocity_u double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityU[i]);
	
	fprintf(outfile,"SCALARS velocity_v double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityV[i]);
	
	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS velocity_w double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g\n",velocityW[i]);
	}
	
	fprintf(outfile,"SCALARS object_tag int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",objectTag[i]);
	
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",soundSpeed[i]);

		

	fclose(outfile);
	
	return 0;
}

/*

int VTKParticleViewer::writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle) {
	
	// alias pointers
	const double* positionX = m_pParticleData->getPositionX();
	const double* positionY = m_pParticleData->getPositionY();
	const double* positionZ = m_pParticleData->getPositionZ();

	const double* velocityU = m_pParticleData->getVelocityU();
	const double* velocityV = m_pParticleData->getVelocityV();
	const double* velocityW = m_pParticleData->getVelocityW();

	const double* volume = m_pParticleData->getVolume();
	const double* pressure = m_pParticleData->getPressure();
	const double* soundSpeed = m_pParticleData->getSoundSpeed();

	const int* objectTag = m_pParticleData->getObjectTag();

	// Create an output file the name "filename"
	string filename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".vtk";
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
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], positionZ[i]);
	}
	fprintf(outfile,"POINT_DATA %ld\n",numParticle);
	
	fprintf(outfile,"VECTORS Velocity double\n");
	if(m_pParticleData->getDimension()==2) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], velocityW[i]);
	}

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",pressure[i]);
		
	fprintf(outfile,"SCALARS volume double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",volume[i]);
	
	fprintf(outfile,"SCALARS velocity_u double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityU[i]);
	
	fprintf(outfile,"SCALARS velocity_v double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityV[i]);
	
	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS velocity_w double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g\n",velocityW[i]);
	}
	
	fprintf(outfile,"SCALARS object_tag int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",objectTag[i]);
	
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",soundSpeed[i]);

		

	fclose(outfile);
	
	return 0;
	
	
}
*/

////////////////////////////////////////////////////////////////////////////////////////
// End of VTKParticleViewer
////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////
// Start of TXTParticleViewer1D
////////////////////////////////////////////////////////////////////////////////////////

TXTParticleViewer1D::TXTParticleViewer1D(ParticleData* data, const std::string& particleType, 
const string& outputfileName, int numDigits) {
	m_pParticleData = data;
	m_sOutputfileName = outputfileName;
	m_iNumDigits = numDigits; 
	m_sParticleType = particleType;
}


int TXTParticleViewer1D::writeResult(double time, size_t writeStep) {
	 
	// alias pointers
	const double* positionX  = m_pParticleData->getPositionX();	
	const double* velocityU  = m_pParticleData->getVelocityU();	
	const double* volume     = m_pParticleData->getVolume();
	const double* pressure   = m_pParticleData->getPressure();
	const double* soundSpeed = m_pParticleData->getSoundSpeed();
	const double* dd1        = m_pParticleData->getDD1();
	const double* dd2_left   = m_pParticleData->getDD2Left();
	const double* dd2_right  = m_pParticleData->getDD2Right();
	const double* dd3_left   = m_pParticleData->getDD3Left();
	const double* dd3_right  = m_pParticleData->getDD3Right();
	
	
	// Create an output file with the name "filename"
	string filename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
	ofstream outfile;
	outfile.open(filename.c_str());
	
	size_t startIndex, numParticle;
	if(m_sParticleType=="all") {
		startIndex = 0;
		numParticle = m_pParticleData->getTotalNum();
	}
	else if(m_sParticleType=="fluid") {
		startIndex = 1;
		numParticle = m_pParticleData->getTotalNum()-1;
	}
	size_t endIndex = startIndex + numParticle;


	if(outfile.is_open()) {
		outfile<<"x"<<setw(24)
			   <<"V"<< setw(24) 
			   <<"vel"<< setw(24) 
			   <<"p"<<setw(24)
			   <<"cs"<<setw(24)
			   <<"dd1[i]"<<setw(24)
			   <<"dd2_left[i]"<<setw(24)
			   <<"dd2_right[i]"<<setw(24)
			   <<"dd3_left[i]"<<setw(24)
			   <<"dd3_right[i]"<<setw(24)
			   <<endl;
		for(size_t i=startIndex; i<endIndex; i++) {	
			outfile.precision(15);
			outfile<<left<<setw(24) 
			<<positionX[i]<<setw(24)
			<<volume[i] <<setw(24) 
			<<velocityU[i]<<setw(24) 
			<<pressure[i]<<setw(24)  
			<<soundSpeed[i]<<setw(24)
			<<dd1[i]<<setw(24)
			<<dd2_left[i]<<setw(24)
			<<dd2_right[i]<<setw(24)
			<<dd3_left[i]<<setw(24)
			<<dd3_right[i]<<setw(24)
			<<endl;
		}
		
	}
	else {
		cout<<"Unable to open file: "<<filename<<endl;
		return 1; // ERROR
	}

	outfile.close();	
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of TXTParticleViewer1D
////////////////////////////////////////////////////////////////////////////////////////
