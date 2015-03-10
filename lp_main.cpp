/*! 
 * \author Chen, Hsin-Chiang <morrischen2008@gmail.com> 
 * \date Mon Nov. 03 2014 
 * 
 * \brief 
 *      
 */

#include "initializer.h"
#include "neighbour_searcher.h"
#include "particle_data.h"
#include "particle_viewer.h"
#include "lp_solver.h"
#include "time_controller.h"

#include <iostream>
#include <cassert>
#include <sys/stat.h> // mkdir()
#include <cstring> // strerror()
#include <vector>
using namespace std;

// execute as ./lp -i <inputfile name> -o <output directory name> 

int main(int argc, const char* argv[]) {	
	
	// create a directory for output based on argv[4]	
	if(mkdir(argv[4], 0777) == -1) { //mkdir takes const char*
		cerr <<"ERROR:"<< strerror(errno) <<endl;
		assert(false);
	}		
		
	// read parameters from input file and initialize boundary and fluid geoemtry/state
	const string inputfileName = argv[2];
	Initializer* init = new Initializer(inputfileName);
	
	// load boundary and fluid geoemtry/state into particleData
	ParticleData* pData = new ParticleData(*init);
	
	// create a particle viewer
	const string outputfileNameAll = string(argv[4]) + "/out";		
	const string outputfileNameFluid = string(argv[4]) + "/out_fluid";
	ParticleViewer* pViewerAll = new VTKParticleViewer(pData, "all", outputfileNameAll);	
	ParticleViewer* pViewerFluid = new VTKParticleViewer(pData, "fluid", outputfileNameFluid);
	vector<ParticleViewer*> viewers({pViewerAll,pViewerFluid});

	NeighbourSearcher* neiSearcher = new OctreeSearcher(*init);
	
	LPSolver* lpSolver = new HyperbolicLPSolver(*init, pData, neiSearcher);	

	TimeController* timeControl = new DefaultTimeController(*init, lpSolver, viewers);

	delete init; // initialization finished
	
	return timeControl->solve();

}
