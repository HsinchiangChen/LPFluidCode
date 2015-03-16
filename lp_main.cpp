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
	
	string inputfileName;
	string outputfileNameAll;		
	string outputfileNameFluid;
	string debugfileName = "debug"; // if no -d option is provided the file "debug" will contain only some init info
	bool ifDebug = false;
	for(int i=1; i<argc; i++) { // argv[0] is name of the executable 
		if(!strcmp(argv[i], "-i")) { // input; strcmp returns 0 when 2 c_str is the same
			if(i+1 >= argc || argv[i+1][0]=='-') { // no inputfile name following -i option
				cout<<"ERROR: Needs to provide inputfile name!"<<endl;
				exit(1);
			}
			inputfileName = argv[i+1];
		}
		else if(!strcmp(argv[i], "-o")) { // output
			if(i+1 >= argc || argv[i+1][0]=='-') { // no outputfile name following -o option
				cout<<"ERROR: Needs to provide outputfile name!"<<endl;
				exit(1);
			}
			if(mkdir(argv[i+1], 0777) == -1) { //mkdir takes const char*
				cerr<<"ERROR:"<<strerror(errno)<<endl;
				exit(1);
			}
			outputfileNameAll = string(argv[i+1]) + "/out";
			outputfileNameFluid = string(argv[i+1]) + "/out_fluid";
		}
		else if(!strcmp(argv[i], "-d")) { // debug
			ifDebug = true;
			if(i+1 >= argc || argv[i+1][0]=='-') // no debugfile name following -d option is fine, use default name
				continue;
			debugfileName = argv[i+1];
		}	
	}
	if(inputfileName.empty()) {
		cout<<"ERROR: Needs to provide -i option!"<<endl;
		exit(1);
	}
	if(outputfileNameAll.empty() || outputfileNameFluid.empty()) {
		cout<<"ERROR: Needs to provide -o option!"<<endl;
		exit(1);
	}

	//TODO use exception
	// create a directory for output based on argv[4]; if the directory already exists assert(false)	
	//if(mkdir(argv[4], 0777) == -1) { //mkdir takes const char*
	//	cerr <<"ERROR:"<< strerror(errno) <<endl;
	//	assert(false);
	//}		
		
	// read parameters from input file and initialize boundary and fluid geoemtry/state
	//const string inputfileName = argv[2];
	//const string debugfileName = "debug";
	
	Initializer* init = new Initializer(inputfileName, debugfileName, ifDebug);
	
	// load boundary and fluid geoemtry/state into particleData
	ParticleData* pData = new ParticleData(*init);
	
	//TODO: do this in command line using eg. -all -fluid or in input file as options entered by user
	//const string outputfileNameAll = string(argv[4]) + "/out";		
	//const string outputfileNameFluid = string(argv[4]) + "/out_fluid";
	ParticleViewer* pViewerAll, *pViewerFluid;
	NeighbourSearcher* neiSearcher;
	LPSolver* lpSolver;
	if(init->getDimension() == 1) {
		pViewerAll = new TXTParticleViewer1D(pData, "all", outputfileNameAll);	
		pViewerFluid = new TXTParticleViewer1D(pData, "fluid", outputfileNameFluid);
		lpSolver = new HyperbolicLPSolver1D(*init, pData);	
	}
	else {	
		pViewerAll = new VTKParticleViewer(pData, "all", outputfileNameAll);	
		pViewerFluid = new VTKParticleViewer(pData, "fluid", outputfileNameFluid);
		neiSearcher = new OctreeSearcher(*init); // neighbour is used only in 2D and 3D
		lpSolver = new HyperbolicLPSolver(*init, pData, neiSearcher);
	}
	vector<ParticleViewer*> viewers({pViewerAll,pViewerFluid});	

	TimeController* timeControl = new DefaultTimeController(*init, lpSolver, viewers);

	delete init; // initialization finished
	
	return timeControl->solve();

}
