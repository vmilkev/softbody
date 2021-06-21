#include "cs_sample.hpp"

int main(int argc, char* argv[]) {

	Sample S;

	if (argc < 2) {
		S.ReportExit(1);
		exit(EXIT_FAILURE);
	}
	else if (argc == 2){
		S.type = "sample";
	}
	else if (argc == 3){
		S.type = "simulation";
	}
	else {
		S.ReportExit(2);
		exit(EXIT_FAILURE);
	}

	std::string fileNameOutput;

	if ( S.CheckProjType() ) {		
		
		#ifdef _WIN64
			fileNameOutput = S.type + "\\" + "sample.smpl";
		#endif

		#ifdef linux
			fileNameOutput = S.type + "/" + "sample.smpl";
		#endif
	
		S.Init( "sample_runtime.log", argv[1] );
		S.SetSample(1); // start with one cell
		S.RunSimulation( S.samplTime );
	}
	else {

		#ifdef _WIN64
			fileNameOutput = S.type + "\\" + "simulation.smpl";
		#endif

		#ifdef linux
			fileNameOutput = S.type + "/" + "simulation.smpl";
		#endif
	
		S.Init("simulation_runtime.log", argv[1]);
		S.RestoreSample(argv[2]);
		S.RunGrowthSimulation( S.runTime );
	}

	S.Finalize(fileNameOutput);

	return 0;
}
