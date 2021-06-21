#include "Parameters.h"
#include "Functions.h"

//#include "../../graphics/src/visualizer.hpp"

int main(int argc, char* argv[]) {

	RuntimeLog(runTimeFile, false, "", 1.0, true, false);

	if (argc < 2) {
		RuntimeLog(runTimeFile, true, "You forget to enter the name of file with a list of parameters !", 0.0, false, false);
		RuntimeLog(runTimeFile, false, "Example:  c:>'path'/dem.exe 'path'/parameters.dat [optional input] ", 0.0, false, false, true);
		RuntimeLog(runTimeFile, true, "EXIT with error !", 0.0, false, true);
		exit(EXIT_FAILURE);
	}
	else if (argc == 2)
		type = "sample";
	else if (argc == 3)
		type = "simulation";
	else {
		RuntimeLog(runTimeFile, true, "Too much input parameters ! It should not exceed 3.", 0.0, false, false);
		RuntimeLog(runTimeFile, true, "EXIT with error !", 0.0, false, true);
		exit(EXIT_FAILURE);
	}

	srand((unsigned)time(NULL));

	RuntimeLog(runTimeFile, true, "READING INITIAL DATA (PARAMETERS) FROM FILE", 0.0, false, false);
	
	ReadParametersFromFile(argv[1]);

	RuntimeLog(runTimeFile, true, "RUN PROJECT MANAGER", 0.0, false, false);

	ProjectManager();
	
	double _sizeYstr;

	if (CheckProjType(type)) {

		std::string fileNameOutput;
		
		#ifdef _WIN64
			fileNameOutput = type + "\\" + "sample.smpl";
		#endif

		#ifdef linux
			fileNameOutput = type + "/" + "sample.smpl";
		#endif

		RuntimeLog(runTimeFile, true, "STARTING 'SAMPLE' PROJECT. INITIATE PARAMETERS.", 0.0, false, false);

		InitiateParameters(_rMin, _rMax);

		_sizeYstr = _sizeY*sampleStretch;

		// initiate the load piston structure, showing no piston load
		wLoadPiston.movementFlag = 0;
		wLoadPiston.loadTime = 0.;
		wLoadPiston.loadRate = 0.;

		_timeInterval = samplTime;
		
		RuntimeLog(runTimeFile, true, "INITIAL SAMPLE CONSTRUCTION:", 0.0, false, false);
		RuntimeLog(runTimeFile, false, "fabricating a (stretched) sample", 0.0, false, false, true);

		SetSample(_rMin, _rMax, _sizeX, _sizeYstr, _nodesRes, dim, grain, stdDev);

		RuntimeLog(runTimeFile, true, "Fabrication has been completed .", 0.0, false, false);

		SetWalls(_sizeX, _sizeYstr, dim, wLoadPiston, grain);
		
		if (recordOutput)
			relaxationFlag = 1;
		else
			relaxationFlag = 0;

		snapshotCounter = snapshotCounterP = snapshotCounterG = snapshotCounterF = snapshotCounterB = 0;
		grnCount = pstCount = wlCount = fCount = bCount = 0;
		grnTop = pstTop = wlTop = fTop = bTop = true;
		snapshotBuffer = 200;
		bSmf = grnSmf = pstSmf = wlSmf = fSmf = true;

 		omp_set_nested(1);
#pragma omp parallel sections
		{
#pragma omp section
			{
				relaxationFlag = DynamicRelaxation(_timeInterval, dim, grain, wLoadPiston, false, b, recordOutput, type, stopMechRatio);
			}
#pragma omp section
			{
				if (uint_recordOutput == 2){
					std::cout<<"visual postprocessing"<<std::endl;
					RunVisualPostProcessing(false, dim);
				}
				else{
					std::cout<<"normal postprocessing"<<std::endl;
					RunPostProcessing(false, dim);
				}
			}
		}		



//exit(1);		
		/* record final state if no record output option */
		if (!recordOutput)
			PostProcessing(false, dim, _pst, wl);
		
		RuntimeLog(runTimeFile, true, "COMPLETED.", 0.0, false, false);
		RuntimeLog(runTimeFile, true, "IMPLEMENT POROSITY:", 0.0, false, false);
		
		GetPorouseSample(_sizeX, _sizeY, dim, grain, porosity);

		ClearForces(grain);

		SetBondsMatrix(b, grain);

		FindBonds(dim, b, grain, shape, scale);

		if (recordOutput)
			relaxationFlag = 1;
		else
			relaxationFlag = 0;

		grnCount = pstCount = wlCount = fCount = bCount = 0;
		grnTop = pstTop = wlTop = fTop = bTop = true;
		snapshotBuffer = 200;
		bSmf = grnSmf = pstSmf = wlSmf = fSmf = true;
		
		omp_set_nested(1);
#pragma omp parallel sections
		{
#pragma omp section
			{
				relaxationFlag = DynamicRelaxation(_timeInterval, dim, grain, wLoadPiston, isCohesion, b, recordOutput, type, stopMechRatio);
			}
#pragma omp section
			{
				if (uint_recordOutput == 2)
					RunVisualPostProcessing(isCohesion, dim);
				else
					RunPostProcessing(isCohesion, dim);
			}
		}
		
		/* record final state if no record output option */
		if (!recordOutput)
			PostProcessing(isCohesion, dim, _pst, wl);
		
		RuntimeLog(runTimeFile, true, "COMPLETED.", 0.0, false, false);
		
		RuntimeLog(runTimeFile, true, "CORRECT THE SIZE OF THE SAMPLE (REMOVE THE EXCESS OF PARTICLES).", 0.0, false, false);
		
		/* After the correction all the grains belong to the new instance '_grain' */
		GetCorrectedSample(_rMax, _sizeX, _sizeY, dim, grain, _grain);

		grainNum = _grain.size(); /* Temporal solution; '_grain' needs to be passed to the post processing functions. */

		SetWalls(_sizeX, _sizeYstr, dim, wLoadPiston, _grain);

		ClearForces(_grain);

		SetBondsMatrix(b, _grain);

		FindBonds(dim, b, _grain, shape, scale);

		if (recordOutput)
			relaxationFlag = 1;
		else
			relaxationFlag = 0;

		grnCount = pstCount = wlCount = fCount = bCount = 0;
		grnTop = pstTop = wlTop = fTop = bTop = true;
		snapshotBuffer = 200;
		bSmf = grnSmf = pstSmf = wlSmf = fSmf = true;

		omp_set_nested(1);
#pragma omp parallel sections
		{
#pragma omp section
			{
				relaxationFlag = DynamicRelaxation(_timeInterval, dim, _grain, wLoadPiston, isCohesion, b, recordOutput, type, stopMechRatio);
			}
#pragma omp section
			{
				if (uint_recordOutput == 2)
					RunVisualPostProcessing(isCohesion, dim);
				else
					RunPostProcessing(isCohesion, dim);
			}
		}

		/* record final state if no record output option */
		if (!recordOutput)
			PostProcessing(isCohesion, dim, _pst, wl);

		RuntimeLog(runTimeFile, true, "COMPLETED.", 0.0, false, false);

		RuntimeLog(runTimeFile, true, "FINALIZING THE PROJECT. SAVING *.SMPL FILE", 0.0, false, false);

		FinalizeSampleProject(fileNameOutput, _grain, b, wLoadPiston, wl);

		RuntimeLog(runTimeFile, true, "SAMPLE PROJECT COMPLETED.", 0.0, false, true);

		return 0;
	}
	else {

		RuntimeLog(runTimeFile, true, "STARTING 'SIMULATION' PROJECT. INITIATE PARAMETERS.", 0.0, false, false);

		InitiateParameters(_rMin, _rMax);

		RuntimeLog(runTimeFile, true, "LOADING FABRICATED SAMPLE.", 0.0, false, false);

		InitializeSimulationProject(argv[2], grain, b, wLoadPiston, wl);

		grainNum = _grain.size(); /* Temporal solution; '_grain' needs to be passed to the post processing functions. */

		_sizeYstr = _sizeY*1.0; /* extend a bit walls in order to avoid grains move outside the sample domain */

		_timeInterval = runTime;

		wLoadPiston.movementFlag = 1;
		wLoadPiston.loadTime = runTime;
		wLoadPiston.loadRate = -lRate;
		wLoadPiston.sizeR = pSize/2;

		SetWalls(_sizeX, _sizeYstr, dim, wLoadPiston, grain);
		
		ClearForces(grain);
		
		SetBondsMatrix(b, grain);
		
		FindBonds(dim, b, grain, shape, scale);
		
		if (recordOutput)
			relaxationFlag = 1;
		else
			relaxationFlag = 0;

		snapshotCounter = snapshotCounterP = snapshotCounterG = snapshotCounterF = snapshotCounterB = 0;
		grnCount = pstCount = wlCount = fCount = bCount = 0;
		grnTop = pstTop = wlTop = fTop = bTop = true;
		snapshotBuffer = 200;
		bSmf = grnSmf = pstSmf = wlSmf = fSmf = true;
		
		omp_set_nested(1);
#pragma omp parallel sections
		{
#pragma omp section
			{
				relaxationFlag = DynamicRelaxation(_timeInterval, dim, grain, wLoadPiston, isCohesion, b, recordOutput, type, stopMechRatio);
			}
#pragma omp section
			{
				if (uint_recordOutput == 2)
					RunVisualPostProcessing(isCohesion, dim);
				else
					RunPostProcessing(isCohesion, dim);
			}
		}

		/* record final state if no record output option */
		if (!recordOutput)
			PostProcessing(isCohesion, dim, _pst, wl);

		RuntimeLog(runTimeFile, true, "COMPLETED.", 0.0, false, false);

		return 0;
	}

	return 0;
}
