#include "cs_sample.hpp"
#include "../../graphics/src/visualizer.hpp"

#include <time.h> 
#include <locale>
#include <utility>
#include <codecvt>

//---------------------------------------------------------------------------------------------------------------------------------

Sample::Sample()
{
}

//---------------------------------------------------------------------------------------------------------------------------------

Sample::~Sample()
{
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::Init(const char* runtimeFileName, const char* paramFile)
{
    runTimeFile = runtimeFileName;
    //paramFile = parFile;

    RuntimeLog(runTimeFile, false, "", 1.0, true, false);

    srand((unsigned)time(NULL));

	RuntimeLog(runTimeFile, true, "INSTANTIATE SAMPLE OBJECT", 0.0, false, false);
	RuntimeLog(runTimeFile, true, "AND READING INITIAL DATA (PARAMETERS) FROM FILE", 0.0, false, false);

    ReadParametersFromFile( paramFile );

    ProjectManager();

    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::SetSample( int cellsNum )
{	

    RuntimeLog(runTimeFile, true, "SEED CELLS, cellsNum = ", cellsNum, false, false, true);

	SeedCells(_rMin, _rMax, _sizeX, _sizeY, cellsNum, aggregates, stdDev);

	RuntimeLog(runTimeFile, true, "SET WALLS.", 0.0, false, false);

	SetWalls(_sizeX, _sizeY);

    RuntimeLog(runTimeFile, true, "CLEAR FORCES.", 0.0, false, false);

	ClearForces(aggregates);

    RuntimeLog(runTimeFile, true, "SET BONDS.", 0.0, false, false);

	SetBondsMatrix(b, aggregates);

	FindBonds(b, aggregates, shape, scale);

    RuntimeLog(runTimeFile, true, "SAMPLE INITIATION HAS BEEN COMPLETED.", 0.0, false, false);
	
    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::RunSimulation(double time)
{
	if (recordOutput)
		relaxationFlag = 1;
	else
		relaxationFlag = 0;

	snapshotCounter = snapshotCounterG = snapshotCounterF = snapshotCounterB = 0;
	grnCount = wlCount = fCount = bCount = 0;
	grnTop = wlTop = fTop = bTop = true;
	snapshotBuffer = 200;
	bSmf = grnSmf = wlSmf = fSmf = true;

	if (recordOutput){
		omp_set_nested(1);
	#pragma omp parallel sections
		{
	#pragma omp section
			{
				relaxationFlag = DynamicRelaxation(time, aggregates, b, recordOutput, type);
			}
	#pragma omp section
			{
				if (uint_recordOutput == 2){
					RunVisualPostProcessing();
				}
				else{
					RunPostProcessing();
				}
			}
		}
	}
	else{
		relaxationFlag = DynamicRelaxation(time, aggregates, b, recordOutput, type);
	}
			

    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::SimulateGrowth2( std::vector<Cell> &cells, double time, bool record )
{
	RuntimeLog(runTimeFile, true, "Initiated growth in:  ", 0.0, false, false);

    double dt = 0.05;

	//size_t stepNum = (size_t)ceil(time / dt);
	size_t step = 0;
	size_t outputPeriod = 10;

	bool newCell = false;

    for ( size_t t = 0; t < time; t++ ){

        size_t sizeCluster = cells.size();

		//std::cout<<"growth, at "<<sizeCluster<<std::endl;

        for ( size_t i = 0; i < sizeCluster; i++ ){
            cells[i].growth(cells,dt);
        }
		
		//std::cout<<"new cell, at "<<sizeCluster<<std::endl;

        for ( size_t i = 0; i < sizeCluster; i++ ){
            if (cells[i].isElongating && !cells[i].isLinked){
                cells.push_back( cells[i].addcell() );
				newCell = true;
			}
        }

/* 		if (newCell){
			RuntimeLog(runTimeFile, false, "-------------- New iteration (time): ", t, false, false, true);
			RuntimeLog(runTimeFile, false, "Total cells: ", cells.size(), false, false, true);
			for ( size_t i = 0; i < sizeCluster; i++ ){
				//if (cells[i].isElongating && !cells[i].isLinked){
					RuntimeLog(runTimeFile, false, " cell number: ", i+1, false, false, true);
					RuntimeLog(runTimeFile, false, "            cell ID: ", cells[i].name, false, false, true);
					RuntimeLog(runTimeFile, false, "          conf. num: ", cells[i].conf_num, false, false, true);
					RuntimeLog(runTimeFile, false, "         elongation: ", cells[i].currentElongation, false, false, true);
					RuntimeLog(runTimeFile, false, "         is virtual: ", cells[i].isVirtual, false, false, true);
				//}
			}
			newCell = false;
		}
 */
		/* record the results of simulation */
		if (step%outputPeriod == 0 && record)
			GetCellsSnapshots(cells);

		step++;
    }

	RuntimeLog(runTimeFile, true, "End of growth.", 0.0, false, false);
	
/* 	RuntimeLog(runTimeFile, true, "Confinement values:  ", 0.0, false, false);
	for ( size_t i = 0; i < cells.size(); i++ ){
		RuntimeLog(runTimeFile, false, "cell no. ", cells[i].conf_degree, false, false, true);
	}
 */
    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::RunGrowthSimulation(double time)
{	
	RuntimeLog(runTimeFile, true, "Growth simulation.", 0.0, false, false);

	//std::cout<<"recordOutput = "<<recordOutput<<"; uint_recordOutput = "<<uint_recordOutput<<std::endl; exit(1);

	if (recordOutput)
		relaxationFlag = 1;
	else
		relaxationFlag = 0;

	snapshotCounter = snapshotCounterG = snapshotCounterF = snapshotCounterB = 0;
	grnCount = wlCount = fCount = bCount = 0;
	grnTop = wlTop = fTop = bTop = true;
	snapshotBuffer = 200;
	bSmf = grnSmf = wlSmf = fSmf = true;

	if (recordOutput){
		omp_set_nested(1);
	#pragma omp parallel sections
		{
	#pragma omp section
			{
				//relaxationFlag = SimulateGrowth2(aggregates, time, recordOutput);
				relaxationFlag = SimulateGrowth(time, aggregates, b, recordOutput, type);
			}
	#pragma omp section
			{
				if (uint_recordOutput == 2){
					RunVisualPostProcessing();
				}
				else{
					RunPostProcessing();
				}
			}
		}
	}
	else{
			//relaxationFlag = SimulateGrowth2(aggregates, time, recordOutput);
			relaxationFlag = SimulateGrowth(time, aggregates, b, recordOutput, type);
	}
			

    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::ReportExit( int exitValue )
{
    RuntimeLog("default_runtime.log", false, "", 1.0, true, false);

    if (exitValue == 1){
        RuntimeLog("default_runtime.log", true, "You forget to enter the name of file with a list of parameters !", 0.0, false, false);
        RuntimeLog("default_runtime.log", false, "Example:  c:>'path'/dem.exe 'path'/parameters.dat [optional input] ", 0.0, false, false, true);
        RuntimeLog("default_runtime.log", true, "EXIT with error !", 0.0, false, true);
    }
    else {
        RuntimeLog("default_runtime.log", true, "Too much input parameters ! It should not exceed 3.", 0.0, false, false);
		RuntimeLog("default_runtime.log", true, "EXIT with error !", 0.0, false, true);
    }
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::ReadParametersFromFile(const char* paramFileName)
{
	std::ifstream rFile;
	rFile.open(paramFileName);
	std::string line;
	std::vector <std::string> lines;

	if (rFile.fail()) {
		RuntimeLog(runTimeFile, true, "Fail in openning parameters file! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
		return EXIT_FAILURE;
	}
	else {
		while (!rFile.eof()) {
			std::getline(rFile, line);
			lines.push_back(line);
		}
	}
	rFile.close();

	_rMin = stod(lines.at(1));
	_rMax = stod(lines.at(3));
	_sizeX = stod(lines.at(5));
	_sizeY = stod(lines.at(7));
	samplTime = stod(lines.at(9));
	runTime = stod(lines.at(11));
	Grav = stod(lines.at(13));
	kn = stod(lines.at(15));
	kt = stod(lines.at(17));
	b_kn = stod(lines.at(19));
	b_kt = stod(lines.at(21));
	sigma = stod(lines.at(23));
	tau = stod(lines.at(25));
	rho = stod(lines.at(27));
	rho_env = stod(lines.at(29));
	stdDev = stod(lines.at(31));
	shape = stod(lines.at(33));
	scale = stod(lines.at(35));
	recordOutput = stoi(lines.at(37));
	uint_recordOutput = stoi(lines.at(37));
	deloneRapture = stoi(lines.at(39));

	// initialise some parameters
	mu = 0.3; // Sliding friction
	mur = 0.01; // Rolling friction
	nur = kt*_rMin; // Rolling "damping"
	Adh = 0.;// Adhesion
	grainNum = 0;
	beta = 0.75; /* damping coefficient */

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

#ifdef _WIN64

	int Sample::ProjectManager()
	{
		std::string projDirVtk;
		std::string projDirState;
		TCHAR NPath[MAX_PATH];

		int success = GetCurrentDirectory(MAX_PATH, NPath);
		if (success == 0) {
			RuntimeLog(runTimeFile, true, "Cannot obtain the current directory path! EXIT with error !", 0.0, false, true);
			return EXIT_FAILURE;
		}

		projDirMain = NPath;
		std::string _type(type.length(), L' ');
		std::copy(type.begin(), type.end(), _type.begin());
		projDirMain += "\\";
		projDirMain += _type;
		projDirVtk = projDirMain;
		projDirVtk += "\\vtk_files";

		if ( !DirExists(projDirMain) ) {
			success = CreateDirectory(projDirMain.c_str(), NULL);
			if (success == 0) {
				RuntimeLog(runTimeFile, true, "Failed to create the main project directory! EXIT with error !", 0.0, false, true);
				return EXIT_FAILURE;
			}
		}

		if ( !DirExists(projDirVtk) ) {
			success = CreateDirectory(projDirVtk.c_str(), NULL);
			if (success == 0) {
				RuntimeLog(runTimeFile, true, "Failed to create the main project vtk sub-directory! EXIT with error !", 0.0, false, true);
				return EXIT_FAILURE;
			}
		}

		return 0;
	}

//---------------------------------------------------------------------------------------------------------------------------------

	bool Sample::DirExists(const std::string& dirName_in)
	{
		DWORD ftyp = GetFileAttributesA(dirName_in.c_str());
		if (ftyp == INVALID_FILE_ATTRIBUTES)
			return false;

		if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
			return true;

		return false;
	}

	std::string ws2s(const std::wstring& wstr)
	{
		using convert_typeX = std::codecvt_utf8<wchar_t>;
		std::wstring_convert<convert_typeX, wchar_t> converterX;

		return converterX.to_bytes(wstr);
	}

#endif

//---------------------------------------------------------------------------------------------------------------------------------

#ifdef linux

	int Sample::ProjectManager(){

		std::string projDirVtk;
		std::string projDirState;

		char NPath[PATH_MAX];
		if (getcwd(NPath, sizeof(NPath)) != NULL) {
			//
		} else {
			RuntimeLog(runTimeFile, true, "Cannot obtain the current directory path! EXIT with error !", 0.0, false, true);
			return EXIT_FAILURE;
		}

		projDirMain = NPath;
		std::string _type(type.length(), L' ');
		std::copy(type.begin(), type.end(), _type.begin());
		projDirMain += "/";
		projDirMain += _type;
		projDirVtk = projDirMain;
		projDirVtk += "/vtk_files";

		int success;

		if ( !DirExists(projDirMain) ) {
			success = mkdir(projDirMain.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if (success != 0) {
				RuntimeLog(runTimeFile, true, "Failed to create the main project directory! EXIT with error !", 0.0, false, true);
				return EXIT_FAILURE;
			}
		}

		if ( !DirExists(projDirVtk) ) {
			success = mkdir(projDirVtk.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if (success != 0) {
				RuntimeLog(runTimeFile, true, "Failed to create the main project vtk sub-directory! EXIT with error !", 0.0, false, true);
				return EXIT_FAILURE;
			}
		}

		return 0;
	}

//---------------------------------------------------------------------------------------------------------------------------------

	bool Sample::DirExists(const std::string& dirName_in)
	{
		struct stat sb;
		
		if (stat(dirName_in.c_str(), &sb) == -1) {
			RuntimeLog(runTimeFile, true, "STAT error while checking directory! Directory does not exists.", 0.0, false, false);
			//RuntimeLog(runTimeFile, true, "STAT error while checking directory!", 0.0, false, true);
			//exit(EXIT_FAILURE);
		}
		else{
			if(sb.st_mode & S_IFDIR != 0)
				return true;
		}

		return false;
	}

#endif

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::Finalize( std::string fileNameOutput )
{
    /* record final state if no record output option */
    if (!recordOutput)
        PostProcessing();
    
    RuntimeLog(runTimeFile, true, "COMPLETED.", 0.0, false, false);		
    RuntimeLog(runTimeFile, true, "FINALIZING THE PROJECT. SAVING *.SMPL FILE", 0.0, false, false);

    FinalizeSampleProject(fileNameOutput, aggregates, b, wl);

    RuntimeLog(runTimeFile, true, "SAMPLE PROJECT COMPLETED.", 0.0, false, true);

    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::FinalizeSampleProject(std::string paramFileName, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, std::vector<WALLS> &wl)
{
	RuntimeLog(runTimeFile, true, "Resulting cells' posotions: ", 0.0, false, false);
	for (size_t i = 0; i < gr.size(); i++){
		RuntimeLog(runTimeFile, false, "cell no ", i, false, false, true);
		RuntimeLog(runTimeFile, false, "R = ", gr[i].R, false, false, true);
		RuntimeLog(runTimeFile, false, "x = ", gr[i].x, false, false, true);
		RuntimeLog(runTimeFile, false, "y = ", gr[i].x, false, false, true);
		RuntimeLog(runTimeFile, false, "z = ", gr[i].x, false, false, true);
	}

	std::ofstream out(paramFileName, std::ios::out | std::ofstream::binary);
	if (!out) {
		RuntimeLog(runTimeFile, true, "Problem with *.SMPL file !EXIT with error !", 0.0, false, true);
		return EXIT_FAILURE;
	}

	size_t numOfGrains = 0;
	try {
		SmallGrain *a; a = (SmallGrain *)malloc(sizeof(SmallGrain));
		if (a == NULL) throw 40;
		for (auto const &e : gr) {
			a->name = e.name;
			a->I = e.I;
			a->m = e.m;
			a->R = e.R;
			a->rhoInd = e.rhoInd;
			a->x = e.x;
			a->y = e.y;
			a->z = e.z;
			out.write(reinterpret_cast<char*> (a), sizeof(SmallGrain));
			numOfGrains++;
		}
		free(a);
		RuntimeLog(runTimeFile, false, "Number of recorded cells at *.SMPL file", numOfGrains, false, false, true);
	}
	catch (int ex) {
		RuntimeLog(runTimeFile, true, "Problem with allocating memory for *.SMPL file !EXIT with error !", 0.0, false, true);
		exit(EXIT_FAILURE);
	}
	catch (std::exception const& e) {
		RuntimeLog(runTimeFile, true, e.what(), 0.0, false, true);
		exit(EXIT_FAILURE);
	}

	out.close();

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::RestoreSample(const char* paramFileName)
{
	RuntimeLog(runTimeFile, true, "LOADING FABRICATED SAMPLE.", 0.0, false, false);

	InitializeSimulationProject(paramFileName, aggregates, b, wl);

    grainNum = aggregates.size();

    SetWalls(_sizeX, _sizeY);
    
    ClearForces(aggregates);
    
    SetBondsMatrix(b, aggregates);
    
    FindBonds(b, aggregates, shape, scale);

    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::InitializeSimulationProject(const char* paramFileName, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, std::vector<WALLS> &wl)
{
	std::ifstream in(paramFileName, std::ios::in | std::ios::binary);
	if (!in) {
		RuntimeLog(runTimeFile, true, "Problem with *.SMPL file !EXIT with error !", 0.0, false, true);
		return EXIT_FAILURE;
	}

	in.seekg(0, in.end);
	size_t length = in.tellg() / sizeof(SmallGrain);
	in.seekg(0, in.beg);

	RuntimeLog(runTimeFile, false, "Number of recorded cells at *.SMPL file", length, false, false, true);
	
/* 	for (size_t i = 0; i < length; i++)
		gr.push_back(Cell());	
 */
	try {
		SmallGrain *a; a = (SmallGrain *)malloc(sizeof(SmallGrain));
		if (a == NULL) throw 40;
		for (size_t i = 0; i < length; i++) {
			in.read(reinterpret_cast<char*> (a), sizeof(SmallGrain));
			glm::vec3 origin(a->x,a->y,a->z);
			gr.push_back(Cell(origin,a->R));
			gr[i].name = a->name;
			gr[i].I = a->I;
			gr[i].m = a->m;
			//gr[i].R = a->R;
			gr[i].rhoInd = a->rhoInd;
			gr[i].x = a->x;
			gr[i].y = a->y;
			gr[i].z = a->z;
		}
		free(a);

		RuntimeLog(runTimeFile, true, "Resulting cells' posotions: ", 0.0, false, false);
		for (size_t i = 0; i < gr.size(); i++){
			RuntimeLog(runTimeFile, false, "cell no ", i, false, false, true);
			RuntimeLog(runTimeFile, false, "R = ", gr[i].R, false, false, true);
			RuntimeLog(runTimeFile, false, "x = ", gr[i].x, false, false, true);
			RuntimeLog(runTimeFile, false, "y = ", gr[i].x, false, false, true);
			RuntimeLog(runTimeFile, false, "z = ", gr[i].x, false, false, true);
		}

	}
	catch (int ex) {
		RuntimeLog(runTimeFile, true, "Problem with allocating memory for reading from *.SMPL file !EXIT with error !", 0.0, false, true);
		exit(EXIT_FAILURE);
	}
	catch (std::exception const& e) {
		RuntimeLog(runTimeFile, true, e.what(), 0.0, false, true);
		exit(EXIT_FAILURE);
	}

	in.close();

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

double Sample::GetSampleHeight(double sizeY, std::vector <Cell> &gr)
{
	double maxHeight = 0;

	for (size_t i = 0; i < gr.size(); i++) {
		if ((gr[i].x + gr[i].R) >= wTop.lowX && (gr[i].x - gr[i].R) <= wTop.highX && (gr[i].z + gr[i].R) >= wTop.lowZ && (gr[i].z - gr[i].R) <= wTop.highZ) {
			if ((gr[i].y + gr[i].R) > maxHeight)
				maxHeight = gr[i].y + gr[i].R;
		}
	}

	return maxHeight;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::SetWalls(double sizeX, double sizeY)
{
	double sizeZ = sizeX;

	// bottom
	wBottom.wF = 0.; wBottom.wV = 0.; wBottom.wA = 0.; wBottom.wM = 1.;

	wBottom.lowX = 0.; wBottom.highX = sizeX;
	wBottom.lowY = 0.; wBottom.highY = 0.;
	wBottom.lowZ = 0.; wBottom.highZ = sizeZ;

	// top
	wTop.wF = 0.; wTop.wV = 0.; wTop.wA = 0.; wTop.wM = 1.;

	wTop.lowX = 0.; wTop.highX = sizeX;
	wTop.lowY = sizeY; wTop.highY = sizeY;
	wTop.lowZ = 0.; wTop.highZ = sizeZ;

	// left
	wLeft.wF = 0.; wLeft.wV = 0.; wLeft.wA = 0.; wLeft.wM = 1.;

	wLeft.lowX = 0.; wLeft.highX = 0.;
	wLeft.lowY = 0.; wLeft.highY = sizeY;
	wLeft.lowZ = 0.; wLeft.highZ = sizeZ;

	//right
	wRight.wF = 0.; wRight.wV = 0.; wRight.wA = 0.; wRight.wM = 1.;

	wRight.lowX = sizeX; wRight.highX = sizeX;
	wRight.lowY = 0.; wRight.highY = sizeY;
	wRight.lowZ = 0.; wRight.highZ = sizeZ;

	// front
	wFront.wF = 0.; wFront.wV = 0.; wFront.wA = 0.; wFront.wM = 1.;

	wFront.lowX = 0.; wFront.highX = sizeX;
	wFront.lowY = 0.; wFront.highY = sizeY;
	wFront.lowZ = 0.; wFront.highZ = 0.;

	// back
	wBack.wF = 0.; wBack.wV = 0.; wBack.wA = 0.; wBack.wM = 1.;

	wBack.lowX = 0.; wBack.highX = sizeX;
	wBack.lowY = 0.; wBack.highY = sizeY;
	wBack.lowZ = sizeZ; wBack.highZ = sizeZ;

}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::ForcesGrain3D(int i, int j, std::vector <Cell> &gr, double _dt)
{
	double rXji, rYji, rZji, rNorm, dn;
	// Compute distance between cells
	rXji = gr[i].x - gr[j].x;
	rYji = gr[i].y - gr[j].y;
	rZji = gr[i].z - gr[j].z;
	rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
	dn = rNorm - (gr[i].R + gr[j].R);
	double r = std::min(gr[i].R, gr[j].R);
	double Reff = 2.*gr[i].R*gr[j].R / (gr[i].R + gr[j].R); // effective radius
	//double nu = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;
	double nu_n = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;
	//double nu_s = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kt)*beta;

	if (dn <= 0.) // Contact
	{
		double Vxji, Vyji, Vzji, Vthji, Xn, Yn, Zn, Xt, Yt, Zt, Vn, Vt, Fn, Ft;
		// Local axes
		Xn = rXji / rNorm;
		Yn = rYji / rNorm;
		Zn = rZji / rNorm;
		Xt = -Yn;
		Yt = Xn;
		Zt = -1.0;

		// Compute the velocity at contact
		Vxji = gr[i].Vx - gr[j].Vx;
		Vyji = gr[i].Vy - gr[j].Vy;
		Vzji = gr[i].Vz - gr[j].Vz;
		Vn = Vxji*Xn + Vyji*Yn + Vzji*Zn;
		Vt = Vxji*Xt + Vyji*Yt + Vzji*Zt - ((gr[i].R - dn / 2)*gr[i].Vth + (gr[j].R - dn / 2)*gr[j].Vth);
		Vthji = (gr[i].Vth - gr[j].Vth)*Reff; // Compute rolling velocity (from Radjai & Dubois, 2011)

		// Compute force in local axes		
		Fn = -kn * dn - nu_n * Vn;
		if (Fn < 0.) Fn = 0.;

		Ft = fabs(kt*Vt);
		if (fabs(Ft) > mu*fabs(Fn)) Ft = mu*fabs(Fn);
		if (Vt > 0.) Ft = -Ft;

		// Compute rolling momen
		double Moment = -Vthji*nur;
		double threshold = mur*Reff*(Fn + Adh);

		if (threshold < 0.)
			threshold = 0.;

		if (fabs(Moment) > threshold){
			if (Vthji > 0.)
				Moment = -threshold;
			else
				Moment = threshold;
		}

		// Calculate sum of forces on i and j in global coordinates
		gr[i].Fx += Fn*Xn + Ft*Xt;
		gr[i].Fy += Fn*Yn + Ft*Yt;
		gr[i].Fz += Fn*Zn + Ft*Zt;
		gr[j].Fx += -Fn*Xn - Ft*Xt;
		gr[j].Fy += -Fn*Yn - Ft*Yt;
		gr[j].Fz += -Fn*Zn - Ft*Zt;
		gr[i].Fth += -Ft*gr[i].R - Moment;
		gr[j].Fth += -Ft*gr[j].R - Moment;

		// Add fn to pressure on cells i and j
		gr[i].p += Fn;
		gr[j].p += Fn;

	}
	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::ForcesWall3D(int i, int wall, std::vector <Cell> &gr)
{
	double nu = 0.;
	double r = gr[i].R;
	nu = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;

	double dn;
	//bottom wall, identifier == 1
	if (wall == 1)
	{
		dn = gr[i].y - gr[i].R - wBottom.lowY;
		if (dn < 0.)
		{
			double vn, fn;
			vn = gr[i].Vy;
			fn = -kn * dn - nu * fabs(vn);
			if (fn < 0.) fn = 0.;
			if (vn < 0.) fn = -kn * dn;
			// Update sum of forces on cells
			gr[i].Fy = gr[i].Fy + fn;
			// Add fn to pressure on cells i
			gr[i].p += fn;
		}
	}
	//top wall, identifier == 2
	else if (wall == 2)
	{
		if (gr[i].y > wTop.highY)
		{			
			dn = wTop.highY - (gr[i].y + gr[i].R);
			if (dn < 0.)
			{
				double vn, fn;
				vn = gr[i].Vy;
				fn = -kn * dn;
				if (vn < 0.) fn = -kn * dn - nu * fabs(vn);
				else if (fn < 0.) fn = 0.;
				gr[i].Fy = gr[i].Fy - fn;
				gr[i].p += fn;
			}			
		}
		if (gr[i].y <= wTop.highY)
		{
			dn = wTop.highY - (gr[i].y + gr[i].R);
			if (dn < 0.)
			{
				double vn, fn;
				vn = gr[i].Vy;
				fn = -kn * dn;
				if (vn < 0.) fn = -kn * dn - nu * fabs(vn);
				else if (fn < 0.) fn = 0.;
				gr[i].Fy = gr[i].Fy - fn;
				gr[i].p += fn;
			}
		}

	}
	//left wall, identifier == 3
	else if (wall == 3)
	{
		dn = gr[i].x - gr[i].R - wLeft.lowX;
		if (dn < 0.)
		{
			double vn, fn;
			vn = gr[i].Vx;
			fn = -kn * dn - nu * fabs(vn);
			if (fn < 0.) fn = 0.;
			if (vn < 0.) fn = -kn * dn;
			// Update sum of forces on cells
			gr[i].Fx = gr[i].Fx + fn;
			// Add fn to pressure on cells i
			gr[i].p += fn;
		}
	}
	//right wall, identifier == 4
	else if (wall == 4)
	{
		dn = wRight.highX - (gr[i].x + gr[i].R);
		if (dn < 0.)
		{
			double vn, fn;
			vn = gr[i].Vx;
			fn = -kn * dn;
			if (vn < 0.) fn = -kn * dn - nu * fabs(vn);
			else if (fn < 0.) fn = 0.;
			// Update sum of forces on cells
			gr[i].Fx = gr[i].Fx - fn;
			// Add fn to pressure on cells i
			gr[i].p += fn;
		}
	}
	//front wall, identifier == 5
	else if (wall == 5)
	{
		dn = gr[i].z - gr[i].R - wFront.lowZ;
		if (dn < 0.)
		{
			double vn, fn;
			// velocity
			vn = gr[i].Vz;
			fn = -kn * dn - nu * fabs(vn);
			if (fn < 0.) fn = 0.;
			if (vn < 0.) fn = -kn * dn;
			// Update sum of forces on cells
			gr[i].Fz = gr[i].Fz + fn;
			// Add fn to pressure on cells i
			gr[i].p += fn;
		}
	}
	//back wall, identifier == 6
	else if (wall == 6)
	{
		dn = wBack.highZ - (gr[i].z + gr[i].R);
		if (dn < 0.)
		{
			double vn, fn;
			vn = gr[i].Vz;
			fn = -kn * dn;
			if (vn < 0.) fn = -kn * dn - nu * fabs(vn);
			else if (fn < 0.) fn = 0.;
			// Update sum of forces on cells
			gr[i].Fz = gr[i].Fz - fn;
			// Add fn to pressure on cells i
			gr[i].p += fn;
		}
	}
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::ForcesGrainCohesive3D(int i, int j, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, double _dt)
{
	double rXji, rYji, rZji, rNorm, dn;
	double Vxji, Vyji, Vzji, Vth, Xn, Yn, Zn, Xt, Yt, Zt, Vn, Vt;

	double r = std::min(gr[i].R, gr[j].R);
	double nu = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;
	double nuBond = 2 * sqrt(b[i][j].M * b[i][j].kn)*beta;

	// Compute distance between cells
	rXji = gr[i].x - gr[j].x;
	rYji = gr[i].y - gr[j].y;
	rZji = gr[i].z - gr[j].z;
	rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
	dn = rNorm - (gr[i].R + gr[j].R);

	double dnBond;
	double Fn = 0.;
	double Ft = 0.;

	if (b[i][j].isActive == true) {
		
		dnBond = rNorm - (gr[i].R + gr[j].R + b[i][j].gap); // Bond overlap

		Vxji = gr[i].Vx - gr[j].Vx;  // sign convention: v<0 --> particles approach each other
		Vyji = gr[i].Vy - gr[j].Vy; // sign convention: v<0 --> particles approach each other
		Vzji = gr[i].Vz - gr[j].Vz;
		Vth = gr[i].Vth - gr[j].Vth;

		Xn = rXji / rNorm;
		Yn = rYji / rNorm;
		Zn = rZji / rNorm;

		Xt = -Yn;
		Yt = Xn;
		Zt = -1.0;

		Vn = Vxji*Xn + Vyji*Yn + Vzji*Zn; // sign convention: vn<0 --> normal distance is decreasing
		Vt = Vxji*Xt + Vyji*Yt + Vzji*Zt - ((gr[i].R - dn / 2)*gr[i].Vth + (gr[j].R - dn / 2)*gr[j].Vth);

		double anglei = gr[i].th + b[i][j].alpha0;
		double anglej = gr[j].th + b[j][i].alpha0;

		double angleiXY = gr[i].th + b[i][j].alpha0XY;
		double anglejXY = gr[j].th + b[j][i].alpha0XY;

		double angleiYZ = gr[i].th + b[i][j].alpha0YZ;
		double anglejYZ = gr[j].th + b[j][i].alpha0YZ;

		// Bond elongation
		double I1I2X = -gr[i].R * cos(anglei) + gr[j].R * cos(anglej) - rXji;
		double I1I2Y = -gr[i].R * sin(anglei) + gr[j].R * sin(anglej) - rYji;
		double I1I2Z = -gr[i].R * cos(angleiYZ) + gr[j].R * cos(anglejYZ) - rZji;

		// Update global position of interaction points
		b[i][j].xI = gr[i].x + gr[i].R * cos(angleiXY);
		b[i][j].yI = gr[i].y + gr[i].R * sin(angleiXY);
		b[i][j].zI = gr[i].z + gr[i].R * cos(angleiYZ);

		b[j][i].xI = gr[j].x + gr[j].R * cos(anglejXY);
		b[j][i].yI = gr[j].y + gr[j].R * sin(anglejXY);
		b[i][j].zI = gr[j].z + gr[j].R * cos(angleiYZ);

		// Bond shearing distance
		double I1I2T = I1I2X*Xt + I1I2Y*Yt + I1I2Z*Zt;

		double gammaij = gr[i].th - gr[j].th;

		if ( b[i][j].rigid /* || (gr[i].conf_num == gr[j].conf_num) */ ) {
			Fn = 0.0;
			Ft = 0.0;
		}
		else {
			if (dn >= 0.) {
				Fn = -b[i][j].kn*dnBond - nuBond*Vn;
				if ( Fn > 0.0 ) Fn = 0.0;
			}
 			else if ( dn < 0.0 && dn >= -0.5*gr[i].R){
				Fn = -kn*dn - nu*Vn;
				if ( Fn < 0.0 ) Fn = 0.0;
			}
			else {
				// we are at the cell-cell interface:
 				double inter_const = 1.0*gr[i].R; // normal interface distance (between origins)
				double inter_d = std::fabs(dn) - inter_const;
				double inter_F = -b[i][j].kn*inter_d - nuBond*Vn*1.0;

				if (inter_d >= 0.0){
					if ( inter_F > 0.0 ) inter_F = 0.0;
				}
				else {
					if ( inter_F < 0.0 ) inter_F = 0.0;
				}
				// double f1 = -kn*dn - nu*Vn /* + b[i][j].kn*dn */;
				// if ( f1 < 0.0 ) f1 = 0.0;
				// f1 = 0.0; //debug
				double f2;
				if (std::fabs( gr[i].conf_degree - gr[j].conf_degree ) > 0.8){
					f2 = -kn*dn * ( gr[i].conf_degree - gr[j].conf_degree )*( gr[i].conf_degree - gr[j].conf_degree ) - nu*Vn*10.0;
 				}
				else
				{
					f2 = 0.0;
				}
				
				if ( f2 < 0.0 ) f2 = 0.0; // avoid attraction so far
				//f2 = 0.0;//debug
/*   				if (std::fabs( gr[i].conf_degree - gr[j].conf_degree ) > 0.8){
					RuntimeLog(runTimeFile, false, "cell i name:", gr[i].name, false, false, true);
					RuntimeLog(runTimeFile, false, "cell j name:", gr[j].name, false, false, true);
					RuntimeLog(runTimeFile, false, "delta:", std::fabs( gr[i].conf_degree - gr[j].conf_degree ), false, false, true);
					RuntimeLog(runTimeFile, false, "--------------------------", 0.0, false, false, true);
					RuntimeLog(runTimeFile, false, "cell name:", gr[i].name, false, false, true);
					RuntimeLog(runTimeFile, false, "age:", gr[i].age, false, false, true);
					RuntimeLog(runTimeFile, false, "elongation:", gr[i].pastElongation, false, false, true);
					RuntimeLog(runTimeFile, false, "isLinked:", gr[i].isLinked, false, false, true);
					RuntimeLog(runTimeFile, false, "linkedID:", gr[i].linkedID, false, false, true);
					RuntimeLog(runTimeFile, false, "isElongating:", gr[i].isElongating, false, false, true);
					RuntimeLog(runTimeFile, false, "growth rate:", gr[i].growthRate, false, false, true);
					RuntimeLog(runTimeFile, false, "conf degree:", gr[i].conf_degree, false, false, true);
					RuntimeLog(runTimeFile, false, "conf degree_upd:", gr[i].conf_degree_upd, false, false, true);
					RuntimeLog(runTimeFile, false, "conf num:", gr[i].conf_num, false, false, true);
					RuntimeLog(runTimeFile, false, "--------------------------", 0.0, false, false, true);
					RuntimeLog(runTimeFile, false, "cell name:", gr[j].name, false, false, true);
					RuntimeLog(runTimeFile, false, "age:", gr[j].age, false, false, true);
					RuntimeLog(runTimeFile, false, "elongation:", gr[j].pastElongation, false, false, true);
					RuntimeLog(runTimeFile, false, "isLinked:", gr[j].isLinked, false, false, true);
					RuntimeLog(runTimeFile, false, "linkedID:", gr[j].linkedID, false, false, true);
					RuntimeLog(runTimeFile, false, "isElongating:", gr[j].isElongating, false, false, true);
					RuntimeLog(runTimeFile, false, "growth rate:", gr[j].growthRate, false, false, true);
					RuntimeLog(runTimeFile, false, "conf degree:", gr[j].conf_degree, false, false, true);
					RuntimeLog(runTimeFile, false, "conf degree_upd:", gr[j].conf_degree_upd, false, false, true);
					RuntimeLog(runTimeFile, false, "conf num:", gr[j].conf_num, false, false, true);
					RuntimeLog(runTimeFile, false, "--------------------------", 0.0, false, false, true);
				}
 */
				/* debug */
/*  				RuntimeLog(runTimeFile, false, "Number of cells:", gr.size(), false, false, true);
				for (size_t i = 0; i < gr.size(); i++){
					RuntimeLog(runTimeFile, false, "cell name:", gr[i].name, false, false, true);
					RuntimeLog(runTimeFile, false, "age:", gr[i].age, false, false, true);
					RuntimeLog(runTimeFile, false, "elongation:", gr[i].pastElongation, false, false, true);
					RuntimeLog(runTimeFile, false, "isLinked:", gr[i].isLinked, false, false, true);
					RuntimeLog(runTimeFile, false, "linkedID:", gr[i].linkedID, false, false, true);
					RuntimeLog(runTimeFile, false, "isElongating:", gr[i].isElongating, false, false, true);
					RuntimeLog(runTimeFile, false, "growth rate:", gr[i].growthRate, false, false, true);
					RuntimeLog(runTimeFile, false, "conf degree:", gr[i].conf_degree, false, false, true);
					RuntimeLog(runTimeFile, false, "conf degree_upd:", gr[i].conf_degree, false, false, true);
					RuntimeLog(runTimeFile, false, "conf num:", gr[i].conf_num, false, false, true);
					RuntimeLog(runTimeFile, false, "--------------------------", 0.0, false, false, true);
				}
 */				

				Fn = f2+inter_F;
			}

			/*Ft = -b[i][j].kt*I1I2T;*/ //this variant does not works, needs debugging!
			Ft = fabs(b[i][j].kt*Vt);
			if (fabs(Ft) > mu*fabs(Fn)) Ft = mu*fabs(Fn);
			if (Vt > 0.) Ft = -Ft;
			Ft = 0.0;// debug
		}

		//Fn = Ft = 0.0; // debug

		double couple = -b[i][j].M*gammaij;

		b[i][j].contactForce = fabs(Ft) + fabs(Fn);

		gr[i].Fx += Fn*Xn + Ft*Xt;
		gr[i].Fy += Fn*Yn + Ft*Yt;
		gr[i].Fz += Fn*Zn + Ft*Zt;
		gr[j].Fx += - Fn*Xn - Ft*Xt;
		gr[j].Fy += - Fn*Yn - Ft*Yt;
		gr[j].Fz += - Fn*Zn - Ft*Zt;

		gr[i].Fth += -Ft*gr[i].R - couple;
		gr[j].Fth += -Ft*gr[j].R + couple;

		double sigmaStress = fabs(Fn);
		double tauStress = fabs(Ft);

		// Check yield criterion
		double rupture = -1;

		if (deloneRapture)
			rupture = (Ft / b[i][j].tau)*(Ft / b[i][j].tau) + (couple / b[i][j].YM)*(couple / b[i][j].YM) + (Fn / b[i][j].sigma) - 1;
		else {
			if (sigmaStress >= b[i][j].sigma || tauStress >= b[i][j].tau)
				rupture = 0.;
		}

		if (rupture >= 0.) {
			b[i][j].isActive = b[j][i].isActive = false;
		}
		
		return;
	}
 	else
		ForcesGrain3D(i, j, gr, _dt);

	return;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::IdentifyCohesiveBonds(int i, int j, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, double *randValues)
{
	double xOiOj = gr[j].x - gr[i].x;
	double yOiOj = gr[j].y - gr[i].y;
	double zOiOj = gr[j].z - gr[i].z;
	double OiOj = 0.;
	
	OiOj = sqrt(xOiOj*xOiOj + yOiOj*yOiOj + zOiOj*zOiOj);
	
	double dn = OiOj - gr[i].R - gr[j].R;
	double alphaXY;
	double alphaYZ;
	//double alpha0;

	double gap = 0.1*std::min(gr[i].R, gr[j].R);

	if ((dn - gap) <= 0.0f) {

		if (xOiOj == 0.) {
			if (yOiOj > 0) alphaXY = M_PI / 2.;
			else alphaXY = -M_PI / 2.;
		}
		else {
			alphaXY = atan2(yOiOj, xOiOj);
		}

		if (zOiOj == 0.) {
			if (yOiOj > 0) alphaYZ = M_PI / 2.;
			else alphaYZ = -M_PI / 2.;
		}
		else {
			alphaYZ = atan2(yOiOj, xOiOj);
		}

		b[i][j].alpha0 = alphaXY - gr[i].th;
		b[j][i].alpha0 = alphaXY + M_PI - gr[j].th;

		b[i][j].alpha0YZ = alphaYZ - gr[i].th;
		b[j][i].alpha0YZ = alphaYZ + M_PI - gr[j].th;

		b[i][j].isActive = true;
		b[i][j].gap = 0.0;

		b[i][j].rigid = false;
		if ( gr[i].isLinked && gr[j].isLinked ){
			if (gr[i].name == gr[j].linkedID)
				b[i][j].rigid = true; // means the cell is growing (elongating)
		}		

		b[i][j].kn = b_kn * randValues[0];  // Tensile
		b[i][j].kt = b_kt * randValues[1];
		b[i][j].M = rho*gap*M_PI*(gap/2)*(gap / 2);

		b[i][j].sigma = sigma;
		b[i][j].tau = tau;
		//b[i][j].cohesion = cohesion;
		//b[i][j].frictionAngle = frictionAngle;

		b[i][j].xI = gr[i].x + gr[i].R*cos(alphaXY);  // x-coordinate of interaction point on g[i]
		b[i][j].yI = gr[i].y + gr[i].R*sin(alphaXY);  // y-coordinate of interaction point on g[i]
		b[i][j].zI = gr[i].z + gr[i].R*cos(alphaYZ);

		b[j][i].xI = gr[j].x + gr[j].R*cos(alphaXY + M_PI);  // x-coordinate of interaction point on g[j]
		b[j][i].yI = gr[j].y + gr[j].R*sin(alphaXY + M_PI);  // y-coordinate of interaction point on g[j]
		b[j][i].zI = gr[i].z + gr[i].R*cos(alphaYZ + M_PI);
	}
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::SetBondsMatrix(std::vector<std::vector<Bond>> &b, std::vector <Cell> &gr)
{
	b.clear();
	for (size_t i = 0; i < gr.size(); i++) {
		std::vector<Bond> vectb;
		for (size_t j = 0; j < gr.size(); j++) {
			Bond b0;
			vectb.push_back(b0);
		}
		b.push_back(vectb);
	}
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::FindBonds(std::vector<std::vector<Bond>> &b, std::vector <Cell> &gr, double shape, double scale)
{
	double randValues[2];
	bool notWeibull = false;
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	if (shape == 0.0 || scale == 0.0) {
		shape = 1.0;
		scale = 1.0;
		notWeibull = true;
		randValues[0] = 1.0;
		randValues[1] = 1.0;
	}
	std::weibull_distribution<double> distribution(shape, scale);

	for (size_t i = 0; i < gr.size(); i++) {
		for (size_t j = i + 1; j < gr.size(); j++) {
			b[i][j].isActive = false;
			if (!notWeibull) {
				randValues[0] = distribution(generator);
				randValues[1] = distribution(generator);
			}
			IdentifyCohesiveBonds(i, j, gr, b, randValues);
		}
	}
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::IsBallInDomain(std::vector <Cell> &gr)
{
	std::vector <Cell>::iterator p;
	bool cond;

	for (size_t i = 0; i < gr.size(); i++) {
		cond = (gr[i].y > _sizeY * 3) || (gr[i].x < wLeft.lowX) || (gr[i].x > wRight.lowX) || (gr[i].z < wFront.lowZ) || (gr[i].z > wBack.lowZ);
		if (cond) {
			p = gr.begin();
			p += i;
			gr.erase(p);
		}
	}
	grainNum = gr.size();
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::ClearForces(std::vector <Cell> &gr)
{
	for (size_t i = 0; i < gr.size(); i++)
	{
		gr[i].Ath = 0.;
		gr[i].Ax = 0.;
		gr[i].Ay = 0.;
		gr[i].Az = 0.;
		gr[i].Fth = 0.;
		gr[i].Fx = 0.;
		gr[i].Fy = 0.;
		gr[i].Fz = 0.;
		gr[i].p = 0.;
		gr[i].th = 0.;
		gr[i].Vth = 0.;
		gr[i].Vx = 0.;
		gr[i].Vy = 0.;
		gr[i].Vz = 0.;
		gr[i].normOfSumF = 0.;
		gr[i].sumOfNormF = 0.;
	}
}

//---------------------------------------------------------------------------------------------------------------------------------

double Sample::Get_dt()
{
	double _kn = std::max(std::max(kn, b_kn), std::max(kt, b_kt));

	return sqrt(rho*M_PI*_rMin*_rMin*_rMin / _kn) *0.001;
}

//---------------------------------------------------------------------------------------------------------------------------------

bool Sample::CheckProjType()
{
	std::string smpl("sample");
	std::string sml("simulation");

	return sml.compare(type);
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::SimulateGrowth( double timeInterval, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, bool record, std::string prType )
{
	//RuntimeLog(runTimeFile, true, "Initiated growth in:  ", 0.0, false, false);
	//RuntimeLog(runTimeFile, true, "End of growth.", 0.0, false, false);

	double dt = Get_dt(); /* time interval for integration, or in other words - time increment. */
	double mechRatio = 1.0e10;

	size_t wallsNum = 6; // number of walls

	size_t stepNum = (size_t)ceil(timeInterval / dt);

	int mRatTimeIntrv = (int)ceil(0.03 / dt);
	int updContLst = (int)ceil(0.5 / sqrt(dt));
	
	int outputPeriod = updContLst * 1; /* record output at the same time as updating the list of cells contacts. */

	size_t step = 0;

	std::string fileLocation;

	#ifdef _WIN64
	fileLocation = prType + "\\" + "Equilibrium_history.txt";
	#endif

	#ifdef linux
	fileLocation = prType + "/" + "Equilibrium_history.txt";
	#endif

	std::ofstream wrtMRatio;

	std::vector <double> _mechRatio_a; // temporal storage to collect mechRatio values.
	double _mechRatio_i = 0.0; /*mech ratio at previous step*/

	RuntimeLog(runTimeFile, true, "Starting the dynamic relaxation:", 0, false, false);
	RuntimeLog(runTimeFile, false, "Number of cells:", gr.size(), false, false, true);
	RuntimeLog(runTimeFile, false, "Critical timestep:", dt, false, false, true);
	RuntimeLog(runTimeFile, false, "Interval for contact list updates:", updContLst, false, false, true);
	
	if ( CheckProjType() ) {
		wrtMRatio.open(fileLocation, std::ios::out);
	}

	RuntimeLog(runTimeFile, false, "Interval for error value updates:", mRatTimeIntrv, false, false, true);
	RuntimeLog(runTimeFile, false, "Allowed maximum value of relaxation cycles:", stepNum, false, false, true);
	RuntimeLog(runTimeFile, true, "CALCULATING ...  ", 0, false, false);

	ClearForces(gr);

	size_t cells = gr.size();
	//grainNum = gr.size();

	while (step < stepNum)
	{
		/*-------- begining of the growth stage --------*/
		size_t sizeCluster = gr.size();

        for ( size_t i = 0; i < sizeCluster; i++ )
            gr[i].growth(gr,dt);
		
        for ( size_t i = 0; i < sizeCluster; i++ ){
            if (gr[i].isElongating && !gr[i].isLinked){
                gr.push_back( gr[i].addcell() );
			}        
		}
		/*-------- end of the growth stage -----------*/

		cells = gr.size();

		/* debug output */
/*  		if ( step%size_t(stepNum/timeInterval) == 0 ) {
			RuntimeLog(runTimeFile, true, "Current Time:", step*dt, false, false, true);
			RuntimeLog(runTimeFile, false, "Number of cells:", gr.size(), false, false, true);
			for (size_t i = 0; i < gr.size(); i++){
				RuntimeLog(runTimeFile, false, "cell name:", gr[i].name, false, false, true);
				RuntimeLog(runTimeFile, false, "age:", gr[i].age, false, false, true);
				RuntimeLog(runTimeFile, false, "elongation:", gr[i].pastElongation, false, false, true);
				RuntimeLog(runTimeFile, false, "isLinked:", gr[i].isLinked, false, false, true);
				RuntimeLog(runTimeFile, false, "linkedID:", gr[i].linkedID, false, false, true);
				RuntimeLog(runTimeFile, false, "isElongating:", gr[i].isElongating, false, false, true);
				RuntimeLog(runTimeFile, false, "growth rate:", gr[i].growthRate, false, false, true);
				RuntimeLog(runTimeFile, false, "conf degree:", gr[i].conf_degree, false, false, true);
				RuntimeLog(runTimeFile, false, "conf degree_upd:", gr[i].conf_degree_upd, false, false, true);
				RuntimeLog(runTimeFile, false, "conf num:", gr[i].conf_num, false, false, true);
				RuntimeLog(runTimeFile, false, "--------------------------", 0.0, false, false, true);
			}
		}
 */	

		SetBondsMatrix(b, gr);
		FindBonds(b, gr, shape, scale);

		//if (step % updContLst == 0)
			GetContactsList(gr, b);

		/* Velocity Verlet's prediction Step */
		for (size_t i = 0; i < cells; i++) {

			// if ( (gr[i].Ath != 0.0 || gr[i].Ax != 0.0 || gr[i].Ay != 0.0 || gr[i].Az != 0.0) || (gr[i].Vth != 0.0 || gr[i].Vx != 0.0 || gr[i].Vy != 0.0 || gr[i].Vz != 0.0) || (gr[i].Fth != 0.0 || gr[i].Fx != 0.0 || gr[i].Fy != 0.0 || gr[i].Fz != 0.0) ){
			//  	std::cout<<"not zero for cell no = "<< i<<std::endl;
			//  	std::cout<< "Ath= " <<gr[i].Ath <<"; Ax= "<< gr[i].Ax<<"; Ay= " << gr[i].Ay <<"; Az= " << gr[i].Az <<"; Vth= " << gr[i].Vth <<"; Vx= " << gr[i].Vx <<"; Vy= " << gr[i].Vy <<"; Vz= " << gr[i].Vz <<"; Fth= " << gr[i].Fth <<"; Fx = " << gr[i].Fx <<"; Fy= "<< gr[i].Fy <<"; Fz= "<< gr[i].Fz << std::endl;
			// }
			
			// update positions
			gr[i].x = gr[i].x + dt * gr[i].Vx + 0.5*dt*dt*gr[i].Ax;
			gr[i].y = gr[i].y + dt * gr[i].Vy + 0.5*dt*dt*gr[i].Ay;
			gr[i].z = gr[i].z + dt * gr[i].Vz + 0.5*dt*dt*gr[i].Az;
			gr[i].th = gr[i].th + dt * gr[i].Vth + 0.5*dt*dt*gr[i].Ath;

			/* correct cells' origins */
			gr[i].correct_origin();

			// update velocities
			gr[i].Vx = gr[i].Vx + 0.5 * dt * gr[i].Ax;
			gr[i].Vy = gr[i].Vy + 0.5 * dt * gr[i].Ay;
			gr[i].Vz = gr[i].Vz + 0.5 * dt * gr[i].Az;
			gr[i].Vth = gr[i].Vth + 0.5 * dt * gr[i].Ath;
		}

		wBottom.lowY = wBottom.lowY + dt * wBottom.wV + 0.5 * dt * dt * wBottom.wA;
		wTop.highY = wTop.highY + dt * wTop.wV + 0.5 * dt * dt * wTop.wA;
		wLeft.lowX = wLeft.lowX + dt * wLeft.wV + 0.5 * dt * dt * wLeft.wA;
		wRight.highX = wRight.highX + dt * wRight.wV + 0.5 * dt * dt * wRight.wA;

		wBottom.wV = wBottom.wV + 0.5 * dt * wBottom.wA;
		wTop.wV = /*wTop.wV + */0.5 * dt * wTop.wA; /* constant loading without acceleration*/
		wLeft.wV = wLeft.wV + 0.5 * dt * wLeft.wA;
		wRight.wV = wRight.wV + 0.5 * dt * wRight.wA;

		wFront.lowZ = wFront.lowZ + dt * wFront.wV + 0.5 * dt * dt * wFront.wA;
		wBack.highZ = wBack.highZ + dt * wBack.wV + 0.5 * dt * dt * wBack.wA;
		wFront.wV = wFront.wV + 0.5 * dt * wFront.wA;
		wBack.wV = wBack.wV + 0.5 * dt * wBack.wA;

		/* calculate mech ratio, only for 'sample' project */
		if (CheckProjType()) {

			double _sumOfNormF = 1.0;

			for (size_t i = 0; i < cells; i++) {
				gr[i].sumOfNormF = sqrt(gr[i].Fx*gr[i].Fx + gr[i].Fy*gr[i].Fy + gr[i].Fz*gr[i].Fz) + sqrt(gr[i].Fth*gr[i].Fth);
				_sumOfNormF += gr[i].sumOfNormF;
			}
			_mechRatio_a.push_back(_sumOfNormF);

			if (step % mRatTimeIntrv == 0) {
				double _mechRatio = accumulate(_mechRatio_a.begin(), _mechRatio_a.end(), 0.0) / _mechRatio_a.size();
				mechRatio = (abs(_mechRatio - _mechRatio_i) / _mechRatio);
				_mechRatio_i = _mechRatio;
				_mechRatio_a.clear();
				if (wrtMRatio)
					wrtMRatio << step << "  " << mechRatio << std::endl;
			}
		}

		/*Set sum of forces acting on cell to 0*/
		for (size_t i = 0; i < cells; i++) {
			gr[i].Fx = 0.;
			gr[i].Fy = 0.;
			gr[i].Fz = 0.;
			gr[i].Fth = 0.;
			gr[i].p = 0.;
		}

		/*Set confining pressure*/
		wBottom.wF = 0;
		wTop.wF = 0;
		wLeft.wF = 0;
		wRight.wF = 0;
		wFront.wF = 0;
		wBack.wF = 0;

		/*Force computation between cells*/
//#pragma omp parallel for 
		for (size_t i = 0; i < cells; i++) {
			for (auto j = gr[i].contactList.begin(); j != gr[i].contactList.end(); ++j) {
				ForcesGrainCohesive3D(i, *j, gr, b, dt);
			}
		}

		/*Force computation between cells and walls*/
//#pragma omp parallel for
/*  		for (size_t i = 0; i < cells; i++)
			for (size_t w = 1; w <= wallsNum; w++)
				ForcesWall3D(i, w, gr);
 */
		/*Calculate corrected acceleration*/
		for (size_t i = 0; i < cells; i++) {
			gr[i].Ax = gr[i].Fx / gr[i].m;
			gr[i].Ay = gr[i].Fy / gr[i].m - Grav + Grav * (rho_env/rho);
			gr[i].Az = gr[i].Fz / gr[i].m;
			gr[i].Ath = gr[i].Fth / gr[i].I;
		}

		wBottom.wA = wBottom.wF / wBottom.wM;
		wTop.wA = wTop.wF / wTop.wM;
		wLeft.wA = wLeft.wF / wLeft.wM;
		wRight.wA = wRight.wF / wRight.wM;


		wFront.wA = wFront.wF / wFront.wM;
		wBack.wA = wBack.wF / wBack.wM;

		/* Correction step for velocities */
		for (size_t i = 0; i < cells; i++) {
			gr[i].Vx = gr[i].Vx + 0.5 * dt * gr[i].Ax;
			gr[i].Vy = gr[i].Vy + 0.5 * dt * gr[i].Ay;
			gr[i].Vz = gr[i].Vz + 0.5 * dt * gr[i].Az;
			gr[i].Vth = gr[i].Vth + 0.5 * dt * gr[i].Ath;
		}

		wBottom.wV = wBottom.wV + 0.5 * dt*wBottom.wA;
		wTop.wV = /*wTop.wV +*/ 0.5 * dt*wTop.wA;/*0.5 * wTop.wA;*/ //constant loading with no acceleration
		wLeft.wV = wLeft.wV + 0.5 * dt*wLeft.wA;
		wRight.wV = wRight.wV + 0.5 * dt*wRight.wA;

		wFront.wV = wFront.wV + 0.5 * dt * wFront.wA;
		wBack.wV = wBack.wV + 0.5 * dt * wBack.wA;

		allWalls._wBack = wBack;
		allWalls._wBottom = wBottom;
		allWalls._wFront = wFront;
		allWalls._wLeft = wLeft;
		allWalls._wRight = wRight;
		allWalls._wTop = wTop;

		/* record the results of simulation */
		if (step%outputPeriod == 0 && record)
			GetCellsSnapshots(gr);
			//GetSnapshots(gr, b);

		if (step % outputPeriod == 0) {
			if (CheckProjType())
				RuntimeLog(runTimeFile, true, "Sample equilibrium error  ", mechRatio, false, false,true);
			//else
				//RuntimeLog(runTimeFile, true, "Completed, %", (int)((step * dt * 100) / timeInterval)+1, false, false, true);
		}

		/* writes the last snapshot of the sample if the project has no record option. */
		if (!record && step == (stepNum - 1))
			GetCellsSnapshots(gr);
			//GetSnapshots(gr, b);

		step++;
	}

    return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::DynamicRelaxation(double timeInterval, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, bool record, std::string prType) {

	//auto n_threads = std::thread::hardware_concurrency();

	double dt = Get_dt(); /* time interval for integration, or in other words - time increment. */
	double mechRatio = 1.0e10;

	size_t wallsNum = 6; // number of walls

	size_t stepNum = (size_t)ceil(timeInterval / dt);

	size_t mRatTimeIntrv = (int)ceil(0.03 / dt);
	size_t updContLst = (int)ceil(0.2 / sqrt(dt));
	
	size_t outputPeriod = updContLst * 1; /* record output at the same time as updating the list of cells contacts. */

	size_t step = 0;

	std::string fileLocation;

	#ifdef _WIN64
	fileLocation = prType + "\\" + "Equilibrium_history.txt";
	#endif

	#ifdef linux
	fileLocation = prType + "/" + "Equilibrium_history.txt";
	#endif

	std::ofstream wrtMRatio;

	std::vector <double> _mechRatio_a; // temporal storage to collect mechRatio values.
	double _mechRatio_i = 0.0; /*mech ratio at previous step*/

	RuntimeLog(runTimeFile, true, "Starting the dynamic relaxation:", 0, false, false);
	RuntimeLog(runTimeFile, false, "Number of cells:", gr.size(), false, false, true);
	RuntimeLog(runTimeFile, false, "Critical timestep:", dt, false, false, true);
	RuntimeLog(runTimeFile, false, "Total timesteps:", stepNum, false, false, true);
	RuntimeLog(runTimeFile, false, "Interval for contact list updates:", updContLst, false, false, true);
	
	if ( CheckProjType() ) {
		wrtMRatio.open(fileLocation, std::ios::out);
	}

	RuntimeLog(runTimeFile, false, "Interval for error value updates:", mRatTimeIntrv, false, false, true);
	RuntimeLog(runTimeFile, false, "Allowed maximum value of relaxation cycles:", stepNum, false, false, true);
	RuntimeLog(runTimeFile, true, "CALCULATING ...  ", 0, false, false);

	ClearForces(gr);

	size_t cells = gr.size();
	grainNum = gr.size();

#ifdef DEBUG
	//------------------------------------
	auto t1_start = chrono::steady_clock::now();
	auto t1_end = chrono::steady_clock::now();
	auto diff_1 = t1_end - t1_start;
	double t1 = chrono::duration <double, milli>(diff_1).count() / 1000; /* time in seconds*/
	t1 = 0.0;
	auto t2_start = chrono::steady_clock::now();
	auto t2_end = chrono::steady_clock::now();
	auto diff_2 = t2_end - t2_start;
	double t2 = chrono::duration <double, milli>(diff_2).count() / 1000; /* time in seconds*/
	t2 = 0.0;
	auto t3_start = chrono::steady_clock::now();
	auto t3_end = chrono::steady_clock::now();
	auto diff_3 = t3_end - t3_start;
	double t3 = chrono::duration <double, milli>(diff_3).count() / 1000; /* time in seconds*/
	t3 = 0.0;
	auto t4_start = chrono::steady_clock::now();
	auto t4_end = chrono::steady_clock::now();
	auto diff_4 = t4_end - t4_start;
	double t4 = chrono::duration <double, milli>(diff_4).count() / 1000; /* time in seconds*/
	t4 = 0.0;
	auto t5_start = chrono::steady_clock::now();
	auto t5_end = chrono::steady_clock::now();
	auto diff_5 = t5_end - t5_start;
	double t5 = chrono::duration <double, milli>(diff_5).count() / 1000; /* time in seconds*/
	t5 = 0.0;
	auto t6_start = chrono::steady_clock::now();
	auto t6_end = chrono::steady_clock::now();
	auto diff_6 = t6_end - t6_start;
	double t6 = chrono::duration <double, milli>(diff_6).count() / 1000; /* time in seconds*/
	t6 = 0.0;
	auto t7_start = chrono::steady_clock::now();
	auto t7_end = chrono::steady_clock::now();
	auto diff_7 = t7_end - t7_start;
	double t7 = chrono::duration <double, milli>(diff_7).count() / 1000; /* time in seconds*/
	t7 = 0.0;
	//------------------------------------
#endif

	while (step < stepNum)
	{

#ifdef DEBUG
		t1_start = chrono::steady_clock::now();
#endif

		if (step % updContLst == 0)
			GetContactsList(gr, b);

#ifdef DEBUG
		t1_end = chrono::steady_clock::now();
		diff_1 = t1_end - t1_start;
		t1 = t1 + chrono::duration <double, milli>(diff_1).count() / 1000; /* time in seconds*/
		t2_start = chrono::steady_clock::now();
#endif

#ifdef DEBUG
		t2_end = chrono::steady_clock::now();
		diff_2 = t2_end - t2_start;
		t2 = t2 + chrono::duration <double, milli>(diff_2).count() / 1000; /* time in seconds*/
		t3_start = chrono::steady_clock::now();
#endif

		/* Velocity Verlet's prediction Step */
		for (size_t i = 0; i < cells; i++) {
			
			// update positions
			gr[i].x = gr[i].x + dt * gr[i].Vx + 0.5*dt*dt*gr[i].Ax;
			gr[i].y = gr[i].y + dt * gr[i].Vy + 0.5*dt*dt*gr[i].Ay;
			gr[i].z = gr[i].z + dt * gr[i].Vz + 0.5*dt*dt*gr[i].Az;
			gr[i].th = gr[i].th + dt * gr[i].Vth + 0.5*dt*dt*gr[i].Ath;

			// update velocities
			gr[i].Vx = gr[i].Vx + 0.5 * dt * gr[i].Ax;
			gr[i].Vy = gr[i].Vy + 0.5 * dt * gr[i].Ay;
			gr[i].Vz = gr[i].Vz + 0.5 * dt * gr[i].Az;
			gr[i].Vth = gr[i].Vth + 0.5 * dt * gr[i].Ath;
		}

		wBottom.lowY = wBottom.lowY + dt * wBottom.wV + 0.5 * dt * dt * wBottom.wA;
		wTop.highY = wTop.highY + dt * wTop.wV + 0.5 * dt * dt * wTop.wA;
		wLeft.lowX = wLeft.lowX + dt * wLeft.wV + 0.5 * dt * dt * wLeft.wA;
		wRight.highX = wRight.highX + dt * wRight.wV + 0.5 * dt * dt * wRight.wA;

		wBottom.wV = wBottom.wV + 0.5 * dt * wBottom.wA;
		wTop.wV = /*wTop.wV + */0.5 * dt * wTop.wA; /* constant loading without acceleration*/
		wLeft.wV = wLeft.wV + 0.5 * dt * wLeft.wA;
		wRight.wV = wRight.wV + 0.5 * dt * wRight.wA;

		wFront.lowZ = wFront.lowZ + dt * wFront.wV + 0.5 * dt * dt * wFront.wA;
		wBack.highZ = wBack.highZ + dt * wBack.wV + 0.5 * dt * dt * wBack.wA;
		wFront.wV = wFront.wV + 0.5 * dt * wFront.wA;
		wBack.wV = wBack.wV + 0.5 * dt * wBack.wA;

		/* calculate mech ratio, only for 'sample' project */
		if (CheckProjType()) {

			double _sumOfNormF = 1.0;

			for (size_t i = 0; i < cells; i++) {
				gr[i].sumOfNormF = sqrt(gr[i].Fx*gr[i].Fx + gr[i].Fy*gr[i].Fy + gr[i].Fz*gr[i].Fz) + sqrt(gr[i].Fth*gr[i].Fth);
				_sumOfNormF += gr[i].sumOfNormF;
			}
			_mechRatio_a.push_back(_sumOfNormF);

			if (step % mRatTimeIntrv == 0) {
				double _mechRatio = accumulate(_mechRatio_a.begin(), _mechRatio_a.end(), 0.0) / _mechRatio_a.size();
				mechRatio = (abs(_mechRatio - _mechRatio_i) / _mechRatio);
				_mechRatio_i = _mechRatio;
				_mechRatio_a.clear();
				if (wrtMRatio)
					wrtMRatio << step << "  " << mechRatio << std::endl;
			}
		}

		/*Set sum of forces acting on cell to 0*/
		for (size_t i = 0; i < cells; i++) {
			gr[i].Fx = 0.;
			gr[i].Fy = 0.;
			gr[i].Fz = 0.;
			gr[i].Fth = 0.;
			gr[i].p = 0.;
		}

		/*Set confining pressure*/
		wBottom.wF = 0;
		wTop.wF = 0;
		wLeft.wF = 0;
		wRight.wF = 0;
		wFront.wF = 0;
		wBack.wF = 0;

#ifdef DEBUG
		t3_end = chrono::steady_clock::now();
		diff_3 = t3_end - t3_start;
		t3 = t3 + chrono::duration <double, milli>(diff_3).count() / 1000; /* time in seconds*/
		t4_start = chrono::steady_clock::now();
#endif

		/*Force computation between cells*/
#pragma omp parallel for 
		for (size_t i = 0; i < cells; i++) {
			for (auto j = gr[i].contactList.begin(); j != gr[i].contactList.end(); ++j) {
				ForcesGrainCohesive3D(i, *j, gr, b, dt);
			}
		}

#ifdef DEBUG
		t4_end = chrono::steady_clock::now();
		diff_4 = t4_end - t4_start;
		t4 = t4 + chrono::duration <double, milli>(diff_4).count() / 1000; /* time in seconds*/
		t5_start = chrono::steady_clock::now();
#endif

		/*Force computation between cells and walls*/
#pragma omp parallel for
		for (size_t i = 0; i < cells; i++)
			for (size_t w = 1; w <= wallsNum; w++)
				ForcesWall3D(i, w, gr);

#ifdef DEBUG
		t5_end = chrono::steady_clock::now();
		diff_5 = t5_end - t5_start;
		t5 = t5 + chrono::duration <double, milli>(diff_5).count() / 1000; /* time in seconds*/
		t6_start = chrono::steady_clock::now();
#endif

		/*Calculate corrected acceleration*/
		for (size_t i = 0; i < cells; i++) {
			gr[i].Ax = gr[i].Fx / gr[i].m;
			gr[i].Ay = gr[i].Fy / gr[i].m - Grav + Grav * (rho_env/rho);
			gr[i].Az = gr[i].Fz / gr[i].m;
			gr[i].Ath = gr[i].Fth / gr[i].I;
		}

		wBottom.wA = wBottom.wF / wBottom.wM;
		wTop.wA = wTop.wF / wTop.wM;
		wLeft.wA = wLeft.wF / wLeft.wM;
		wRight.wA = wRight.wF / wRight.wM;


		wFront.wA = wFront.wF / wFront.wM;
		wBack.wA = wBack.wF / wBack.wM;

		/* Correction step for velocities */
		for (size_t i = 0; i < cells; i++) {
			gr[i].Vx = gr[i].Vx + 0.5 * dt * gr[i].Ax;
			gr[i].Vy = gr[i].Vy + 0.5 * dt * gr[i].Ay;
			gr[i].Vz = gr[i].Vz + 0.5 * dt * gr[i].Az;
			gr[i].Vth = gr[i].Vth + 0.5 * dt * gr[i].Ath;
		}

		wBottom.wV = wBottom.wV + 0.5 * dt*wBottom.wA;
		wTop.wV = /*wTop.wV +*/ 0.5 * dt*wTop.wA;/*0.5 * wTop.wA;*/ //constant loading with no acceleration
		wLeft.wV = wLeft.wV + 0.5 * dt*wLeft.wA;
		wRight.wV = wRight.wV + 0.5 * dt*wRight.wA;

		wFront.wV = wFront.wV + 0.5 * dt * wFront.wA;
		wBack.wV = wBack.wV + 0.5 * dt * wBack.wA;

		allWalls._wBack = wBack;
		allWalls._wBottom = wBottom;
		allWalls._wFront = wFront;
		allWalls._wLeft = wLeft;
		allWalls._wRight = wRight;
		allWalls._wTop = wTop;

#ifdef DEBUG
		t6_end = chrono::steady_clock::now();
		diff_6 = t6_end - t6_start;
		t6 = t6 + chrono::duration <double, milli>(diff_6).count() / 1000; /* time in seconds*/
		t7_start = chrono::steady_clock::now();
#endif

		/* correct cells' origins */
		for (size_t i = 0; i < cells; i++)
			gr[i].correct_origin();

		/* record the results of simulation */
		if (step%outputPeriod == 0 && record)
			GetSnapshots(gr, b);

		if (step % outputPeriod == 0) {
			if (CheckProjType())
				RuntimeLog(runTimeFile, true, "Sample equilibrium error  ", mechRatio, false, false,true);
			else
				RuntimeLog(runTimeFile, true, "Completed, %", (int)((step * dt * 100) / timeInterval)+1, false, false, true);
		}

		/* writes the last snapshot of the sample if the project has no record option. */
		if (!record && step == (stepNum - 1))
			GetSnapshots(gr, b);

#ifdef DEBUG
		t7_end = chrono::steady_clock::now();
		diff_7 = t7_end - t7_start;
		t7 = t7 + chrono::duration <double, milli>(diff_7).count() / 1000; /* time in seconds*/
#endif

		step++;
	}

#ifdef DEBUG
	RuntimeLog(runTimeFile, false, "t1:", t1, false, false, true);
	RuntimeLog(runTimeFile, false, "t2:", t2, false, false, true);
	RuntimeLog(runTimeFile, false, "t3:", t3, false, false, true);
	RuntimeLog(runTimeFile, false, "t4:", t4, false, false, true);
	RuntimeLog(runTimeFile, false, "t5:", t5, false, false, true);
	RuntimeLog(runTimeFile, false, "t6:", t6, false, false, true);
	RuntimeLog(runTimeFile, false, "t7:", t7, false, false, true);
#endif

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

std::vector<shortBond> Sample::GetSnapshotB(std::vector<std::vector<Bond>> &b1)
{
	std::vector<shortBond> allBnds;
	for (size_t k = 0; k < grainNum; k++) {
		for (size_t l = k + 1; l < grainNum; l++) {
			if (b1[k][l].isActive) {
				shortBond sB;
				sB.i = k;
				sB.j = l;
				sB.f = b1[k][l].contactForce;
				allBnds.push_back(sB);
			}
		}
	}

	return allBnds;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::GetSnapshots(std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b1)
{
	/*omp_set_num_threads(1);*/
#pragma omp parallel sections
	{
#pragma omp section
		{
			wlSmf = false;
			wl.push_back(allWalls);
			long int psz = wlCount - 1;
			if (psz >= snapshotBuffer)
			{
				if (wlCount <= wl.size())
				{
					wl.erase(wl.begin(), wl.begin() + psz);
					wlTop = true;
				}
			}
			wlSmf = true;
		}
#pragma omp section
		{
			grnSmf = false;
			grainSnapshot.push_back(gr);
			long int sz = grnCount - 1;
			if (sz >= snapshotBuffer)
			{
				if (grnCount <= grainSnapshot.size())
				{
					grainSnapshot.erase(grainSnapshot.begin(), grainSnapshot.begin() + sz);
					grnTop = true;
				}
			}
			grnSmf = true;

			fSmf = false;
			grainSnapshotF.push_back(gr);
			sz = fCount - 1;
			if (sz >= snapshotBuffer)
			{
				if (fCount <= grainSnapshotF.size())
				{
					grainSnapshotF.erase(grainSnapshotF.begin(), grainSnapshotF.begin() + sz);
					fTop = true;
				}
			}
			fSmf = true;
		}
#pragma omp section
		{
			bSmf = false;
			grainSnapshotB.push_back(gr);
			allBondsSnapshot.push_back(GetSnapshotB(b1));
			long int sz = bCount - 1;
			if (sz >= snapshotBuffer)
			{
				if (bCount <= grainSnapshotB.size())
				{
					grainSnapshotB.erase(grainSnapshotB.begin(), grainSnapshotB.begin() + sz);
					allBondsSnapshot.erase(allBondsSnapshot.begin(), allBondsSnapshot.begin() + sz);
					bTop = true;
				}
			}
			bSmf = true;			
		}
	}

	snapshotCounter++;

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::GetCellsSnapshots(std::vector <Cell> &gr)
{
	grnSmf = false;
	grainSnapshot.push_back(gr);
	long int sz = grnCount - 1;
	if (sz >= snapshotBuffer)
	{
		if (grnCount <= grainSnapshot.size())
		{
			grainSnapshot.erase(grainSnapshot.begin(), grainSnapshot.begin() + sz);
			grnTop = true;
		}
	}
	grnSmf = true;

	fSmf = false;
	grainSnapshotF.push_back(gr);
	sz = fCount - 1;
	if (sz >= snapshotBuffer)
	{
		if (fCount <= grainSnapshotF.size())
		{
			grainSnapshotF.erase(grainSnapshotF.begin(), grainSnapshotF.begin() + sz);
			fTop = true;
		}
	}
	fSmf = true;

	snapshotCounter++;

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::PostProcGrain(int numfile)
{

	if (grnTop) {
		grnCount = 0;
		grnTop = false;
	}

	size_t g_sz = grainSnapshot.size() - 1;

	if (grnSmf && grnCount <= g_sz) {
		std::vector <Cell> grn;
		try
		{
			grn = grainSnapshot.at(grnCount);
		}
		catch (const std::out_of_range& e) {
			RuntimeLog(runTimeFile, true, "Out of range error -> Cell! Call exit(EXIT_FAILURE) from 'PostProcGrain' !):", 0.0, false, true);
		}
		grnCount++;

		size_t grainNum = grn.size();
		size_t i;
		char fnameG[256]; // file name

		#ifdef _WIN64
		if (CheckProjType())
			snprintf(fnameG, 200, "sample\\vtk_files\\cells%04d.vtk", numfile);
		else
			snprintf(fnameG, 200, "simulation\\vtk_files\\cells%04d.vtk", numfile);
		#endif

		#ifdef linux
		if (CheckProjType())
			snprintf(fnameG, 200, "sample/vtk_files/cells%04d.vtk", numfile);
		else
			snprintf(fnameG, 200, "simulation/vtk_files/cells%04d.vtk", numfile);
		#endif

		std::ofstream wrtToFileG(fnameG, std::ios::out);

		if (!wrtToFileG.is_open()) {
			RuntimeLog(runTimeFile, true, "Error when opening file for GRAIN solution! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
			exit(1);
		}

		if (wrtToFileG)
		{
			wrtToFileG.precision(5); wrtToFileG << std::scientific;
			wrtToFileG << "# vtk DataFile Version 3.0" << std::endl;
			wrtToFileG << "My cells" << std::endl;
			wrtToFileG << "ASCII" << std::endl;
			wrtToFileG << "DATASET UNSTRUCTURED_GRID" << std::endl;
			wrtToFileG << "POINTS " << grainNum << " float" << std::endl;
			for (i = 0; i < grainNum; i++)
				wrtToFileG << grn[i].x << " " << grn[i].y << " " << grn[i].z << std::endl;
			wrtToFileG << "POINT_DATA " << grainNum << std::endl;
			wrtToFileG << "VECTORS Radius float" << std::endl;
			for (i = 0; i < grainNum; i++)
				wrtToFileG << grn[i].R << " 0.0 0.0" << std::endl;

			wrtToFileG << "SCALARS Pressure float" << std::endl;
			wrtToFileG << "LOOKUP_TABLE default" << std::endl;
			for (i = 0; i < grainNum; i++)
				wrtToFileG << grn[i].p << std::endl;

			wrtToFileG << "SCALARS Rotation float" << std::endl;
			wrtToFileG << "LOOKUP_TABLE default" << std::endl;
			for (i = 0; i < grainNum; i++)
				wrtToFileG << grn[i].th << std::endl;
		}

		wrtToFileG.close();

		snapshotCounterG++;

		if (!relaxationFlag) {
			int sz = grnCount - 1;
			if (sz >= snapshotBuffer) {
				grainSnapshot.erase(grainSnapshot.begin(), grainSnapshot.begin() + sz);
				grnTop = true;
			}
		}
	}
	else if (!relaxationFlag)
		grainSnapshot.clear();

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::PostProcGrainVisual(std::vector <std::vector <float> > &data)
{

	std::vector <std::vector <float> > t_data;

	if (grnTop) {
		grnCount = 0;
		grnTop = false;
	}

	size_t g_sz = grainSnapshot.size() /* - (size_t) 1 */;

	if ( grnSmf && (grnCount </* = */ g_sz) ) {
		std::vector <Cell> grn;
		try
		{
			grn = grainSnapshot.at(grnCount);
		}
		catch (const std::out_of_range& e) {
			RuntimeLog(runTimeFile, true, "Out of range error -> Cell! from 'PostProcGrainVisual' !:", 0.0, false, false);

			RuntimeLog(runTimeFile, false, "grnSmf = ", grnSmf, false, false, true);
			RuntimeLog(runTimeFile, false, "grnCount = ", grnCount, false, false, true);
			RuntimeLog(runTimeFile, false, "g_sz = ", g_sz, false, false, true);
			RuntimeLog(runTimeFile, false, "grainSnapshot.size() = ", grainSnapshot.size(), false, false, true);
			RuntimeLog(runTimeFile, true, "Call exit(EXIT_FAILURE) from 'PostProcGrainVisual' !):", 0.0, false, true);
			exit(EXIT_FAILURE);
		}
		grnCount++;

		size_t grainNum = grn.size();
		size_t i;

		std::vector <float> tvect(6, 0.0f);

		// get sacaler for pressure
		double pscaler = 0.0f;
		double avoidzerro = 0.1f;
		for (i = 0; i < grainNum; i++){
			double pmax = grn[i].p;
			if ( pmax > pscaler )
				pscaler = pmax;
		}
		pscaler = pscaler + avoidzerro;
		
		double scaler = std::min(_sizeX, _sizeY) * 0.5f;//grn[0].R;

		for (i = 0; i < grainNum; i++){
			tvect[0] = static_cast <float> ( grn[i].x / scaler );
			tvect[1] = static_cast <float> ( grn[i].y / scaler );
			tvect[2] = static_cast <float> ( grn[i].z / scaler );
			tvect[3] = static_cast <float> ( grn[i].R / scaler * 2.0f );
			tvect[4] = static_cast <float> ( (grn[i].p + avoidzerro) / pscaler );
			tvect[5] = static_cast <float> (grn[i].th);

			t_data.push_back(tvect);
		}

		std::swap(t_data, data);

		snapshotCounterG++;

		if (!relaxationFlag) {
			int sz = grnCount - 1;
			if (sz >= snapshotBuffer) {
				grainSnapshot.erase(grainSnapshot.begin(), grainSnapshot.begin() + sz);
				grnTop = true;
			}
		}
	}
	else if (!relaxationFlag)
		grainSnapshot.clear();

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::PostProcForce(int numfile)
{
	if (fTop) {
		fCount = 0;
		fTop = false;
	}
	size_t f_sz = grainSnapshotF.size() - 1;
	if (fSmf && fCount <= f_sz)
	{
		std::vector <Cell> grn = grainSnapshotF.at(fCount);
		fCount++;

		size_t grainNum = grn.size();

		size_t maxForceNum = grainNum * 12;
 		// double *forceFn; forceFn = new double[maxForceNum];
		// size_t *force_i; force_i = new size_t[maxForceNum];
		// size_t *force_j; force_j = new size_t[maxForceNum];
		std::vector<double> forceFn;
		std::vector<size_t> force_i;
		std::vector<size_t> force_j;
		size_t forceNum = 0; // number of forces
		char fnameF[256]; // file name

		std::string smpl("sample");

		#ifdef _WIN64
		if (!smpl.compare(type)) snprintf(fnameF, 200, "sample\\vtk_files\\forces%04d.vtk", numfile);
		else snprintf(fnameF, 200, "simulation\\vtk_files\\forces%04d.vtk", numfile);
		#endif

		#ifdef linux
		if (!smpl.compare(type)) snprintf(fnameF, 200, "sample/vtk_files/forces%04d.vtk", numfile);
		else snprintf(fnameF, 200, "simulation/vtk_files/forces%04d.vtk", numfile);
		#endif

		std::ofstream wrtToFileF(fnameF, std::ios::out);

		if (!wrtToFileF.is_open())
		{
			RuntimeLog(runTimeFile, true, "Error when opening file for FORCE solution! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
			exit(1);
		}

		if (wrtToFileF)
		{
			// Count and store forces in vectors
			for (size_t i = 0; i < grainNum; i++)
			{
				for (size_t j = i + 1; j < grainNum; j++)
				{
					double rXji, rYji, rZji, dn;
					double rNorm = 0.;
					// Compute distance between cells
					rXji = grn[i].x - grn[j].x;
					rYji = grn[i].y - grn[j].y;
					rZji = grn[i].z - grn[j].z;
					rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
					dn = rNorm - (grn[i].R + grn[j].R);
					if (dn < 0.) // Contact
					{
						// forceFn[forceNum] = -kn * dn;
						// force_i[forceNum] = i;
						// force_j[forceNum] = j;
						forceFn.push_back(-kn * dn);
						force_i.push_back(i);
						force_j.push_back(j);

						forceNum++;
					}
				}
			}
			wrtToFileF.precision(5); wrtToFileF << std::scientific;
			wrtToFileF << "# vtk DataFile Version 3.0" << std::endl;
			wrtToFileF << "My forces" << std::endl;
			wrtToFileF << "ASCII" << std::endl;
			wrtToFileF << "DATASET POLYDATA" << std::endl;
			wrtToFileF << "POINTS " << 2 * forceNum << " float" << std::endl;
			for (size_t k = 0; k < forceNum; k++)
			{
				wrtToFileF << grn[force_i[k]].x << " " << grn[force_i[k]].y << " " << grn[force_i[k]].z << std::endl;
				wrtToFileF << grn[force_j[k]].x << " " << grn[force_j[k]].y << " " << grn[force_j[k]].z << std::endl;
			}
			wrtToFileF << "LINES " << forceNum << " " << 3 * forceNum << std::endl;
			for (size_t k = 0; k < forceNum; k++)
			{
				wrtToFileF << "2 " << 2 * k << " " << 2 * k + 1 << std::endl;
			}
			wrtToFileF << "POINT_DATA " << 2 * forceNum << std::endl;
			wrtToFileF << "SCALARS Force float" << std::endl;
			wrtToFileF << "LOOKUP_TABLE default" << std::endl;
			for (size_t k = 0; k < forceNum; k++)
			{
				wrtToFileF << forceFn[k] << std::endl;
				wrtToFileF << forceFn[k] << std::endl;
			}
		}

		wrtToFileF.close();

		snapshotCounterF++;

		// delete[] force_i;
		// delete[] force_j;
		// delete[] forceFn;
		force_i.clear();
		force_j.clear();
		forceFn.clear();
		forceFn.shrink_to_fit();
		force_i.shrink_to_fit();
		force_j.shrink_to_fit();


		if (!relaxationFlag)
		{
			int sz = fCount - 1;
			if (sz >= snapshotBuffer)
			{
				grainSnapshotF.erase(grainSnapshotF.begin(), grainSnapshotF.begin() + sz);
				fTop = true;
			}
		}
	}
	else if (!relaxationFlag) grainSnapshotF.clear();

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::PostProcBond(int numfile)
{
	if (bTop) {
		bCount = 0;
		bTop = false;
	}

	size_t b_sz = grainSnapshotB.size() - 1;

	if (bSmf && bCount <= b_sz) {

		std::vector <Cell> grn;
		std::vector<shortBond> _bond;

		try {
			grn = grainSnapshotB.at(bCount);
			_bond = allBondsSnapshot.at(bCount);
		}
		catch (const std::out_of_range& e) {
			RuntimeLog(runTimeFile, true, "Out of range error -> Bond! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
			exit(1);
		}
		
		bCount++;

		auto bondsNum = _bond.size();
		//size_t _grainNum = grn.size();
		std::vector <double> forceFn;
		std::vector <size_t> force_i;
		std::vector <size_t> force_j;
		size_t forceNum = 0; // number of forces
		char fnameB[256]; // file name

		std::string smpl("sample");

		#ifdef _WIN64
		if (!smpl.compare(type)) snprintf(fnameB, 200, "sample\\vtk_files\\bonds%04d.vtk", numfile);
		else snprintf(fnameB, 200, "simulation\\vtk_files\\bonds%04d.vtk", numfile);
		#endif

		#ifdef linux
		if (!smpl.compare(type)) snprintf(fnameB, 200, "sample/vtk_files/bonds%04d.vtk", numfile);
		else snprintf(fnameB, 200, "simulation/vtk_files/bonds%04d.vtk", numfile);
		#endif

		std::ofstream wrtToFileB(fnameB, std::ios::out);

		if (!wrtToFileB.is_open()) {
			RuntimeLog(runTimeFile, true, "Error when opening file for solution -> Bond! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
			exit(1);
		}
		
		if (wrtToFileB) {

			for (size_t k = 0; k < bondsNum; k++) {

				try {

					forceFn.push_back(_bond[k].f);
					force_i.push_back(_bond[k].i);
					force_j.push_back(_bond[k].j);
					forceNum++;
				}
				catch (const std::bad_alloc& e) {
					RuntimeLog(runTimeFile, true, e.what(), 0.0, false, true);
				}
			}

			wrtToFileB.precision(5); wrtToFileB << std::scientific;
			wrtToFileB << "# vtk DataFile Version 3.0" << std::endl;
			wrtToFileB << "My forces" << std::endl;
			wrtToFileB << "ASCII" << std::endl;
			wrtToFileB << "DATASET POLYDATA" << std::endl;
			wrtToFileB << "POINTS " << 2 * forceNum << " float" << std::endl;
			
			for (size_t k = 0; k < forceNum; k++) {
				wrtToFileB << grn[force_i[k]].x << " " << grn[force_i[k]].y << " " << grn[force_i[k]].z << std::endl;
				wrtToFileB << grn[force_j[k]].x << " " << grn[force_j[k]].y << " " << grn[force_j[k]].z << std::endl;
			}

			wrtToFileB << "LINES " << forceNum << " " << 3 * forceNum << std::endl;
			
			for (size_t k = 0; k < forceNum; k++)
				wrtToFileB << "2 " << 2 * k << " " << 2 * k + 1 << std::endl;

			wrtToFileB << "POINT_DATA " << 2 * forceNum << std::endl;
			wrtToFileB << "SCALARS Force float" << std::endl;
			wrtToFileB << "LOOKUP_TABLE default" << std::endl;
			
			for (size_t k = 0; k < forceNum; k++) {

				wrtToFileB << forceFn[k] << std::endl;
				wrtToFileB << forceFn[k] << std::endl;
			}

		}

		wrtToFileB.close();

		snapshotCounterB++;

		forceFn.clear();
		force_i.clear();
		force_j.clear();
		forceFn.shrink_to_fit();
		force_i.shrink_to_fit();
		force_j.shrink_to_fit();
		_bond.clear();
		_bond.shrink_to_fit();

		if (!relaxationFlag) {
			int sz = bCount - 1;
			if (sz >= snapshotBuffer) {
				grainSnapshotB.erase(grainSnapshotB.begin(), grainSnapshotB.begin() + sz);
				allBondsSnapshot.erase(allBondsSnapshot.begin(), allBondsSnapshot.begin() + sz);
				bTop = true;
			}
		}
	}
	else if (!relaxationFlag) {
		grainSnapshotB.clear();
		allBondsSnapshot.clear();
	}

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::PostProcWall(int numfile, std::vector<WALLS> &wl)
{
	if (wlTop) {
		wlCount = 0;
		wlTop = false;
	}
	size_t w_sz = wl.size() - 1;
	if (wlSmf && wlCount <= w_sz)
	{
		//WALLS _allWalls = wl.front();
		WALLS _allWalls = wl.at(wlCount);
		wlCount++;

/* 		Wall wBottom = _allWalls._wBottom;
		Wall wTop = _allWalls._wTop;
		Wall wLeft = _allWalls._wLeft;
		Wall wRight = _allWalls._wRight;
 */		Wall wFront = _allWalls._wFront;
		Wall wBack = _allWalls._wBack;

		char fnameW[256]; // file name

		std::string smpl("sample");

		#ifdef _WIN64
		if (!smpl.compare(type)) snprintf(fnameW, 200, "sample\\vtk_files\\walls%04d.vtk", numfile);
		else snprintf(fnameW, 200, "simulation\\vtk_files\\walls%04d.vtk", numfile);
		#endif

		#ifdef linux
		if (!smpl.compare(type)) snprintf(fnameW, 200, "sample/vtk_files/walls%04d.vtk", numfile);
		else snprintf(fnameW, 200, "simulation/vtk_files/walls%04d.vtk", numfile);
		#endif

		std::ofstream wrtToFileW(fnameW, std::ios::out);

		if (!wrtToFileW.is_open())
		{
			//std::cout << "Error when opening file for solution... Exit." << std::endl;
			RuntimeLog(runTimeFile, true, "Error when opening file for solution (Wall)... Exit.", 0.0, false, true);
			exit(1);
		}

		if (wrtToFileW)
		{
			wrtToFileW.precision(3); wrtToFileW << std::scientific;
			wrtToFileW << "# vtk DataFile Version 3.0" << std::endl;
			wrtToFileW << "My walls" << std::endl;
			wrtToFileW << "ASCII" << std::endl;
			wrtToFileW << "DATASET POLYDATA" << std::endl;
			wrtToFileW << "POINTS 8 float" << std::endl;

			// front wall
			wrtToFileW << wFront.lowX << " " << wFront.lowY << " " << wFront.lowZ << std::endl;
			wrtToFileW << wFront.highX << " " << wFront.lowY << " " << wFront.lowZ << std::endl;
			wrtToFileW << wFront.highX << " " << wFront.highY << " " << wFront.lowZ << std::endl;
			wrtToFileW << wFront.lowX << " " << wFront.highY << " " << wFront.lowZ << std::endl;
			// back wall
			wrtToFileW << wBack.lowX << " " << wBack.lowY << " " << wBack.highZ << std::endl;
			wrtToFileW << wBack.highX << " " << wBack.lowY << " " << wBack.highZ << std::endl;
			wrtToFileW << wBack.highX << " " << wBack.highY << " " << wBack.highZ << std::endl;
			wrtToFileW << wBack.lowX << " " << wBack.highY << " " << wBack.highZ << std::endl;
			
			wrtToFileW << "POLYGONS 6 30" << std::endl;

			wrtToFileW << "4 0 1 2 3" << std::endl;
			wrtToFileW << "4 4 5 6 7" << std::endl;
			wrtToFileW << "4 0 1 5 4" << std::endl;
			wrtToFileW << "4 2 3 7 6" << std::endl;
			wrtToFileW << "4 0 4 7 3" << std::endl;
			wrtToFileW << "4 1 2 6 5" << std::endl;
		}

		wrtToFileW.close();

		snapshotCounterW++;

		if (!relaxationFlag)
		{
			int psz = wlCount - 1;
			if (psz >= snapshotBuffer)
			{
				wl.erase(wl.begin(), wl.begin() + psz);
				wlTop = true;
			}
		}

	}
	else if (!relaxationFlag) wl.clear();

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::PostProcessing()
{
	/*omp_set_num_threads(1);*/
#pragma omp parallel sections
	{
		
#pragma omp section
		{
			if (!grainSnapshot.empty())
				PostProcGrain(snapshotCounterG);
		}
		
#pragma omp section
		{
			if (!grainSnapshotF.empty())
				PostProcForce(snapshotCounterF);
		}
		
#pragma omp section
		{
			if (!allBondsSnapshot.empty() )
				PostProcBond(snapshotCounterB);
		}
		
#pragma omp section
		{
			if (!wl.empty())
				PostProcWall(snapshotCounterW, wl);
		}
	}
	
	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::RunPostProcessing()
{//std::<<"relaxationFlag = "<<relaxationFlag<<std::endl;
	do {
		PostProcessing();
	} while (!grainSnapshot.empty() || relaxationFlag);

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::RunVisualPostProcessing()
{
	std::vector <std::vector <float> > grainData;
	bool accessLock = true;
	bool escSignal = false;
	bool renderSignal = true;


omp_set_num_threads(2);
#pragma omp parallel sections //shared(accessLock, renderSignal, escSignal)
	{
#pragma omp section
		{
			do {
					if ( accessLock && renderSignal ){
						accessLock = false;
						PostProcGrainVisual(grainData);
						accessLock = true;
						renderSignal = false;
					}
			} while (!grainSnapshot.empty() || relaxationFlag);			
			escSignal = true;
		}
#pragma omp section
		{
			oglu::visualizer scene;
			scene.show(escSignal, accessLock, renderSignal, grainData);
		}
	}

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

int Sample::GetContactsList(std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b)
{
//#pragma omp parallel for
	for (size_t i = 0; i < gr.size(); i++) {
		gr[i].contactList.clear();
		for (size_t j = i + 1; j < gr.size(); j++) {
			double rXji, rYji, rZji, dn;
			double rNorm = 0.;
			// Compute distance between cells
			rXji = gr[i].x - gr[j].x;
			rYji = gr[i].y - gr[j].y;
			rZji = gr[i].z - gr[j].z;
			rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
			dn = rNorm - (gr[i].R + gr[j].R);

			if (b[i][j].isActive == true) {
				gr[i].contactList.push_back(j);
				continue;
			}
			else if (dn < 0.)
				gr[i].contactList.push_back(j);

		}
	}
	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Sample::RuntimeLog(std::string fileName, bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose, bool showValue) {

	//char ctime_no_newline[26];
	//char* _ctime_no_newline;
	//char delim[] = "\n";
	//char *next_token;

	/* Initial message to the log file */
	if (makeFileOpen) {
		runtimelogStr.open(fileName, std::fstream::out);

		auto timeNow = std::chrono::system_clock::now();
		auto outTime = std::chrono::system_clock::to_time_t(timeNow);
		//ctime_s(ctime_no_newline, sizeof ctime_no_newline, &outTime);
		//_ctime_no_newline = strtok_s(ctime_no_newline, delim, &next_token);

		if (runtimelogStr) {
			//runtimelogStr << "DEM PROGRAM RUNTIME LOG FILE   " << _ctime_no_newline << std::endl;
			runtimelogStr << "DEM PROGRAM RUNTIME LOG FILE   " << ctime(&outTime) << std::endl;
			runtimelogStr << std::endl;
		}
	}

	auto timeNow = std::chrono::system_clock::now();
	auto outTime = std::chrono::system_clock::to_time_t(timeNow);
	//ctime_s(ctime_no_newline, sizeof ctime_no_newline, &outTime);
	//_ctime_no_newline = strtok_s(ctime_no_newline, delim, &next_token);

	/* Conditional messages to the log file */

	if (wrTime && !showValue) {/* Shows only the time and the message */
		if (runtimelogStr) {
			runtimelogStr << std::endl;
			//runtimelogStr << _ctime_no_newline << "   " << message << std::endl;
			runtimelogStr << ctime(&outTime) << "   " << message << std::endl;
			logFilePos_end = runtimelogStr.tellg();
		}
	}

	if (!wrTime && showValue) {/* Shows only the value and the message */
		if (runtimelogStr) {
			runtimelogStr << "                           " << message << ".........." << someValue << std::endl;
			logFilePos_end = runtimelogStr.tellg();
		}
	}
	else if (wrTime && showValue) {/* Shows the value, time and the message */
		if (runtimelogStr) {
			runtimelogStr.seekg(logFilePos_end, runtimelogStr.beg);
			//runtimelogStr << _ctime_no_newline << "             " << message << ".........." << someValue << std::endl;
			runtimelogStr << ctime(&outTime) << "             " << message << ".........." << someValue << std::endl;
		}
	}

	if (makeFileClose)
		runtimelogStr.close();

}

//-------------------------------------------------------------------------------------------------------------------------

void Sample::SeedCells(double rMin, double rMax, double sizeD, double sizeH, double cellsNum, std::vector <Cell> &gr, double stdDev) {

	for ( size_t i = 0; i < cellsNum; i++) {

		double R, _R;
		double x, y, z;

		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::default_random_engine generator1(seed + seed);
		std::default_random_engine generator2(seed + seed + seed);
		
		bool noGauss = false;
		if (stdDev == 0.0) {
			//stdDev = 0.000001;
			noGauss = true;
		}

		_R = (rMax + rMin) / 2;

		std::normal_distribution<double> distributionR( _R, stdDev );
		std::uniform_real_distribution<double> distributionX( 0.0f, sizeD );
		std::uniform_real_distribution<double> distributionY( 0.0f, sizeH );
		std::uniform_real_distribution<double> distributionZ( 0.0f, sizeD );
		
		if (!noGauss) {
			do
				R = distributionR(generator);
			while ( !(R >= rMin - stdDev && R <= rMax + stdDev) );
		}
		else {
			R = _R;
		}

		x = distributionX(generator);
		y = distributionY(generator1);
		z = distributionZ(generator2);

		glm::dvec3 origin(x,y,z);

/* 		std::cout<<"sizeD = "<<sizeD<<"; sizeH = "<<sizeH<<std::endl; 
		std::cout<<"sizeD-_R = "<<sizeD-_R<<"; sizeH-_R = "<<sizeH-_R<<std::endl;
		std::cout<<"x = "<<origin.x<<"; y = "<<origin.y<<"; z = "<<origin.z<<std::endl;
*/
		gr.push_back( Cell(origin,R) );

		gr[i].m = (4. / 3.) * M_PI*rho*R*R*R;
		gr[i].I = (2. / 5.)*gr[i].m*R*R;
		gr[i].Ath = 0.; gr[i].Ax = 0.; gr[i].Ay = 0.; gr[i].Az = 0.;
		gr[i].Vth = 0.; gr[i].Vx = 0.; gr[i].Vy = 0.; gr[i].Vz = 0.; gr[i].th = 0.;
		gr[i].isActive = true;
	}

}

//---------------------------------------------------------------------------------------------------------------------------------
