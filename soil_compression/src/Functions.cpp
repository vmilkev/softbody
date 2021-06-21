#include "Functions.h"
#include <time.h> 
#include <locale>
#include <utility>
#include <codecvt>

#include "../../graphics/src/visualizer.hpp"

void InitiateParameters(double rMin, double rMax)
{
	b_YM = 0.5;
	mu = 0.3; // Sliding friction
	mur = 0.01; // Rolling friction
	nur = kt*rMin; // Rolling "damping"
	Adh = 0.;// Adhesion
	grainNum = 0;
	beta = 0.75; /* damping coefficient */
	_nodesRes = 1.0;
	sampleStretch = 1.8;
}

int ReadParametersFromFile(const char* paramFileName)
{
	extern int dim;
	extern bool isCohesion;
	extern double _rMin;
	extern double _rMax;
	extern double _sizeX;
	extern double _sizeY;
	extern double _timeInterval;	
	extern double porosity;
	extern double lRate;
	extern double pSize;
	extern double runTime, samplTime;
	extern double Grav;
	extern double stdDev;
	extern double shape;
	extern double scale;
	extern uint uint_recordOutput;
	extern bool recordOutput;
	extern std::string type;
	extern double stopMechRatio;
	extern bool deloneRapture;

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

	dim = stoi(lines.at(1));
	std::string tru("true");
	std::string fls("false");
	std::string act = lines.at(3);

	if ( tru.compare(act) == 0 )
		isCohesion = true;
	else if ( fls.compare(act) == 0 )
		isCohesion = false;
	else{
		RuntimeLog(runTimeFile, true, "Fail in reading isCohesion data from file! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
	}

	_rMin = stod(lines.at(5));
	_rMax = stod(lines.at(7));
	_sizeX = stod(lines.at(9));
	_sizeY = stod(lines.at(11));
	samplTime = stod(lines.at(13));
	porosity = stod(lines.at(15));
	lRate = stod(lines.at(17));
	pSize = stod(lines.at(19));
	runTime = stod(lines.at(21));
	Grav = stod(lines.at(23));
	kn = stod(lines.at(25));
	kt = stod(lines.at(27));
	b_kn = stod(lines.at(29));
	b_kt = stod(lines.at(31));
	sigma = stod(lines.at(33));
	tau = stod(lines.at(35));
	rho = stod(lines.at(37));
	stdDev = stod(lines.at(39));
	shape = stod(lines.at(41));
	scale = stod(lines.at(43));
	recordOutput = stoi(lines.at(45));
	uint_recordOutput = stoi(lines.at(45));
	stopMechRatio = stod(lines.at(47));
	deloneRapture = stoi(lines.at(49));

	return 0;
}

#ifdef _WIN64

	int ProjectManager()
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

	bool DirExists(const std::string& dirName_in)
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

#ifdef linux

	int ProjectManager(){

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

	bool DirExists(const std::string& dirName_in)
	{
		struct stat sb;
		
		if (stat(dirName_in.c_str(), &sb) == -1) {
			RuntimeLog(runTimeFile, true, "STAT error while checking directory!", 0.0, false, false);
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


int FinalizeSampleProject(std::string paramFileName, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, struct Piston pist, std::vector<WALLS> &wl)
{
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
			a->grainId = e.grainId;
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
		RuntimeLog(runTimeFile, false, "Number of recorded grains at *.SMPL file", numOfGrains, false, false, true);
	}
	catch (int ex) {
		RuntimeLog(runTimeFile, true, "Problem with allocating memory for *.SMPL file !EXIT with error !", 0.0, false, true);
		return EXIT_FAILURE;
	}
	catch (std::exception const& e) {
		RuntimeLog(runTimeFile, true, e.what(), 0.0, false, true);
		return EXIT_FAILURE;
	}

	out.close();

	return 0;
}

int InitializeSimulationProject(const char* paramFileName, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, struct Piston pist, std::vector<WALLS> &wl)
{
	std::ifstream in(paramFileName, std::ios::in | std::ios::binary);
	if (!in) {
		RuntimeLog(runTimeFile, true, "Problem with *.SMPL file !EXIT with error !", 0.0, false, true);
		return EXIT_FAILURE;
	}

	in.seekg(0, in.end);
	size_t length = in.tellg() / sizeof(SmallGrain);
	in.seekg(0, in.beg);

	RuntimeLog(runTimeFile, false, "Number of recorded grains at *.SMPL file", length, false, false, true);
	
	for (size_t i = 0; i < length; i++)
		gr.push_back(Grain());	

	try {
		SmallGrain *a; a = (SmallGrain *)malloc(sizeof(SmallGrain));
		if (a == NULL) throw 40;
		for (size_t i = 0; i < length; i++) {
			in.read(reinterpret_cast<char*> (a), sizeof(SmallGrain));
			gr[i].grainId = a->grainId;
			gr[i].I = a->I;
			gr[i].m = a->m;
			gr[i].R = a->R;
			gr[i].rhoInd = a->rhoInd;
			gr[i].x = a->x;
			gr[i].y = a->y;
			gr[i].z = a->z;
		}
		free(a);
	}
	catch (int ex) {
		RuntimeLog(runTimeFile, true, "Problem with allocating memory for reading from *.SMPL file !EXIT with error !", 0.0, false, true);
		return EXIT_FAILURE;
	}
	catch (std::exception const& e) {
		RuntimeLog(runTimeFile, true, e.what(), 0.0, false, true);
		return EXIT_FAILURE;
	}

	in.close();

	return 0;
}

void SetSample2D(double rMin, double rMax, double sizeX, double sizeY, double nodesRes, std::vector <Grain> &gr, double stdDev) {

	int sampleGridNodesX = (int)ceil(nodesRes*sizeX / (rMax + rMin));
	int sampleGridNodesY = (int)ceil(nodesRes*sizeY / (rMax + rMin));

	double gridStep = std::min(sizeX / sampleGridNodesX, sizeY / sampleGridNodesY);

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	bool noGauss = false;
	if (stdDev == 0.0) {
		stdDev = 1;
		noGauss = true;
	}
	std::normal_distribution<double> distribution((rMax + rMin) / 2, stdDev);

	int iStep = (int)floor((0.5*(rMax + rMin)) / gridStep);  // if the grid not is not ocupied then define how much nodes can be ocupied by the grain gr[grainNum] with radius R

	for (int i = sampleGridNodesY - iStep - 1; i > iStep; --i) {
		for (int j = iStep + 1; j < sampleGridNodesX - iStep; j++) { // select the grid node
			gr.push_back(Grain());
			if (!noGauss) {
				do {
					gr[grainNum].R = distribution(generator);
				} while (!(gr[grainNum].R >= rMin - stdDev && gr[grainNum].R <= rMax + stdDev));
			}
			else {
				gr[grainNum].R = (double)rand() / (RAND_MAX + 1L)*(rMax - rMin) + rMin;
			}
			gr[grainNum].m = M_PI*rho*gr[grainNum].R*gr[grainNum].R;
			gr[grainNum].I = 0.5*gr[grainNum].m*gr[grainNum].R*gr[grainNum].R;
			gr[grainNum].Ath = 0.; gr[grainNum].Ax = 0.; gr[grainNum].Ay = 0.;
			gr[grainNum].Vth = 0.; gr[grainNum].Vx = 0.; gr[grainNum].Vy = 0.; gr[grainNum].th = 0.;
			gr[grainNum].x = j*gridStep;
			gr[grainNum].y = i*gridStep;
			gr[grainNum].grainId = grainNum;
			gr[grainNum].isActive = true;
			grainNum++;
		}
	}

}

void SetSample3D(double rMin, double rMax, double sizeD, double sizeH, double nodesRes, std::vector <Grain> &gr, double stdDev) {

	double sizeZ = sizeD;

	int sampleGridNodesX = (int)ceil(nodesRes*sizeD / (rMax + rMin));
	int sampleGridNodesY = (int)ceil(nodesRes*sizeH / (rMax + rMin));
	int sampleGridNodesZ = (int)ceil(nodesRes*sizeD / (rMax + rMin));

	double gridStep = std::min(sizeD / sampleGridNodesX, sizeH / sampleGridNodesY);

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	bool noGauss = false;
	if (stdDev == 0.0) {
		stdDev = 1;
		noGauss = true;
	}
	std::normal_distribution<double> distribution((rMax + rMin) / 2, stdDev);

	int iStep = (int)floor(0.5*(rMax + rMin) / gridStep);

	for (int i = sampleGridNodesY - iStep - 1; i > iStep; --i) {
		for (int j = iStep + 1; j < sampleGridNodesX - iStep; j++) {
			for (int k = iStep + 1; k < sampleGridNodesZ - iStep; k++) {// select the grid node
				gr.push_back(Grain());
				if (!noGauss) {
					do
						gr[grainNum].R = distribution(generator);
					while (!(gr[grainNum].R >= rMin - stdDev && gr[grainNum].R <= rMax + stdDev));
				}
				else {
					gr[grainNum].R = (double)rand() / (RAND_MAX + 1L)*(rMax - rMin) + rMin;
				}
				gr[grainNum].m = (4. / 3.) * M_PI*rho*gr[grainNum].R*gr[grainNum].R*gr[grainNum].R;
				gr[grainNum].I = (2. / 5.)*gr[grainNum].m*gr[grainNum].R*gr[grainNum].R;
				gr[grainNum].Ath = 0.; gr[grainNum].Ax = 0.; gr[grainNum].Ay = 0.; gr[grainNum].Az = 0.;
				gr[grainNum].Vth = 0.; gr[grainNum].Vx = 0.; gr[grainNum].Vy = 0.; gr[grainNum].Vz = 0.; gr[grainNum].th = 0.;
				gr[grainNum].x = j*gridStep;
				gr[grainNum].y = i*gridStep;
				gr[grainNum].z = k*gridStep;
				gr[grainNum].grainId = grainNum;
				gr[grainNum].isActive = true;
				grainNum++;
			}
		}
	}

}

void SetSample(double rMin, double rMax, double sizeX, double sizeY, double nodesRes, int dim, std::vector <Grain> &gr, double stdDev)
{
	double sizeH = sizeY;
	double sizeD = sizeX;

	if (stdDev >= rMin) {
		RuntimeLog(runTimeFile, true, "Standart deviation of Normal distribution is grater than minimal grain radius !", 0.0, false, false);
		RuntimeLog(runTimeFile, true, "EXIT with error !", 0.0, false, true);
		exit(EXIT_FAILURE);
	}

	if (dim == 2)
		SetSample2D(rMin, rMax, sizeX, sizeY, nodesRes, gr, stdDev);
	else
		SetSample3D(rMin, rMax, sizeD, sizeH, nodesRes, gr, stdDev);
}

double GetSampleHeight(double sizeY, int dim/*, struct Piston pist*/, std::vector <Grain> &gr)
{
	double maxHeight = 0;

	for (size_t i = 0; i < gr.size(); i++) {
		if (dim == 2) {
			if ((gr[i].x + gr[i].R) >= wTop.lowX && (gr[i].x - gr[i].R) <= wTop.highX) {
				if ((gr[i].y + gr[i].R) > maxHeight)
					maxHeight = gr[i].y + gr[i].R;
			}
		}
		else {
			if ((gr[i].x + gr[i].R) >= wTop.lowX && (gr[i].x - gr[i].R) <= wTop.highX && (gr[i].z + gr[i].R) >= wTop.lowZ && (gr[i].z - gr[i].R) <= wTop.highZ) {
				if ((gr[i].y + gr[i].R) > maxHeight)
					maxHeight = gr[i].y + gr[i].R;
			}
		}
	}

	return maxHeight;
}

void SetWalls(double sizeX, double sizeY, int dim, struct Piston pist, std::vector <Grain> &gr)
{
	double sizeZ = sizeX;

	// bottom
	wBottom.wF = 0.; wBottom.wV = 0.; wBottom.wA = 0.; wBottom.wM = 1.;

	wBottom.lowX = 0.; wBottom.highX = sizeX;
	wBottom.lowY = 0.; wBottom.highY = 0.;
	wBottom.lowZ = 0.; wBottom.highZ = sizeZ;

	// top
	wTop.wF = 0.; wTop.wV = 0.; wTop.wA = 0.; wTop.wM = 1.;

	if (!pist.movementFlag) {
		wTop.lowX = 0.; wTop.highX = sizeX;
		wTop.lowY = sizeY; wTop.highY = sizeY;
		wTop.lowZ = 0.; wTop.highZ = sizeZ;
	}
	else {
		double midPointX = sizeX / 2.;
		double midPointZ = sizeZ / 2.;
		wTop.lowX = midPointX - pist.sizeR; wTop.highX = midPointX + pist.sizeR;
		wTop.lowZ = midPointZ - pist.sizeR; wTop.highZ = midPointZ + pist.sizeR;
		double new_sizeY = GetSampleHeight(sizeY, dim/*, pist*/, gr);
		wTop.lowY = new_sizeY; wTop.highY = new_sizeY;
	}

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

	if (dim == 3)
	{
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
}

int ForcesGrain2D(int i, int j, std::vector <Grain> &gr, double _dt)
{
	double rXji, rYji, rNorm, dn;
	// Compute distance between grains
	rXji = gr[i].x - gr[j].x;
	rYji = gr[i].y - gr[j].y;
	rNorm = sqrt(rXji*rXji + rYji*rYji);
	dn = rNorm - (gr[i].R + gr[j].R);
	double r = std::min(gr[i].R, gr[j].R);
	double Reff = 2.*gr[i].R*gr[j].R / (gr[i].R + gr[j].R); // effective radius
	double nu_n = 2 * sqrt(rho*M_PI*r*r*kn)*beta;
	double nu_s = 2 * sqrt(rho*M_PI*r*r*kt)*beta;

	if (dn < 0.) // Contact
	{
		double Vxji, Vyji, Xn, Yn, Xt, Yt, Vn, Vt, Fn, Ft, Vthji;
		// Local axes
		Xn = rXji / rNorm;
		Yn = rYji / rNorm;
		Xt = -Yn;
		Yt = Xn;
		// Compute the velocity at contact
		Vxji = gr[i].Vx - gr[j].Vx;
		Vyji = gr[i].Vy - gr[j].Vy;
		Vthji = (gr[i].Vth - gr[j].Vth)*Reff; // Compute rolling velocity (from Radjai & Dubois, 2011)
		Vn = Vxji*Xn + Vyji*Yn;
		Vt = Vxji*Xt + Vyji*Yt - ((gr[i].R - dn / 2)*gr[i].Vth + (gr[j].R - dn / 2)*gr[j].Vth);

		// Compute force in local axes
		Fn = -kn * dn - nu_n * (Vn);
		if (Fn < 0.) Fn = 0.;

		Ft = fabs(kt*Vt) /**_dt+ nu_s * fabs(Vt)*/;
		if (fabs(Ft) > mu*fabs(Fn)) Ft = mu*fabs(Fn);
		if (Vt > 0.) Ft = -Ft;

		// Compute rolling momen
		double Moment = -Vthji*nur;
		double threshold = mur*Reff*(Fn + Adh);

		if (threshold < 0.)
			threshold = 0.;

		if (fabs(Moment) > threshold) {
			if (Vthji > 0.)
				Moment = -threshold;
			else
				Moment = threshold;
		}

		// Calculate sum of forces on i and j in global coordinates		
		gr[i].Fx += (Fn*Xn + Ft*Xt);
		gr[i].Fy += (Fn*Yn + Ft*Yt);
		gr[j].Fx += -(Fn*Xn + Ft*Xt);
		gr[j].Fy += -(Fn*Yn + Ft*Yt);
		gr[i].Fth += (-Ft*gr[i].R - Moment);
		gr[j].Fth += (-Ft*gr[j].R - Moment);

		// Add fn to pressure on grains i and j
		gr[i].p += Fn;
		gr[j].p += Fn;

	}
	return 0;
}

int ForcesGrain3D(int i, int j, std::vector <Grain> &gr, double _dt)
{
	double rXji, rYji, rZji, rNorm, dn;
	// Compute distance between grains
	rXji = gr[i].x - gr[j].x;
	rYji = gr[i].y - gr[j].y;
	rZji = gr[i].z - gr[j].z;
	rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
	dn = rNorm - (gr[i].R + gr[j].R);
	double r = std::min(gr[i].R, gr[j].R);
	double Reff = 2.*gr[i].R*gr[j].R / (gr[i].R + gr[j].R); // effective radius
	double nu = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;
	double nu_n = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;
	double nu_s = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kt)*beta;

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

		if (fabs(Moment) > threshold)
			if (Vthji > 0.)
				Moment = -threshold;
			else
				Moment = threshold;

		// Calculate sum of forces on i and j in global coordinates
		gr[i].Fx += Fn*Xn + Ft*Xt;
		gr[i].Fy += Fn*Yn + Ft*Yt;
		gr[i].Fz += Fn*Zn + Ft*Zt;
		gr[j].Fx += -Fn*Xn - Ft*Xt;
		gr[j].Fy += -Fn*Yn - Ft*Yt;
		gr[j].Fz += -Fn*Zn - Ft*Zt;
		gr[i].Fth += -Ft*gr[i].R - Moment;
		gr[j].Fth += -Ft*gr[j].R - Moment;

		// Add fn to pressure on grains i and j
		gr[i].p += Fn;
		gr[j].p += Fn;

	}
	return 0;
}

void ForcesGrain(int i, int j, int dim, std::vector <Grain> &gr, double _dt)
{
	if (dim == 2)
		ForcesGrain2D(i, j, gr, _dt);
	else
		ForcesGrain3D(i, j, gr, _dt);
}

void ForcesWall(int i, int dim, int wall, std::vector <Grain> &gr, struct Piston pist)
{
	if (dim == 2)
		ForcesWall2D(i, dim, wall, gr, pist);
	else
		ForcesWall3D(i, dim, wall, gr, pist);
}


void ForcesWall2D(int i, int dim, int wall, std::vector <Grain> &gr, struct Piston pist)
{
	double nu = 0.;
	double r = gr[i].R;
	nu = 2 * sqrt(rho*M_PI*r*r*kn)*beta;

	double dn;
	//bottom wall, identifier == 1
	if (wall == 1)
	{
		dn = gr[i].y - gr[i].R - wBottom.lowY;
		if (dn < 0.)
		{
			double vn, fn;
			// velocity
			vn = gr[i].Vy;
			fn = -kn * dn - nu * fabs(vn);
			if (fn < 0.) fn = 0.;
			if (vn < 0.) fn = -kn * dn;

			// Update sum of forces on grains
			gr[i].Fy = gr[i].Fy + fn;
			// Add fn to pressure on grains i
			gr[i].p += fn;
		}
	}
	//top wall, identifier == 2
	else if (wall == 2)
	{
		if (gr[i].y > wTop.highY)
		{
			if (!pist.movementFlag)
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
			if (pist.movementFlag)
			{
				/* we are below of the loading piston */
				if ((gr[i].y - gr[i].R) - wTop.highY < 0.)
				{
					bool tempFlag = false;

					//-----------------------------------------------------------
					if (gr[i].x < wTop.lowX)
					{
						if ( gr[i].x > wTop.lowX - gr[i].R)
						{
							double locX = wTop.lowX - gr[i].x;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = (gr[i].y - locY) - wTop.highY;
							tempFlag = true;
						}
					}
					if ( gr[i].x > wTop.highX)
					{
						if (gr[i].x < wTop.highX + gr[i].R)
						{
							double locX = gr[i].x - wTop.highX;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = (gr[i].y - locY) - wTop.highY;
							tempFlag = true;
						}
					}
					if ((gr[i].x >= wTop.lowX) && (gr[i].x <= wTop.highX))
					{
							dn = (gr[i].y - gr[i].R) - wTop.highY;
							tempFlag = true;
					}
					
					//-----------------------------------------------------------
					if (!tempFlag)
						dn = 0;

					if (dn < 0.)
					{
						double vn, fn;
						vn = gr[i].Vy;
						fn = -kn * dn - nu * fabs(vn);
						if (fn < 0.) fn = 0.;
						if (vn < 0.) fn = -kn * dn;
						gr[i].Fy = gr[i].Fy + fn;
						gr[i].p += fn;
					}
				}
			}
		}
		if (gr[i].y <= wTop.highY)
		{
			if (!pist.movementFlag)
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
			if (pist.movementFlag)
			{
				if (wTop.highY - (gr[i].y + gr[i].R) < 0.)
				{
					/* we are above of the loading piston */

					bool tempFlag = false;

					//-------------------------------------------------------------------------
					if (gr[i].x < wTop.lowX)
					{
						if (gr[i].x > wTop.lowX - gr[i].R)
						{
							double locX = wTop.lowX - gr[i].x;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = wTop.highY - (gr[i].y + locY);
							tempFlag = true;
						}
					}
					if (gr[i].x > wTop.highX)
					{
						if (gr[i].x < wTop.highX + gr[i].R)
						{
							double locX = gr[i].x - wTop.highX;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = wTop.highY - (gr[i].y + locY);
							tempFlag = true;
						}
					}
					if ((gr[i].x >= wTop.lowX) && (gr[i].x <= wTop.highX))
					{
							dn = wTop.highY - (gr[i].y + gr[i].R);
							tempFlag = true;
					}

					//-------------------------------------------------------------------------

					if (!tempFlag) dn = 0;

					if (dn < 0.)
					{
						double vn, fn;
						vn = gr[i].Vy + wTop.wV;
						fn = -kn * dn;
						if (vn < 0.) fn = -kn * dn - nu * fabs(vn);
						else if (fn < 0.) fn = 0.;
						gr[i].Fy = gr[i].Fy - fn;
						gr[i].p += fn;
					}
				}
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
			// velocity
			vn = gr[i].Vx;
			fn = -kn * dn - nu * fabs(vn);
			if (fn < 0.) fn = 0.;
			if (vn < 0.) fn = -kn * dn;
			// Update sum of forces on grains
			gr[i].Fx = gr[i].Fx + fn;
			// Add fn to pressure on grains i
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
			// Update sum of forces on grains
			gr[i].Fx = gr[i].Fx - fn;
			// Add fn to pressure on grains i
			gr[i].p += fn;
		}
	}
}

void ForcesWall3D(int i, int dim, int wall, std::vector <Grain> &gr, struct Piston pist)
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
			// Update sum of forces on grains
			gr[i].Fy = gr[i].Fy + fn;
			// Add fn to pressure on grains i
			gr[i].p += fn;
		}
	}
	//top wall, identifier == 2
	else if (wall == 2)
	{
		if (gr[i].y > wTop.highY)
		{
			if (!pist.movementFlag)
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
			if (pist.movementFlag)
			{
				/* we are below of the loading piston */
				if ((gr[i].y - gr[i].R) - wTop.highY < 0.)
				{
					bool tempFlag = false;

					//-----------------------------------------------------------
					if ((gr[i].x < wTop.lowX) && (gr[i].z > wTop.lowZ) && (gr[i].z < wTop.highZ))
					{
						if (gr[i].x > wTop.lowX - gr[i].R)
						{
							double locX = wTop.lowX - gr[i].x;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = (gr[i].y - locY) - wTop.highY;
							tempFlag = true;
						}
					}
					if ((gr[i].x > wTop.highX) && (gr[i].z < wTop.highZ) && (gr[i].z > wTop.lowZ))
					{
						if (gr[i].x < wTop.highX + gr[i].R)
						{
							double locX = gr[i].x - wTop.highX;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = (gr[i].y - locY) - wTop.highY;
							tempFlag = true;
						}
					}
					if ((gr[i].x >= wTop.lowX) && (gr[i].x <= wTop.highX))
					{
						if ((gr[i].z <= wTop.highZ) && (gr[i].z >= wTop.lowZ))
						{
							dn = (gr[i].y - gr[i].R) - wTop.highY;
							tempFlag = true;
						}
					}

					if ((gr[i].z < wTop.lowZ) && (gr[i].x > wTop.lowX) && (gr[i].x < wTop.highX))
					{
						if (gr[i].z > wTop.lowZ - gr[i].R)
						{
							double locZ = wTop.lowZ - gr[i].z;
							double locY = sqrt((gr[i].R - locZ) * (2.0 * gr[i].R - (gr[i].R - locZ)));
							dn = (gr[i].y - locY) - wTop.highY;
							tempFlag = true;
						}
					}
					if ((gr[i].z > wTop.highZ) && (gr[i].x > wTop.lowX) && (gr[i].x < wTop.highX))
					{
						if (gr[i].z < wTop.highZ + gr[i].R)
						{
							double locZ = gr[i].z - wTop.highZ;
							double locY = sqrt((gr[i].R - locZ) * (2.0 * gr[i].R - (gr[i].R - locZ)));
							dn = (gr[i].y - locY) - wTop.highY;
							tempFlag = true;
						}
					}
					//-----------------------------------------------------------
					if (!tempFlag)
						dn = 0;

					if (dn < 0.)
					{
						double vn, fn;
						vn = gr[i].Vy;
						fn = -kn * dn - nu * fabs(vn);
						if (fn < 0.) fn = 0.;
						if (vn < 0.) fn = -kn * dn;
						gr[i].Fy = gr[i].Fy + fn;
						gr[i].p += fn;
					}
				}
			}
		}
		if (gr[i].y <= wTop.highY)
		{
			if (!pist.movementFlag)
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
			if (pist.movementFlag)
			{
				if (wTop.highY - (gr[i].y + gr[i].R) < 0.)
				{
					/* we are above of the loading piston */

					bool tempFlag = false;

					//-------------------------------------------------------------------------
					if ((gr[i].x < wTop.lowX) && (gr[i].z > wTop.lowZ) && (gr[i].z < wTop.highZ))
					{
						if (gr[i].x > wTop.lowX - gr[i].R)
						{
							double locX = wTop.lowX - gr[i].x;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = wTop.highY - (gr[i].y + locY);
							tempFlag = true;
						}
					}
					if ((gr[i].x > wTop.highX) && (gr[i].z < wTop.highZ) && (gr[i].z > wTop.lowZ))
					{
						if (gr[i].x < wTop.highX + gr[i].R)
						{
							double locX = gr[i].x - wTop.highX;
							double locY = sqrt((gr[i].R - locX) * (2.0 * gr[i].R - (gr[i].R - locX)));
							dn = wTop.highY - (gr[i].y + locY);
							tempFlag = true;
						}
					}
					if ((gr[i].x >= wTop.lowX) && (gr[i].x <= wTop.highX))
					{
						if ((gr[i].z <= wTop.highZ) && (gr[i].z >= wTop.lowZ))
						{
							dn = wTop.highY - (gr[i].y + gr[i].R);
							tempFlag = true;
						}
					}

					if ((gr[i].z < wTop.lowZ) && (gr[i].x > wTop.lowX) && (gr[i].x < wTop.highX))
					{
						if (gr[i].z > wTop.lowZ - gr[i].R)
						{
							double locZ = wTop.lowZ - gr[i].z;
							double locY = sqrt((gr[i].R - locZ) * (2.0 * gr[i].R - (gr[i].R - locZ)));
							dn = wTop.highY - (gr[i].y + locY);
							tempFlag = true;
						}
					}
					if ((gr[i].z > wTop.highZ) && (gr[i].x > wTop.lowX) && (gr[i].x < wTop.highX))
					{
						if (gr[i].z < wTop.highZ + gr[i].R)
						{
							double locZ = gr[i].z - wTop.highZ;
							double locY = sqrt((gr[i].R - locZ) * (2.0 * gr[i].R - (gr[i].R - locZ)));
							dn = wTop.highY - (gr[i].y + locY);
							tempFlag = true;
						}
					}
					//-------------------------------------------------------------------------

					if (!tempFlag) dn = 0;

					if (dn < 0.)
					{
						double vn, fn;
						vn = gr[i].Vy + wTop.wV;
						fn = -kn * dn;
						if (vn < 0.) fn = -kn * dn - nu * fabs(vn);
						else if (fn < 0.) fn = 0.;
						gr[i].Fy = gr[i].Fy - fn;
						gr[i].p += fn;
					}
				}
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
			// Update sum of forces on grains
			gr[i].Fx = gr[i].Fx + fn;
			// Add fn to pressure on grains i
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
			// Update sum of forces on grains
			gr[i].Fx = gr[i].Fx - fn;
			// Add fn to pressure on grains i
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
			// Update sum of forces on grains
			gr[i].Fz = gr[i].Fz + fn;
			// Add fn to pressure on grains i
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
			// Update sum of forces on grains
			gr[i].Fz = gr[i].Fz - fn;
			// Add fn to pressure on grains i
			gr[i].p += fn;
		}
	}
}

int GetCorrectedSample(double rMax, double sizeX, double sizeY, int dim, std::vector <Grain> &gr, std::vector <Grain> &_gr)
{
	size_t j = 0;

	for (size_t i = 0; i < gr.size(); i++) {
		if (gr[i].y - gr[i].R < sizeY) {
			_gr.push_back(Grain());
			_gr[j] = gr[i];
			j++;
		}
	}
	grainNum = _gr.size();

	return 0;
}

// New functions, in case of cohesion

void ForcesGrainCohesive3D(int i, int j, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double _dt)
{
	double rXji, rYji, rZji, rNorm, dn;
	double Vxji, Vyji, Vzji, Vth, Xn, Yn, Zn, Xt, Yt, Zt, Vn, Vt;

	double r = std::min(gr[i].R, gr[j].R);
	double nu = 2 * sqrt((4. / 3.)*rho*M_PI*r*r*r*kn)*beta;
	double nuBond = 2 * sqrt(b[i][j].M * b[i][j].kn)*beta;

	// Compute distance between grains
	rXji = gr[i].x - gr[j].x;
	rYji = gr[i].y - gr[j].y;
	rZji = gr[i].z - gr[j].z;
	rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
	dn = rNorm - (gr[i].R + gr[j].R);

	double dnBond;
	double Fn = 0.;
	double Ft = 0.;

	if (b[i][j].isActive == true)
	{
		dnBond = rNorm - (gr[i].R + gr[j].R + b[i][j].gap); // bond overlap

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

		if (dn >= 0.) {
			Fn = -b[i][j].kn*dnBond - nuBond*Vn;
		}
		else {
			Fn = -kn*dn - nu*Vn;
			if (Fn < 0) Fn = 0.;
		}

		/*Ft = -b[i][j].kt*I1I2T;*/ //this variant does not works, needs debugging!
		Ft = fabs(b[i][j].kt*Vt);
		if (fabs(Ft) > mu*fabs(Fn)) Ft = mu*fabs(Fn);
		if (Vt > 0.) Ft = -Ft;/**/


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
		b[i][j].contactForce = sigmaStress + tauStress;

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

void ForcesGrainCohesive2D(int i, int j, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double _dt)
{
	double rXji, rYji, rNorm, dn;

	double r = std::min(gr[i].R, gr[j].R);
	double nu = 2 * sqrt(rho*M_PI*r*r*kn)*beta;
	// Compute distance between grains
	rXji = gr[i].x - gr[j].x;
	rYji = gr[i].y - gr[j].y;
	rNorm = sqrt(rXji*rXji + rYji*rYji);
	dn = rNorm - (gr[i].R + gr[j].R); // overlap

	double dnBond;
	double vxji, vyji, vthji, xn, yn, xt, yt, vn, vt;
	double fn = 0.;
	double ft = 0.;

	if (b[i][j].isActive == true)
	{
		double nuBond = 2 * sqrt(b[i][j].M * b[i][j].kn)*beta;

		dnBond = rNorm - (gr[i].R + gr[j].R + b[i][j].gap); // bond overlap

		vxji = gr[j].Vx - gr[i].Vx;  // sign convention: v<0 --> particles approach each other
		vyji = gr[j].Vy - gr[i].Vy;  // sign convention: v<0 --> particles approach each other
		vthji = gr[j].Vth - gr[i].Vth;

		// unit vectors (normal and tangential)
		xn = -rXji / rNorm;
		yn = -rYji / rNorm;
		xt = -yn;
		yt = xn;

		vn = vxji*xn + vyji*yn; // sign convention: vn<0 --> normal distance is decreasing
		vt = vxji*xt + vyji*yt - ((gr[i].R)*gr[i].Vth + (gr[j].R)*gr[j].Vth);

		double anglei = gr[i].th + b[i][j].alpha0;
		double anglej = gr[j].th + b[j][i].alpha0;

		// Bond elongation
		double I1I2X = -gr[i].R * cos(anglei) + gr[j].R * cos(anglej) - rXji;
		double I1I2Y = -gr[i].R * sin(anglei) + gr[j].R * sin(anglej) - rYji;

		// Update global position of interaction points
		b[i][j].xI = gr[i].x + gr[i].R * cos(anglei);
		b[i][j].yI = gr[i].y + gr[i].R * sin(anglei);

		b[j][i].xI = gr[j].x + gr[j].R * cos(anglej);
		b[j][i].yI = gr[j].y + gr[j].R * sin(anglej);

		// Bond shearing distance
		double I1I2T = I1I2X*xt + I1I2Y*yt;

		double gammaij = gr[i].th - gr[j].th;

		if (dn >= 0.) {//if traction or compression depends on sign of dnBond)
			fn = -b[i][j].kn*dnBond - nuBond*(vn);
		}
		else {
			fn = -kn*dn - nu*(vn);
			if (fn < 0.) fn = 0.;
		}

		//ft = -b[i][j].kt*(I1I2T);
		ft = fabs(b[i][j].kt*vt);
		if (fabs(ft) > mu*fabs(fn)) ft = mu*fabs(fn);
		if (vt > 0.) ft = -ft;/**/

		b[i][j].contactForce = fabs(ft) + fabs(fn);

		double couple = -b[i][j].M*(gammaij)*1.0;

		gr[i].Fx += -fn*xn - ft*xt;
		gr[i].Fy += -fn*yn - ft*yt;
		gr[j].Fx += (fn*xn + ft*xt);
		gr[j].Fy += (fn*yn + ft*yt);

		gr[i].Fth += -ft*gr[i].R + couple;
		gr[j].Fth += -ft*gr[j].R - couple;

		double sigmaStress = fabs(fn);
		double tauStress = fabs(ft);
		b[i][j].contactForce = sigmaStress + tauStress;

		// Check yield criterion
		double rupture = -1;
		
		if (deloneRapture)
			rupture = (ft / b[i][j].tau)*(ft / b[i][j].tau) + (couple / b[i][j].YM)*(couple / b[i][j].YM) + (fn / b[i][j].sigma) - 1;
		else {
			if (sigmaStress >= b[i][j].sigma || tauStress >= b[i][j].tau)
				rupture = 0.;
		}

		if (rupture >= 0.) {
			b[i][j].isActive = b[j][i].isActive = false;
		}

		return;

	}
	else ForcesGrain2D(i, j, gr, _dt);

	return;
}

void ForcesGrainCohesive(int i, int j, int dim, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double _dt)
{
	if (dim == 2)
		ForcesGrainCohesive2D(i, j, gr, b, _dt);
	else
		ForcesGrainCohesive3D(i, j, gr, b, _dt);
}

void IdentifyCohesiveBonds(int i, int j, int dim, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double *randValues)
{
	double xOiOj = gr[j].x - gr[i].x;
	double yOiOj = gr[j].y - gr[i].y;
	double zOiOj = gr[j].z - gr[i].z;
	double OiOj = 0.;
	
	if (dim == 2)
		OiOj = sqrt(xOiOj*xOiOj + yOiOj*yOiOj);
	else
		OiOj = sqrt(xOiOj*xOiOj + yOiOj*yOiOj + zOiOj*zOiOj);
	
	double dn = OiOj - gr[i].R - gr[j].R;
	double alphaXY;
	double alphaYZ;
	double alpha0;

	double gap = 0.2*std::min(gr[i].R, gr[j].R);

	//if (dn <= gap)
	if ((dn - gap) < 0.)
	{
		if (xOiOj == 0.) {
			if (yOiOj > 0) alphaXY = M_PI / 2.;
			else alphaXY = -M_PI / 2.;
		}
		else {
			alphaXY = atan2(yOiOj, xOiOj);
		}

		if (zOiOj == 0.)
		{
			if (yOiOj > 0) alphaYZ = M_PI / 2.;
			else alphaYZ = -M_PI / 2.;
		}
		else
		{
			alphaYZ = atan2(yOiOj, xOiOj);
		}

		b[i][j].alpha0 = alphaXY - gr[i].th;
		b[j][i].alpha0 = alphaXY + M_PI - gr[j].th;

		if (dim == 3)
		{
			b[i][j].alpha0YZ = alphaYZ - gr[i].th;
			b[j][i].alpha0YZ = alphaYZ + M_PI - gr[j].th;
		}

		b[i][j].isActive = true;
		/*if (dn > 0.) b[i][j].gap = dn;
		else */b[i][j].gap = 0.0;

		b[i][j].kn = b_kn * randValues[0];  // Tensile
		b[i][j].kt = b_kt * randValues[1];
		b[i][j].M = rho*gap*M_PI*(gap/2)*(gap / 2);
		//b[i][j].Yn = b_Yn;
		//b[i][j].Yt = b_Yt;
		b[i][j].YM = b_YM;

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

void SetBondsMatrix(std::vector<std::vector<bond>> &b, std::vector <Grain> &gr)
{
	b.clear();
	size_t grains = gr.size();

	for (auto i = 0; i < grains; i++) {
		std::vector<bond> vectb;
		for (auto j = 0; j < grains; j++) {
			bond b0;
			vectb.push_back(b0);
		}
		b.push_back(vectb);
	}
}

void FindBonds(int dim, std::vector<std::vector<bond>> &b, std::vector <Grain> &gr, double shape, double scale)
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

	size_t grains = gr.size();

	for (auto i = 0; i < grains; i++) {
		for (auto j = i + 1; j < grains; j++) {
			b[i][j].isActive = false;
			if (!notWeibull) {
				randValues[0] = distribution(generator);
				randValues[1] = distribution(generator);
			}
			IdentifyCohesiveBonds(i, j, dim, gr, b, randValues);
		}
	}
}

void IsBallInDomain(int dim, std::vector <Grain> &gr)
{
	std::vector <Grain>::iterator p;
	bool cond;

	for (size_t i = 0; i < gr.size(); i++) {
		cond = (gr[i].y > _sizeY * 3) || (gr[i].x < wLeft.lowX) || (gr[i].x > wRight.lowX);
		if (dim == 3)
			cond = (gr[i].y > _sizeY * 3) || (gr[i].x < wLeft.lowX) || (gr[i].x > wRight.lowX) || (gr[i].z < wFront.lowZ) || (gr[i].z > wBack.lowZ);
		if (cond) {
			p = gr.begin();
			p += i;
			gr.erase(p);
		}
	}
	grainNum = gr.size();
}

void GetPorouseSample(double sizeX, double sizeY, int dim, std::vector <Grain> &gr, double poros)
{
	double sampleVolume = 0.;
	double grainsVolume = 0.;
	double grainVolumeAver = 0.;
	size_t reqiredToRemoveGrainNum;
	std::vector <Grain>::iterator p;

	if (dim == 2)
		sampleVolume = sizeX * sizeY;
	else
		sampleVolume = sizeX * sizeX * sizeY;

	for (size_t i = 0; i < gr.size(); i++)
		grainsVolume += gr[i].m / rho;

	grainVolumeAver = grainsVolume / gr.size();

	reqiredToRemoveGrainNum = (size_t)ceil(gr.size() * poros);

	for (size_t i = 0; i < reqiredToRemoveGrainNum; i++) {
		p = gr.begin();
		p += rand() % gr.size();
		gr.erase(p);
	}
	grainNum = gr.size();
}

void ClearForces(std::vector <Grain> &gr)
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

double Get_dt(int dim, bool cohesion)
{
	double dt = 0.;
	extern double _rMin, _rMax;
	double _kn;

	if (cohesion)
		_kn = std::max(std::max(kn, b_kn), std::max(kt, b_kt));
	else
		_kn = std::max(kn, kt);

	double avR = (_rMin + _rMax) / 2;

	if (dim == 2)
		dt = sqrt(rho*M_PI*_rMin*_rMin / _kn)*0.001;
	else
		dt = sqrt(rho*M_PI*_rMin*_rMin*_rMin / _kn)*0.01;

	return dt;
}

bool CheckProjType(std::string prType)
{
	std::string smpl("sample");
	std::string sml("simulation");

	return sml.compare(prType);
}

int DynamicRelaxation(double timeInterval, int dim, std::vector <Grain> &gr, struct Piston pist, bool activeCohesion, std::vector<std::vector<bond>> &b, bool record, std::string prType, double _mechRat) {

	auto n_threads = std::thread::hardware_concurrency();

	double dt = Get_dt(dim, activeCohesion); /* time interval for integration, or in other words - time increment. */
	double mechRatio = 1.0e10;
	//if (_mechRat > dt*0.05) _mechRat = dt*0.05;

	int wallsNum; // number of walls
	if (dim == 2)
		wallsNum = 4;
	else
		wallsNum = 6;

	size_t stepNum = (size_t)ceil(timeInterval / dt);
	size_t loadSteps = (size_t)ceil(pist.loadTime / dt);
	/*double pistonLoad = pist.loadRate / loadSteps;*/
	double pistonLoad = pist.loadRate / dt;

	int mRatTimeIntrv = (int)ceil(0.03 / dt);
	int updContLst = (int)ceil(0.2 / sqrt(dt));
	
	int outputPeriod = updContLst * 10; /* record output at the same time as updating the list of grains contacts. */

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
	RuntimeLog(runTimeFile, false, "Number of grains:", gr.size(), false, false, true);
	RuntimeLog(runTimeFile, false, "Critical timestep:", dt, false, false, true);
	RuntimeLog(runTimeFile, false, "Interval for contact list updates:", updContLst, false, false, true);
	
	if (CheckProjType(prType)) {
		wrtMRatio.open(fileLocation, std::ios::out);
		RuntimeLog(runTimeFile, false, "Allowed error value of sample equilibrium:", _mechRat, false, false, true);
	}
	else
		RuntimeLog(runTimeFile, false, "Constant loading piston acceleration at each relaxation cycle:", pistonLoad, false, false, true);

	RuntimeLog(runTimeFile, false, "Interval for error value updates:", mRatTimeIntrv, false, false, true);
	RuntimeLog(runTimeFile, false, "Allowed maximum value of relaxation cycles:", stepNum, false, false, true);
	RuntimeLog(runTimeFile, true, "CALCULATING ...  ", 0, false, false);

	ClearForces(gr);

	size_t grains = gr.size();
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
			GetContactsList(dim, gr, b, activeCohesion);

#ifdef DEBUG
		t1_end = chrono::steady_clock::now();
		diff_1 = t1_end - t1_start;
		t1 = t1 + chrono::duration <double, milli>(diff_1).count() / 1000; /* time in seconds*/
		t2_start = chrono::steady_clock::now();
#endif

		if (step % 2*updContLst == 0) {
			if (wLoadPiston.movementFlag)
				IsBallInDomain(dim, gr);
		}

#ifdef DEBUG
		t2_end = chrono::steady_clock::now();
		diff_2 = t2_end - t2_start;
		t2 = t2 + chrono::duration <double, milli>(diff_2).count() / 1000; /* time in seconds*/
		t3_start = chrono::steady_clock::now();
#endif

		/* Velocity Verlet's prediction Step */
		for (auto i = 0; i < grains; i++) {
			
			// update positions
			gr[i].x = gr[i].x + dt * gr[i].Vx + 0.5*dt*dt*gr[i].Ax;
			gr[i].y = gr[i].y + dt * gr[i].Vy + 0.5*dt*dt*gr[i].Ay;
			if (dim == 3) gr[i].z = gr[i].z + dt * gr[i].Vz + 0.5*dt*dt*gr[i].Az;
			gr[i].th = gr[i].th + dt * gr[i].Vth + 0.5*dt*dt*gr[i].Ath;

			// update velocities
			gr[i].Vx = gr[i].Vx + 0.5 * dt * gr[i].Ax;
			gr[i].Vy = gr[i].Vy + 0.5 * dt * gr[i].Ay;
			if (dim == 3) gr[i].Vz = gr[i].Vz + 0.5 * dt * gr[i].Az;
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

		if (dim == 3) {
			wFront.lowZ = wFront.lowZ + dt * wFront.wV + 0.5 * dt * dt * wFront.wA;
			wBack.highZ = wBack.highZ + dt * wBack.wV + 0.5 * dt * dt * wBack.wA;
			wFront.wV = wFront.wV + 0.5 * dt * wFront.wA;
			wBack.wV = wBack.wV + 0.5 * dt * wBack.wA;
		}

		/* calculate mech ratio, only for 'sample' project */
		if (CheckProjType(prType)) {

			double _sumOfNormF = 1.0;

			for (auto i = 0; i < grainNum; i++) {

				if (dim == 3)
					gr[i].sumOfNormF = sqrt(gr[i].Fx*gr[i].Fx + gr[i].Fy*gr[i].Fy + gr[i].Fz*gr[i].Fz) + sqrt(gr[i].Fth*gr[i].Fth);
				else
					gr[i].sumOfNormF = sqrt(gr[i].Fx*gr[i].Fx + gr[i].Fy*gr[i].Fy) + sqrt(gr[i].Fth*gr[i].Fth);

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

		/*Set sum of forces acting on grain to 0*/
		for (auto i = 0; i < grains; i++) {
			gr[i].Fx = 0.;
			gr[i].Fy = 0.;
			if (dim == 3) gr[i].Fz = 0.;
			gr[i].Fth = 0.;
			gr[i].p = 0.;
		}

		/*Set confining pressure*/
		wBottom.wF = 0;
		if (!pist.movementFlag)
			wTop.wF = 0;
		else {
			//int loadSteps = (int) /*round*/ceil(pist.loadTime / dt);
			if (step <= loadSteps) {
				wTop.wF = pistonLoad;
			}
			else {
				wTop.wF = 0;
				wTop.wV = 0; /*std::cout<<"NO LOAD!"<<std::endl;*/
			}
		}

		wLeft.wF = 0;
		wRight.wF = 0;

		if (dim == 3) {
			wFront.wF = 0;
			wBack.wF = 0;
		}

#ifdef DEBUG
		t3_end = chrono::steady_clock::now();
		diff_3 = t3_end - t3_start;
		t3 = t3 + chrono::duration <double, milli>(diff_3).count() / 1000; /* time in seconds*/
		t4_start = chrono::steady_clock::now();
#endif

		/*Force computation between grains*/
#pragma omp parallel for 
		for (auto i = 0; i < grains; i++) {
			for (auto j = gr[i].contactList.begin(); j != gr[i].contactList.end(); ++j) {
				if (activeCohesion)
					ForcesGrainCohesive(i, *j, dim, gr, b, dt);
				else
					ForcesGrain(i, *j, dim, gr, dt);
			}
		}

#ifdef DEBUG
		t4_end = chrono::steady_clock::now();
		diff_4 = t4_end - t4_start;
		t4 = t4 + chrono::duration <double, milli>(diff_4).count() / 1000; /* time in seconds*/
		t5_start = chrono::steady_clock::now();
#endif

		/*Force computation between grains and walls*/
#pragma omp parallel for
		for (auto i = 0; i < grains; i++)
			for (auto w = 1; w <= wallsNum; w++)
				ForcesWall(i, dim, w, gr, pist);

#ifdef DEBUG
		t5_end = chrono::steady_clock::now();
		diff_5 = t5_end - t5_start;
		t5 = t5 + chrono::duration <double, milli>(diff_5).count() / 1000; /* time in seconds*/
		t6_start = chrono::steady_clock::now();
#endif

		/*Calculate corrected acceleration*/
		for (auto i = 0; i < grains; i++) {
			gr[i].Ax = gr[i].Fx / gr[i].m;
			gr[i].Ay = gr[i].Fy / gr[i].m - Grav;
			if (dim == 3) gr[i].Az = gr[i].Fz / gr[i].m;
			gr[i].Ath = gr[i].Fth / gr[i].I;
		}

		wBottom.wA = wBottom.wF / wBottom.wM;
		wTop.wA = wTop.wF / wTop.wM;
		wLeft.wA = wLeft.wF / wLeft.wM;
		wRight.wA = wRight.wF / wRight.wM;

		if (dim == 3) {
			wFront.wA = wFront.wF / wFront.wM;
			wBack.wA = wBack.wF / wBack.wM;
		}

		/* Correction step for velocities */
		for (auto i = 0; i < grains; i++) {
			gr[i].Vx = gr[i].Vx + 0.5 * dt * gr[i].Ax;
			gr[i].Vy = gr[i].Vy + 0.5 * dt * gr[i].Ay;
			if (dim == 3) gr[i].Vz = gr[i].Vz + 0.5 * dt * gr[i].Az;
			gr[i].Vth = gr[i].Vth + 0.5 * dt * gr[i].Ath;
		}

		wBottom.wV = wBottom.wV + 0.5 * dt*wBottom.wA;
		wTop.wV = /*wTop.wV +*/ 0.5 * dt*wTop.wA;/*0.5 * wTop.wA;*/ //constant loading with no acceleration
		wLeft.wV = wLeft.wV + 0.5 * dt*wLeft.wA;
		wRight.wV = wRight.wV + 0.5 * dt*wRight.wA;

		if (dim == 3) {
			wFront.wV = wFront.wV + 0.5 * dt * wFront.wA;
			wBack.wV = wBack.wV + 0.5 * dt * wBack.wA;
		}

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

		/* record the results of simulation */
		if (step%outputPeriod == 0 && record)
			GetSnapshots(activeCohesion, gr, pist, b);

		if (step % outputPeriod == 0) {
			if (CheckProjType(prType))
				RuntimeLog(runTimeFile, true, "Sample equilibrium error  ", mechRatio, false, false,true);
			else
				RuntimeLog(runTimeFile, true, "Completed, %", (int)((step * dt * 100) / timeInterval)+1, false, false, true);
		}

		/* This part is only for the 'sample' project.
		If the producing sample is in equilibrium: exit from relaxation part. */
		if (CheckProjType(prType) && step % mRatTimeIntrv == 0) {
			if (mechRatio < _mechRat) {
				if (!record)
					GetSnapshots(activeCohesion, gr, pist, b);
				
				break;
			}
		}

		/* writes the last snapshot of the sample if the project has no record option. */
		if (!record && step == (stepNum - 1))
			GetSnapshots(activeCohesion, gr, pist, b);

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

std::vector<shortBond> GetSnapshotB(std::vector<std::vector<bond>> &b1)
{
	std::vector<shortBond> allBnds;
	for (auto k = 0; k < grainNum; k++) {
		for (auto l = k + 1; l < grainNum; l++) {
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

int GetSnapshots(bool _activeCohesion, std::vector <Grain> &gr, struct Piston pist, std::vector<std::vector<bond>> &b1)
{
	/*omp_set_num_threads(1);*/
#pragma omp parallel sections
	{
#pragma omp section
		{
			/*pstSmf = false;
			pst.push_back(pist);
			long int psz = pstCount - 1;
			if (psz >= snapshotBuffer)
			{
				if (pstCount <= pst.size())
				{
					pst.erase(pst.begin(), pst.begin() + psz);
					pstTop = true;
				}
			}
			pstSmf = true;*/
			
			pstSmf = false;
			_pst.push_back(allWalls);
			long int psz = pstCount - 1;
			if (psz >= snapshotBuffer)
			{
				if (pstCount <= _pst.size())
				{
					_pst.erase(_pst.begin(), _pst.begin() + psz);
					pstTop = true;
				}
			}
			pstSmf = true;

			wlSmf = false;
			wl.push_back(allWalls);
			psz = wlCount - 1;
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
			if (_activeCohesion)
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
	}

	snapshotCounter++;

	return 0;
}

int PostProcPiston(int numfile, int dim, std::vector<WALLS> &_pst/*std::vector<Piston> &pst*/)
{
	if (pstTop) {
		pstCount = 0;
		pstTop = false;
	}
	int p_sz = _pst.size() - 1;
	if (pstSmf && pstCount <= p_sz)
	{
		/*Piston _wTop = pst.at(pstCount);
		pstCount++;*/
		WALLS _allWalls = _pst.at(pstCount);
		pstCount++;

		Wall wBottom = _allWalls._wBottom;
		Wall wTop = _allWalls._wTop;
		Wall wLeft = _allWalls._wLeft;
		Wall wRight = _allWalls._wRight;
		Wall wFront = _allWalls._wFront;
		Wall wBack = _allWalls._wBack;


		char fnameP[256]; // file name

		extern std::string type;
		std::string smpl("sample");

		#ifdef _WIN64
		if (!smpl.compare(type)) snprintf(fnameP, 200, "sample\\vtk_files\\piston%04d.vtk", numfile);
		else snprintf(fnameP, 200, "simulation\\vtk_files\\piston%04d.vtk", numfile);
		#endif

		#ifdef linux
		if (!smpl.compare(type)) snprintf(fnameP, 200, "sample/vtk_files/piston%04d.vtk", numfile);
		else snprintf(fnameP, 200, "simulation/vtk_files/piston%04d.vtk", numfile);
		#endif

		std::ofstream wrtToFileP(fnameP, std::ios::out);

		if (!wrtToFileP.is_open()) {
			RuntimeLog(runTimeFile, true, "Error when opening file for solution (in PostProcPiston())! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);			
			exit(1);
		}

		if (wrtToFileP)
		{
			if (dim == 2) {

				wrtToFileP.precision(3); wrtToFileP << std::scientific;
				wrtToFileP << "# vtk DataFile Version 3.0" << std::endl;
				wrtToFileP << "My walls" << std::endl;
				wrtToFileP << "ASCII" << std::endl;
				wrtToFileP << "DATASET POLYDATA" << std::endl;
				wrtToFileP << "POINTS 2 float" << std::endl;

				// upper wall
				wrtToFileP << wTop.lowX << " " << wTop.highY << " 0" << std::endl;
				wrtToFileP << wTop.highX << " " << wTop.highY << " 0" << std::endl;

				wrtToFileP << "LINES 1 3" << std::endl;
				wrtToFileP << "2 0 1" << std::endl;

			}
			if (dim == 3) {

				wrtToFileP.precision(3); wrtToFileP << std::scientific;
				wrtToFileP << "# vtk DataFile Version 3.0" << std::endl;
				wrtToFileP << "My walls" << std::endl;
				wrtToFileP << "ASCII" << std::endl;
				wrtToFileP << "DATASET POLYDATA" << std::endl;
				wrtToFileP << "POINTS 4 float" << std::endl;

				// Top wall
				wrtToFileP << wTop.lowX << " " << wTop.highY << " " << wTop.lowZ << std::endl;
				wrtToFileP << wTop.highX << " " << wTop.highY << " " << wTop.lowZ << std::endl;
				wrtToFileP << wTop.highX << " " << wTop.highY << " " << wTop.highZ << std::endl;
				wrtToFileP << wTop.lowX << " " << wTop.highY << " " << wTop.highZ << std::endl;

				wrtToFileP << "POLYGONS 1 5" << std::endl;

				wrtToFileP << "4 0 1 2 3" << std::endl;
			}
		}

		wrtToFileP.close();

		snapshotCounterP++;

		if (!relaxationFlag)
		{
			int psz = pstCount - 1;
			if (psz >= snapshotBuffer)
			{
				_pst.erase(_pst.begin(), _pst.begin() + psz);
				pstTop = true;
			}
		}
	}
	else if (!relaxationFlag) pst.clear();

	return 0;
}

int PostProcGrain(int numfile, int dim)
{

	if (grnTop) {
		grnCount = 0;
		grnTop = false;
	}

	int g_sz = grainSnapshot.size() - 1;

	if (grnSmf && grnCount <= g_sz) {
		std::vector <Grain> grn;
		try
		{
			grn = grainSnapshot.at(grnCount);
		}
		catch (const std::out_of_range& e) {
			RuntimeLog(runTimeFile, true, "Out of range error -> Grain! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
		}
		grnCount++;

		size_t grainNum = grn.size();
		size_t i;
		char fnameG[256]; // file name

		#ifdef _WIN64
		if (CheckProjType(type))
			snprintf(fnameG, 200, "sample\\vtk_files\\grains%04d.vtk", numfile);
		else
			snprintf(fnameG, 200, "simulation\\vtk_files\\grains%04d.vtk", numfile);
		#endif

		#ifdef linux
		if (CheckProjType(type))
			snprintf(fnameG, 200, "sample/vtk_files/grains%04d.vtk", numfile);
		else
			snprintf(fnameG, 200, "simulation/vtk_files/grains%04d.vtk", numfile);
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
			wrtToFileG << "My grains" << std::endl;
			wrtToFileG << "ASCII" << std::endl;
			wrtToFileG << "DATASET UNSTRUCTURED_GRID" << std::endl;
			wrtToFileG << "POINTS " << grainNum << " float" << std::endl;
			if (dim == 2) {
				for (i = 0; i < grainNum; i++)
					wrtToFileG << grn[i].x << " " << grn[i].y << " 0.0" << std::endl;
			}
			if (dim == 3) {
				for (i = 0; i < grainNum; i++)
					wrtToFileG << grn[i].x << " " << grn[i].y << " " << grn[i].z << std::endl;
			}
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

int PostProcGrainVisual(int dim, std::vector <std::vector <float> > &data)
{
	std::vector <std::vector <float> > t_data;

	if (grnTop) {
		grnCount = 0;
		grnTop = false;
	}

	int g_sz = grainSnapshot.size() - 1;

	if (grnSmf && grnCount <= g_sz) {
		std::vector <Grain> grn;
		try
		{
			grn = grainSnapshot.at(grnCount);
		}
		catch (const std::out_of_range& e) {
			RuntimeLog(runTimeFile, true, "Out of range error -> Grain! Call exit(EXIT_FAILURE) function !):", 0.0, false, true);
		}
		grnCount++;

		size_t grainNum = grn.size();
		size_t i;

		std::vector <float> tvect(6, 0.0f);

		// get sacaler for pressure
		double pscaler = 0.0;
		double avoidzerro = 0.1;
		for (i = 0; i < grainNum; i++){
			double pmax = grn[i].p;
			if ( pmax > pscaler )
				pscaler = pmax;
		}
		pscaler = pscaler + avoidzerro;
		
		double scaler = 0.083;

		for (i = 0; i < grainNum; i++){
			tvect[0] = static_cast <float> ( grn[i].x / scaler );
			tvect[1] = static_cast <float> ( grn[i].y / scaler );
			if (dim == 2)
				tvect[1] = 0.0;
			else
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

int PostProcForce(int numfile, int dim)
{
	if (fTop) {
		fCount = 0;
		fTop = false;
	}
	int f_sz = grainSnapshotF.size() - 1;
	if (fSmf && fCount <= f_sz)
	{
		std::vector <Grain> grn = grainSnapshotF.at(fCount);
		fCount++;

		size_t grainNum = grn.size();

		size_t maxForceNum = grainNum * 12;
		double *forceFn; forceFn = new double[maxForceNum];
		size_t *force_i; force_i = new size_t[maxForceNum];
		size_t *force_j; force_j = new size_t[maxForceNum];
		size_t forceNum = 0; // number of forces
		char fnameF[256]; // file name

		extern std::string type;
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
			for (auto i = 0; i < grainNum; i++)
			{
				for (auto j = i + 1; j < grainNum; j++)
				{
					double rXji, rYji, rZji, dn;
					double rNorm = 0.;
					// Compute distance between grains
					rXji = grn[i].x - grn[j].x;
					rYji = grn[i].y - grn[j].y;
					if (dim == 2) rNorm = sqrt(rXji*rXji + rYji*rYji);
					if (dim == 3) { rZji = grn[i].z - grn[j].z; rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji); }
					dn = rNorm - (grn[i].R + grn[j].R);
					if (dn < 0.) // Contact
					{
						forceFn[forceNum] = -kn * dn;
						force_i[forceNum] = i;
						force_j[forceNum] = j;
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
			for (auto k = 0; k < forceNum; k++)
			{
				if (dim == 2) {
					wrtToFileF << grn[force_i[k]].x << " " << grn[force_i[k]].y << " 0" << std::endl;
					wrtToFileF << grn[force_j[k]].x << " " << grn[force_j[k]].y << " 0" << std::endl;
				}
				if (dim == 3) {
					wrtToFileF << grn[force_i[k]].x << " " << grn[force_i[k]].y << " " << grn[force_i[k]].z << std::endl;
					wrtToFileF << grn[force_j[k]].x << " " << grn[force_j[k]].y << " " << grn[force_j[k]].z << std::endl;
				}
			}
			wrtToFileF << "LINES " << forceNum << " " << 3 * forceNum << std::endl;
			for (auto k = 0; k < forceNum; k++)
			{
				wrtToFileF << "2 " << 2 * k << " " << 2 * k + 1 << std::endl;
			}
			wrtToFileF << "POINT_DATA " << 2 * forceNum << std::endl;
			wrtToFileF << "SCALARS Force float" << std::endl;
			wrtToFileF << "LOOKUP_TABLE default" << std::endl;
			for (auto k = 0; k < forceNum; k++)
			{
				wrtToFileF << forceFn[k] << std::endl;
				wrtToFileF << forceFn[k] << std::endl;
			}
		}

		wrtToFileF.close();

		snapshotCounterF++;

		delete[] force_i;
		delete[] force_j;
		delete[] forceFn;

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
int PostProcBond(int numfile, int dim)
{
	if (bTop) {
		bCount = 0;
		bTop = false;
	}

	int b_sz = grainSnapshotB.size() - 1;

	if (bSmf && bCount <= b_sz) {

		std::vector <Grain> grn;
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
		size_t _grainNum = grn.size();
		std::vector <double> forceFn;
		std::vector <size_t> force_i;
		std::vector <size_t> force_j;
		size_t forceNum = 0; // number of forces
		char fnameB[256]; // file name

		extern std::string type;
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

			for (auto k = 0; k < bondsNum; k++) {

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
			
			for (auto k = 0; k < forceNum; k++) {

				if (dim == 2) {
					wrtToFileB << grn[force_i[k]].x << " " << grn[force_i[k]].y << " 0" << std::endl;
					wrtToFileB << grn[force_j[k]].x << " " << grn[force_j[k]].y << " 0" << std::endl;
				}
				else {
					wrtToFileB << grn[force_i[k]].x << " " << grn[force_i[k]].y << " " << grn[force_i[k]].z << std::endl;
					wrtToFileB << grn[force_j[k]].x << " " << grn[force_j[k]].y << " " << grn[force_j[k]].z << std::endl;
				}
			}

			wrtToFileB << "LINES " << forceNum << " " << 3 * forceNum << std::endl;
			
			for (auto k = 0; k < forceNum; k++)
				wrtToFileB << "2 " << 2 * k << " " << 2 * k + 1 << std::endl;

			wrtToFileB << "POINT_DATA " << 2 * forceNum << std::endl;
			wrtToFileB << "SCALARS Force float" << std::endl;
			wrtToFileB << "LOOKUP_TABLE default" << std::endl;
			
			for (auto k = 0; k < forceNum; k++) {

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

int PostProcWall(int numfile, int dim, std::vector<WALLS> &wl)
{
	if (wlTop) {
		wlCount = 0;
		wlTop = false;
	}
	int w_sz = wl.size() - 1;
	if (wlSmf && wlCount <= w_sz)
	{
		//WALLS _allWalls = wl.front();
		WALLS _allWalls = wl.at(wlCount);
		wlCount++;

		Wall wBottom = _allWalls._wBottom;
		Wall wTop = _allWalls._wTop;
		Wall wLeft = _allWalls._wLeft;
		Wall wRight = _allWalls._wRight;
		Wall wFront = _allWalls._wFront;
		Wall wBack = _allWalls._wBack;

		char fnameW[256]; // file name

		extern std::string type;
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

			if (dim == 2) {

				wrtToFileW.precision(3); wrtToFileW << std::scientific;
				wrtToFileW << "# vtk DataFile Version 3.0" << std::endl;
				wrtToFileW << "My walls" << std::endl;
				wrtToFileW << "ASCII" << std::endl;
				wrtToFileW << "DATASET POLYDATA" << std::endl;


				wrtToFileW << "POINTS 8 float" << std::endl;
				// lower wall
				wrtToFileW << wBottom.lowX << " " << wBottom.lowY << " 0" << std::endl;
				wrtToFileW << wBottom.highX << " " << wBottom.lowY << " 0" << std::endl;
				// upper wall
				wrtToFileW << wTop.lowX << " " << wTop.highY << " 0" << std::endl;
				wrtToFileW << wTop.highX << " " << wTop.highY << " 0" << std::endl;
				// left wall
				wrtToFileW << wLeft.lowX << " " << wLeft.lowY << " 0" << std::endl;
				wrtToFileW << wLeft.lowX << " " << wLeft.highY << " 0" << std::endl;
				// right wall
				wrtToFileW << wRight.highX << " " << wRight.lowY << " 0" << std::endl;
				wrtToFileW << wRight.highX << " " << wRight.highY << " 0" << std::endl;

				wrtToFileW << "LINES 4 12" << std::endl;
				wrtToFileW << "2 0 1" << std::endl;
				wrtToFileW << "2 2 3" << std::endl;
				wrtToFileW << "2 4 5" << std::endl;
				wrtToFileW << "2 6 7" << std::endl;

				wrtToFileW << "POINT_DATA 8" << std::endl;
				wrtToFileW << "SCALARS Rayon float" << std::endl;
				wrtToFileW << "LOOKUP_TABLE default" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
				wrtToFileW << "1" << std::endl;
			}
			if (dim == 3) {

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

int PostProcessing(bool _activeCohesion, int dim, std::vector<WALLS> &_pst/*std::vector<Piston> &pst*/, std::vector<WALLS> &wl)
{
	/*omp_set_num_threads(1);*/
#pragma omp parallel sections
	{
#pragma omp section
		{
			if (!_pst.empty())
				PostProcPiston(snapshotCounterP, dim, _pst);
		}
		
#pragma omp section
		{
			if (!grainSnapshot.empty())
				PostProcGrain(snapshotCounterG, dim);
		}
		
#pragma omp section
		{
			if (!grainSnapshotF.empty())
				PostProcForce(snapshotCounterF, dim);
		}
		
#pragma omp section
		{
			if (_activeCohesion)
				if (!allBondsSnapshot.empty() )
					PostProcBond(snapshotCounterB, dim);
		}
		
#pragma omp section
		{
			if (!wl.empty())
				PostProcWall(snapshotCounterW, dim, wl);
		}
	}
	
	return 0;
}

int RunPostProcessing(bool _activeCohesion, int dim)
{
	do {
		PostProcessing(_activeCohesion, dim, _pst, wl);
	} while (!grainSnapshot.empty() || relaxationFlag);

	return 0;
}

int RunVisualPostProcessing(bool _activeCohesion, int dim)
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
				//bool emptyFlag = grainSnapshot.empty(); std::cout<<"grainSnapshot.size() = "<<grainSnapshot.size()<<std::endl;
				//std::cout<<"grainSnapshot.size() = "<<grainSnapshot.size()<<"; renderSignal = "<<renderSignal<<"; accessLock = "<<accessLock<<std::endl;
				
				//if ( !grainSnapshot.empty() && renderSignal ){
					//std::cout<<"here 1"<<std::endl;
					if ( accessLock && renderSignal ){
						accessLock = false;
						PostProcGrainVisual(dim, grainData);
						accessLock = true;
						renderSignal = false;
					}

					//std::cout<<"grainData.size() = "<<grainData.size()<<std::endl;
				//}
				//std::cout<<"here 3"<<std::endl;
				//std::cout<<"passing throuugh! accessLock = "<<accessLock<<std::endl;
			} while (!grainSnapshot.empty() || relaxationFlag);
			
			escSignal = true;
			//std::cout<<"escSignal = "<<escSignal<<std::endl;
		}
#pragma omp section
		{
			oglu::visualizer scene/*(grainData)*/;
			scene.show(escSignal, accessLock, renderSignal, grainData);
		}
	}

	return 0;
}

int GetContactsList(int dim, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, bool _activeCohesion)
{
#pragma omp parallel for
	for (auto i = 0; i < gr.size(); i++) {
		gr[i].contactList.clear();
		for (auto j = i + 1; j < gr.size(); j++) {
			double rXji, rYji, rZji, dn;
			double rNorm = 0.;
			// Compute distance between grains
			rXji = gr[i].x - gr[j].x;
			rYji = gr[i].y - gr[j].y;
			if (dim == 2)
				rNorm = sqrt(rXji*rXji + rYji*rYji);
			else {
				rZji = gr[i].z - gr[j].z;
				rNorm = sqrt(rXji*rXji + rYji*rYji + rZji*rZji);
			}
			dn = rNorm - (gr[i].R + gr[j].R);
			if (_activeCohesion) {
				if (b[i][j].isActive == true) {
					gr[i].contactList.push_back(j);
				}
				else if (dn < 0.)
					gr[i].contactList.push_back(j);
			}
			else if (dn < 0.)
				gr[i].contactList.push_back(j);
		}
	}
	return 0;
}

void RuntimeLog(std::string fileName, bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose, bool showValue) {

	char ctime_no_newline[26];
	char* _ctime_no_newline;
	char delim[] = "\n";
	char *next_token;

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