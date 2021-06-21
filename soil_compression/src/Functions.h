#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "allincludes.h"
#include "Structures.h"
#include "ExternalDeclarations.h"


int ReadParametersFromFile(const char* paramFileName);

int ProjectManager();

bool DirExists(const std::string& dirName_in);

std::string ws2s(const std::wstring& wstr);

int FinalizeSampleProject(std::string paramFileName, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, struct Piston pist, std::vector<WALLS> &wl);

int InitializeSimulationProject(const char* paramFileName, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, struct Piston pist, std::vector<WALLS> &wl);

void InitiateParameters(double rMin, double rMax);

void SetSample2D(double rMin, double rMax, double sizeX, double sizeY, double nodesRes, std::vector <Grain> &gr, double stdDev);

void SetSample3D(double rMin, double rMax, double sizeD, double sizeH, double nodesRes, std::vector <Grain> &gr, double stdDev);

void SetSample(double rMin, double rMax, double sizeX, double sizeY, double nodesRes, int dim, std::vector <Grain> &gr, double stdDev);

double GetSampleHeight(double sizeY, int dim, std::vector <Grain> &gr);

void SetWalls(double sizeX, double sizeY, int dim, struct Piston pist, std::vector <Grain> &gr);

int ForcesGrain2D(int i, int j, std::vector <Grain> &gr, double _dt);

int ForcesGrain3D(int i, int j, std::vector <Grain> &gr, double _dt);

void ForcesGrain(int i, int j, int dim, std::vector <Grain> &gr, double _dt);

void ForcesWall(int i, int dim, int wall, std::vector <Grain> &gr, struct Piston pist);

void ForcesWall2D(int i, int dim, int wall, std::vector <Grain> &gr, struct Piston pist);

void ForcesWall3D(int i, int dim, int wall, std::vector <Grain> &gr, struct Piston pist);

int GetCorrectedSample(double rMax, double sizeX, double sizeY, int dim, std::vector <Grain> &gr, std::vector <Grain> &_gr);

void ForcesGrainCohesive3D(int i, int j, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double _dt);

void ForcesGrainCohesive2D(int i, int j, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double _dt);

void ForcesGrainCohesive(int i, int j, int dim, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double _dt);

void IdentifyCohesiveBonds(int i, int j, int dim, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, double *randValues);

void SetBondsMatrix(std::vector<std::vector<bond>> &b, std::vector <Grain> &gr);

void FindBonds(int dim, std::vector<std::vector<bond>> &b, std::vector <Grain> &gr, double shape, double scale);

void IsBallInDomain(int dim, std::vector <Grain> &gr);

void GetPorouseSample(double sizeX, double sizeY, int dim, std::vector <Grain> &gr, double poros);

void ClearForces(std::vector <Grain> &gr);

double Get_dt(int dim, bool cohesion);

bool CheckProjType(std::string prType);

int DynamicRelaxation(double timeInterval, int dim, std::vector <Grain> &gr, struct Piston pist, bool activeCohesion, std::vector<std::vector<bond>> &b, bool record, std::string prType, double _mechRat);

std::vector<shortBond> GetSnapshotB(std::vector<std::vector<bond>> &b1);

int GetSnapshots(bool _activeCohesion, std::vector <Grain> &gr, struct Piston pist, std::vector<std::vector<bond>> &b1/*, std::vector<Piston> &pst, std::vector<std::vector <Grain>> &grainSnapshot, std::vector<std::vector<std::vector<bond> >> &bondSnapshot, std::vector<WALLS> &wl, std::vector<std::vector <Grain>> &grainSnapshotF, std::vector<std::vector <Grain>> &grainSnapshotB*/);

int PostProcPiston(int numfile, int dim, std::vector<WALLS> &_pst);

int PostProcGrain(int numfile, int dim);

int PostProcGrainVisual(int dim, std::vector <std::vector <float> > &data);

int PostProcForce(int numfile, int dim);

int PostProcBond(int numfile, int dim);

int PostProcWall(int numfile, int dim, std::vector<WALLS> &wl);

int PostProcessing(bool _activeCohesion, int dim, std::vector<WALLS> &_pst, std::vector<WALLS> &wl);

int RunPostProcessing(bool _activeCohesion, int dim);

int RunVisualPostProcessing(bool _activeCohesion, int dim);

int GetContactsList(int dim, std::vector <Grain> &gr, std::vector<std::vector<bond>> &b, bool _activeCohesion);

void RuntimeLog(std::string fileName, bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose, bool showValue = false);

#endif