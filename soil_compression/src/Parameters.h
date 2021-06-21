#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include <string>
#include "Structures.h"
#include <fstream>

double Grav;

// Properties of grains: all grains posses the same properties
double kn; // Normal stiffness
double kt; // Tangential stiffness
double b_kn;  // Tensile
double b_kt; //b_kn;
double b_M;
double b_YM;
bool deloneRapture;
double sigma;
double tau;
double rho; // Density of the grains
double mu; // Sliding friction
double mur; // Rolling friction
double nur; // Rolling "damping"
double Adh;// Adhesion
size_t grainNum;
double beta; // Damping ratio

int dim;
bool isCohesion;
double _rMin;
double _rMax;
double _sizeX;
double _sizeY;
double _nodesRes;
double _timeInterval;
double sampleStretch;
double porosity;
double lRate;
double pSize;
double runTime, loadTime, samplTime;
double stopMechRatio;

double stdDev; // standart deviation in Gauss distribution
double shape; // alpha - shape parameter in Weibull distribution
double scale; // betha - scale parameter in Weibull distribution

int relaxationFlag;
uint uint_recordOutput;
bool recordOutput;

std::string type;
std::string projDirMain;

size_t snapshotCounter;
size_t snapshotCounterP;
size_t snapshotCounterG;
size_t snapshotCounterF;
size_t snapshotCounterB;
size_t snapshotCounterW;

//-----------------------------------------------------------------------------------------------
Wall wBottom, wTop, wLeft, wRight, wFront, wBack;
WALLS allWalls;
Piston wLoadPiston;
std::vector<std::vector<bond> > b, _b;
std::vector <Grain> grain;
std::vector <Grain> _grain;
std::vector<std::vector<std::vector<int> > > cells;
//----------------------------------------------------------------------------------------------
std::vector<WALLS> wl;
int wlCount;
bool wlTop;
bool wlSmf;
std::vector<Piston> pst;
std::vector<WALLS> _pst; /* Associate the piston to the top wall. This wall is going to be reduced to the piston size if LOADING. */
int pstCount;
bool pstTop;
bool pstSmf;
std::vector<std::vector <Grain>> grainSnapshot;
int grnCount;
bool grnTop;
bool grnSmf;
std::vector<std::vector <Grain>> grainSnapshotF;
int fCount;
bool fTop;
bool fSmf;
std::vector<std::vector <Grain>> grainSnapshotB;
int bCount; // access index in snapshot container
bool bTop; //flag (some kind of semaphore), if TRUE - nullify bCount
bool bSmf; //mutex, blocking access to snapshot container during its modification
std::vector<std::vector<std::vector<bond> >> bondSnapshot;
std::vector<shortBond> allBonds;
std::vector<std::vector<shortBond>> allBondsSnapshot;
int snapshotBuffer;

std::fstream runtimelogStr;
std::string runTimeFile("DEM_RUNTIME.txt");
size_t logFilePos_bgn;
size_t logFilePos_end;

#endif
