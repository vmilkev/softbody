#ifndef EXTERNALDECLARATIONS_H
#define EXTERNALDECLARATIONS_H

#include <vector>

extern double Grav;

extern double kn; // Normal stiffness
extern double kt; // Tangential stiffness

extern double b_kn;  // Tensile
extern double b_kt; //b_kn;
extern double b_M;
extern double b_YM;

extern bool deloneRapture;

extern double sigma;
extern double tau;

extern double rho; // Density of the grains
extern double mu; // Sliding friction

extern double mur; // Rolling friction
extern double nur; // Rolling "damping"

extern double Adh;// Adhesion

extern size_t grainNum;

extern double _nodesRes;
extern double sampleStretch;

extern double beta; // Damping ratio

extern int relaxationFlag;
extern double _sizeY;

extern size_t snapshotCounter;
extern size_t snapshotCounterP;
extern size_t snapshotCounterG;
extern size_t snapshotCounterF;
extern size_t snapshotCounterB;
extern size_t snapshotCounterW;

extern struct Wall wBottom, wTop, wLeft, wRight, wFront, wBack;
extern struct WALLS allWalls;
extern struct Piston wLoadPiston;
extern std::vector<std::vector< struct bond> > b, _b;
extern std::vector < struct Grain> grain;
extern std::vector < struct Grain> _grain;
extern std::vector<std::vector<std::vector<int> > > cells;

extern std::vector<WALLS> wl;
extern int wlCount;
extern bool wlTop;
extern bool wlSmf;
extern std::vector<Piston> pst;
extern std::vector<WALLS> _pst;
extern int pstCount;
extern bool pstTop;
extern bool pstSmf;
extern std::vector<std::vector <Grain>> grainSnapshot;
extern int grnCount;
extern bool grnTop;
extern bool grnSmf;
extern std::vector<std::vector <Grain>> grainSnapshotF;
extern int fCount;
extern bool fTop;
extern bool fSmf;
extern std::vector<std::vector <Grain>> grainSnapshotB;
extern int bCount;
extern bool bTop;
extern bool bSmf;
extern std::vector<std::vector<std::vector<bond> >> bondSnapshot;
extern std::vector<shortBond> allBonds;
extern std::vector<std::vector<shortBond>> allBondsSnapshot;
extern int snapshotBuffer;

extern std::fstream runtimelogStr;
extern std::string runTimeFile;
extern size_t logFilePos_bgn;
extern size_t logFilePos_end;

extern std::string type;
extern std::string projDirMain;

#endif