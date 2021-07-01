#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include "allincludes.h"
#include "Structures.h"
#include "cs_cell.hpp"
#include "cs_bond.hpp"

class Sample
{
private:

    /* data */
    double Grav;
    double kn; // Normal stiffness
    double kt; // Tangential stiffness
    double b_kn;  // Tensile
    double b_kt; //b_kn;
    bool deloneRapture;
    double sigma;
    double tau;
    double rho; // Density of the grains
    double rho_env; // Density of the environment
    double mu; // Sliding friction
    double mur; // Rolling friction
    double nur; // Rolling "damping"
    double Adh;// Adhesion
    size_t grainNum;
    double beta; // Damping ratio
    double _rMin;
    double _rMax;
    double _sizeX;
    double _sizeY;
    double stdDev; // standart deviation in Gauss distribution
    double shape; // alpha - shape parameter in Weibull distribution
    double scale; // betha - scale parameter in Weibull distribution
    int relaxationFlag;
    uint uint_recordOutput;
    bool recordOutput;
    std::vector<std::vector<Bond> > b;
    std::vector <Cell> aggregates;
    std::vector<WALLS> wl;
    std::string projDirMain;
    size_t snapshotCounter;
    size_t snapshotCounterG;
    size_t snapshotCounterF;
    size_t snapshotCounterB;
    size_t snapshotCounterW;
    Wall wBottom, wTop, wLeft, wRight, wFront, wBack;
    WALLS allWalls;
    size_t wlCount;
    bool wlTop;
    bool wlSmf;
    std::vector<std::vector <Cell>> grainSnapshot;
    size_t grnCount;
    bool grnTop;
    bool grnSmf;
    std::vector<std::vector <Cell>> grainSnapshotF;
    size_t fCount;
    bool fTop;
    bool fSmf;
    std::vector<std::vector <Cell>> grainSnapshotB;
    size_t bCount; // access index in snapshot container
    bool bTop; //flag (some kind of semaphore), if TRUE - nullify bCount
    bool bSmf; //mutex, blocking access to snapshot container during its modification
    std::vector<std::vector<std::vector<Bond> >> bondSnapshot;
    std::vector<shortBond> allBonds;
    std::vector<std::vector<shortBond>> allBondsSnapshot;
    int snapshotBuffer;
    std::fstream runtimelogStr;
    const char*  runTimeFile;
    //const char*  paramFile;
    //static constexpr const char* runTimeFile = "DEM_RUNTIME2.txt";
    size_t logFilePos_bgn;
    size_t logFilePos_end;

public:

    std::string type;
    double runTime;
    double samplTime;

public:

    Sample();
    ~Sample();

    int Init(const char* runtimeFileName, const char* paramFile);
    int SetSample( int cellsNum );
    void ReportExit(int exitValue);
    bool CheckProjType();
    int RunSimulation(double time);
    int Finalize( std::string fileNameOutput );
    int RestoreSample(const char* paramFileName);
    int RunGrowthSimulation(double grwtime, double simtime);

private:

    int ReadParametersFromFile(const char* paramFileName);
    int ProjectManager();
    bool DirExists(const std::string& dirName_in);
    std::string ws2s(const std::wstring& wstr);
    int InitializeSimulationProject(const char* paramFileName, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, std::vector<WALLS> &wl);
    double GetSampleHeight(double sizeY, std::vector <Cell> &gr);
    void SetWalls(double sizeX, double sizeY);
    int ForcesGrain3D(int i, int j, std::vector <Cell> &gr, double _dt);
    void ForcesWall3D(int i, int wall, std::vector <Cell> &gr);
    void ForcesGrainCohesive3D(int i, int j, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, double _dt);
    void IdentifyCohesiveBonds(int i, int j, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, double *randValues);
    void SetBondsMatrix(std::vector<std::vector<Bond>> &b, std::vector <Cell> &gr);
    void FindBonds(std::vector<std::vector<Bond>> &b, std::vector <Cell> &gr, double shape, double scale);
    void IsBallInDomain(std::vector <Cell> &gr);
    void ClearForces(std::vector <Cell> &gr);
    double Get_dt();
    int DynamicRelaxation(double timeInterval, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, bool record, std::string prType);
    std::vector<shortBond> GetSnapshotB(std::vector<std::vector<Bond>> &b1);
    int GetSnapshots(std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b1);
    int GetCellsSnapshots(std::vector <Cell> &gr);
    int PostProcGrain(int numfile);
    int PostProcGrainVisual(std::vector <std::vector <float> > &data);
    int PostProcForce(int numfile);
    int PostProcBond(int numfile);
    int PostProcWall(int numfile, std::vector<WALLS> &wl);
    int FinalizeSampleProject(std::string paramFileName, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, std::vector<WALLS> &wl);
    int PostProcessing();
    int RunPostProcessing();
    int RunVisualPostProcessing();
    int GetContactsList(std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b);
    void RuntimeLog(std::string fileName, bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose, bool showValue = false);
    void SeedCells(double rMin, double rMax, double sizeD, double sizeH, double cellsNum, std::vector <Cell> &gr, double stdDev);
    int SimulateGrowth( double timeInterval, std::vector <Cell> &gr, std::vector<std::vector<Bond>> &b, bool record, std::string prType );
    int SimulateGrowth2( double grwtime, double timeInterval, std::vector<Cell> &cells, std::vector<std::vector<Bond>> &b, bool record, std::string prType );

};

#endif