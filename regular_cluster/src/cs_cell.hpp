#ifndef CELL_HPP
#define CELL_HPP

#include "allincludes.h"
#include "glmath.hpp"

class Cell
{
public:

    /* oled data */
	double m; // mass
	double rhoInd; //individual density
	double R; // radius
	double I; // inertia
	double x, y, z, th; // position
	double Vx, Vy, Vz, Vth; // velocity
	double Ax, Ay, Az, Ath; // acceleration
	double Fx, Fy, Fz, Fth; // sum of forces acting on the grain
	double p; // pressure on grain
	double normOfSumF;
	double sumOfNormF;
	size_t name;
	bool isActive;
	std::vector <size_t> contactList;

	/* new stuff */
	glm::dvec3 origin;
	glm::dvec3 growthdir;
	glm::dvec3 divisiondir;
	bool isElongating;
	double growthRate;
	bool isLinked;
	size_t linkedID;
	double conf_degree;
	double conf_degree_upd;
	bool isVirtual;
	
	size_t conf_num;
	double age;
	double currentElongation;
	double pastElongation;

private:

	//double age;
	//size_t conf_num;
	std::vector < glm::dvec3 > confinements;
	glm::dvec3 origin_nov;
	double growthRate_nov;
	//double currentElongation;

	double initElongationTime = 5.0;
	double growthRateConstant = 0.0005;
	double maxElongation;
	size_t maxID = 1000000;
	double growthThreshold = 0.9;
	double confLimit = 1.0;

public:

    Cell( glm::dvec3 orig, double radius );
    ~Cell();

	Cell addcell();
	void growth( std::vector<Cell> &cells, double time );
	void correct_origin();

private:

	glm::dvec3 getNormale( glm::dvec3 p0, glm::dvec3 p1, glm::dvec3 p2 );
	glm::dvec3 getCross( glm::dvec3 p1, glm::dvec3 p2 );
	glm::dvec3 getOrthVect( glm::dvec3 p );
	void getConfinements( std::vector<Cell> &cells );
	double updConfDegree( std::vector<Cell> &cells );
	void getGrowthDirections( std::vector<Cell> &cells );
	glm::dvec3 getNewOrigin( glm::dvec3 currentorigin, glm::dvec3 growthdirection, double cellsshift );
};

#endif