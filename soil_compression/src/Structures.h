#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>

struct bond
{
	bool isActive;
	double alpha0;
	double alpha0XY;
	double alpha0YZ;
	double kn; // Tensile
	double kt;
	double M;
	double Yn; // Yield thresholds
	double Yt;
	double YM;
	double xI, yI, zI;  // coordinates of contact point
	double sigma;
	double tau;
	double frictionAngle;
	double cohesion;
	double gap;
	double contactForce;
};

struct shortBond
{
	int i;
	int j;
	double f;
};

struct Grain
{
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
	int grainId;
	bool isActive;
	std::vector <size_t> contactList;
};

struct SmallGrain
{
	double x, y, z; // position
	double m; // mass
	double R; // radius
	double I; // inertia
	double rhoInd; //individual density
	int grainId;
};

struct Wall
{
	double wF;
	double wV;
	double wA;
	double lowX, lowY, lowZ, highX, highY, highZ;
	double wM;
};

struct WALLS
{
	Wall _wBottom, _wTop, _wLeft, _wRight, _wFront, _wBack;
};

struct Piston
{
	double lowX, lowY, lowZ, highX, highY, highZ;
	int movementFlag;
	double loadTime;
	double loadRate;
	double sizeR;
};

#endif
