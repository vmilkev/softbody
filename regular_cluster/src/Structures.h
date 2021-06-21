#ifndef STRUCTURES_H
#define STRUCTURES_H

struct shortBond
{
	int i;
	int j;
	double f;
};

struct SmallGrain
{
	double x, y, z; // position
	double m; // mass
	double R; // radius
	double I; // inertia
	double rhoInd; //individual density
	size_t name;
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

#endif
