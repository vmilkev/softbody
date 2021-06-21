#ifndef BOND_HPP
#define BOND_HPP

class Bond
{
public:

    /* data */
	bool isActive;
	bool rigid;
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

    /* methods */

public:
    Bond();
    ~Bond();
};

#endif