// Simulation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;

// Compute compressed radius for object 1

double getC1(double c2, double d)
{
	return d - c2;
}

// Compute compressed radius for object 2

double getC2(double r1, double r2, double d)
{
	return (d / 2) - (pow(r1, 2) - pow(r2, 2)) / (2 * d);
}

// ToDo: Merge getK1 and getK2

double getK1(double c1, double r1)
{
	return c1 / r1;
}

double getK2(double c2, double r2)
{
	return c2 / r2;
}

// ToDo: Merge getF1 and getF2

double getF1(double F, double k1)
{
	return F * k1;
}

double getF2(double F, double k2)
{
	return F * k2;
}

// ToDo: Merge getFx and getFy

double getFx(double x1, double x2, double fTotal, double d)
{
	return ((x2 - x1) * fTotal) / d;
}

double getFy(double y1, double y2, double fTotal, double d)
{
	return ((y2 - y1) * fTotal) / d;
}

// ToDo: Merge getCMx and getCMy

double getCMx(double m1, double m2, double x1, double x2)
{
	return (x1*m1 + x2*m2) / (m1 + m2);
}

double getCMy(double m1, double m2, double y1, double y2)
{
	return (y1*m1 + y2*m2) / (m1 + m2);
}

// Compute the acceleration

double getAccel(double x, double CM, double F, double m)
{
	if (x < CM)
	{
		return (F / m) * -1;
	}
	else
	{
		return (F / m);
	}
}

// Compute the velocity

double getVel(double v0, double a, double t)
{
	return v0 + a*t;
}

// Compute the position

double getPos(double x0, double v, double t)
{
	return x0 + v*t;
}

// Compute the straight-line distance

double getDist(double x1, double x2, double y1, double y2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

// Main method

int main() {

	double m1 = 1;				// Mass of object 1 (kg)
	double m2 = 1;				// Mass of object 2 (kg)

	double r1 = 0.5;				// Radius of object 1 (m)
	double r2 = 0.5;				// Radius of object 2 (m)

	double rTotal = r1 + r2;	// Combined radii of objects 1 and 2 (m)

	double x1 = -0.4;			// Object 1, x position (m)
	double y1 = -0.3;			// Object 1, y position (m)
	double x2 = 2.4;			// Object 2, x position (m)
	double y2 = 0.6;			// Object 2, y position (m)

	double vx1 = 1.2;			// Object 1, velocity in x (ms-1)
	double vy1 = 0.6;			// Object 1, velocity in y (ms-1)
	double vx2 = -1.5;			// Object 2, velocity in x (ms-1)
	double vy2 = 0.3;			// Object 2, velocity in y (ms-1)

	double F = 1960;			// Experimentally found

	double d = getDist(x1, x2, y1, y2);	// The objects' distance from one another

	double deltaT = 0.0001;	// The time intervals (should this be user - defined?)

	// ----------------------------------------
	// Everything after here is to be computed
	// ----------------------------------------

	double ax1;					// Object 1, acceleration in x (ms-2)
	double ay1;					// Object 1, acceleration in y (ms-2)

	double ax2;					// Object 2, acceleration in x (ms-2)
	double ay2;					// Object 2, acceleration in y (ms-2)

	double c1;					// Object 1, compressed radius (m)
	double c2;					// Object 2, compressed radius (m)

	double k1;					// Object 1, ratio of compression 
	double k2;					// Object 2, ratio of compression

	double fTotal;				// Net compression force within the system

	double fx;					// The force that is applied in x (N)
	double fy;					// The force that is applied in y (N)

	double CMx;					// The center of mass of the system in x
	double CMy;					// The center of mass of the system in y

	double t = 0;				// The time variable, essentially just a counter

	bool haveCollided = false;	// Have the balls collided yet?
	bool backApart = false;		// Have the balls bounced back apart yet?

	// toDo: Have user input and also options to solve for certain things (i.e. only the final velocities, only the final positions after time, etc.)

	while (!backApart)
	{

		if (d > rTotal)				// NonContact case
		{

			if (haveCollided)		// If they have already collided
			{
				backApart = true;	// They are therefore now back apart
			}

			else
			{

				x1 = getPos(x1, vx1, deltaT);	// Compute the position of object 1 in x
				y1 = getPos(y1, vy1, deltaT);	// Compute the position of object 1 in y

				x2 = getPos(x2, vx2, deltaT);	// Compute the position of object 2 in x
				y2 = getPos(y2, vy2, deltaT);	// Compute the position of object 2 in y

				d = getDist(x1, x2, y1, y2);	// Compute the straight-line distance between the centers of the two objects

				t += deltaT;					// Compute the current time 

				cout << "At time t = " << t << "s..." << endl;
				cout << "x1 = " << x1 << endl;
				cout << "x2 = " << x2 << endl;
				cout << "vx1 = " << vx1 << endl;
				cout << "vy1 = " << vy1 << endl;
				cout << "vx2 = " << vx2 << endl;
				cout << "vy2 = " << vy2 << endl;
				cout << "d = " << d << endl;				
			}
		}

		else // Contact case
		{

			cout << "BOOM" << endl;					// Debug statement, will remove

			haveCollided = true;					// We now know that the two balls have collided

			c2 = getC2(r1, r2, d);					// Compute the compressed radii of object 1
			c1 = getC1(c2, d);						// Compute the compressed radii of object 2

			k1 = getK1(c1, r1);						// Compute the ratio of compression for object 1
			k2 = getK2(c2, r2);						// Compute the ratio of compression for object 2

			fTotal = getF1(F, k1) + getF2(F, k2);	// Compute the system's net force

			fx = getFx(x1, x2, fTotal, d);			// Compute the force applied in x
			fy = getFy(y1, y2, fTotal, d);			// Compute the force applied in y

			CMx = getCMx(m1, m2, x1, x2);			// Compute the center of mass in x
			CMy = getCMy(m1, m2, y1, y2);			// Compute the center of mass in y

			ax1 = getAccel(x1, CMx, fx, m1);		// Compute the acceleration of object 1 in x
			ay1 = getAccel(y1, CMy, fy, m1);		// Compute the acceleration of object 1 in y

			ax2 = getAccel(x2, CMx, fx, m2);		// Compute the acceleration of object 2 in x
			ay2 = getAccel(y2, CMy, fy, m2);		// Compute the acceleration of object 2 in y

			vx1 = getVel(vx1, ax1, deltaT);			// Compute the velocity of object 1 in x
			vy1 = getVel(vy1, ay1, deltaT);			// Compute the velocity of object 1 in y

			vx2 = getVel(vx2, ax2, deltaT);			// Compute the velocity of object 2 in x
			vy2 = getVel(vy2, ay2, deltaT);			// Compute the velocity of object 2 in y

			x1 = getPos(x1, vx1, deltaT);			// Compute the position of object 1 in x
			y1 = getPos(y1, vy1, deltaT);			// Compute the position of object 1 in y

			x2 = getPos(x2, vx2, deltaT);			// Compute the position of object 2 in x
			y2 = getPos(y2, vy2, deltaT);			// Compute the position of object 2 in y

			d = getDist(x1, x2, y1, y2);			// Compute the distance between the centers of the two objects

			t += deltaT;							// Compute the current time

			// Output results

			cout << "At time t = " << t << "s..." << endl;
			cout << "x1 = " << x1 << endl;
			cout << "x2 = " << x2 << endl;
			cout << "vx1 = " << vx1 << endl;
			cout << "vy1 = " << vy1 << endl;
			cout << "vx2 = " << vx2 << endl;
			cout << "vy2 = " << vy2 << endl;
			cout << "d = " << d << endl;

		}
	}

	system("pause");

}
