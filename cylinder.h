#ifndef CYLINDER_H
#define CYLINDER_H

class Cylinder
{
	public:
		Cylinder(double radius, double height);
		
		bool contains(double x[]);
		bool insideHeight(double x[]);
		bool insideRadius(double x[]);

		bool reflect(double v[], double x[], bool state[2]);

	private:
		void reflectHeight(double v[]);
		void reflectRadius(double v[], double x[]);

		double fRadius; ///< radius of cylinder
		double fRSquared; ///< squared radius of cylinder
		double fHeight; ///< height of cylinder
};

#endif // CYLINDER_H
