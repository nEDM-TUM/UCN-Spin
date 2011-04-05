#include <iostream>
#include <cmath>
#include "../roots.h"
#include "../polynom.h"

using namespace std;

int main() {
	try {
		clog << scientific;
		clog.precision(15);

		double left = -10;
		double right = 10;

		Polynom f(1., 0, -2., 2.);
		Polynom df = f.derivative();

		clog << "Function values:" << endl;
		const char *tab = "\t\t\t";
		cout << "#x" << tab << f.toString() << tab << df.toString() << endl;
		cout.precision(15);
		for (double x = left; x <= right; x += 0.1)
			cout << x << tab << f(x) << tab << df(x) << endl;
		clog << endl << endl;

		double result = Roots::safeNewton(f, df, left, right, 1e-14);

		clog << endl << "Result:" << endl << endl;
		clog << "f(" << result << ") = " << f(result) << endl;
	}
	catch (char const *s) {
		cerr << endl << "ERROR: " << s << endl;
	}
}
