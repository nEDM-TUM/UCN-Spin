#include <iostream>
#include "../polynom.h"

using namespace std;

int main() {
	Polynom p1(1,2,3);
	Polynom p2(1,2);

	cout << "(" << p1.toString() << ") * (" << p2.toString() << ") = " << (p1*p2).toString() << endl;

	return 0;
}
