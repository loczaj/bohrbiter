#ifndef COHEN_HPP
#define COHEN_HPP

#include <map>
#include <tuple>

using namespace std;

struct CohenConfiguration {

	typedef tuple<double, double, double, double, double, double, bool> CohenOrbit;

	map<string, CohenOrbit> orbits;

	CohenConfiguration() {
		orbits["1:1s1"] = CohenOrbit(1.0000, 0.0000, 0.0000, 0.9535, 0.0000, 0.0000, true);

		orbits["2:1s1"] = CohenOrbit(0.5714, 0.0000, 0.0000, 1.6686, 0.0000, 0.0000, true);
		orbits["2:1s2"] = CohenOrbit(0.5714, 3.1416, 0.0000, 1.6686, 3.1416, 0.0000, false);
	}

	vector3D position(Element const element, string orbit) {
		return vector3D();
	}

	vector3D velocity(Element const element, string orbit) {
		return vector3D();
	}

	bool spin(Element const element, string orbit) {
		return true;
	}

};

#endif /* COHEN_HPP */
