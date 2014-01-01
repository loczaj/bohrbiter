#ifndef COHEN_HPP
#define COHEN_HPP

#include <map>
#include <tuple>

#include "elements.hpp"

using namespace std;

struct CohenConfiguration {

	typedef tuple<double, double, double, double, double, double, bool> CohenOrbit;

	map<string, CohenOrbit> orbits;

	CohenConfiguration() {
		orbits["1:1s1"] = CohenOrbit(1.0000, 0.00, 0.0000, 0.9535, 0.00, 0.0000, true);

		orbits["2:1s1"] = CohenOrbit(0.5714, 0.00, 0.0000, 1.6686, 0.00, 0.0000, true);
		orbits["2:1s2"] = CohenOrbit(0.5714, M_PI, 0.0000, 1.6686, M_PI, 0.0000, false);
	}

	vector3D position(const Element &element, const string &orbit) {
		CohenOrbit cohenOrbit = orbits.at(key(element, orbit));
		vector3D position(get<0>(cohenOrbit), get<1>(cohenOrbit), get<2>(cohenOrbit));
		return position.convertFromSphericalToCartesian();
	}

	vector3D momentum(const Element &element, const string &orbit) {
		CohenOrbit cohenOrbit = orbits.at(key(element, orbit));
		vector3D momentum(get<3>(cohenOrbit), get<4>(cohenOrbit), get<5>(cohenOrbit));
		return momentum.convertFromSphericalToCartesian();
	}

	bool spin(const Element &element, const string &orbit) {
		return get<6>(orbits.at(key(element, orbit)));
	}

	string key(const Element &element, const string &orbit) {
		return to_string(PeriodicTable::atomicNumber(element)) + ":" + orbit;
	}

};

#endif /* COHEN_HPP */
