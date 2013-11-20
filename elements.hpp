#ifndef ORBITALS_HPP
#define ORBITALS_HPP

#include <map>
#include <string>
#include <vector>

using namespace std;

// TODO finish declaring elements over z=18
enum class Element
	: int {
		H = 1, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar
};

int atomicNumber(Element const value) {
	return static_cast<unsigned int>(value);
}

struct orbitals: public map<Element, vector<string>> {

	orbitals() {
		at(Element::H) = {"1s1"};
		at(Element::He) = {"1s1", "1s2"};

		at(Element::Li) = {"1s1", "1s2", "2s1"};
		at(Element::Be) = {"1s1", "1s2", "2s1", "2s2"};

		at(Element::B) = {"1s1", "1s2", "2s1", "2s2", "2p1"};
		at(Element::C) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2"};
		at(Element::N) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3"};
		at(Element::O) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4"};
		at(Element::F) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5"};
		at(Element::Ne) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6"};

		at(Element::Na) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1"};
		at(Element::Mg) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2"};

		at(Element::Al) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1"};
		at(Element::Si) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2"};
		at(Element::P) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3"};
		at(Element::S) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3", "3p4"};
		at(Element::Cl) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3", "3p4", "3p5"};
		at(Element::Ar) = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3", "3p4", "3p5", "3p6"};
	}
};

#endif /* ORBITALS_HPP */
