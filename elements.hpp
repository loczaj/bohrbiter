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

struct orbitals {
	map<Element, vector<string>> names;

	orbitals() {
		names[Element::H] = {"1s1"};
		names[Element::He] = {"1s1", "1s2"};

		names[Element::Li] = {"1s1", "1s2", "2s1"};
		names[Element::Be] = {"1s1", "1s2", "2s1", "2s2"};

		names[Element::B] = {"1s1", "1s2", "2s1", "2s2", "2p1"};
		names[Element::C] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2"};
		names[Element::N] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3"};
		names[Element::O] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4"};
		names[Element::F] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5"};
		names[Element::Ne] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6"};

		names[Element::Na] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1"};
		names[Element::Mg] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2"};

		names[Element::Al] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1"};
		names[Element::Si] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2"};
		names[Element::P] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3"};
		names[Element::S] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3", "3p4"};
		names[Element::Cl] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3", "3p4", "3p5"};
		names[Element::Ar] = {"1s1", "1s2", "2s1", "2s2", "2p1", "2p2", "2p3", "2p4", "2p5", "2p6", "3s1", "3s2", "3p1", "3p2", "3p3", "3p4", "3p5", "3p6"};
	}
};

#endif /* ORBITALS_HPP */
