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

template<class T>
vector<T> append(const vector<T> &a, const vector<T> &b) {
	vector<T> sum(a);
	sum.insert(sum.end(), b.begin(), b.end());
	return sum;
}

struct PeriodicTable {
	map<Element, vector<string>> orbitals;

	static int atomicNumber(Element const element) {
		return static_cast<int>(element);
	}

	vector<string> atomicOrbitals(Element const element) {
		return orbitals.at(element);
	}

	PeriodicTable() {
		orbitals[Element::H] = {"1s1"};
		orbitals[Element::He] = {"1s1", "1s2"};

		orbitals[Element::Li] = {"1s1", "1s2", "2s1"};
		orbitals[Element::Be] = {"1s1", "1s2", "2s1", "2s2"};

		auto Be = orbitals.at(Element::Be);
		orbitals[Element::B] = append(Be, {"2p1"});
		orbitals[Element::C] = append(Be, {"2p1", "2p2"});
		orbitals[Element::N] = append(Be, {"2p1", "2p2", "2p3"});
		orbitals[Element::O] = append(Be, {"2p1", "2p2", "2p3", "2p4"});
		orbitals[Element::F] = append(Be, {"2p1", "2p2", "2p3", "2p4", "2p5"});
		orbitals[Element::Ne] = append(Be, {"2p1", "2p2", "2p3", "2p4", "2p5", "2p6"});

		auto Ne = orbitals.at(Element::Ne);
		orbitals[Element::Na] = append(Ne, {"3s1"});
		orbitals[Element::Mg] = append(Ne, {"3s1", "3s2"});

		auto Mg = orbitals.at(Element::Mg);
		orbitals[Element::Al] = append(Mg, {"3p1"});
		orbitals[Element::Si] = append(Mg, {"3p1", "3p2"});
		orbitals[Element::P] = append(Mg, {"3p1", "3p2", "3p3"});
		orbitals[Element::S] = append(Mg, {"3p1", "3p2", "3p3", "3p4"});
		orbitals[Element::Cl] = append(Mg, {"3p1", "3p2", "3p3", "3p4", "3p5"});
		orbitals[Element::Ar] = append(Mg, {"3p1", "3p2", "3p3", "3p4", "3p5", "3p6"});
	}

	static double nucleusMassInAU(Element const element, double atomicMass) {
		return atomicMass / 5.485799095 * 10000 - 1.0 * atomicNumber(element);
	}
};

#endif /* ORBITALS_HPP */
