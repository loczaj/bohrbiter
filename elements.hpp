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

inline int atomicNumber(Element const value) {
	return static_cast<int>(value);
}

template<class T>
vector<T> append(const vector<T> &a, const vector<T> &b) {
	vector<T> sum(a);
	sum.insert(sum.end(), b.begin(), b.end());
	return sum;
}

struct orbitals: public map<Element, vector<string>> {

	orbitals() {
		operator[](Element::H) = {"1s1"};
		operator[](Element::He) = {"1s1", "1s2"};

		operator[](Element::Li) = {"1s1", "1s2", "2s1"};
		operator[](Element::Be) = {"1s1", "1s2", "2s1", "2s2"};

		operator[](Element::B) = append(at(Element::Be), {"2p1"});
		operator[](Element::C) = append(at(Element::Be), {"2p1", "2p2"});
		operator[](Element::N) = append(at(Element::Be), {"2p1", "2p2", "2p3"});
		operator[](Element::O) = append(at(Element::Be), {"2p1", "2p2", "2p3", "2p4"});
		operator[](Element::F) = append(at(Element::Be), {"2p1", "2p2", "2p3", "2p4", "2p5"});
		operator[](Element::Ne) = append(at(Element::Be), {"2p1", "2p2", "2p3", "2p4", "2p5", "2p6"});

		operator[](Element::Na) = append(at(Element::Ne), {"3s1"});
		operator[](Element::Mg) = append(at(Element::Ne), {"3s1", "3s2"});

		operator[](Element::Al) = append(at(Element::Mg), {"3p1"});
		operator[](Element::Si) = append(at(Element::Mg), {"3p1", "3p2"});
		operator[](Element::P) = append(at(Element::Mg), {"3p1", "3p2", "3p3"});
		operator[](Element::S) = append(at(Element::Mg), {"3p1", "3p2", "3p3", "3p4"});
		operator[](Element::Cl) = append(at(Element::Mg), {"3p1", "3p2", "3p3", "3p4", "3p5"});
		operator[](Element::Ar) = append(at(Element::Mg), {"3p1", "3p2", "3p3", "3p4", "3p5", "3p6"});
	}
};

#endif /* ORBITALS_HPP */
