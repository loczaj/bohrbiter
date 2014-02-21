#ifndef ABRINES_PERCIVAL_HPP
#define ABRINES_PERCIVAL_HPP

#include <random>
#include <simulbody/simulator.hpp>

#include "atom.hpp"

using namespace simulbody;

class AbrinesPercivalAtom: public Atom {

public:
	AbrinesPercivalAtom(System* system, Element element, double atomicMass);
	AbrinesPercivalAtom(System* system, Element electronConfig, Element nucleusElement, double atomicMass);

	virtual void install() override;
	virtual void randomize(std::mt19937_64 &randomEngine) override;
	virtual void createInteractions() override;

private:
	double solveKeplerEquation(double thetaN, double epsilon, double tolerance,
			std::mt19937_64 &randomEngine);
};

#endif /* ABRINES_PERCIVAL_HPP */
