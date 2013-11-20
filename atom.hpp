#ifndef ATOM_HPP
#define ATOM_HPP

#include <vector>
#include <simulbody/simulator.hpp>

#include "elements.hpp"

using namespace simulbody;

class Atom {
	identifier nucleus;
	std::map<std::string, identifier> electrons;

	std::vector<Interaction> interactions;

public:
	Atom(Element element, unsigned int massNumber);
	Atom(Element element, unsigned int massNumber, unsigned int electronNumber);

	virtual void moveTo(System &system, vector3D position);
	virtual void speedUp(System &system, vector3D velocity);

	virtual void install(System &system) = 0;
	virtual void setInteractions() = 0;

	virtual double getOrbitalEnergy(System &system, identifier electron) const;
	virtual double getAngularMomentum(System &system, identifier electron) const;

	virtual ~Atom();
};

#endif /* ATOM_HPP */
