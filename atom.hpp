#ifndef ATOM_HPP
#define ATOM_HPP

#include <vector>
#include <simulbody/simulator.hpp>

#include "elements.hpp"

using namespace simulbody;

class Atom {
protected:
	Element element;
	Element electronConfiguration;
	unsigned int massNumber;

	double nucleusMass;
	double nucleusCharge;

	identifier nucleus;
	std::map<std::string, identifier> electrons;
	std::vector<Interaction> interactions;

public:
	Atom(System &system, Element element, unsigned int massNumber);
	Atom(System &system, Element element, unsigned int massNumber, Element electronConfiguration);

	identifier getNucleus();
	std::vector<identifier> getElectrons();
	identifier getElectron(std::size_t orbit);
	identifier getElectron(std::string orbitName);

	virtual void moveTo(System &system, vector3D position);
	virtual void speedUp(System &system, vector3D velocity);

	virtual void install(System &system) = 0;
	virtual void setInteractions() = 0;

	virtual double getOrbitalEnergy(System &system, identifier electron) const;
	virtual double getAngularMomentum(System &system, identifier electron) const;

	virtual ~Atom();

	static constexpr double electronMass = 1.0;
	static orbitals orbits;
};

#endif /* ATOM_HPP */
