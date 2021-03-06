#ifndef ATOM_HPP
#define ATOM_HPP

#include <random>
#include <vector>
#include <simulbody/simulator.hpp>

#include "elements.hpp"

using namespace simulbody;

class Atom {
protected:
	System* system;

	Element electronConfiguration;
	Element nucleusElement;

	double reducedMass;
	double nucleusMass;
	double nucleusCharge;

	identifier nucleus;
	std::vector<std::string> orbitNames;
	std::map<std::string, identifier> electrons;
	std::vector<Interaction*> interactions;

public:

	Atom(System* system, Element element, double atomicMass);
	Atom(System* system, Element electronConfiguration, Element nucleusElement, double atomicMass);

	identifier getNucleus() const;
	std::vector<identifier> getElectrons() const;
	std::vector<identifier> getBodies() const;
	identifier getElectron(std::size_t orbit) const;
	identifier getElectron(std::string orbitName) const;
	std::vector<Interaction*> getInteractions() const;

	vector3D getPosition() const;
	vector3D getVelocity() const;
	vector3D getImpulse() const;

	double getMass() const;
	double getReducedMass() const;
	double getNucleusMass() const;
	double getNucleusCharge() const;

	Element getNucleusElement() const;
	Element getElectronConfiguration() const;
	std::vector<string> getOrbitNames() const;

	void setPosition(vector3D position);
	void setVelocity(vector3D velocity);

	virtual void install() = 0;
	virtual void randomize(std::mt19937_64 &randomEngine) = 0;
	virtual void createInteractions() = 0;

	virtual double getEnergy() const;
	virtual double getIonizationEnergy(std::string orbit) const;
	virtual double getOrbitalEnergy(std::string orbit) const;
	virtual vector3D getOrbitalAngularMomentum(std::string orbit) const;

	virtual ~Atom();

	static constexpr double electronMass = 1.0;
	static constexpr double protonMass = 1836.1527;
};

#endif /* ATOM_HPP */
