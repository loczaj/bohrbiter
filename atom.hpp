#ifndef ATOM_HPP
#define ATOM_HPP

#include <random>
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
	std::vector<Interaction*> interactions;

public:
	std::vector<std::string> orbitNames;

	Atom(System &system, Element element, unsigned int massNumber);
	Atom(System &system, Element element, unsigned int massNumber, Element electronConfiguration);

	identifier getNucleus() const;
	std::vector<identifier> getElectrons() const;
	std::vector<identifier> getBodies() const;
	identifier getElectron(std::size_t orbit) const;
	identifier getElectron(std::string orbitName) const;

	vector3D getPosition(System &system) const;
	vector3D getVelocity(System &system) const;
	vector3D getImpulse(System &system) const;
	double getMass(System &system) const;

	void setPosition(System &system, vector3D position);
	void setVelocity(System &system, vector3D velocity);

	virtual void install(System &system) = 0;
	virtual void randomize(System &system, std::mt19937_64 &randomGenerator) = 0;
	virtual void createInteractions(System &system) = 0;

	virtual double getEnergy(System &system) const;
	virtual double getOrbitalEnergy(System &system, std::string orbit) const;
	virtual vector3D getOrbitalAngularMomentum(System &system, std::string orbit) const;

	virtual ~Atom();

	static constexpr double electronMass = 1.0;
};

#endif /* ATOM_HPP */
