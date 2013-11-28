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

	Atom(System* system, Element element, unsigned int massNumber);
	Atom(System* system, Element element, unsigned int massNumber, Element electronConfiguration);

	identifier getNucleus() const;
	std::vector<identifier> getElectrons() const;
	std::vector<identifier> getBodies() const;
	identifier getElectron(std::size_t orbit) const;
	identifier getElectron(std::string orbitName) const;

	vector3D getPosition() const;
	vector3D getVelocity() const;
	vector3D getImpulse() const;
	double getMass() const;

	void setPosition(vector3D position);
	void setVelocity(vector3D velocity);

	virtual void install() = 0;
	virtual void randomize(std::mt19937_64 &randomGenerator) = 0;
	virtual void createInteractions() = 0;

	virtual double getEnergy() const;
	virtual double getOrbitalEnergy(std::string orbit) const;
	virtual vector3D getOrbitalAngularMomentum(std::string orbit) const;

	virtual ~Atom();

	static constexpr double electronMass = 1.0;
	static constexpr double protonMass = 1836.1527;
};

#endif /* ATOM_HPP */
