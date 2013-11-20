#include "atom.hpp"

using namespace simulbody;
using namespace std;

Atom::Atom(System &system, Element element, unsigned int massNumber)
		: Atom(system, element, massNumber, element) {
}

Atom::Atom(System &system, Element element, unsigned int massNumber, Element electronConfiguration)
		: element(element), massNumber(massNumber), electronConfiguration(electronConfiguration) {

	nucleusMass = 1836.0 * (double) massNumber;
	nucleusCharge = (double) atomicNumber(element);
	nucleus = system.createBody(nucleusMass);

	orbitNames = vector<string>(orbitals()[electronConfiguration]);
	for (string orbitName : orbitNames) {
		electrons[orbitName] = system.createBody(electronMass);
	}

	setInteractions();
}

identifier Atom::getNucleus() const {
	return nucleus;
}

vector<identifier> Atom::getElectrons() const {
	vector<identifier> result;
	for (string orbitName : orbitNames) {
		result.push_back(getElectron(orbitName));
	}
	return result;
}

std::vector<identifier> Atom::getBodies() const {
	vector<identifier> bodies = getElectrons();
	bodies.push_back(nucleus);
	return bodies;
}

identifier Atom::getElectron(size_t orbit) const {
	return getElectron(orbitNames.at(orbit));
}

identifier Atom::getElectron(string orbitName) const {
	return electrons[orbitName];
}

vector3D Atom::getPosition(System &system) const {
	return system.getCenterOfMass(getBodies());
}

vector3D Atom::getImpulse(System &system) const {
	return system.getImpulse(getBodies());
}

vector3D Atom::getVelocity(System &system) const {
	return getImpulse(system) / getMass(system);
}

double Atom::getMass(System &system) const {
	return system.getMass(getBodies());
}

void Atom::setPosition(System &system, vector3D position) {
	vector3D delta = position - getPosition(system);
	for (identifier body : getBodies()) {
		system.setBodyPosition(body, system.getBodyPosition(body) + delta);
	}
}

void Atom::setVelocity(System &system, vector3D velocity) {
	vector3D delta = velocity - getVelocity(system);
	for (identifier body : getBodies()) {
		system.setBodyVelocity(body, system.getBodyVelocity(body) + delta);
	}
}

double Atom::getEnergy(System &system) const {
	double energy = 0.0;

	for (identifier e1 : getElectrons()) {
		for (identifier e2 : getElectrons()) {
			if (e1 < e2)
				energy += system.getPairPotentialEnergy(e1, e2);
		}
		energy += system.getPairPotentialEnergy(e1, nucleus);
	}

	for (identifier body : getBodies()) {
		energy += system.getBodyKineticEnergyReferenced(body, getVelocity(system));
	}

	return energy;
}

double Atom::getOrbitalEnergy(System &system, std::string orbit) const {
	double energy = 0.0;

	identifier e1 = getElectron(orbit);
	for (identifier e2 : getElectrons()) {
		if (e1 != e2)
			energy += system.getPairPotentialEnergy(e1, e2);
	}

	energy += system.getPairPotentialEnergy(e1, nucleus);
	energy += system.getBodyKineticEnergyReferenced(e1, getVelocity(system));

	return energy;
}

vector3D Atom::getAngularMomentum(System &system, std::string orbit) const {
	identifier electron = getElectron(orbit);
	return system.getBodyAngularMomentum(electron, nucleus);
}
