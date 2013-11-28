#include "atom.hpp"

using namespace simulbody;
using namespace std;

Atom::Atom(System* system, Element element, unsigned int massNumber)
		: Atom(system, element, massNumber, element) {
}

Atom::Atom(System* system, Element element, unsigned int massNumber, Element electronConfiguration)
		: system(system), element(element), electronConfiguration(electronConfiguration), massNumber(
				massNumber) {

	nucleusMass = 1836.1527 * (double) massNumber;
	nucleusCharge = (double) PeriodicTable::atomicNumber(element);
	nucleus = system->createBody(nucleusMass);

	orbitNames = PeriodicTable().atomicOrbitals(electronConfiguration);
	for (string orbitName : orbitNames) {
		electrons[orbitName] = system->createBody(electronMass);
	}
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
	vector<identifier> bodies = append( { nucleus }, getElectrons());
	return bodies;
}

identifier Atom::getElectron(size_t orbit) const {
	return getElectron(orbitNames.at(orbit));
}

identifier Atom::getElectron(string orbitName) const {
	return electrons.at(orbitName);
}

vector3D Atom::getPosition() const {
	return system->getGroupCenterOfMass(getBodies());
}

vector3D Atom::getImpulse() const {
	return system->getGroupImpulse(getBodies());
}

vector3D Atom::getVelocity() const {
	return getImpulse() / getMass();
}

double Atom::getMass() const {
	return system->getGroupMass(getBodies());
}

void Atom::setPosition(vector3D position) {
	vector3D delta = position - getPosition();
	for (identifier body : getBodies()) {
		system->setBodyPosition(body, system->getBodyPosition(body) + delta);
	}
}

void Atom::setVelocity(vector3D velocity) {
	vector3D delta = velocity - getVelocity();
	for (identifier body : getBodies()) {
		system->setBodyVelocity(body, system->getBodyVelocity(body) + delta);
	}
}

double Atom::getEnergy() const {
	double energy = 0.0;

	for (identifier e1 : getElectrons()) {
		for (identifier e2 : getElectrons()) {
			if (e1 < e2)
				energy += system->getPairPotentialEnergy(e1, e2);
		}
		energy += system->getPairPotentialEnergy(e1, nucleus);
	}

	for (identifier body : getBodies()) {
		energy += system->getBodyKineticEnergyReferenced(body, getVelocity());
	}

	return energy;
}

double Atom::getOrbitalEnergy(std::string orbit) const {
	double energy = 0.0;

	identifier e1 = getElectron(orbit);
	for (identifier e2 : getElectrons()) {
		if (e1 != e2)
			energy += system->getPairPotentialEnergy(e1, e2);
	}

	energy += system->getPairPotentialEnergy(e1, nucleus);
	energy += system->getBodyKineticEnergyReferenced(e1, nucleus);

	return energy;
}

vector3D Atom::getOrbitalAngularMomentum(std::string orbit) const {
	identifier electron = getElectron(orbit);
	return system->getBodyAngularMomentum(electron, nucleus);
}

Atom::~Atom() {
}
