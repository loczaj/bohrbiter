#include "atom.hpp"

using namespace simulbody;
using namespace std;

Atom::Atom(System* system, Element element, double atomicMass)
		: Atom(system, element, element, atomicMass) {
}

Atom::Atom(System* system, Element electronConfiguration, Element nucleusElement, double atomicMass)
		: system(system), electronConfiguration(electronConfiguration), nucleusElement(nucleusElement) {

	nucleusMass = PeriodicTable::nucleusMassInAU(nucleusElement, atomicMass);
	nucleusCharge = (double) PeriodicTable::atomicNumber(nucleusElement);
	nucleus = system->createBody(nucleusMass);
	reducedMass = (electronMass * nucleusMass) / (electronMass + nucleusMass);

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

std::vector<Interaction*> Atom::getInteractions() const {
	return interactions;
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

double Atom::getReducedMass() const {
	return reducedMass;
}

double Atom::getNucleusMass() const {
	return nucleusMass;
}

double Atom::getNucleusCharge() const {
	return nucleusCharge;
}

Element Atom::getNucleusElement() const {
	return nucleusElement;
}

Element Atom::getElectronConfiguration() const {
	return electronConfiguration;
}

std::vector<string> Atom::getOrbitNames() const {
	return orbitNames;
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
