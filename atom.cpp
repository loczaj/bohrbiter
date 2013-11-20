#include "atom.hpp"

using namespace simulbody;

Atom::Atom(System &system, Element element, unsigned int massNumber)
		: Atom(system, element, massNumber, element) {
}

Atom::Atom(System &system, Element element, unsigned int massNumber, Element electronConfiguration)
		: element(element), massNumber(massNumber), electronConfiguration(electronConfiguration) {

	nucleusMass = 1836.0 * (double) massNumber;
	nucleusCharge = (double) atomicNumber(element);

	nucleus = system.createBody(nucleusMass);

	for (std::string orbitName : orbits.names[electronConfiguration]) {
		electrons[orbitName] = system.createBody(electronMass);
	}

	setInteractions();
}

identifier Atom::getNucleus() {
	return nucleus;
}

std::vector<identifier> Atom::getElectrons() {
	std::vector<identifier> result;
	for (std::string orbitName : orbits.names[electronConfiguration]) {
		result.push_back(getElectron(orbitName));
	}
	return result;
}

identifier Atom::getElectron(std::size_t orbit) {
	return getElectron(orbits.names[electronConfiguration].at(orbit));
}

identifier Atom::getElectron(std::string orbitName) {
	return electrons[orbitName];
}

void Atom::moveTo(System &system, vector3D position) {
	std::vector<identifier> bodies = getElectrons();
	bodies.push_back(nucleus);

	vector3D delta = position - system.getCenterOfMass(bodies);
	for (identifier body : bodies) {
		system.setBodyPosition(body, system.getBodyPosition(body) + delta);
	}
}

void Atom::speedUp(System &system, vector3D velocity) {
	std::vector<identifier> bodies = getElectrons();
	bodies.push_back(nucleus);

	vector3D delta = velocity - system.getImpulse(bodies) / system.getMass(bodies);
	for (identifier body : bodies) {
		system.setBodyVelocity(body, system.getBodyVelocity(body) + delta);
	}
}
