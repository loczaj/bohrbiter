#include <stdexcept>
#include <cmath>

#include "abrines-percival.hpp"
#include <simulbody/interactions/coulomb.hpp>


using namespace simulbody;

AbrinesPercivalAtom::AbrinesPercivalAtom(System* system, Element element, double atomicMass)
		: AbrinesPercivalAtom(system, element, element, atomicMass) {
}

AbrinesPercivalAtom::AbrinesPercivalAtom(System* system, Element electronConfig, Element nucleusElement,
		double atomicMass)
		: Atom(system, electronConfig, nucleusElement, atomicMass) {

	if (electronConfiguration != Element::H && electronConfiguration != Element::He) {
		throw std::invalid_argument("AbrinesPercivalAtom supports only Hydrogen and Helium.");
	}

	createInteractions();
}

void AbrinesPercivalAtom::install() {
	vector3D r(0, 1.0 / (reducedMass * nucleusCharge), 0);
	vector3D v(0, 0, nucleusCharge);

	system->setBodyPosition(getElectron("1s1"), system->getBodyPosition(nucleus) + r);
	system->setBodyVelocity(getElectron("1s1"), system->getBodyVelocity(nucleus) + v);

	if (electronConfiguration == Element::He) {
		system->setBodyPosition(getElectron("1s2"), system->getBodyPosition(nucleus) - r);
		system->setBodyVelocity(getElectron("1s2"), system->getBodyVelocity(nucleus) - v);
	}

	this->setPosition(vector3D(0, 0, 0));
	this->setVelocity(vector3D(0, 0, 0));
}

void AbrinesPercivalAtom::randomize(std::mt19937_64 &randomEngine) {

	std::uniform_real_distribution<double> distMinusPiPi(-M_PI, M_PI);
	std::uniform_real_distribution<double> distMinusOneOne(-1, 1);
	std::uniform_real_distribution<double> distNull2Pi(0, 2 * M_PI);
	std::uniform_real_distribution<double> distNullOne(0, 1);

	double phi = distMinusPiPi(randomEngine);
	double eta = distMinusPiPi(randomEngine);
	double theta = acos(distMinusOneOne(randomEngine));
	double epsilon = sqrt(distNullOne(randomEngine));
	double thetaN = distNull2Pi(randomEngine);

	double u = solveKeplerEquation(thetaN, epsilon, 1e-10, randomEngine);

	double a = nucleusCharge / (2.0 * 0.5 * reducedMass * pow(nucleusCharge, 2));
	double b = sqrt(2.0 * 0.5 * reducedMass * reducedMass * pow(nucleusCharge, 2));

	vector3D C00(0, a * sqrt(1 - epsilon * epsilon) * sin(u), a * (cos(u) - epsilon));
	vector3D P00(0, b * sqrt(1 - epsilon * epsilon) * cos(u) / (1 - epsilon * cos(u)),
			-b * sin(u) / (1 - epsilon * cos(u)));

	vector3D C0 = C00.eulerRotation(phi, theta, eta);
	vector3D P0 = P00.eulerRotation(phi, theta, eta);

	system->setBodyPosition(getElectron("1s1"), system->getBodyPosition(nucleus) + C0);
	system->setBodyVelocity(getElectron("1s1"), system->getBodyVelocity(nucleus) + P0 / reducedMass);

	if (electronConfiguration == Element::He) {
		system->setBodyPosition(getElectron("1s2"), system->getBodyPosition(nucleus) - C0);
		system->setBodyVelocity(getElectron("1s2"), system->getBodyVelocity(nucleus) - P0 / reducedMass);
	}
}

void AbrinesPercivalAtom::createInteractions() {
	interactions.clear();

	for (identifier e1 : getElectrons()) {
		for (identifier e2 : getElectrons()) {
			if (e1 < e2)
				interactions.push_back(new CoulombInteraction(1.0, e1, e2));
		}
		interactions.push_back(new CoulombInteraction(-1.0 * nucleusCharge, e1, nucleus));
	}

	for (Interaction* interaction : interactions) {
		system->addInteraction(interaction);
	}
}

double AbrinesPercivalAtom::solveKeplerEquation(double thetaN, double epsilon, double tolerance,
		std::mt19937_64 &randomEngine) {

	std::uniform_real_distribution<double> dist(0.1, 6.0);
	double u0 = dist(randomEngine);
	double u = u0;
	int rounds = 0;

	do {
		u0 = u;
		u = u - (thetaN + epsilon * sin(u) - u) / (epsilon * cos(u) - 1);
		rounds++;
	} while ((abs(u0 - u) > tolerance) && (rounds < 100));

	if (abs(u0 - u) > tolerance)
		return solveKeplerEquation(thetaN, epsilon, tolerance, randomEngine);
	else
		return u;
}
