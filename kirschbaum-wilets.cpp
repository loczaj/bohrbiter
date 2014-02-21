#include <cmath>

#include "kirschbaum-wilets.hpp"

using namespace simulbody;

HeisenbergInteraction::HeisenbergInteraction(double alpha, double xi, identifier nucleus, identifier electron)
		: alpha(alpha), xi(xi) {
	xi2 = xi * xi;
	xi4 = xi2 * xi2;
	this->setBodies(nucleus, electron);
}

void HeisenbergInteraction::setBodyMasses(double earthMass, double moonMass) {
	Interaction::setBodyMasses(earthMass, moonMass);

	muSlashXi2 = reducedMass / xi2;
	xi2DotMu = xi2 * reducedMass;
	xi2SlashAlphaSlashMuSlash2 = xi2 / alpha / reducedMass / 2;
}

void HeisenbergInteraction::apply(const Phase &x, Phase &dxdt, const double t) {

	calculateRelativePosition(x);
	r2 = relativePosition.scalarProduct(relativePosition);
	r4 = r2 * r2;

	calculateRelativeVelocity(x);
	p = relativeVelocity * reducedMass;
	p2 = p.scalarProduct(p);
	p4 = p2 * p2;

	exponent = exp(alpha * (1 - r4 * p4 / xi4));
	forceFactor = xi2SlashAlphaSlashMuSlash2 * exponent / r4 + p4 * exponent / xi2DotMu;
	actingForce = relativePosition * forceFactor;

	applyForceOnMoon(dxdt, actingForce);
	applyForceOnEarth(dxdt, -actingForce);

	// A Here we assume moonMass=1 (Heisenberg-interaction applies to electrons)
	velocityFactor = -p2 * muSlashXi2 * r2 * exponent;
	actingVelocity = relativeVelocity * velocityFactor;

	addVelocityOnMoon(dxdt, actingVelocity);
	addVelocityOnEarth(dxdt, -actingVelocity / earthMass);
}

double HeisenbergInteraction::getEnergy(const Phase &phase) {
	calculateRelativePosition(phase);
	r2 = relativePosition.scalarProduct(relativePosition);
	r4 = r2 * r2;

	calculateRelativeVelocity(phase);
	p = relativeVelocity * reducedMass;
	p2 = p.scalarProduct(p);
	p4 = p2 * p2;

	return xi2 * exp(alpha * (1 - r4 * p4 / xi4)) / (4 * alpha * r2 * reducedMass);
}

HeisenbergInteraction::~HeisenbergInteraction() {
}

KirschbaumWiletsAtom::KirschbaumWiletsAtom(System* system, Element element, double atomicMass)
		: KirschbaumWiletsAtom(system, element, element, atomicMass) {
}

KirschbaumWiletsAtom::KirschbaumWiletsAtom(System* system, Element electronConfig, Element nucleusElement,
		double atomicMass)
		: Atom(system, electronConfig, nucleusElement, atomicMass) {

	createInteractions();
}

void KirschbaumWiletsAtom::install() {

	CohenConfiguration configuration;

	system->setBodyPosition(nucleus, vector3D(0.0, 0.0, 0.0));
	system->setBodyVelocity(nucleus, vector3D(0.0, 0.0, 0.0));

	for (string orbit : orbitNames) {
		system->setBodyPosition(getElectron(orbit), configuration.position(electronConfiguration, orbit));
		system->setBodyVelocity(getElectron(orbit), configuration.momentum(electronConfiguration, orbit));
	}
}

void KirschbaumWiletsAtom::randomize(std::mt19937_64 &randomEngine) {
	std::uniform_real_distribution<double> distMinusPiPi(-M_PI, M_PI);
	std::uniform_real_distribution<double> distMinusOneOne(-1, 1);

	double phi = distMinusPiPi(randomEngine);
	double theta = acos(distMinusOneOne(randomEngine));
	double eta = distMinusPiPi(randomEngine);

	install();

	for (identifier electron : getElectrons()) {
		system->setBodyPosition(electron, system->getBodyPosition(electron).eulerRotation(phi, theta, eta));
		system->setBodyVelocity(electron, system->getBodyVelocity(electron).eulerRotation(phi, theta, eta));
	}
}

void KirschbaumWiletsAtom::createInteractions() {
	interactions.clear();

	for (identifier e1 : getElectrons()) {
		for (identifier e2 : getElectrons()) {
			if (e1 < e2)
				interactions.push_back(new CoulombInteraction(1.0, e1, e2));
		}
		interactions.push_back(new CoulombInteraction(-1.0 * nucleusCharge, nucleus, e1));
		interactions.push_back(new HeisenbergInteraction(45.0, 1.257, nucleus, e1));
	}

	for (Interaction* interaction : interactions) {
		system->addInteraction(interaction);
	}
}
