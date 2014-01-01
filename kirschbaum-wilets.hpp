#ifndef KIRSCHBAUM_WILETS_HPP
#define KIRSCHBAUM_WILETS_HPP

#include <cmath>
#include <random>
#include <simulbody/simulator.hpp>

#include "atom.hpp"
#include "cohen.hpp"

using namespace simulbody;

class HeisenbergInteraction: public Interaction {

private:
	double alpha;
	double mu;
	double xi, xi2, xi4;
	double nucleusMass;

	vector3D cv;

	double r2 = 0, r4 = 0;
	double p2 = 0, p4 = 0;
	double exponent = 0, factor = 0, cvfactor = 0;

public:
	HeisenbergInteraction(double alpha, double mu, double xi, identifier nucleus, identifier electron,
			double nucleusMass)
			: alpha(alpha), mu(mu), xi(xi), nucleusMass(nucleusMass) {
		xi2 = xi * xi;
		xi4 = xi2 * xi2;
		this->setBodies(nucleus, electron);
	}

	virtual void apply(const Phase &x, Phase &dxdt, const double t) override {

		calculateRelativePositionR(x);
		r2 = r.scalarProduct(r);
		r4 = r2 * r2;

		calculateRelativeVelocityV(x);
		p2 = v.scalarProduct(v); // electron mass = 1.0
		p4 = p2 * p2;

		exponent = exp(alpha * (1 - r4 * p4 / xi4));
		factor = xi2 * exponent / (2 * alpha * mu * r4) + p4 * exponent / (xi2 * mu);
		F = r * factor;

		applyForceOnMoon(dxdt, F);
		applyForceOnEarth(dxdt, -F);

		cvfactor = -p2 * r2 * exponent / (xi2 * mu);
		cv = v * cvfactor;

		addCollateralVelocityOnMoon(dxdt, cv);
		addCollateralVelocityOnEarth(dxdt, -cv / nucleusMass);
	}

	virtual double getEnergy(const Phase &phase) override {
		calculateRelativePositionR(phase);
		r2 = r.scalarProduct(r);
		r4 = r2 * r2;

		calculateRelativeVelocityV(phase);
		p2 = v.scalarProduct(v);
		p4 = p2 * p2;

		return xi2 * exp(alpha * (1 - r4 * p4 / xi4)) / (4 * alpha * r2 * mu);
	}

	virtual ~HeisenbergInteraction() override {
	}

};

class KirschbaumWiletsAtom: public Atom {

public:
	KirschbaumWiletsAtom(System* system, Element element, double atomicMass)
			: KirschbaumWiletsAtom(system, element, element, atomicMass) {
	}

	KirschbaumWiletsAtom(System* system, Element electronConfig, Element nucleusElement, double atomicMass)
			: Atom(system, electronConfig, nucleusElement, atomicMass) {

		createInteractions();
	}

	virtual void install() {

		CohenConfiguration configuration;

		system->setBodyPosition(nucleus, vector3D(0.0, 0.0, 0.0));
		system->setBodyVelocity(nucleus, vector3D(0.0, 0.0, 0.0));

		for (string orbit : orbitNames) {
			system->setBodyPosition(getElectron(orbit), configuration.position(electronConfiguration, orbit));
			system->setBodyVelocity(getElectron(orbit), configuration.momentum(electronConfiguration, orbit));
		}
	}

	virtual void randomize(std::mt19937_64 &randomEngine) {
		std::uniform_real_distribution<double> distMinusPiPi(-M_PI, M_PI);
		std::uniform_real_distribution<double> distMinusOneOne(-1, 1);

		double phi = distMinusPiPi(randomEngine);
		double eta = distMinusPiPi(randomEngine);
		double theta = acos(distMinusOneOne(randomEngine));

		install();

		for (string orbit : orbitNames) {
			identifier electron = getElectron(orbit);
			system->setBodyPosition(electron,
					system->getBodyPosition(electron).eulerRotation(phi, theta, eta));
			system->setBodyVelocity(electron,
					system->getBodyVelocity(electron).eulerRotation(phi, theta, eta));
		}
	}

	virtual void createInteractions() {
		interactions.clear();

		for (identifier e1 : getElectrons()) {
			for (identifier e2 : getElectrons()) {
				if (e1 < e2)
					interactions.push_back(new CoulombInteraction(1.0, e1, e2));
			}
			interactions.push_back(new CoulombInteraction(-1.0 * nucleusCharge, nucleus, e1));
			interactions.push_back(
					new HeisenbergInteraction(5.0, reducedMass, 0.9535, nucleus, e1, nucleusMass));
		}

		for (Interaction* interaction : interactions) {
			system->addInteraction(interaction);
		}
	}

};

#endif /* KIRSCHBAUM_WILETS_HPP */
