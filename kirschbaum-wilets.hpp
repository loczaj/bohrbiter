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
	double xi, xi2, xi4;

	vector3D p;

	double r2 = 0, r4 = 0;
	double p2 = 0, p4 = 0;
	double exponent = 0, factor = 0, vcollfactor = 0;

	double muSlashXi2 = 0;
	double xi2DotMu = 0;
	double xi2SlashAlphaSlashMuSlash2 = 0;

public:
	HeisenbergInteraction(double alpha, double xi, identifier nucleus, identifier electron)
			: alpha(alpha), xi(xi) {
		xi2 = xi * xi;
		xi4 = xi2 * xi2;
		this->setBodies(nucleus, electron);
	}

	virtual void setBodyMasses(double earthMass, double moonMass) override {
		Interaction::setBodyMasses(earthMass, moonMass);

		muSlashXi2 = reducedMass / xi2;
		xi2DotMu = xi2 * reducedMass;
		xi2SlashAlphaSlashMuSlash2 = xi2 / alpha / reducedMass / 2;
	}

	virtual void apply(const Phase &x, Phase &dxdt, const double t) override {

		calculateRelativePositionR(x);
		r2 = r.scalarProduct(r);
		r4 = r2 * r2;

		calculateRelativeVelocityV(x);
		p = v * reducedMass;
		p2 = p.scalarProduct(p);
		p4 = p2 * p2;

		exponent = exp(alpha * (1 - r4 * p4 / xi4));
		factor = xi2SlashAlphaSlashMuSlash2 * exponent / r4 + p4 * exponent / xi2DotMu;
		F = r * factor;

		applyForceOnMoon(dxdt, F);
		applyForceOnEarth(dxdt, -F);

		// A Here we assume moonMass=1 (Heisenberg-interaction applies to electrons)
		vcollfactor = -p2 * muSlashXi2 * r2 * exponent;
		vcoll = v * vcollfactor;

		addCollateralVelocityOnMoon(dxdt, vcoll);
		addCollateralVelocityOnEarth(dxdt, -vcoll / earthMass);
	}

	virtual double getEnergy(const Phase &phase) override {
		calculateRelativePositionR(phase);
		r2 = r.scalarProduct(r);
		r4 = r2 * r2;

		calculateRelativeVelocityV(phase);
		p = v * reducedMass;
		p2 = p.scalarProduct(p);
		p4 = p2 * p2;

		return xi2 * exp(alpha * (1 - r4 * p4 / xi4)) / (4 * alpha * r2 * reducedMass);
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
		double theta = distMinusPiPi(randomEngine); //acos(distMinusOneOne(randomEngine));

		install();

		for (identifier electron : getElectrons()) {
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
			interactions.push_back(new HeisenbergInteraction(5.0, 1.257, nucleus, e1));
		}

		for (Interaction* interaction : interactions) {
			system->addInteraction(interaction);
		}
	}

};

#endif /* KIRSCHBAUM_WILETS_HPP */
