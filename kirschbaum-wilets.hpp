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

	double earthMass;
	double moonMass;

	double r2 = 0, r4 = 0;
	double p2 = 0, p4 = 0;
	double exponent = 0, factor = 0;

	vector3D p, pEarth;

	void calculateP4(const Phase& phase) {
		p = phase.getBodyVelocity(moon);
		pEarth = phase.getBodyVelocity(earth);
		p *= moonMass;
		pEarth *= earthMass;
		p -= pEarth;
		p2 = p.scalarProduct(p);
		p4 = p2 * p2;
	}

public:
	HeisenbergInteraction(double alpha, double mu, double xi, double earthMass, double moonMass,
			identifier earth, identifier moon)
			: alpha(alpha), mu(mu), xi(xi), earthMass(earthMass), moonMass(moonMass) {
		xi2 = xi * xi;
		xi4 = xi2 * xi2;
		this->setBodies(earth, moon);
	}

	virtual void apply(const Phase &phase, const double t) override {

		calculateR(phase);
		r2 = r.scalarProduct(r);
		r4 = r2 * r2;

		calculateP4(phase);

		exponent = exp(alpha * (1 - r4 * p4 / xi4));
		factor = xi2 * exponent / (2 * alpha * mu * r4) + p4 * exponent / (xi2 * mu);

		F = r;
		F *= factor;

		applyFOnMoon(phase);
		applyFOnEarth(phase);
	}

	virtual double getEnergy(const Phase &phase) override {
		calculateR(phase);
		r2 = r.scalarProduct(r);
		r4 = r2 * r2;

		calculateP4(phase);

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
	}

	virtual void randomize(std::mt19937_64 &randomEngine) {
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
					new HeisenbergInteraction(5, reducedMass, 1.0488, nucleusMass, electronMass, nucleus,
							e1));
		}

		for (Interaction* interaction : interactions) {
			system->addInteraction(interaction);
		}
	}

};

#endif /* KIRSCHBAUM_WILETS_HPP */
