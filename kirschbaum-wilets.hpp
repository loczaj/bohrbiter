#ifndef KIRSCHBAUM_WILETS_HPP
#define KIRSCHBAUM_WILETS_HPP

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
	double exponent = 0, forceFactor = 0, velocityFactor = 0;

	double muSlashXi2 = 0;
	double xi2DotMu = 0;
	double xi2SlashAlphaSlashMuSlash2 = 0;

public:
	HeisenbergInteraction(double alpha, double xi, identifier nucleus, identifier electron);

	virtual void setBodyMasses(double earthMass, double moonMass) override;
	virtual void apply(const Phase &x, Phase &dxdt, const double t) override;
	virtual double getEnergy(const Phase &phase) override;

	virtual ~HeisenbergInteraction();
};

class KirschbaumWiletsAtom: public Atom {

public:
	KirschbaumWiletsAtom(System* system, Element element, double atomicMass);
	KirschbaumWiletsAtom(System* system, Element electronConfig, Element nucleusElement, double atomicMass);

	virtual void install() override;
	virtual void randomize(std::mt19937_64 &randomEngine) override;
	virtual void createInteractions() override;
};

#endif /* KIRSCHBAUM_WILETS_HPP */
