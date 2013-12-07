#ifndef AP_HYDROGEN_HPP
#define AP_HYDROGEN_HPP

#include <random>
#include <simulbody/simulator.hpp>
#include <cmath>

#include "atom.hpp"

using namespace simulbody;

class AbrinesPercivalHydrogen: public Atom {
	double reducedMass;

public:
	AbrinesPercivalHydrogen(System* system)
			: AbrinesPercivalHydrogen(system, Element::H, 1.00782503207) {
	}

	AbrinesPercivalHydrogen(System* system, Element nucleus, double atomicMass)
			: Atom(system, nucleus, atomicMass, Element::H) {
		reducedMass = (electronMass * nucleusMass) / (electronMass + nucleusMass);
		createInteractions();
		install();
	}

	virtual void install() {
		system->setBodyPosition(getNucleus(), vector3D(0, 0, 0));
		system->setBodyVelocity(getNucleus(), vector3D(0, 0, 0));

		system->setBodyPosition(getElectron("1s1"), vector3D(0, (1.0 / (reducedMass * nucleusCharge)), 0));
		system->setBodyVelocity(getElectron("1s1"), vector3D(0, 0, nucleusCharge));

		this->setPosition(vector3D(0, 0, 0));
		this->setVelocity(vector3D(0, 0, 0));
	}

	virtual void randomize(std::mt19937_64 &randomGenerator) {
		std::uniform_real_distribution<double> distMinusPiPi(-M_PI, M_PI);
		std::uniform_real_distribution<double> distMinusOneOne(-1, 1);
		std::uniform_real_distribution<double> distNull2Pi(0, 2 * M_PI);
		std::uniform_real_distribution<double> distNullOne(0, 1);

		double phi = distMinusPiPi(randomGenerator);
		double eta = distMinusPiPi(randomGenerator);
		double theta = acos(distMinusOneOne(randomGenerator));
		double epsilon = sqrt(distNullOne(randomGenerator));
		double thetaN = distNull2Pi(randomGenerator);

		double u = solveKeplerEquation(thetaN, epsilon, 1e-10, randomGenerator);

		double a = nucleusCharge / (2.0 * 0.5 * reducedMass * pow(nucleusCharge, 2));
		double b = sqrt(2.0 * 0.5 * reducedMass * reducedMass * pow(nucleusCharge, 2));

		vector3D C00(0, a * sqrt(1 - epsilon * epsilon) * sin(u), a * (cos(u) - epsilon));
		vector3D P00(0, b * sqrt(1 - epsilon * epsilon) * cos(u) / (1 - epsilon * cos(u)),
				-b * sin(u) / (1 - epsilon * cos(u)));

		vector3D C0 = C00.eulerRotation(phi, theta, eta);
		vector3D P0 = P00.eulerRotation(phi, theta, eta);

		system->setBodyPosition(getElectron("1s1"), system->getBodyPosition(nucleus) + C0);
		system->setBodyVelocity(getElectron("1s1"), system->getBodyVelocity(nucleus) + P0 / reducedMass);
	}

	virtual void createInteractions() {
		interactions.resize(1);
		interactions[0] = new CoulombInteraction(-1.0 * nucleusCharge, getNucleus(), getElectron("1s1"));
		system->addInteraction(interactions[0]);
	}

private:
	double solveKeplerEquation(double thetaN, double epsilon, double tolerance,
			std::mt19937_64 &randomGenerator) {

		std::uniform_real_distribution<double> dist(0.1, 6.0);
		double u0 = dist(randomGenerator);
		double u = u0;
		int rounds = 0;

		do {
			u0 = u;
			u = u - (thetaN + epsilon * sin(u) - u) / (epsilon * cos(u) - 1);
			rounds++;
		} while ((abs(u0 - u) > tolerance) && (rounds < 100));

		if (abs(u0 - u) > tolerance)
			return solveKeplerEquation(thetaN, epsilon, tolerance, randomGenerator);
		else
			return u;
	}
};

#endif /* AP_HYDROGEN_HPP */
