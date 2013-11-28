#ifndef AP_HYDROGEN_HPP
#define AP_HYDROGEN_HPP

#include <random>
#include <simulbody/simulator.hpp>
#include <cmath>

#include "atom.hpp"

using namespace simulbody;

class AbrinesPercivalHydrogen: public Atom {
public:
	AbrinesPercivalHydrogen(System &system)
			: Atom(system, Element::H, 1) {
		createInteractions(system);
		install(system);
	}

	virtual void install(System &system) {
		system.setBodyPosition(getNucleus(), vector3D(0, 0, 0));
		system.setBodyVelocity(getNucleus(), vector3D(0, 0, 0));

		system.setBodyPosition(getElectron("1s1"), vector3D(0, 1, 0));
		system.setBodyVelocity(getElectron("1s1"), vector3D(0, 0, 1));

		this->setPosition(system, vector3D(0, 0, 0));
		this->setVelocity(system, vector3D(0, 0, 0));
	}

	virtual void randomize(System &system, std::mt19937_64 &randomGenerator) {
		std::uniform_real_distribution<double> distMinusPiPi(-M_PI, M_PI);
		std::uniform_real_distribution<double> distMinusOneOne(-1, 1);
		std::uniform_real_distribution<double> distNull2Pi(0, 2 * M_PI);
		std::uniform_real_distribution<double> distNullOne(0, 1);

		double phi = distMinusPiPi(randomGenerator);
		double eta = distMinusPiPi(randomGenerator);
		double theta = acos(distMinusOneOne(randomGenerator));
		double epsilon = sqrt(distNullOne(randomGenerator));
		double thetaN = distNull2Pi(randomGenerator);

		double u = solveKeplerEquation(thetaN, epsilon, 0.5, 1e-10);
		double a = nucleusCharge / (2.0 * 0.5);
		double b = sqrt(2.0 * 0.5);

		vector3D c00(0, a * sqrt(1 - epsilon * epsilon) * sin(u), a * (cos(u) - epsilon));
		vector3D p00(0, b * sqrt(1 - epsilon * epsilon) * cos(u) / (1 - epsilon * cos(u)),
				-b * sin(u) / (1 - epsilon * cos(u)));

		vector3D c0 = c00.eulerRotation(phi, theta, eta);
		vector3D p0 = p00.eulerRotation(phi, theta, eta);

		system.setBodyPosition(getElectron("1s1"), system.getBodyPosition(nucleus) + c0);
		system.setBodyVelocity(getElectron("1s1"), system.getBodyVelocity(nucleus) + p0);
	}

	virtual void createInteractions(System &system) {
		interactions.resize(1);
		interactions[0] = new CoulombInteraction(-1.0, getNucleus(), getElectron("1s1"));
		system.addInteraction(interactions[0]);
	}

private:
	double solveKeplerEquation(double thetaN, double epsilon, double u0, double tolerance) {
		double u = u0;
		do {
			u0 = u;
			u = u - (thetaN + epsilon * sin(u) - u) / (epsilon * cos(u) - 1);
		} while (abs(u0 - u) > tolerance);

		return u;
	}
};

#endif /* AP_HYDROGEN_HPP */
