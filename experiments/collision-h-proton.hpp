#ifndef COLLISION_H_PROTON_HPP
#define COLLISION_H_PROTON_HPP

#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../abrines-percival.hpp"
#include "../experiment.hpp"

using namespace std;

class CollisionAbrinesPercivalHydrogenWithProton: public Experiment {

	ofstream stream;

	System bbsystem;
	Printer* printer;
	PositionPrintField printField;

	identifier projectile;
	AbrinesPercivalAtom* hydrogen;
	DistanceCondition* condition;
	Interaction* coulombProjectileElectron;
	Interaction* coulombProjectileNucleus;

	double b2max;
	double projectileVelocity;
	double absoluteStepperError;
	double relativeStepperError;
	double relativeEnergyError;

	int ionization = 0;
	int ecapture = 0;
	int extended = 0;

public:

	CollisionAbrinesPercivalHydrogenWithProton(double impact2max, double energykeV,
			double absoluteStepperError, double relativeStepperError, double relativeEnergyError)
			: b2max(impact2max), absoluteStepperError(absoluteStepperError), relativeStepperError(
					relativeStepperError), relativeEnergyError(relativeEnergyError) {

		projectileVelocity = Utils::calculateAcceleratedVelocityInAU(Atom::protonMass, 1.0, energykeV);

		projectile = bbsystem.createBody(Atom::protonMass);
		hydrogen = new AbrinesPercivalAtom(&bbsystem, Element::H, 1.00782503207);
		condition = new DistanceCondition(projectile, hydrogen->getNucleus(), 51.0);
		printer = nullptr;

		coulombProjectileElectron = new CoulombInteraction(-1.0, projectile, hydrogen->getElectron("1s1"));
		coulombProjectileNucleus = new CoulombInteraction(hydrogen->getNucleusCharge(), projectile,
				hydrogen->getNucleus());

		bbsystem.addInteraction(coulombProjectileElectron);
		bbsystem.addInteraction(coulombProjectileNucleus);
	}

	int open(int numberOfRounds, bool seedRandom) {
		stream.open("result.csv");
		stream.precision(10);

		return 0;
	}

	int run(int round, bool tracking, bool skipUntracked) {
		runge_kutta_dopri5<Phase> stepper;
		auto ctrdStepper = make_controlled(absoluteStepperError, relativeStepperError, stepper);
		Simulator<decltype(ctrdStepper)> simulator(ctrdStepper, &bbsystem);

		uniform_real_distribution<double> distributionNullB2Max(0, b2max);
		double b = sqrt(distributionNullB2Max(randomEngine));

		hydrogen->randomize(randomEngine);
		hydrogen->setPosition(vector3D(0, 0, 0));
		hydrogen->setVelocity(vector3D(0, 0, 0));

		bbsystem.setBodyPosition(projectile, vector3D(0, b, -50.0));
		bbsystem.setBodyVelocity(projectile, vector3D(0, 0, projectileVelocity));

		if (skipUntracked && !tracking) {
			return 0;
		}

		stream << bbsystem.phase;
		stream.flush();

		if (tracking) {
			printer = new Printer(to_string(round) + ".csv");
			printer->addField(&printField);
			simulator.setObserver(*printer);
		}

		if (condition->evaluate(bbsystem.phase, 0)) {
			stream << "\t" << "Initial condition error" << endl;
			return -3;
		}

		double energy = bbsystem.getSystemEnergy();
		double time = simulator.simulate(0.0, 1.0, 0.0001, *condition, 100);
		bool eBoundToTarget;
		bool eBoundToProjec;

		while (true) {

			if (time < 0.0) {
				stream << "\t" << "Distance not reached error" << endl;
				return -1;
			}

			if (abs((energy - bbsystem.getSystemEnergy()) / energy) > relativeEnergyError) {
				stream << "\t" << "Energy error: " << energy << " vs. " << bbsystem.getSystemEnergy() << endl;
				return -2;
			}

			eBoundToTarget = Utils::isBound(bbsystem, hydrogen->getElectron("1s1"), hydrogen->getNucleus());
			eBoundToProjec = Utils::isBound(bbsystem, hydrogen->getElectron("1s1"), projectile);

			if (!eBoundToTarget || !eBoundToProjec) {
				break;
			}

			stream << " " << round << " Extend Run";
			simulator.simulate(time, time + 1.0, 0.0001);
			time += 1.0;
			extended++;
		}

		if (eBoundToTarget && !eBoundToProjec) {
			stream << endl;
		}

		if (!eBoundToTarget && eBoundToProjec) {
			stream << "\t" << round << " --> Electron Capture" << endl;
			ecapture++;
		}

		if (!eBoundToTarget && !eBoundToProjec) {
			stream << "\t" << round << " --> Ionization" << endl;
			ionization++;
		}

		if (printer != nullptr)
			delete printer;

		stream.flush();
		return 0;
	}

	int close(int successfulRounds) {
		double ionizationRate = ((double) ionization) / ((double) successfulRounds);
		double ecaptureRate = ((double) ecapture) / ((double) successfulRounds);
		double extendedRate = ((double) extended) / ((double) successfulRounds);
		cout << "Ionization: " << ionization << " (" << ionizationRate * 100.0 << " %)" << endl;
		cout << "El.Capture: " << ecapture << " (" << ecaptureRate * 100.0 << " %)" << endl;
		cout << "Extended run: " << extended << " (" << extendedRate * 100.0 << " %)" << endl << endl;
		cout << "Cross sections:" << endl;

		cout << "\t Ionization: " << ionizationRate * M_PI * b2max << endl;
		cout << "\t El.Capture: " << ecaptureRate * M_PI * b2max << endl << endl;

		stream.close();
		return 0;
	}
};

#endif /* COLLISION_H_PROTON_HPP */
