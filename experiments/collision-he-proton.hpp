#ifndef COLLISION_HE_PROTON_HPP
#define COLLISION_HE_PROTON_HPP

#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../kirschbaum-wilets.hpp"
#include "../experiment.hpp"

using namespace std;

class CollisionKirschbaumWiletsHeliumWithProton: public Experiment {

	ofstream stream;

	System bbsystem;
	Printer* printer;
	PositionPrintField printField;

	identifier projectile;
	KirschbaumWiletsAtom* helium;
	DistanceCondition* condition;
	Interaction* coulombProjectile1s1;
	Interaction* coulombProjectile1s2;
	Interaction* coulombProjectileNucleus;
	Interaction* heisenbergProjectile1s1;
	Interaction* heisenbergProjectile1s2;

	double b2max;
	double projectileVelocity;
	double absoluteStepperError;
	double relativeStepperError;
	double relativeEnergyError;

	int ionization1 = 0;
	int ionization2 = 0;
	int ecapture1 = 0;
	int ecapture2 = 0;
	int ionizationAndCapture = 0;
	int extended = 0;

public:

	CollisionKirschbaumWiletsHeliumWithProton(double impact2max, double energykeV,
			double absoluteStepperError, double relativeStepperError, double relativeEnergyError)
			: b2max(impact2max), absoluteStepperError(absoluteStepperError), relativeStepperError(
					relativeStepperError), relativeEnergyError(relativeEnergyError) {

		projectileVelocity = Utils::calculateAcceleratedVelocityInAU(Atom::protonMass, 1.0, energykeV);

		projectile = bbsystem.createBody(Atom::protonMass);
		helium = new KirschbaumWiletsAtom(&bbsystem, Element::He, 4.00260325);
		condition = new DistanceCondition(projectile, helium->getNucleus(), 51.0);
		printer = nullptr;

		coulombProjectile1s1 = new CoulombInteraction(-1.0, projectile, helium->getElectron("1s1"));
		coulombProjectile1s2 = new CoulombInteraction(-1.0, projectile, helium->getElectron("1s2"));
		coulombProjectileNucleus = new CoulombInteraction(helium->getNucleusCharge(), projectile,
				helium->getNucleus());

		heisenbergProjectile1s1 = new HeisenbergInteraction(5.0, 1.0, projectile, helium->getElectron("1s1"));
		heisenbergProjectile1s2 = new HeisenbergInteraction(5.0, 1.0, projectile, helium->getElectron("1s2"));

		bbsystem.addInteraction(coulombProjectile1s1);
		bbsystem.addInteraction(coulombProjectile1s2);
		bbsystem.addInteraction(coulombProjectileNucleus);
		bbsystem.addInteraction(heisenbergProjectile1s1);
		bbsystem.addInteraction(heisenbergProjectile1s2);
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

		helium->randomize(randomEngine);
		helium->setPosition(vector3D(0, 0, 0));
		helium->setVelocity(vector3D(0, 0, 0));

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
			simulator.setPrinter(*printer);
		}

		if (condition->evaluate(bbsystem.phase, 0)) {
			stream << "\t" << "Initial condition error" << endl;
			return -3;
		}

		double energy = bbsystem.getSystemEnergy();
		double time = simulator.simulate(0.0, 1.0, 0.0001, *condition, 100);
		bool e1s1BoundToTarget, e1s2BoundToTarget;
		bool e1s1BoundToProjec, e1s2BoundToProjec;

		while (true) {

			if (time < 0.0) {
				stream << "\t" << "Distance not reached error" << endl;
				return -1;
			}

			if (abs((energy - bbsystem.getSystemEnergy()) / energy) > relativeEnergyError) {
				stream << "\t" << "Energy error: " << energy << " vs. " << bbsystem.getSystemEnergy() << endl;
				return -2;
			}

			e1s1BoundToTarget = Utils::isBound(bbsystem, helium->getElectron("1s1"), helium->getNucleus());
			e1s2BoundToTarget = Utils::isBound(bbsystem, helium->getElectron("1s2"), helium->getNucleus());
			e1s1BoundToProjec = Utils::isBound(bbsystem, helium->getElectron("1s1"), projectile);
			e1s2BoundToProjec = Utils::isBound(bbsystem, helium->getElectron("1s2"), projectile);

			if ((!e1s1BoundToTarget || !e1s1BoundToProjec) && (!e1s2BoundToTarget || !e1s2BoundToProjec)) {
				break;
			}

			stream << " " << round << " Extend Run";
			simulator.simulate(time, time + 1.0, 0.0001);
			time += 1.0;
			extended++;
		}

		string bindings;
		for (bool bound : { e1s1BoundToTarget, e1s2BoundToTarget, e1s1BoundToProjec, e1s2BoundToProjec }) {
			if (bound)
				bindings.append("-");
			else
				bindings.append("+");
		}

		switch (Utils::hash(bindings.c_str())) {
		case Utils::hash("--++"):
			stream << endl;
			break;
		case Utils::hash("+-++"):
		case Utils::hash("-+++"):
			stream << "\t" << round << " --> Single Ionization" << endl;
			ionization1++;
			break;
		case Utils::hash("++++"):
			stream << "\t" << round << " --> Dual Ionization" << endl;
			ionization2++;
			break;
		case Utils::hash("+--+"):
		case Utils::hash("-++-"):
			stream << "\t" << round << " --> Single electron Capture" << endl;
			ecapture1++;
			break;
		case Utils::hash("++--"):
			stream << "\t" << round << " --> Dual electron Capture" << endl;
			ecapture2++;
			break;
		case Utils::hash("++-+"):
		case Utils::hash("+++-"):
			stream << "\t" << round << " --> Ionization And Capture" << endl;
			ionizationAndCapture++;
			break;
		default:
			throw new std::logic_error("Unhandled energy configuration.");
		}

		if (printer != nullptr)
			delete printer;

		stream.flush();
		return 0;
	}

	int close(int successfulRounds) {
		double ionizationRate1 = ((double) ionization1) / ((double) successfulRounds);
		double ionizationRate2 = ((double) ionization2) / ((double) successfulRounds);
		double ecaptureRate1 = ((double) ecapture1) / ((double) successfulRounds);
		double ecaptureRate2 = ((double) ecapture2) / ((double) successfulRounds);
		double ionAndEcaptRate = ((double) ionizationAndCapture) / ((double) successfulRounds);
		double extendedRate = ((double) extended) / ((double) successfulRounds);
		cout << "Single ionization: " << ionization1 << " (" << ionizationRate1 * 100.0 << " %)" << endl;
		cout << "Dual ionization  : " << ionization2 << " (" << ionizationRate2 * 100.0 << " %)" << endl;
		cout << "Single el.Capture: " << ecapture1 << " (" << ecaptureRate1 * 100.0 << " %)" << endl;
		cout << "Dual el.Capture  : " << ecapture2 << " (" << ecaptureRate2 * 100.0 << " %)" << endl;
		cout << "Ionizat & Capture: " << ionizationAndCapture << " (" << ionAndEcaptRate * 100.0 << " %)"
				<< endl;
		cout << "Extended run: " << extended << " (" << extendedRate * 100.0 << " %)" << endl << endl;
		cout << "Cross sections:" << endl;

		cout << "\t Single ionization: " << ionizationRate1 * M_PI * b2max * 0.28003 << endl;
		cout << "\t Dual ionization  : " << ionizationRate2 * M_PI * b2max * 0.28003 << endl;
		cout << "\t Single el.Capture: " << ecaptureRate1 * M_PI * b2max * 0.28003 << endl;
		cout << "\t Dual el.Capture  : " << ecaptureRate2 * M_PI * b2max * 0.28003 << endl;
		cout << "\t Ionizat & Capture: " << ionAndEcaptRate * M_PI * b2max * 0.28003 << endl << endl;

		stream.close();
		return 0;
	}
};

#endif /* COLLISION_HE_PROTON_HPP */
