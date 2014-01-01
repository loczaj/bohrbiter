#ifndef COLLISION_HE_PROTON_HPP
#define COLLISION_HE_PROTON_HPP

#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../kirschbaum-wilets.hpp"
#include "../experiment.hpp"

class CollisionHeliumWithProton: public Experiment {

	std::ofstream stream;

	System bbsystem;
	PositionPrintField printField;

	KirschbaumWiletsAtom *helium = nullptr;

public:

	int open(int numberOfRounds, bool seedRandom) override {
		helium = new KirschbaumWiletsAtom(&bbsystem, Element::He, 4.00260325);
		return 0;
	}

	int run(int round, bool tracking, bool skipUntracked) override {
		runge_kutta_dopri5<Phase> stepper;
		auto ctrstepper = make_controlled(1e-6, 1e-6, stepper);
		Simulator<decltype(ctrstepper)> simulator(ctrstepper, &bbsystem);

		//helium->randomize(randomEngine);
		helium->install();

		Printer printer("helium-" + std::to_string(round) + ".csv");
		printer.addField(new PositionPrintField());
		simulator.setPrinter(printer);

		cout << "Electron config: " << static_cast<int>(helium->getElectronConfiguration()) << endl;
		cout << "Nucleus element: " << static_cast<int>(helium->getNucleusElement()) << endl;
		cout << "Nucleus charge: " << helium->getNucleusCharge() << endl;
		cout << "Nucleus mass: " << helium->getNucleusMass() << endl;
		cout << "Reduced mass: " << helium->getReducedMass() << endl;
		cout << "Number of bodies: " << helium->getBodies().size() << endl;
		cout << "Number of inters: " << helium->getInteractions().size() << endl;
		cout << "Orbits: ";
		for (string orbit : helium->getOrbitNames()) {
			cout << orbit << " ";
		}
		cout << endl << endl;

		cout << "Energy: " << bbsystem.getSystemEnergy() << endl;
		cout << "  1s1 : " << helium->getOrbitalEnergy("1s1") << endl;
		cout << "  1s2 : " << helium->getOrbitalEnergy("1s2") << endl << endl;

		double energy = bbsystem.getSystemEnergy();
		double time = simulator.simulate(0.0, 10000.0, 0.0001);

		cout << "Energy: " << bbsystem.getSystemEnergy() << endl;
		cout << "  1s1 : " << helium->getOrbitalEnergy("1s1") << endl;
		cout << "  1s2 : " << helium->getOrbitalEnergy("1s2") << endl;

		if (time < 0.0)
			return -1;

		if (abs((energy - bbsystem.getSystemEnergy()) / energy) > 1e-5)
			return -2;

		return 0;
	}

	int close(int successfulRounds) override {
		return 0;
	}

};

#endif /* COLLISION_HE_PROTON_HPP */

