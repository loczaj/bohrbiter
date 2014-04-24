#ifndef HELIUM_KW_HPP
#define HELIUM_KW_HPP

#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../kirschbaum-wilets.hpp"
#include "../experiment.hpp"

class KirschbaumWiletsHeliumExperiment: public Experiment {

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
		auto ctrstepper = make_controlled(1e-8, 1e-8, stepper);
		Simulator<decltype(ctrstepper)> simulator(ctrstepper, &bbsystem);

		//helium->randomize(randomEngine);
		helium->install();

		Printer printer("helium-" + std::to_string(round) + ".csv");
		printer.addField(new PositionPrintField());
//		printer.addField(new TimePrintField());
//		printer.addField(new BodyPrintField(1, { Coord::x, Coord::y, Coord::z, Coord::vx, Coord::vy, Coord::vz }));
//		printer.addField(new InteractionPrintField(helium->getInteractions()[1], InteractionAttribute::actingForce));
//		printer.addField(new InteractionPrintField(helium->getInteractions()[2], InteractionAttribute::actingForce));
//		printer.addField(new InteractionPrintField(helium->getInteractions()[2], InteractionAttribute::actingVelocity));
		simulator.setObserver(printer);

		cout << "Electron config: " << static_cast<int>(helium->getElectronConfiguration()) << endl;
		cout << "Nucleus element: " << static_cast<int>(helium->getNucleusElement()) << endl;
		cout << "Nucleus charge: " << helium->getNucleusCharge() << endl;
		cout << "Nucleus mass: " << helium->getNucleusMass() << endl;
		cout << "Reduced mass: " << helium->getReducedMass() << endl;
		cout << "Energy: " << helium->getEnergy() << endl;
		cout << "Number of bodies: " << helium->getBodies().size() << endl;
		cout << "Number of inters: " << helium->getInteractions().size() << endl;
		cout << "Orbits: " << endl;
		for (string orbit : helium->getOrbitNames()) {
			cout << orbit << " (ion. en.: " << helium->getIonizationEnergy(orbit) << ")" << endl;
		}
		cout << endl;

		cout << "Energy: " << bbsystem.getSystemEnergy() << endl;
		cout << "  1s1 : " << helium->getOrbitalEnergy("1s1") << endl;
		cout << "  1s2 : " << helium->getOrbitalEnergy("1s2") << endl << endl;

		cout << "Mass: " << endl;
		cout << "  nuc : " << bbsystem.getBodyMass(helium->getNucleus()) << endl;
		cout << "  1s1 : " << bbsystem.getBodyMass(helium->getElectron("1s1")) << endl;
		cout << "  1s2 : " << bbsystem.getBodyMass(helium->getElectron("1s2")) << endl << endl;

		double energy = bbsystem.getSystemEnergy();
		double time = simulator.simulate(0.0, 150.0, 0.0001);

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

#endif /* HELIUM_KW_HPP */

