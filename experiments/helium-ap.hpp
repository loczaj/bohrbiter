#ifndef HELIUM_AP_HPP
#define HELIUM_AP_HPP

#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>
#include <simulbody/printer.hpp>


#include "../abrines-percival.hpp"
#include "../experiment.hpp"

class AbrinesPercivalHeliumExperiment: public Experiment {

	std::ofstream stream;

	System bbsystem;
	PositionPrintField printField;

	AbrinesPercivalAtom *helium = nullptr;

public:

	int open(int numberOfRounds, bool seedRandom) {
		helium = new AbrinesPercivalAtom(&bbsystem, Element::He, 4.00260325);
		return 0;
	}

	int run(int round, bool tracking, bool skipUntracked) {
		runge_kutta_dopri5<Phase> stepper;
		auto ctrdStepper = make_controlled(1e-10, 1e-10, stepper);
		Simulator<decltype(ctrdStepper)> simulator(ctrdStepper, &bbsystem);

		helium->randomize(randomEngine);
		//helium->install();

		Printer printer("helium-" + std::to_string(round) + ".csv");
		printer.addField(new PositionPrintField());
		simulator.setObserver(printer);

		double energy = bbsystem.getSystemEnergy();
		double time = simulator.simulate(0.0, 200.0, 0.0001);

		if (time < 0.0)
			return -1;

		if (abs((energy - bbsystem.getSystemEnergy()) / energy) > 1e-6)
			return -2;

		return 0;
	}

	int close(int successfulRounds) {
		return 0;
	}

};

#endif /* HELIUM_AP_HPP */

