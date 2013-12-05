#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <cmath>
#include <iostream>
#include <simulbody/simulator.hpp>

using namespace simulbody;

class Experiment {
public:

	virtual int open(int quantity) = 0;
	virtual int run(int index) = 0;
	virtual int close() = 0;

	static int carryOut(Experiment &experiment, int rounds) {
		int result = experiment.open(rounds);
		if (result != 0) {
			experiment.close();
			std::cout << "Failed to open experiment. (" << result << ")" << std::endl;
			return result;
		}

		int displayed = 0;
		int star = 1;
		for (int index = 0; index < rounds; index++) {
			while (displayed < 100 * (index + 1)) {
				if (star % 10 == 0)
					std::cout << star / 10;
				else
					std::cout << "*";

				std::cout.flush();
				displayed += rounds;
				star++;
			}

			result = experiment.run(index);
			if (result != 0) {
				experiment.close();
				std::cout << std::endl << "Failed to run experiment. (" << result << ")" << std::endl;
				return result;
			}
		}

		result = experiment.close();
		if (result != 0) {
			std::cout << std::endl << "Failed to close experiment. (" << result << ")" << std::endl;
		}

		std::cout << std::endl << "Experiment completed." << std::endl;
		return result;
	}

	virtual ~Experiment() {
	}
};

struct Utils {
	static double calculateAcceleratedVelocityInAU(double massAU, double chargeAU, double voltageKV) {
		double voltageAU = 1000 * voltageKV / 27.211383;
		return sqrt(2 * chargeAU * voltageAU / massAU);
	}

	static double isBound(System system, identifier body, identifier reference) {
		double energy = system.getBodyKineticEnergyReferenced(body, reference);
		energy += system.getPairPotentialEnergy(body, reference);
		return energy < 0.0;
	}
};

#endif /* EXPERIMENT_HPP */
