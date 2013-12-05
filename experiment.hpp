#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <cmath>
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
			return result;
		}

		for (int index = 0; index < rounds; index++) {
			int result = experiment.run(index);
			if (result != 0) {
				experiment.close();
				return result;
			}
		}

		return experiment.close();
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
