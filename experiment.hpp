#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <cmath>

class Experiment {
public:

	virtual int open(int quantity) = 0;
	virtual int run(int index) = 0;
	virtual int close() = 0;

	static int carryOut(Experiment &experiment, int quantity) {
		int result = experiment.open(quantity);
		if (result != 0)
			return result;

		for (int index = 0; index < quantity; index++) {
			int result = experiment.run(index);
			if (result != 0)
				return result;
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
};

#endif /* EXPERIMENT_HPP */