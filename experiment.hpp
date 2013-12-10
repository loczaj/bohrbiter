#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <cmath>
#include <iostream>
#include <simulbody/simulator.hpp>

using namespace simulbody;

class Experiment {
public:

	virtual int open(int numberOfRounds) = 0;
	virtual int run(int round, bool tracking, bool skipUntracked) = 0;
	virtual int close(int successfulRounds) = 0;

	static int track(Experiment &experiment, std::initializer_list<int> roundsToTrack) {
		int rounds = *std::max_element(roundsToTrack.begin(), roundsToTrack.end());
		return carryOut(experiment, rounds, roundsToTrack, true);
	}

	static int carryOut(Experiment &experiment, int numberOfRounds, std::initializer_list<int> roundsToTrack =
			{ }, bool skipUntracked = false) {
		int result = experiment.open(numberOfRounds);
		if (result != 0) {
			experiment.close(0);
			std::cout << "Failed to open experiment. (" << result << ")" << std::endl;
			return result;
		}

		int round = 0;
		int successfulRounds = 0;
		int displayed = 0;
		int star = 1;

		for (round = 0; round < numberOfRounds; round++) {
			while (displayed < 100 * (round + 1)) {
				if (star % 10 == 0)
					std::cout << star / 10;
				else
					std::cout << "*";

				std::cout.flush();
				displayed += numberOfRounds;
				star++;
			}

			bool tracking = true;
			if (std::find(roundsToTrack.begin(), roundsToTrack.end(), (round + 1)) == roundsToTrack.end())
				tracking = false;

			result = experiment.run(round + 1, tracking, skipUntracked);
			if (result != 0) {
				std::cout << std::endl << "Round " << (round + 1) << " failed with: " << result << " ";
			} else {
				successfulRounds++;
			}
		}

		std::cout << std::endl;
		result = experiment.close(successfulRounds);

		std::cout << successfulRounds << " rounds passed." << std::endl;

		if (numberOfRounds != successfulRounds)
			std::cout << numberOfRounds - successfulRounds << " rounds failed." << std::endl;

		if (result != 0)
			std::cout << "Failed to close experiment. (" << result << ")" << std::endl;
		else
			std::cout << "Experiment completed." << std::endl;

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

	static bool isBound(System system, identifier body, identifier reference) {
		double energy = system.getBodyKineticEnergyReferenced(body, reference);
		energy += system.getPairPotentialEnergy(body, reference);
		return (energy < 0.0);
	}
};

#endif /* EXPERIMENT_HPP */
