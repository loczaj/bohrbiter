#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <cmath>
#include <iostream>
#include <random>
#include <simulbody/simulator.hpp>

using namespace simulbody;
using namespace std;

class Experiment {
protected:

	mt19937_64 randomEngine;

public:

	virtual int open(int numberOfRounds, bool seedRandom) = 0;
	virtual int run(int round, bool tracking, bool skipUntracked) = 0;
	virtual int close(int successfulRounds) = 0;

	static int track(Experiment &experiment, initializer_list<int> roundsToTrack) {
		int rounds = *max_element(roundsToTrack.begin(), roundsToTrack.end());
		return carryOut(experiment, rounds, false, roundsToTrack, true);
	}

	static int carryOut(Experiment &experiment, int numberOfRounds, bool seedRandom = false,
			initializer_list<int> roundsToTrack = { }, bool skipUntracked = false) {

		if (seedRandom) {
			random_device rdev { };
			experiment.randomEngine.seed(rdev());
		} else {
			experiment.randomEngine.seed(mt19937_64::default_seed);
		}

		int result = experiment.open(numberOfRounds, seedRandom);
		if (result != 0) {
			experiment.close(0);
			cout << "Failed to open experiment. (" << result << ")" << endl;
			return result;
		}

		int round = 0;
		int successfulRounds = 0;
		int displayed = 0;
		int star = 1;

		for (round = 0; round < numberOfRounds; round++) {

			bool tracking = true;
			if (find(roundsToTrack.begin(), roundsToTrack.end(), (round + 1)) == roundsToTrack.end())
				tracking = false;

			result = experiment.run(round + 1, tracking, skipUntracked);
			if (result != 0) {
				cout << endl << "Round " << (round + 1) << " failed with: " << result << " ";
			} else {
				successfulRounds++;
			}

			while (displayed < 100 * (round + 1)) {
				if (star % 10 == 0)
					cout << star / 10;
				else
					cout << "*";

				cout.flush();
				displayed += numberOfRounds;
				star++;
			}
		}

		cout << endl;
		result = experiment.close(successfulRounds);

		cout << successfulRounds << " rounds passed." << endl;

		if (numberOfRounds != successfulRounds)
			cout << numberOfRounds - successfulRounds << " rounds failed." << endl;

		if (result != 0)
			cout << "Failed to close experiment. (" << result << ")" << endl;
		else
			cout << "Experiment completed." << endl;

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
