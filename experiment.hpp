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

	int track(vector<int> roundsToTrack) {
		int rounds = *max_element(roundsToTrack.begin(), roundsToTrack.end());
		return carryOut(rounds, false, roundsToTrack, true);
	}

	int carryOut(int numberOfRounds, bool seedRandom = false, vector<int> roundsToTrack = { },
			bool skipUntracked = false) {

		if (seedRandom) {
			random_device rdev { };
			this->randomEngine.seed(rdev());
		} else {
			this->randomEngine.seed(mt19937_64::default_seed);
		}

		int result = this->open(numberOfRounds, seedRandom);
		if (result != 0) {
			this->close(0);
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

			result = this->run(round + 1, tracking, skipUntracked);
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
		result = this->close(successfulRounds);

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

	static constexpr unsigned int hash(const char *str, int offset = 0) {
		return !str[offset] ? 5381 : (hash(str, offset + 1) * 33) ^ str[offset];
	}
};

#endif /* EXPERIMENT_HPP */
