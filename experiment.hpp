#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

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

	int track(vector<int> roundsToTrack);

	int carryOut(int numberOfRounds, bool seedRandom = false, vector<int> roundsToTrack = { },
			bool skipUntracked = false);

	virtual ~Experiment();
};

struct Utils {

	static double calculateAcceleratedVelocityInAU(double massAU, double chargeAU, double voltageKV);

	static bool isBound(System system, identifier body, identifier reference);

	static constexpr unsigned int hash(const char *str, int offset = 0) {
		return !str[offset] ? 5381 : (hash(str, offset + 1) * 33) ^ str[offset];
	}
};

#endif /* EXPERIMENT_HPP */
