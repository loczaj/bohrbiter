#include <boost/program_options.hpp>
#include <vector>

#include "experiment.hpp"
#include "experiments/collision-h-proton.hpp"
#include "experiments/collision-he-proton.hpp"
#include "experiments/helium-ap.hpp"
#include "experiments/helium-kw.hpp"
#include "experiments/sandbox.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[]) {

	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help,h", "Produce this help message")
	    ("name,n", po::value<std::string>(), "Experiment to carry out")
	    ("random,r", "Use real random numbers")
	    ("track,t", po::value<std::vector<int>>(), "MC iteration to track")
	    ("iterations,i", po::value<int>(), "Number of MC iterations to do")
	    ("b2max,b", po::value<double>(), "Maximal impact parameter square [au]")
	    ("energy,e", po::value<double>(), "Projectile energy [keV]")
	;

	po::positional_options_description p;
	p.add("name", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}

	int iterations = 1;
	if (vm.count("iterations"))
		iterations = vm["iterations"].as<int>();

	double b2max = 0;
	if (vm.count("b2max"))
		b2max = vm["b2max"].as<double>();

	double energy = 0;
	if (vm.count("energy"))
		energy = vm["energy"].as<double>();

	Experiment* experiment = nullptr;
	if (vm.count("name")) {

		if (vm["name"].as<string>() == "sb") {
			std::cout << "Carry out sandbox experiment." << std::endl;
			experiment = new SandboxExperiment();

		} else if (vm["name"].as<string>() == "p+H") {
			std::cout << "Carry out proton + hidrogen collision experiment." << std::endl;
			experiment = new CollisionAbrinesPercivalHydrogenWithProton(b2max, energy, 1e-9, 1e-9, 1e-6);

		} else if (vm["name"].as<string>() == "p+He") {
			std::cout << "Carry out proton + helium collision experiment." << std::endl;
			experiment = new CollisionKirschbaumWiletsHeliumWithProton(b2max, energy, 1e-8, 1e-8, 1e-6);

		} else if (vm["name"].as<string>() == "apHe") {
			std::cout << "Carry out Abrines-Percival helium experiment." << std::endl;
			experiment = new AbrinesPercivalHeliumExperiment();

		} else if (vm["name"].as<string>() == "kwHe") {
			std::cout << "Carry out Kirschbaum-Wilets helium experiment." << std::endl;
			experiment = new KirschbaumWiletsHeliumExperiment();
		}
	}

	if (experiment == nullptr) {
		std::cout << "No experiment chosen." << std::endl;
		return 1;
	}

	if (vm.count("track")) {
		return experiment->track(vm["track"].as<std::vector<int>>());
	} else {
		return experiment->carryOut(iterations, vm.count("random"));
	}
}

int Experiment::track(vector<int> roundsToTrack) {
	int rounds = *max_element(roundsToTrack.begin(), roundsToTrack.end());
	return carryOut(rounds, false, roundsToTrack, true);
}

int Experiment::carryOut(int numberOfRounds, bool seedRandom, vector<int> roundsToTrack, bool skipUntracked) {

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

Experiment::~Experiment() {
}

// static
double Utils::calculateAcceleratedVelocityInAU(double massAU, double chargeAU, double voltageKV) {
	double voltageAU = 1000 * voltageKV / 27.211383;
	return sqrt(2 * chargeAU * voltageAU / massAU);
}

// static
bool Utils::isBound(System system, identifier body, identifier reference) {
	double energy = system.getBodyKineticEnergyReferenced(body, reference);
	energy += system.getPairPotentialEnergy(body, reference);
	return (energy < 0.0);
}
