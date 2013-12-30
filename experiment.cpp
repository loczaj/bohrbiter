#include <boost/program_options.hpp>
#include <vector>

#include "experiment.hpp"
#include "experiments/collision-h-proton.hpp"
#include "experiments/collision-he-proton.hpp"
#include "experiments/helium.hpp"
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
	    ("voltage,v", po::value<double>(), "Accelerator voltage [kV]")
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

	double voltage = 0;
	if (vm.count("voltage"))
		voltage = vm["voltage"].as<double>();

	Experiment* experiment = nullptr;
	if (vm.count("name")) {

		if (vm["name"].as<string>() == "sb") {
			std::cout << "Carry out sandbox experiment." << std::endl;
			experiment = new SandboxExperiment();

		} else if (vm["name"].as<string>() == "p+H") {
			std::cout << "Carry out proton + hidrogen collision experiment." << std::endl;
			experiment = new CollisionAbrinesPercivalHydrogenWithProton(b2max, voltage, 1e-9, 1e-9, 1e-6);

		} else if (vm["name"].as<string>() == "He") {
			std::cout << "Carry out helium experiment." << std::endl;
			experiment = new HeliumExperiment();

		} else if (vm["name"].as<string>() == "p+He") {
			std::cout << "Carry out proton + helium collision experiment." << std::endl;
			experiment = new CollisionHeliumWithProton();
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
