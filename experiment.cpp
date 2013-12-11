#include "experiment.hpp"
#include "experiments/collision-h-proton.hpp"
#include "experiments/sandbox.hpp"

int main(int argc, char* argv[]) {

	// SandboxExperiment experiment;
	CollisionAbrinesPercivalHydrogenWithProton experiment(15, 50, 1e-10, 1e-10, 1e-7);
	// HeliumExperiment experiment;
	return experiment.carryOut(experiment, 10000);
}
