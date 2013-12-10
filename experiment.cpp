#include "experiment.hpp"
#include "experiments/collision-aph-proton.hpp"
#include "experiments/sandbox.hpp"

int main(int argc, char* argv[]) {

	// SandboxExperiment experiment;
	CollisionAbrinesPercivalHydrogenWithProton experiment(15, 50, 1e-10, 1e-10, 1e-8);
	return experiment.carryOut(experiment, 10000);
}
