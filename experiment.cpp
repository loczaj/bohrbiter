#include "experiment.hpp"
#include "experiments/collision-h-proton.hpp"
#include "experiments/collision-he-proton.hpp"
#include "experiments/sandbox.hpp"

int main(int argc, char* argv[]) {

	// SandboxExperiment experiment;
	CollisionAbrinesPercivalHydrogenWithProton experiment(15, 50, 1e-9, 1e-9, 1e-7);
	// CollisionHeliumWithProton experiment;
	return experiment.carryOut(experiment, 5000);
}
