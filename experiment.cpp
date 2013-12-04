#include "experiment.hpp"
#include "experiments/sandbox.hpp"

int main(int argc, char* argv[]) {

	SandboxExperiment experiment;
	return experiment.carryOut(experiment, 1);
}
