#include "experiment.hpp"
#include "experiments/collision-aph-proton.hpp"

int main(int argc, char* argv[]) {

	CollisionAbrinesPercivalHydrogenWithProton experiment(0.5, 50);
	return experiment.carryOut(experiment, 100);
}
