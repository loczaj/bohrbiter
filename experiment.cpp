#include "experiment.hpp"
#include "experiments/collision-aph-proton.hpp"

int main(int argc, char* argv[]) {

	CollisionAbrinesPercivalHydrogenWithPoton experiment;
	return experiment.carryOut(experiment, 1);
}
