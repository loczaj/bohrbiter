#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../ap-hydrogen.hpp"

using namespace simulbody;
using namespace std;

int main(int argc, char* argv[]) {
	System bbsystem;
	AbrinesPercivalHydrogen hydrogen(bbsystem);

	std::mt19937_64 randomGenerator;
	randomGenerator.seed(std::mt19937_64::default_seed);

	cout.precision(10);
	cout << hydrogen.getNucleus() << endl;
	cout << hydrogen.getElectron("1s1") << endl;
	cout << hydrogen.getEnergy(bbsystem) << endl;
	cout << hydrogen.getMass(bbsystem) << endl;
	cout << hydrogen.getPosition(bbsystem) << endl;
	cout << hydrogen.getVelocity(bbsystem) << endl;
	cout << hydrogen.getImpulse(bbsystem) << endl;
	cout << hydrogen.getOrbitalEnergy(bbsystem, "1s1") << endl;
	cout << hydrogen.getOrbitalAngularMomentum(bbsystem, "1s1") << endl << endl;

	cout << hydrogen.getOrbitalEnergy(bbsystem, "1s1") << "\t" << hydrogen.getEnergy(bbsystem) << endl;

	for (int i = 0; i < 10; i++) {
		hydrogen.randomize(bbsystem, randomGenerator);
		cout << hydrogen.getOrbitalEnergy(bbsystem, "1s1") << "\t" << hydrogen.getEnergy(bbsystem) << endl;
	}

	return 0;
}
