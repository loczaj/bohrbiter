#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../ap-h-atom.hpp"

using namespace simulbody;
using namespace std;

int main(int argc, char* argv[]) {
	System bbsystem;
	AbrinesPercivalHydrogen hydrogen(&bbsystem);

	std::mt19937_64 randomGenerator;
	randomGenerator.seed(std::mt19937_64::default_seed);

	cout.precision(10);
	cout << hydrogen.getNucleus() << endl;
	cout << hydrogen.getElectron("1s1") << endl;
	cout << hydrogen.getEnergy() << endl;
	cout << hydrogen.getMass() << endl;
	cout << hydrogen.getPosition() << endl;
	cout << hydrogen.getVelocity() << endl;
	cout << hydrogen.getImpulse() << endl;
	cout << hydrogen.getOrbitalEnergy("1s1") << endl;
	cout << hydrogen.getOrbitalAngularMomentum("1s1") << endl << endl;

	cout << hydrogen.getOrbitalEnergy("1s1") << "\t" << hydrogen.getEnergy() << endl;

	for (int i = 0; i < 10; i++) {
		hydrogen.randomize(randomGenerator);
		cout << hydrogen.getOrbitalEnergy("1s1") << "\t" << hydrogen.getEnergy() << endl;
	}

	return 0;
}
