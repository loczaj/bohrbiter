#include <fstream>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <simulbody/simulator.hpp>

#include "../aph-atom.hpp"
#include "../experiment.hpp"

class CollisionAbrinesPercivalHydrogenWithProton: public Experiment {

	std::ofstream stream;
	std::mt19937_64 randomGenerator;

	System bbsystem;
	DistanceCondition* condition;

	identifier projectile;
	AbrinesPercivalHydrogen* target;
	Interaction* coulombProjectileElectron;
	Interaction* coulombProjectileNucleus;

	double b2max;
	double projectileVelocity;
	std::uniform_real_distribution<double> distributionNullB2Max;

public:

	CollisionAbrinesPercivalHydrogenWithProton(double impact2max, double voltagekV) {
		b2max = impact2max;
		projectileVelocity = Utils::calculateAcceleratedVelocityInAU(Atom::protonMass, 1.0, voltagekV);
		distributionNullB2Max = std::uniform_real_distribution<double>(0, b2max);

		projectile = bbsystem.createBody(Atom::protonMass);
		target = new AbrinesPercivalHydrogen(&bbsystem, Element::H, 1.00782503207);
		condition = new DistanceCondition(projectile, target->getNucleus(), 35.0);

		coulombProjectileElectron = new CoulombInteraction(-1.0, projectile, target->getElectron("1s1"));
		coulombProjectileNucleus = new CoulombInteraction(target->getNucleusCharge(), projectile,
				target->getNucleus());

		bbsystem.addInteraction(coulombProjectileElectron);
		bbsystem.addInteraction(coulombProjectileNucleus);
	}

	int open(int rounds) {
		stream.open("result.csv");
		stream.precision(10);

		randomGenerator.seed(std::mt19937_64::default_seed);
		return 0;
	}

	int run(int index, bool tracking) {
		runge_kutta_dopri5<Phase> stepper;
		auto ctrdStepper = make_controlled(1e-10, 1e-10, stepper);
		Simulator<decltype(ctrdStepper)> simulator(ctrdStepper, &bbsystem);

		double b = sqrt(distributionNullB2Max(randomGenerator));

		target->randomize(randomGenerator);
		target->setPosition(vector3D(0, 0, 0));
		target->setVelocity(vector3D(0, 0, 0));

		bbsystem.setBodyPosition(projectile, vector3D(0, b, -20));
		bbsystem.setBodyVelocity(projectile, vector3D(0, 0, projectileVelocity));

		stream << bbsystem.phase;
		stream.flush();

		double energy = bbsystem.getSystemEnergy();
		double time = simulator.simulate(0.0, 1.0, 0.0001, *condition, 50);

		if (time < 0.0) {
			stream << "Condition error" << std::endl;
			return -1;
		}

		if (abs(energy - bbsystem.getSystemEnergy()) / energy > 1e-8) {
			stream << "Energy error: " << energy << " vs. " << bbsystem.getSystemEnergy() << std::endl;
			return -2;
		}

		bool eBoundToTarget = Utils::isBound(bbsystem, target->getElectron("1s1"), target->getNucleus());
		bool eBoundToProjec = Utils::isBound(bbsystem, target->getElectron("1s1"), projectile);

		if (eBoundToTarget && eBoundToProjec)
			stream << "\t > Molecule" << std::endl;

		if (eBoundToTarget && !eBoundToProjec)
			stream << std::endl;

		if (!eBoundToTarget && eBoundToProjec)
			stream << "\t > Transfer" << std::endl;

		if (!eBoundToTarget && !eBoundToProjec)
			stream << "\t > Ionization" << std::endl;

		stream.flush();
		return 0;
	}

	int close() {
		stream.close();
		return 0;
	}
};
