#ifndef ATOM_HPP
#define ATOM_HPP

#include <vector>

#include"simulbody/simulator.hpp"

class Atom {
	sizeT nucleus;
	std::vector<sizeT> electrons;
	std::vector<Interaction> interactions;

public:
	virtual Atom(unsigned int atomicNumber, unsigned int electronNumber, unsigned int massNumber);

	virtual void locate(System &system, vector3D position);
	virtual void speedUp(System &system, vector3D velocity);

	virtual void install(System &system) = 0;
	virtual void setInteractions() = 0;

	virtual double getOrbitalEnergy(System &system, sizeT electron) const;
	virtual double getAngularMomentum(System &system, sizeT electron) const;

	virtual ~Atom();
};

#endif /* ATOM_HPP */
