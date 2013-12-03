#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

class Experiment {
public:

	virtual int open(int quantity) = 0;
	virtual int run(int index) = 0;
	virtual int close() = 0;

	static int carryOut(Experiment &experiment, int quantity) {
		int result = experiment.open(quantity);
		if (result != 0)
			return result;

		for (int index = 0; index < quantity; index++) {
			int result = experiment.run(index);
			if (result != 0)
				return result;
		}

		return experiment.close();
	}

	virtual ~Experiment() {
	}
};

#endif /* EXPERIMENT_HPP */
