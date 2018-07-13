#include <random>
#include <ctime>
using namespace std;

#ifndef _UniformRandomNumbers_
#define _UniformRandomNumbers_

#define _DEFAULT_RANDOM_SEED_ 3141592
class UniformRandomNumbers
{
public:

	UniformRandomNumbers()
	{
    seed = _DEFAULT_RANDOM_SEED_;
	randomGenerator.seed(seed);

	// Initialize the distribution to be uniform in the interval [-1,1]

	uniform_real_distribution<double>::param_type distParams(-1.0,1.0);
	distribution.param(distParams);
	}

	void initialize()
	{
	seed = _DEFAULT_RANDOM_SEED_;
	randomGenerator.seed(seed);

	// Initialize the distribution to be uniform in the interval [-1,1]

	uniform_real_distribution<double>::param_type distParams(-1.0,1.0);
	distribution.param(distParams);
	}

	void resetSeed(int seed)
    {
    this->seed = seed;
    randomGenerator.seed(seed);
    }

    void resetWithRandomSeed()
    {
    this->seed = (int)time(0);
    resetSeed(seed);
    }

    double getNextRandomNumber()
    {
    	return distribution(randomGenerator);
    }
	int                                            seed;
    mt19937_64                          randomGenerator;
	uniform_real_distribution<double>      distribution;
};

#undef _DEFAULT_RANDOM_SEED_
#endif
