//
// RandOp2d
//
// This operator initializes the values of a GridFunction2d instance
// with random values that are uniformly distributed in [-1,1].
//
// The 64 bit random number generator used is "Mersenne Twister"
// generator mt19937_64 available in <random>.
//
// Use of this header requires a compiler that supports c++11
// constructs.
//
// Initial Construction : Oct.  7, 2015
// Latest  Revision     : Nov. 26, 2015
/*
#############################################################################
#
# Copyright  2015 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/
#include <random>
using namespace std;

#include "../GridFunctionNd/SCC_GridFunction2d.h"
#include "../DoubleVectorNd/SCC_DoubleVector2d.h"

#ifndef _SCC_RandOp2d_
#define _SCC_RandOp2d_

#define RandomGridFunction2dOp_DEFAULT_SEED 3141592

namespace SCC
{
class RandOp2d
{
    public :

	RandOp2d()
    {
		// Initialize the random generator with the default seed

	    seed = RandomGridFunction2dOp_DEFAULT_SEED;
	    randomGenerator.seed(seed);

	    // Initialize the distribution to be uniform in the interval [-1,1]

	    uniform_real_distribution<double>::param_type distParams(-1.0,1.0);
	    distribution.param(distParams);
    }

    virtual ~RandOp2d()
    {}

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

    virtual void randomize(GridFunction2d& v)
    {
        long i;
        double* valuesPtr = v.getDataPointer();
        for(i = 0; i < v.getDimension(); i++)
        {
        valuesPtr[i] =  distribution(randomGenerator);
        }
    }

    virtual void randomize(DoubleVector2d& v)
    {
        long i;
        double* valuesPtr = v.getDataPointer();
        for(i = 0; i < v.getDimension(); i++)
        {
        valuesPtr[i] =  distribution(randomGenerator);
        }
    }

    int                                            seed;
    mt19937_64                          randomGenerator;
	uniform_real_distribution<double>      distribution;
};
 
}
#endif
