//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Sam Thompson
 * @file RNGController.h
 * @brief Contains a generic random number generator.
 *
 * To the extent possible under law, the author has dedicated all copyright and related and neighboring rights to this
 * software to the public domain worldwide. This software is distributed without any warranty.
 * See <http://creativecommons.org/publicdomain/zero/1.0/>.
 * 
 * The definitions for the constants defined here should not be altered.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef FATTAIL_H
#define FATTAIL_H

#include <cstdio>
#include <string>
#include <iomanip>

#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include "Logging.h"
#include "Xoroshiro256plus.h"

using namespace std;
/**
 * @brief These variables contain the special numbers for random number generation.
 */
#define PARAM_R 3.44428647676

/* tabulated values for the heigt of the Ziggurat levels */
constexpr static double ytab[128] = {
        1, 0.963598623011, 0.936280813353, 0.913041104253,
        0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
        0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
        0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
        0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
        0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
        0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
        0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
        0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
        0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
        0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
        0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
        0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
        0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
        0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
        0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
        0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
        0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
        0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
        0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
        0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
        0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
        0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
        0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
        0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
        0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
        0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
        0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
        0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
        0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
        0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
        0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
constexpr static unsigned long ktab[128] = {
        0, 12590644, 14272653, 14988939,
        15384584, 15635009, 15807561, 15933577,
        16029594, 16105155, 16166147, 16216399,
        16258508, 16294295, 16325078, 16351831,
        16375291, 16396026, 16414479, 16431002,
        16445880, 16459343, 16471578, 16482744,
        16492970, 16502368, 16511031, 16519039,
        16526459, 16533352, 16539769, 16545755,
        16551348, 16556584, 16561493, 16566101,
        16570433, 16574511, 16578353, 16581977,
        16585398, 16588629, 16591685, 16594575,
        16597311, 16599901, 16602354, 16604679,
        16606881, 16608968, 16610945, 16612818,
        16614592, 16616272, 16617861, 16619363,
        16620782, 16622121, 16623383, 16624570,
        16625685, 16626730, 16627708, 16628619,
        16629465, 16630248, 16630969, 16631628,
        16632228, 16632768, 16633248, 16633671,
        16634034, 16634340, 16634586, 16634774,
        16634903, 16634972, 16634980, 16634926,
        16634810, 16634628, 16634381, 16634066,
        16633680, 16633222, 16632688, 16632075,
        16631380, 16630598, 16629726, 16628757,
        16627686, 16626507, 16625212, 16623794,
        16622243, 16620548, 16618698, 16616679,
        16614476, 16612071, 16609444, 16606571,
        16603425, 16599973, 16596178, 16591995,
        16587369, 16582237, 16576520, 16570120,
        16562917, 16554758, 16545450, 16534739,
        16522287, 16507638, 16490152, 16468907,
        16442518, 16408804, 16364095, 16301683,
        16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
constexpr static double wtab[128] = {
        1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
        3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
        3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
        4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
        5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
        5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
        5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
        6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
        6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
        6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
        7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
        7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
        7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
        8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
        8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
        8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
        9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
        9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
        9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
        1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
        1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
        1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
        1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
        1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
        1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
        1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
        1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
        1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
        1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
        1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
        1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
        1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

/**
 * @brief Contains the functions for random number generation, based on the Xoroshiro256+ algorithm.
 */
class RNGController : public virtual Xoroshiro256plus
{

private:
    bool seeded;
    uint64_t seed;
    // for the L value of the dispersal kernel (the width - does not affect the shape).
    double tau;
    // for the sigma value of the dispersal kernel (the variance of a normal distribution).
    double sigma;

    typedef double (RNGController::*fptr)(); // once setup will contain the dispersal function to use for this simulation.
    fptr dispersalFunction;

    // once setup will contain the dispersal function for the minimum dispersal distance.
    typedef double (RNGController::*fptr2)(const double &min_distance);

    fptr2 dispersalFunctionMinDistance;

    // the probability that dispersal comes from the uniform distribution. This is only relevant for uniform dispersals.
    double m_prob;
    // the cutoff for the uniform dispersal function i.e. the maximum value to be drawn from the uniform distribution.
    double cutoff;
public:

    /**
     * @brief Standard constructor.
     */
    RNGController() : Xoroshiro256plus(), seeded(false), seed(0), tau(0.0), sigma(0.0), dispersalFunction(nullptr),
                      dispersalFunctionMinDistance(nullptr), m_prob(0.0), cutoff(0.0)
    {

    }

    /**
     * @brief Sets the seed to the given input.
     * Is only seeded if the seed hasn't already been provided.
     * @param seed the input seed.
     */
    void setSeed(uint64_t seed) override
    {
        if(!seeded)
        {
            Xoroshiro256plus::setSeed(seed);
            this->seed = seed;
            seeded = true;
        }
        else
        {
            throw runtime_error("Trying to set the seed again: this can only be set once.");
        }
    }

    /**
     * @brief Clears the seed, if it has already been set.
     * Keeps other simulation parameters, such as sigma and tau.
     */
    void wipeSeed()
    {
        seeded = false;
    }

    /**
     * @brief Generates a random number uniformly from 0 to the maximum value provided.
     * @param max the maximum number.
     * @return an integer of the produced random number.
     */
    unsigned long i0(unsigned long max)
    {
        return (unsigned long) (d01() * (max + 1));
    }

    /**
     * @brief Generates a normally distributed number
     * Uses the standard normal distribution using the Ziggurat method.
     * @return the random number from a normal distribution.
     */
    double norm()
    {
        unsigned long U, sign, i, j;
        double x, y;

        while(true)
        {
            U = i0(UINT32_MAX);
            i = U & 0x0000007F;        /* 7 bit to choose the step */
            sign = U & 0x00000080;    /* 1 bit for the sign */
            j = U >> 8;            /* 24 bit for the x-value */

            x = j * wtab[i];
            if(j < ktab[i])
            {
                break;
            }

            if(i < 127)
            {
                double y0, y1;
                y0 = ytab[i];
                y1 = ytab[i + 1];
                y = y1 + (y0 - y1) * d01();
            }
            else
            {
                x = PARAM_R - log(1.0 - d01()) / PARAM_R;
                y = exp(-PARAM_R * (x - 0.5 * PARAM_R)) * d01();
            }
            if(y < exp(-0.5 * x * x))
            {
                break;
            }
        }
        return sign ? sigma * x : -sigma * x;
    }

    /**
     * @brief Returns a random distance from a 2 dimensional normal distribution, also called the rayleigh distribution.
     * @return dispersal distance of a rayleigh distribution
     */
    double rayleigh()
    {
        return sigma * pow(-2 * log(d01()), 0.5);
    }

    /**
     * @brief Generates a random distance from a rayleigh distribution, given that the distance is more than some
     * minimum.
     * @param dist the minimum distance to generate
     * @return a random distance greater than the minimum provided
     */
    double rayleighMinDist(const double &dist)
    {
        double min_prob = rayleighCDF(dist);
        double rand_prob = min_prob + (1 - min_prob) * d01();
        double out = sigma * pow(-2 * log(rand_prob), 0.5);
        if(out < dist)
        {
            // This probably means that the rayleigh distribution has a less-than-machine-precision probability of
            // producing this distance.
            // Therefore, we just return the distance
            return dist;
        }
        return out;
    }

    /**
     * @brief Gets the cumulative probability of a distance from the rayleigh distribution
     * @param dist the distance to obtain the probability of
     * @return the probability of producing the given distance
     */
    double rayleighCDF(const double &dist)
    {
        return 1 - exp(-pow(dist, 2.0) / (2.0 * pow(sigma, 2.0)));
    }

    /**
     * @brief Sets the dispersal parameters, avoiding requirement to provide these numbers each function call.
     * This is only relevant for fat-tailed dispersal calls.
     * @param sigmain the fatness of the fat-tailed dispersal kernel.
     * @param tauin the width of the fat-tailed dispersal kernel.
     */
    void setDispersalParams(const double sigmain, const double tauin)
    {
        sigma = sigmain;
        tau = tauin; // used to invert the sign here, doesn't any more.
    }

    /**
     * @brief Call from the fat-tailed dispersal kernel with the provided sigma.
     * @deprecated This is the original version used in J Rosindell's codebase, and has been altered for
     * a version which approximates the gaussian distribution at extreme limits.
     * @param z the desired sigma.
     * @return a random number drawn from the fat-tailed dispersal kernel.
     */
    double fattail(double z)
    {
        double result;
        result = pow((pow(d01(), (1.0 / (1.0 - z))) - 1.0), 0.5);
        return result;
    }

    /**
     * @brief Gets the cumulative probability density of travelling the distance
     * @param distance the distance to obtain the cumulative probability for
     * @return the probability of dispersing less than or equal to distance
     */
    double fattailCDF(const double &distance)
    {
        // TODO remove this
        //        return (1.0 / (2.0 * M_PI * sigma * sigma)) *
        //               pow(1 + (distance * distance / (tau * sigma * sigma)), -(tau + 2.0) / 2.0);
        return (1 - pow((1 + ((distance * distance) / (tau * sigma * sigma))), (-tau / 2)));
    }

    /**
     * @brief Gets a fat-tailed random distance greater than some minimum
     * @param min_distance the minimum distance to return
     * @return a fat-tailed distance greater than the minimum
     */
    double fattailMinDistance(const double &min_distance)
    {
        double prob = fattailCDF(min_distance);
        double random_number = 1 - (prob + d01() * (1 - prob));
        double result = sigma * pow((tau * (pow(random_number, -2.0 / tau)) - 1.0), 0.5);
        // This is an approximation for scenarios when the probability of dispersing very long distances
        // is less than machine precision
        if(result < min_distance)
        {
            result = fattail() + min_distance;
        }
        return result;
    }

    // this new version corrects the 1.0 to 2.0 and doesn't require the values to be passed every time.
    /**
     * @brief Call from fat-tailed dispersal kernel.
     * This function requires setDispersalParams() has already been called.
     * @deprecated deprecated, kept for testing purposes only
     * @return a random number drawn from the fat-tailed dispersal kernel.
     */
    double fattail()
    {
        double result;
        // old function version (kept for reference)
        //		result = (tau * pow((pow(d01(),(2.0/(2.0-sigma)))-1.0),0.5));
        result = (sigma * pow((tau * (pow(d01(), -2.0 / tau)) - 1.0), 0.5));
        return result;
    }

    /**
     * @brief Old version of the function call reparameterised for different nu and sigma.
     * @deprecated Kept only for testing purposes.
     * @return a random number drawn from the fat-tailed dispersal kernel.
     */
    double fattail_old()
    {
        double result;
        result = (sigma * pow((pow(d01(), (2.0 / (2.0 + tau))) - 1.0), 0.5));
        return result;
    }

    /**
     * @brief Generates a direction in radians.
     * @return the direction in radians
     */
    double direction()
    {
        return (d01() * 2 * M_PI);
    }

    /**
     * @brief For a given event probability, returns the probability that the event has occured.
     * @param event_probability the event probability.
     * @return whether or not the event has occured.
     */
    bool event(double event_probability)
    {
        if(event_probability < 0.000001)
        {
            if(d01() <= 0.000001)
            {
                return (event(event_probability * 1000000.0));
            }
            return false;
        }
        if(event_probability > 0.999999)
        {
            return (!(event(1.0 - event_probability)));
        }
        return (d01() <= event_probability);

    }

    /**
     * @brief Normal distribution, with percentage chance to choose a uniform distribution instead.
     * @note This function will not produce the same output as norm() for the same parameters, even with a
     * zero chance of picking from the uniform distribution (due to random number draws).
     * @return normally (or uniformly) distributed number
     */
    double normUniform()
    {
        // Check if the dispersal event comes from the uniform distribution
        if(d01() < m_prob)
        {
            // Then it does come from the uniform distribution
            return uniform();
        }
        return rayleigh();
    }

    /**
     * @brief Generates a random distance from a norm-uniform distribution, given that the distance is more than some
     * minimum.
     * @param dist the minimum distance to generate
     * @return a random distance greater than the minimum provided
     */
    double normUniformMinDistance(const double &min_distance)
    {
        if(d01() < m_prob)
        {
            // Then it does come from the uniform distribution
            return uniformMinDistance(min_distance);
        }
        return rayleighMinDist(min_distance);
    }

    /**
     * @brief Draws a random number from a uniform distribution between 0 and cutoff
     * @return a random number in (0, cutoff)
     */
    double uniform()
    {
        return d01() * cutoff;
    }

    /**
     * @brief Generates a random distance from a uniform distribution, given that the distance is more than some
     * minimum.
     * @param dist the minimum distance to generate
     * @return a random distance greater than the minimum provided
     */
    double uniformMinDistance(const double &min_distance)
    {
        if(min_distance > cutoff)
        {
            // Note this may introduce problems for studies of extremely isolated islands
            // I've left this in to make it much easier to deal with scenarios where the
            // disappearing habitat pixel is further from the nearest habitat pixel than
            // the maximum dispersal distance
            return min_distance;
        }
        return min_distance + d01() * (cutoff - min_distance);
    }

    /**
     * @brief Two uniform distributions, the first between 0 and 0.1*cutoff, and the second between 0.9*cutoff and
     * cutoff. Selects from both distributions equally.
     * @note The mean for this function should be identical to a uniform distribution between 0 and cutoff.
     * @return uniformly distributed number
     */
    double uniformUniform()
    {
        if(d01() < 0.5)
        {
            // Then value comes from the first uniform distribution
            return (uniform() * 0.1);
        }
        // Then the value comes from the second uniform distribution
        return 0.9 * cutoff + (uniform() * 0.1);
    }

    /**
     * @brief Generates a random distance from a uniform-uniform distribution, given that the distance is more than some
     * minimum.
     * @param dist the minimum distance to generate
     * @return a random distance greater than the minimum provided
     */
    double uniformUniformMinDistance(const double &min_distance)
    {
        if(d01() < 0.5)
        {
            // Then value comes from the first uniform distribution
            if(min_distance > cutoff * 0.1)
            {
                return uniformMinDistance(min_distance * 10) * 0.1;
            }
        }
        // Then the value comes from the second uniform distribution
        return uniformMinDistance(max(min_distance, 0.9 * cutoff));
    }

    /**
     * @brief Sets the dispersal method by creating the link between dispersalFunction() and the correct
     * dispersal character
     * @param dispersal_method string containing the dispersal type. Can be one of [normal, fat-tail, norm-uniform]
     * @param m_probin the probability of drawing from the uniform distribution. Only relevant for uniform dispersals.
     * @param cutoffin the maximum value to be drawn from the uniform dispersal. Only relevant for uniform dispersals.
     */
    void setDispersalMethod(const string &dispersal_method, const double &m_probin, const double &cutoffin)
    {
        if(dispersal_method == "normal")
        {
            dispersalFunction = &RNGController::rayleigh;
            dispersalFunctionMinDistance = &RNGController::rayleighMinDist;
            if(sigma < 0)
            {
                throw invalid_argument("Cannot have negative sigma with normal dispersal");
            }
        }
        else if(dispersal_method == "fat-tail" || dispersal_method == "fat-tailed")
        {
            dispersalFunction = &RNGController::fattail;
            dispersalFunctionMinDistance = &RNGController::fattailMinDistance;
            if(tau < 0 || sigma < 0)
            {
                throw invalid_argument("Cannot have negative sigma or tau with fat-tailed dispersal");
            }
        }
        else if(dispersal_method == "norm-uniform")
        {
            dispersalFunction = &RNGController::normUniform;
            dispersalFunctionMinDistance = &RNGController::normUniformMinDistance;
            if(sigma < 0)
            {
                throw invalid_argument("Cannot have negative sigma with normal dispersal");
            }
        }
        else if(dispersal_method == "uniform-uniform")
        {
            // This is just here for testing purposes
            dispersalFunction = &RNGController::uniformUniform;
            dispersalFunctionMinDistance = &RNGController::uniformUniformMinDistance;
        }
            // Also provided the old version of the fat-tailed dispersal kernel
        else if(dispersal_method == "fat-tail-old")
        {
            dispersalFunction = &RNGController::fattail_old;
            dispersalFunctionMinDistance = &RNGController::fattailMinDistance;

            if(tau > -2 || sigma < 0)
            {
                throw invalid_argument(
                        "Cannot have sigma < 0 or tau > -2 with fat-tailed dispersal (old implementation).");
            }
        }
        else
        {
            throw runtime_error("Dispersal method not detected. Check implementation exists");
        }
        m_prob = m_probin;
        cutoff = cutoffin;
    }

    /**
     * @brief Runs the dispersal with the allocated dispersal function.
     *
     * @note This function will never return a value larger than the size of LONG_MAX to avoid issues of converting
     * doubles to integers. For dispersal distance within coalescence simulations, this is seemed a reasonable
     * assumption, but may cause issues if code is re-used in later projects.
     *
     * @return distance the dispersal distance
     */
    double dispersal()
    {
        return min(double(LONG_MAX), (this->*dispersalFunction)());
    }

    /**
     * @brief Get a dispersal distance with some minimum
     * @param min_distance the minimum distance to disperse
     * @return the random dispersal distance greater than or equal to the minimum
     */
    double dispersalMinDistance(const double &min_distance)
    {
        return min(double(LONG_MAX), (this->*dispersalFunctionMinDistance)(min_distance));
    }

    /**
     * @brief Sample from a logarithmic distribution
     *
     * Uses the LK sampling method for generating random numbers from a logarithmic distribution, as described by
     * Kemp (1981).
     *
     * @param alpha alpha parameter for the logarithmic distribution
     * @return the randomly generated logarithmic number
     */
    unsigned long randomLogarithmic(long double alpha)
    {
        double u_2 = d01();
        if(u_2 > alpha)
        {
            return 1;
        }
        long double h = log(1 - alpha);
        long double q = 1 - exp(d01() * h);
        if(u_2 < (q * q))
        {
            return static_cast<unsigned long>(floor(1 + log(u_2) / log(q)));
        }
        else if(u_2 > q)
        {
            return 1;
        }
        else
        {
            return 2;
        }

    }

    /**
     * @brief Outputs the NRrand object to the output stream.
     * Used for saving the object to file.
     * @param os the output stream.
     * @param r the NRrand object to output.
     * @return the output stream.
     */
    friend ostream &operator<<(ostream &os, const RNGController &r)
    {
        os << setprecision(64);
        os << r.seed << "," << r.seeded << ",";
        os << r.tau << "," << r.sigma << "," << r.m_prob << "," << r.cutoff << ",";
        os << static_cast<const Xoroshiro256plus &>(r);
        return os;
    }

    /**
     * @brief Inputs the NRrand object from the input stream.
     * Used for reading the NRrand object from a file.
     * @param is the input stream.
     * @param r the NRrand object to input to.
     * @return the input stream.
     */
    friend istream &operator>>(istream &is, RNGController &r)
    {
        char delim;
        is >> r.seed >> delim >> r.seeded >> delim >> r.tau >> delim >> r.sigma >> delim >> r.m_prob >> delim;
        is >> r.cutoff >> delim >> static_cast<Xoroshiro256plus &>(r);
        return is;
    }
};

#endif
