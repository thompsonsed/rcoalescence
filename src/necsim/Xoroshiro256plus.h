//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Sam Thompson, based on a functionally equivalent example from David Blackman and Sebastiano Vigna
 * @file Xoroshiro256plus.h
 * @brief Contains a generic random number generator based on the Xoroshi256+ method.
 *
 * Modified from code by David Blackman and Sebastiano Vigna (vigna@acm.org)
 *
 * To the extent possible under law, the author has dedicated all copyright and related and neighboring rights to this
 * software to the public domain worldwide. This software is distributed without any warranty.
 * See <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 * Samuel Thompson (thomsonsed@gmail.com).
 *
 * The definitions for the constants defined here should not be altered.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef NECSIM_XOROSHIRO256PLUS_H
#define NECSIM_XOROSHIRO256PLUS_H

#include <iostream>
#include <cstdint>
#include <array>

/**
 * @brief Rotate the value a specified number of times.
 * @param value the value to rotate
 * @param bits the number of bits to shift
 * @return the rotated and shifted value
 */
static inline uint64_t rotl(const uint64_t &value, int bits)
{
    return (value << bits) | (value >> (64 - bits));
}

/**
 * @brief Converts a uniform integer on [0, UINT64_MAX) to a double
 * @param i the integer to convert
 * @return double equivalent of the integer
 */
static inline double intToDouble(const uint64_t &i)
{
    union Conversion
    {
        uint64_t i;
        double d;
    };
    Conversion conversion{};
    conversion.i = UINT64_C(0x3FF) << 52 | i >> 12;
    return conversion.d - 1.0;
}

/**
 * @brief A random number generator using the splitmix64 algorithm - this is provided for generating the shuffle table
 * within the main Xoroshiro256+ algorithm.
 */
class SplitMix64
{
private:
    uint64_t x; /* The state can be seeded with any value. */

public:

    /**
     * @brief Default constructor, taking the RNG seed.
     * @param seed the seed to use
     */
    explicit SplitMix64(uint64_t seed) : x(seed)
    { }

    /**
     * @brief Generates the next random integer
     * @return a random integer in the range of 0 to 2^64
     */
    uint64_t next()
    {
        uint64_t z = (x += UINT64_C(0x9E3779B97F4A7C15));
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
        return z ^ (z >> 31);
    }

    /**
     * @brief Shuffle the random number 8 times.
     */
    void shuffle()
    {
        for(unsigned int i = 0; i < 8; i++)
        {
            next();
        }
    }

    /**
     * @brief Gets a random double in the range [0, 1)
     * @return a random double
     */
    double d01()
    {
        return intToDouble(next());
    }
};

/**
 * @brief A random number generator using the xoroshiro
 *
 * @note From the authors of the algorthim:
 * This is xoshiro256+ 1.0, our best and fastest generator for floating-point numbers. We suggest to use its upper bits
 * for floating-point generation, as it is slightly faster than xoshiro256**. It passes all tests we are aware of except
 * for the lowest three bits, which might fail linearity tests (and just those), so if low linear complexity is not
 * considered an issue (as it is usually the case) it can be used to generate 64-bit outputs, too. We suggest to use a
 * sign test to extract a random Boolean value, and right shifts to extract subsets of bits. The state must be seeded so
 * that it is not everywhere zero. If you have a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 * output to fill the shuffle table.
 */
class Xoroshiro256plus
{
protected:
    std::array<uint64_t, 4> shuffle_table;
public:

    Xoroshiro256plus() : shuffle_table()
    {

    }

    virtual ~Xoroshiro256plus() = default;

    /**
     * @brief Explicit constructor which sets the rng seed.
     * @param seed the random seed
     */
    explicit Xoroshiro256plus(uint64_t seed) : shuffle_table()
    {
        setSeed(seed);
    }

    virtual void setSeed(uint64_t seed)
    {
        SplitMix64 seed_generator(seed);
        // Shuffle the seed generator 8 times
        seed_generator.shuffle();
        for(auto &item : shuffle_table)
        {
            item = seed_generator.next();
        }
    }

    /**
     * @brief Generates the next random integer.
     * @return a random integer from 0 to max of 2^64
     */
    uint64_t next()
    {
        const uint64_t result_plus = shuffle_table[0] + shuffle_table[3];

        const uint64_t t = shuffle_table[1] << 17;

        shuffle_table[2] ^= shuffle_table[0];
        shuffle_table[3] ^= shuffle_table[1];
        shuffle_table[1] ^= shuffle_table[2];
        shuffle_table[0] ^= shuffle_table[3];

        shuffle_table[2] ^= t;

        shuffle_table[3] = rotl(shuffle_table[3], 45);

        return result_plus;

    }

    /**
     * @brief Generates a random number in the range [0, 1)
     * @return a random double
     */
    double d01()
    {
        return intToDouble(next());
    }

    /**
     * @brief Jumps the generator forwards by the equivalent of 2^128 calls of next() - useful for parallel
     * computations where different random number sequences are required.
     */
    void jump()
    {
        static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for(unsigned long i : JUMP)
            for(int b = 0; b < 64; b++)
            {
                if(i & UINT64_C(1) << b)
                {
                    s0 ^= shuffle_table[0];
                    s1 ^= shuffle_table[1];
                    s2 ^= shuffle_table[2];
                    s3 ^= shuffle_table[3];
                }
                next();
            }

        shuffle_table[0] = s0;
        shuffle_table[1] = s1;
        shuffle_table[2] = s2;
        shuffle_table[3] = s3;
    }

    /**
     * @brief Writes the random number generator to an output stream.
     * @param os the output stream to write to
     * @param x the random number generator to save
     * @return
     */
    friend std::ostream &operator<<(std::ostream &os, const Xoroshiro256plus &x)
    {
        for(const auto &item: x.shuffle_table)
        {
            os << item << ",";
        }
        return os;
    }

    /**
     * @brief Inputs the NRrand object from the input stream.
     * Used for reading the NRrand object from a file.
     * @param is the input stream.
     * @param r the NRrand object to input to.
     * @return the input stream.
     */
    friend std::istream &operator>>(std::istream &is, Xoroshiro256plus &x)
    {
        char delim;
        for(auto &i : x.shuffle_table)
        {
            is >> i >> delim;
        }
        return is;
    }
};

#endif //NECSIM_XOROSHIRO256PLUS_H
