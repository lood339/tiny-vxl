// This is core/vnl/vnl_random.h
#ifndef vnl_random_h
#define vnl_random_h

#include <ctime>
#include <cmath>
#include <cassert>

#include "vnl/vnl_export.h"

//:
// \file
// \author Aaron Kotcheff (Manchester)
// \brief A superior random number generator

inline constexpr unsigned int vnl_random_array_size = 37;

//: A superior random number generator.
// Implements a new random number generator that
// recently appeared in the literature. It generates 32 bit
// numbers with a higher degree of randomness than previous
// generators and has a cycle of 10^354 i.e. so huge that in
// practice it never cycles.
// For the mathematics behind it see:
// "A New Class of Random Number Generators" G. Marsaglia and A. Zaman,
// Annals of Applied Probability 1991, Vol. 1, No. 3, 462.
class VNL_EXPORT vnl_random
{
    enum {linear_congruential_multiplier = 1664525, mz_previous1 = 24};
    unsigned long linear_congruential_previous;
    unsigned long mz_seed_array[vnl_random_array_size];
    unsigned long mz_array[vnl_random_array_size];
    unsigned int mz_array_position{0UL};
    int mz_borrow{0};
    unsigned long linear_congruential_lrand32();

    double mz_previous_normal;
    int mz_previous_normal_flag{0};

  public:
    //: Default constructor.
    // Initializes the random number generator non-deterministically.
    // i.e. it will generate a different series of random numbers each
    // time the program is run.
    vnl_random();

    //: Destructor
    ~vnl_random();

    //: Construct with seed.
    //  Initializes the random number generator deterministically
    //  using a single ulong as the 'seed'. A linear congruential
    //  generator is used to generate the 37 ulongs needed
    //  as the real seed. The same seed will produce the
    //  same series of random numbers.
    //
    //  9667566  is a good seed.
    vnl_random(unsigned long seed);

    //: Construct with seed.
    //  Initializes the random number generator deterministically
    //  using 37 ulongs as the 'seed'. The same seed will
    //  produce the same series of random numbers.
    vnl_random(unsigned long seed[vnl_random_array_size]);

    //: Copy constructor.
    //  Initializes/sets the random number generator to exactly
    //  the same state as the argument, i.e. both will generate exactly
    //  the same series of random numbers from then on.
    vnl_random(const vnl_random&);

    //: Copy operator.
    //  Initializes/sets the random number generator to exactly
    //  the same state as the argument, i.e. both will generate exactly
    //  the same series of random numbers from then on.
    vnl_random& operator=(const vnl_random&);

    //: Starts a new non-deterministic sequence from an already declared generator.
    void reseed();

    //: Starts a new deterministic sequence from an already declared generator using the provided seed.
    void reseed(unsigned long);

    //: Starts a new deterministic sequence from an already declared generator using the provided seed.
    void reseed(const unsigned long[vnl_random_array_size]);

    //: This restarts the sequence of random numbers.
    //  Restarts so that it repeats
    //  from the point at which you declared the generator, last
    //  initialized it, or last called a 'reseed'.
    void restart();

    //: Generates a random unsigned 32-bit number.
    unsigned long lrand32();

    //: Generates a random unsigned long in [a,b]
    int lrand32(int a, int b);

    //: Generates a random unsigned long in [0,b]
    int lrand32(int b) {return lrand32(0, b);}

    //: Generates a random unsigned long in [a,b]
    int lrand32(int a, int b, int&);

    //:  Generates a random double in the range a <= x <= b with 32 bit randomness.
    //   drand32(1,0) is random down to about the 10th decimal place.
    double drand32(double a, double b);

    //: Generates a random unsigned integer in [0,n)
    // This function allows the random number generator to be used as
    // a functor, e.g. with std::random_shuffle()
    unsigned long operator()(unsigned n) { return lrand32(0, n-1); }

    //:  Generates a random double in the range 0 <= x <= b with 32 bit randomness.
    //   drand32(1.0) is random down to about the 10th decimal place.
    double drand32(double b) {return drand32(0.0, b);}

    //:  Generates a random double in the range 0 <= x <= 1 with 32 bit randomness.
    //   drand32() is random down to about the 10th decimal place.
    double drand32() {return drand32(0.0, 1.0);}

    //: Generates a random double in the range a <= x <= b with 64 bit randomness.
    //  Completely random down to the accuracy of an IEEE double.
    double drand64(double a, double b);

    //: Generates a random double in the range 0 <= x <= b with 64 bit randomness.
    //  Completely random down to the accuracy of an IEEE double.
    double drand64(double b) {return drand64(0.0, b);}

    //: Generates a random double in the range 0 <= x <= 1 with 64 bit randomness.
    //  Completely random down to the accuracy of an IEEE double.
    double drand64() {return drand64(0.0, 1.0);}

    //: Random value from a unit normal distribution about zero.
    // Uses a drand32() as its underlying generator.
    // Because the function uses a probability transform, the randomness (and
    // quantisation) is non-linearly dependent on the value. The further the
    // sample is from zero, the lower the number of bits on which it is random.
    double normal();

    //: Random value from a unit normal distribution about zero.
    // Uses a drand64() as its underlying generator.
    // Because the function uses a probability transform, the randomness (and
    // quantisation) is non-linearly dependent on the value. The further the
    // sample is from zero, the lower the number of bits on which it is random.
    double normal64();
};

// copy from .cpp
inline unsigned long
vnl_random::linear_congruential_lrand32()
{
    return linear_congruential_previous =
    (linear_congruential_previous * linear_congruential_multiplier + 1) & 0xffffffff;
}

//: Construct with seed
inline vnl_random::vnl_random(unsigned long seed)
: linear_congruential_previous(seed)
, mz_array_position(0UL)
, mz_borrow(0)
, mz_previous_normal_flag(0)
{
    reseed(seed);
}

//: Construct with seed
inline vnl_random::vnl_random(unsigned long seed[vnl_random_array_size])
: mz_array_position(0UL)
, mz_borrow(0)
, mz_previous_normal_flag(0)
{
    reseed(seed);
}

inline vnl_random::vnl_random(const vnl_random & r)
: linear_congruential_previous(r.linear_congruential_previous)
, mz_array_position(r.mz_array_position)
, mz_borrow(r.mz_borrow)
, mz_previous_normal_flag(r.mz_previous_normal_flag)
{
    for (unsigned int i = 0; i < vnl_random_array_size; ++i)
    {
        mz_seed_array[i] = r.mz_seed_array[i];
        mz_array[i] = r.mz_array[i];
    }
}

inline vnl_random &
vnl_random::operator=(const vnl_random & r)
{
    linear_congruential_previous = r.linear_congruential_previous;
    mz_array_position = r.mz_array_position;
    mz_borrow = r.mz_borrow;
    mz_previous_normal_flag = r.mz_previous_normal_flag;
    for (unsigned int i = 0; i < vnl_random_array_size; ++i)
    {
        mz_seed_array[i] = r.mz_seed_array[i];
        mz_array[i] = r.mz_array[i];
    }
    return *this;
}

inline  vnl_random::vnl_random()

{
    reseed();
}

inline vnl_random::~vnl_random()
{
    for (unsigned int i = 0; i < vnl_random_array_size; ++i)
    {
        mz_seed_array[i] = 0;
        mz_array[i] = 0;
    }
}

inline void
vnl_random::reseed()
{
    reseed((unsigned long)std::time(nullptr));
}

inline void
vnl_random::reseed(unsigned long seed)
{
    mz_array_position = 0UL;
    mz_borrow = 0;
    
    linear_congruential_previous = seed;
    // Use the lc generator to fill the array
    for (unsigned int i = 0; i < vnl_random_array_size; ++i)
    {
        mz_seed_array[i] = linear_congruential_lrand32();
        mz_array[i] = mz_seed_array[i];
    }
    
    // Warm up with 1000 randoms
    for (int j = 0; j < 1000; j++)
        lrand32();
}

inline void
vnl_random::reseed(const unsigned long seed[vnl_random_array_size])
{
    mz_array_position = 0UL;
    mz_borrow = 0L;
    
    for (unsigned int i = 0; i < vnl_random_array_size; ++i)
    {
        mz_array[i] = seed[i];
        mz_seed_array[i] = seed[i];
    }
}

inline void
vnl_random::restart()
{
    mz_array_position = 0UL;
    
    for (unsigned int i = 0; i < vnl_random_array_size; ++i)
    {
        mz_array[i] = mz_seed_array[i];
    }
}

inline double
vnl_random::normal()
{
    if (mz_previous_normal_flag)
    {
        mz_previous_normal_flag = 0;
        return mz_previous_normal;
    }
    else
    {
        double x, y, r2;
        do
        {
            x = drand32(-1.0, 1.0);
            y = drand32(-1.0, 1.0);
            r2 = x * x + y * y;
        } while (r2 >= 1.0 || r2 == 0.0);
        double fac = std::sqrt(-2.0 * std::log(r2) / r2);
        mz_previous_normal = x * fac;
        mz_previous_normal_flag = 1;
        return y * fac;
    }
}


//: Random value from a unit normal distribution about zero
// Uses a drand64() as its underlying generator.
// Because the function uses a probability transform, the randomness (and
// quantisation) is non-linearly dependent on the value. The further the sample
// is from zero, the lower the number of bits on which it is random.
inline double
vnl_random::normal64()
{
    if (mz_previous_normal_flag)
    {
        mz_previous_normal_flag = 0;
        return mz_previous_normal;
    }
    else
    {
        double x, y, r2;
        do
        {
            x = drand64(-1.0, 1.0);
            y = drand64(-1.0, 1.0);
            r2 = x * x + y * y;
        } while (r2 >= 1.0 || r2 == 0.0);
        double fac = std::sqrt(-2.0 * std::log(r2) / r2);
        mz_previous_normal = x * fac;
        mz_previous_normal_flag = 1;
        return y * fac;
    }
}

inline unsigned long
vnl_random::lrand32()
{
    unsigned long p1 = mz_array[(vnl_random_array_size + mz_array_position - mz_previous1) % vnl_random_array_size];
    unsigned long p2 = (p1 - mz_array[mz_array_position] - mz_borrow) & 0xffffffff;
    if (p2 < p1)
        mz_borrow = 0;
    if (p2 > p1)
        mz_borrow = 1;
    mz_array[mz_array_position] = p2;
    mz_array_position = (mz_array_position + 1) % vnl_random_array_size;
    return p2;
}

inline int
vnl_random::lrand32(int lower, int upper)
{
    assert(lower <= upper);
    
    // Note: we have to reject some numbers otherwise we get a very slight bias
    // towards the lower part of the range lower - upper. See below
    
    unsigned long range = upper - lower + 1;
    unsigned long denom = 0xffffffff / range;
    unsigned long ran;
    while ((ran = lrand32()) >= denom * range)
        ;
    return lower + int(ran / denom);
}


inline int
vnl_random::lrand32(int lower, int upper, int & count)
{
    assert(lower <= upper);
    
    // Note: we have to reject some numbers otherwise we get a very slight bias
    // towards the lower part of the range lower - upper. Hence this is a "count"
    // version of the above function that returns the number of lrand32()
    // calls made.
    
    unsigned long range = upper - lower + 1;
    unsigned long denom = 0xffffffff / range;
    unsigned long ran;
    count = 1;
    while ((ran = lrand32()) >= denom * range)
        ++count;
    return lower + int(ran / denom);
}

inline double
vnl_random::drand32(double lower, double upper)
{
    assert(lower <= upper);
    return (double(lrand32()) / 0xffffffff) * (upper - lower) + lower;
}

inline double
vnl_random::drand64(double lower, double upper)
{
    assert(lower <= upper);
    return (double(lrand32()) / 0xffffffff + double(lrand32()) / (double(0xffffffff) * double(0xffffffff))) *
    (upper - lower) +
    lower;
}

#endif // vnl_random_h
