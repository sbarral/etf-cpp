#ifndef ZIGGURAT_NORMAL_HPP
#define ZIGGURAT_NORMAL_HPP

#include <cmath>
#include <cstddef>
#include <vector>

#include <etf/random_digits.hpp>

#include "tail_dist.hpp"




// An arbitrary precision Ziggurat implementation.
//
// This is a generic, portable and hopefuly more correct re-implementation of
// the original algorithm (Marsaglia and Tsang, 2000). It is generic over the
// RNG, the floating point type and the number of random bits W to be requested
// from the RNG.
//
// Bit width W should not be greater than that of the longest unsigned integer
// available on the platform and must be able to hold more than the 7 bits
// required by the Ziggurat algorithm to generate the table index.
//
// The RNG passed as an argument to operator() should ideally have min=0 and
// max=2^N-1 (for some arbitrary N), but min=1 and/or max=2^N-2 are tolerated.
// If the raw random numbers generated by the RNG have less than W bits, 2 or
// more raw numbers are generated, as appropriate.
//
// The main differences with the original algorithm are:
//  - a bug in the original algorithm was fixed; the original algorithm applies
//    abs() to a signed integer which is occasionally equal to the smallest
//    possible negative value, the result of which is undefined since the absolute
//    value of the smallest negative integer cannot be represented by a positive
//    value; abs() will then typically return the same _negative_ integer,
//    causing a bug in the acceptance test,
//  - it does not require integer types to have exactly 32 or 64 bits (as
//    this is not guarranteed by the C++ standard),
//  - it does not rely on casting unsigned int values to out-of-range int values
//    as a way to obtain a cheap 2-complement negative numbers, because the C
//    standard says it is not portable.
//  - both the sign and the table index are determined from the upper RNG bits
//    because many RNGs have worse quality lower bits (e.g. xorshift family
//    RNGs).
//
// The overhead in term of execution speed compared to the original seems to be
// very modest.
template<typename RealType, std::size_t W>
class ZigguratNormalDistribution
{
private:
    using IntType = typename etf::integer_traits<W>::int_fast_t;
    using UIntType = typename etf::integer_traits<W>::uint_fast_t;

public:
    ZigguratNormalDistribution();    

    template<class G>
    RealType operator()(G& g)
    {
        using std::abs;
        using std::exp;

        while (true) {
            auto r = etf::generate_random_integer<UIntType, W>(g);
            // Integer u is made of bits 0:(W-8); unlike the original ziggurat,
            // we avoid relying on casting to obtain a potentially negative
            // signed integer from an unsigned integer since the C standard does
            // not guarantee the portability of overflowing casts. Instead, we
            // cast an unsigned integer that is small enough to be representable
            // as a signed integer and then apply a complement to obtain a
            // potentially negative integer.
            // Since the signed integer can hold at least W-1 digits while the
            // actual value has at most W-8 digits, using abs() on a negative
            // value is guarranteed to give a representable positive value
            // (this fixes the bug in the original ziggurat whereby abs() could
            // be potentially applied to the smallest signed integer, with
            // undefined consequences).
            constexpr UIntType m_mask = (UIntType(1) << (W - 7)) - 1;
            constexpr IntType complement =
                static_cast<IntType>((UIntType(1) << (W - 8)) - 1);
            IntType u = complement - IntType(r & m_mask);
            // The table index is made of bits (W-7):(W-1).
            std::size_t i = r >> (W - 7);
            
            if (abs(u)<k_[i])
                return u*w_[i];

            if (i==0)
                return u>0 ? tail_dist_(g) : -tail_dist_(g);

            RealType x = u*w_[i];
            RealType v = etf::generate_random_real<RealType, W>(g);
            if (f_[i] + v*(f_[i-1] - f_[i]) <= exp(RealType(-0.5)*x*x))
                return x;
        }
    }

private:
    std::vector<RealType> w_;
    std::vector<IntType> k_;
    std::vector<RealType> f_;
    NormalTailDistribution<RealType, W> tail_dist_;
};



template<typename RealType, std::size_t W>
ZigguratNormalDistribution<RealType, W>::ZigguratNormalDistribution()
:   tail_dist_(3.442619855899)
{
    const std::size_t n = 128;
    const RealType xt = 3.442619855899;
    const RealType weight = 9.91256303526217e-3;
    const IntType scale = static_cast<IntType>(UIntType(1) << (W-8));

    w_.resize(n);
    k_.resize(n);
    f_.resize(n);

    // Build the tables.
    f_[0] = 1;
    f_[n-1] = std::exp(RealType(-0.5)*xt*xt);
    w_[0] = weight/(f_[n-1]*scale);
    w_[n-1] = xt/scale;
    k_[0] = static_cast<IntType>(xt/w_[0]);
    k_[1] = 0;
    
    RealType xmax_prev = xt;
    for(std::size_t i=n-2; i!=0; --i)
    {
        RealType xmax = std::sqrt(-RealType(2)*std::log(weight/xmax_prev + f_[i+1]));
        k_[i+1] = static_cast<IntType>((xmax/xmax_prev)*scale);
        w_[i] = xmax/scale;
        f_[i] = std::exp(RealType(-0.5)*xmax*xmax);
        xmax_prev = xmax;
    }
}



#endif // ZIGGURAT_NORMAL_HPP
