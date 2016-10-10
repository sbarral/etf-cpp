#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <etf/random_digits.hpp>

#include "tail_dist.hpp"


namespace detail {

// A C++ implementation of the original ziggurat algorithmn.
//
// This is a straightforward implementation of the original algorithm
// (Marsaglia and Tsang, 2000).
// It is neither generic nor portable:
//  - it requires integer types with exactly 32 or 64 bits, which existence is
//    not guarranteed by the C++ standard,
//  - it relies on a cast from unsigned int to int with values that are not
//    within the range of representable int values in order to create negative
//    integers; the result of such cast is platform-dependent according to the
//    C standard.
//
// BEWARE: this original algorithm has a bug because it applies abs() to a
// signed integer which occasionally has the smallest possible negative value,
// the result of which is undefined since the absolute value of the smallest
// negative integer cannot be represented by a positive value and in fact abs()
// will then typically return the same _negative_ integer, causing a bug in the
// acceptance test.
//
// Note that the signed and unsigned integers given as template parameters must
// have exactly the same size. In order to prevent mistakes, two public
// typedefs are provided for 32 and 64 bits.
//
template<typename RealType, typename IntType, typename UIntType>
class OriginalZigguratNormalDistribution
{
public:
    OriginalZigguratNormalDistribution();   

    template<class G>
    RealType operator()(G& g)
    {
        using std::abs;
        using std::exp;

        constexpr int M = std::numeric_limits<IntType>::digits;

        while (true) {
            // The potentially overflowing cast from UIntType to IntType is
            // from the original ziggurat implementation but is definitely
            // not portable; the C standard says that casting out-of-range
            // values is implementation-defined, so we cross our fingers...
            auto u = static_cast<IntType>(
                    etf::generate_random_integer<UIntType, M+1>(g));
            std::size_t i = u & 0b1111111; // take the lower 7 bits
            
            if (abs(u)<k_[i])
                return u*w_[i];

            if (i==0)
                return u>0? tail_dist_(g) : -tail_dist_(g);

            RealType x = u*w_[i];
            RealType v = etf::generate_random_real<RealType, M+1>(g);
            if (f_[i] + v*(f_[i-1] - f_[i]) < exp(RealType(-0.5)*x*x))
                return x;
        }
    }

private:
    std::vector<RealType> w_;
    std::vector<IntType> k_;
    std::vector<RealType> f_;
    NormalTailDistribution<
        RealType, std::numeric_limits<UIntType>::digits> tail_dist_;
};

template<typename RealType, typename IntType, typename UIntType>
OriginalZigguratNormalDistribution<
    RealType, IntType, UIntType>::OriginalZigguratNormalDistribution()
:   tail_dist_(3.442619855899)
{
    constexpr int M = std::numeric_limits<IntType>::digits;
    const std::size_t n = 128;
    const RealType xt = 3.442619855899;
    const RealType weight = 9.91256303526217e-3;
    const RealType scale = RealType(UIntType(1) << M);

    w_.resize(n);
    k_.resize(n);
    f_.resize(n);

    // Build the tables.
    f_[0] = 1;
    f_[n-1] = std::exp(-RealType(0.5)*xt*xt);
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

} // namespace detail

template<typename RealType> using OriginalZigguratNormalDistribution32 =
    detail::OriginalZigguratNormalDistribution<RealType,
                                               std::int32_t, std::uint32_t>;

template<typename RealType> using OriginalZigguratNormalDistribution64 =
    detail::OriginalZigguratNormalDistribution<RealType,
                                               std::int64_t, std::uint64_t>;
