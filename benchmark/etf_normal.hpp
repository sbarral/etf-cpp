#ifndef ETF_LIB_NORMAL_HPP
#define ETF_LIB_NORMAL_HPP

#include <array>
#include <cstddef>
#include <cmath>
#include <limits>
#include <vector>

#include <etf/distribution.hpp>
#include <etf/util.hpp>

#include "tail_dist.hpp"


// ETF-based central normal distribution.
template<typename RealType, std::size_t W, std::size_t N>
class EtfNormalDistribution
    : public etf::central_distribution<RealType, W, N, RealType (*)(RealType),
                                       NormalTailDistribution<RealType, W>>
{
private:
    using Parent =
        etf::central_distribution<RealType, W, N, RealType (*)(RealType),
                                  NormalTailDistribution<RealType, W>>;

public:
    EtfNormalDistribution();

private:
    static RealType pdf(RealType x) {
        return std::exp(RealType(-0.5)*x*x);
    }
};


template<typename RealType, std::size_t W, std::size_t N>
EtfNormalDistribution<RealType, W, N>::EtfNormalDistribution()
{
    const std::size_t n = std::size_t(1) << N;
    
    // The tail is sampled with Marsaglia's algorithm.
    // For high precision (W large), the position of the tail can be chosen
    // rather freely; for low W values, though, the tail position should be
    // chosen such that the area of the tail relatively to the whole area
    // sampled (upper rectangles + tail) is a multiple of 1/2^(W-N-1) so as
    // to avoid excessive rounding errors on the sampling probability.
    // Magical values that are closest to the empirical optimum of the tail
    // position (around 3.25) are tabulated for the common cases N=7 and N=8.
    RealType xtail;
    if (N==7 && W>=11) {
        RealType magic_xtail[] =
            { 1.532095304, 1.859950459, 2.150455371, 2.413614185,
              2.655703474, 2.880953316, 3.092363645, 3.292145211 };
        xtail = magic_xtail[std::min(W - 11, N)];
    }
    else if (N==8 && W>=12) {
        RealType magic_xtail[] =
            { 1.533103263, 1.861331463, 2.152146391, 2.415553089,
              2.657829951, 2.883210552, 3.094702254, 3.294526271 };
        xtail = magic_xtail[std::min(W - 12, N)];
    }
    else {
        // let's hope W is large...
        xtail = 3.25;
    }
    
    // Tail area.
    const RealType sqrt_pi_over_two = 1.2533141373155001;
    const RealType tail_area =
        sqrt_pi_over_two*std::erfc(xtail/std::sqrt(RealType(2.0)));

    // Compute the quantiles.
    const double rel_tol = std::numeric_limits<RealType>::epsilon()
                           * RealType(1e4);

    auto x_guess = etf::trapezoidal_rule_prepartition(pdf, RealType(0.0),
                                                       xtail, n);

    auto dpdf = [](RealType x) { return -x*std::exp(RealType(-0.5)*x*x); };
    auto p = etf::newton_partition_monotonic(
        pdf, dpdf,
        x_guess.begin(), x_guess.end(),
        rel_tol);
    
    *static_cast<Parent*>(this) =
        etf::make_central_distribution<RealType, W, N>(
            p.x.begin(), p.x.end(), p.finf.begin(), p.fsup.begin(),
            &pdf, NormalTailDistribution<RealType, W>(xtail), tail_area);
}

#endif // ETF_LIB_NORMAL_HPP

