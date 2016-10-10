#ifndef TAIL_DIST_HPP
#define TAIL_DIST_HPP

#include <cmath>

#include <etf/random_digits.hpp>


// Marsaglia's algorithm for tail sampling.
template <typename RealType, int W>
class NormalTailDistribution {
public:
    NormalTailDistribution() = default;

    NormalTailDistribution(RealType xt): xt_(xt), inv_xt_(RealType(1)/xt) {}

    
    template<class G>
    RealType operator()(G& g) const
    {
        RealType dx;
        RealType y;
        do {
            dx = std::log(RealType(1.0)
                          - etf::generate_random_real<RealType, W>(g))
                 *inv_xt_;
            y  = std::log(RealType(1.0)
                          - etf::generate_random_real<RealType, W>(g));
        } while (-2*y < dx*dx);
        return xt_ - dx;
    }


private:
    RealType xt_;
    RealType inv_xt_;
};

#endif // TAIL_DIST_HPP

