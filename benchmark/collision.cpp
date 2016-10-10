#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <unordered_set>

#include "ziggurat_normal.hpp"
#include "etf_normal.hpp"


// Right p-value for Poisson distribution.
template<typename RealType, typename UIntType>
RealType right_pvalue(RealType mean, UIntType k) {
    // A brut force but reliable computation based on an evaluation of the CDF.
    // One could also use a regularized Gamma function:
    //   boost::math::gamma_p(k, mean)
    // Results are similar with both methods.
    RealType s = 0;
    RealType p = std::exp(-mean);
    for (UIntType i=1; i<k; ++i) {
        s += std::log(mean/i);
        p += std::exp(-mean + s);
    }

    return RealType(1) - (p<=1 ? p : RealType(1));
}


// Map number [0,1) to urns numbered from 0 to 2^dim - 1
template<typename RealType, typename UIntType>
class UrnMap
{
public:
    UrnMap(int dim) : n_(RealType(UIntType(1) << (dim-1))*2) {}

    // Sort a number in [0, 1) into an urn.
    UIntType operator[](RealType x) {
        RealType i = n_*x;
        assert(i<n_);
        return static_cast<UIntType>(i);
    }
private:
    RealType n_;
};


// Perform the Knuth collision test.
//
// The test simulates randomly throwing n balls into m urns where m=2^dim
// using a uniform distribution in [0,1).
// The number of balls is computed with the m/n ratio:
//  m/n = (2^dim)/n = 264
// and the test is repeated several times.
// Knuth (1981) suggested n=2^14 and m=2^20 and hence m/n=64, but when m>=2^30
// the right p-value estimates computed with the ratio m/n=64 for ideal
// inversion sampling look strangely biased towards 1. In practice though,
// using m/n=64 rather than m/n=256 does not appear to change the thresholds
// at which the different methods give right p-values below the 5% threshold.
template<typename RealType, typename UIntType>
struct Experiment
{
    template<typename U>
    void run(int min_dim, int max_dim, int repeat, U random_real_func) {
        std::cout << "[dimensions | trial | p-value]" << std::endl;
        for (int dim = min_dim; dim<=max_dim; ++dim) {
            UrnMap<RealType, UIntType> urn_map(dim);
            UIntType n = RealType(UIntType(1) << (dim - 1))/RealType(128);
            RealType m = RealType(UIntType(1) << (dim - 1))*2;
            RealType expectation = RealType(n)*RealType(n)/(2*m);
            for (int iter=1; iter<=repeat; ++iter) {
                std::unordered_set<UIntType> filled_urns;
                UIntType collisions = 0;
                for (UIntType k=0; k!=n; ++k) {
                    RealType r = random_real_func();
                    collisions +=
                        (filled_urns.insert(urn_map[r])).second == false;
                }

                RealType pvalue = right_pvalue(expectation, collisions);
                std::cout << dim << " " << iter << " " << pvalue << std::endl;
            }
            std::cout << "\n";
        }   
    }
};


// Cumulative distribution function of the normal distribution.
template<typename RealType>
struct Cdf {
    Cdf() : inv_sqrt2_(std::sqrt(0.5)) {}

    RealType operator()(RealType x) {
        return RealType(0.5) * (1 + std::erf(x*inv_sqrt2_));
    }

public:
    const RealType inv_sqrt2_;
};




int main() {
    using GeneratorType = std::mt19937;
    using RealType = double;
    using UIntType = std::uint_least32_t; // must hold at least max_dim bits

    constexpr UIntType W = 32;
    int min_dim = 26;
    int max_dim = 31;
    int repeat = 10;
        

    GeneratorType rng;

    Cdf<RealType> cdf;
    ZigguratNormalDistribution<RealType, W> ziggurat_dist;
    EtfNormalDistribution<RealType, W, 7> etf_dist;

    Experiment<RealType, UIntType> experiment;
    
    std::cout << "Statistics for inversion sampling (theoretical)." << std::endl;
    experiment.run(min_dim, max_dim, repeat,
        [&]() { return etf::generate_random_real<RealType, W>(rng); });

    std::cout << "Statistics for ziggurat." << std::endl;
    experiment.run(min_dim, max_dim, repeat,
        [&]() { return cdf(ziggurat_dist(rng)); });

    std::cout << "Statistics for ETF." << std::endl;
    experiment.run(min_dim, max_dim, repeat,
        [&]() { return cdf(etf_dist(rng)); });
}
