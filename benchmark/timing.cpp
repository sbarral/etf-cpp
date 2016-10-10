#include <iostream>
#include <random>
#include <vector>

#define NONIUS_RUNNER
#include <nonius_single.h++>

#include "original_ziggurat_normal.hpp"
#include "ziggurat_normal.hpp"
#include "etf_normal.hpp"

#include "etf_chi_squared.hpp"
#include "etf_chi_squared_low_dof.hpp"

#define NB_ITER 10000

#define DIST_SUM(RNG, DIST) (\
[](nonius::chronometer meter) { \
    std::size_t nb_iter = NB_ITER; \
    RNG rng; \
    auto dist = DIST; \
    meter.measure([&] { \
        double s = 0.0; \
        for (std::size_t i=0; i!=nb_iter; ++i) { \
            s += dist(rng); \
        } \
        return s; \
    }); \
})


//NONIUS_BENCHMARK("ETF chi-squared, k=1", DIST_SUM(std::mt19937, (EtfChiSquaredLowDofDistribution<double, 64, 8>(1.0, 1e-4, 10.0))));
//NONIUS_BENCHMARK("standard library chi-squared, k=1", DIST_SUM(std::mt19937, std::chi_squared_distribution<double>(1.0)));
//
//NONIUS_BENCHMARK("ETF chi-squared, k=2", DIST_SUM(std::mt19937, (EtfChiSquaredDistribution<double, 64, 8>(2.0, 10.0))));
//NONIUS_BENCHMARK("standard library chi-squared, k=2", DIST_SUM(std::mt19937, std::chi_squared_distribution<double>(2.0)));
//
//NONIUS_BENCHMARK("ETF chi-squared, k=5", DIST_SUM(std::mt19937, (EtfChiSquaredDistribution<double, 64, 8>(5.0, 16.0))));
//NONIUS_BENCHMARK("standard library chi-squared, k=5", DIST_SUM(std::mt19937, std::chi_squared_distribution<double>(5.0)));
//
//NONIUS_BENCHMARK("ETF chi-squared, k=12", DIST_SUM(std::mt19937, (EtfChiSquaredDistribution<double, 64, 8>(12.0, 28.0))));
//NONIUS_BENCHMARK("standard library chi-squared, k=12", DIST_SUM(std::mt19937, std::chi_squared_distribution<double>(12.0)));


NONIUS_BENCHMARK("original ziggurat normal (32-bit)", DIST_SUM(std::mt19937, OriginalZigguratNormalDistribution32<double>()));
NONIUS_BENCHMARK("ziggurat normal (32-bit)", DIST_SUM(std::mt19937, (ZigguratNormalDistribution<double, 32>())));
NONIUS_BENCHMARK("ETF normal (32-bit)", DIST_SUM(std::mt19937, (EtfNormalDistribution<double, 32, 7>())));

NONIUS_BENCHMARK("original ziggurat (64-bit)", DIST_SUM(std::mt19937_64, OriginalZigguratNormalDistribution64<double>()));
NONIUS_BENCHMARK("ziggurat normal (64-bit)", DIST_SUM(std::mt19937_64, (ZigguratNormalDistribution<double, 64>())));
NONIUS_BENCHMARK("ETF normal (64-bit)", DIST_SUM(std::mt19937_64, (EtfNormalDistribution<double, 64, 7>())));
NONIUS_BENCHMARK("standard library normal (64-bit)", DIST_SUM(std::mt19937_64, std::normal_distribution<double>()));

