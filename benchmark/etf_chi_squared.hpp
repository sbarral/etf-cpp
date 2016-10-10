
#ifndef ETF_CHI_SQUARED_HPP
#define ETF_CHI_SQUARED_HPP

#include <array>
#include <cstddef>
#include <cmath>
#include <limits>
#include <vector>

#include <etf/distribution.hpp>
#include <etf/util.hpp>



#ifndef CHI_SQUARED_PDF
#define CHI_SQUARED_PDF

template<typename RealType>
class ChiSquaredPdf {
public:
    ChiSquaredPdf() = default;

    ChiSquaredPdf(RealType k) : m_(RealType(0.5)*k - RealType(1.0)) {}

    RealType operator()(RealType x) {
        return x==RealType(0.0) ?
            0.0 : std::exp(std::log(x)*m_ - RealType(0.5)*x);
    }

private:
    RealType m_;
};

#endif // CHI_SQUARED_PDF


// ETF-based chi-squared distribution.
template<typename RealType, std::size_t W, std::size_t N>
class EtfChiSquaredDistribution
    : public etf::distribution<RealType, W, N,
                               ChiSquaredPdf<RealType>,
                               etf::weibull_tail_distribution<RealType, W>,
                               etf::weibull_pdf<RealType>>
{
private:
    using Parent =
        etf::distribution<RealType, W, N,
                          ChiSquaredPdf<RealType>,
                          etf::weibull_tail_distribution<RealType, W>,
                          etf::weibull_pdf<RealType>>;



public:
    EtfChiSquaredDistribution() = default;

    EtfChiSquaredDistribution(RealType k, RealType xtail);
};


template<typename RealType, std::size_t W, std::size_t N>
EtfChiSquaredDistribution<RealType, W, N>::EtfChiSquaredDistribution(
    RealType k, RealType xtail)
{
    const std::size_t n = std::size_t(1) << N;
    const RealType m = RealType(0.5)*k - RealType(1.0);
    
    // Tail.
    const RealType b = xtail/(RealType(0.5)*xtail - m);
    const RealType w = b*std::pow(xtail, m)*std::exp(-m);
    etf::weibull_tail_distribution<RealType, W> tail_dist(xtail, 1.0, b, 0.0);
    etf::weibull_pdf<RealType> tail_pdf(1.0, b, 0.0, w);
    const RealType tail_area = tail_pdf.tail_area(xtail);


    // Compute the quantiles.
    ChiSquaredPdf<RealType> pdf(k);
    const double rel_tol = std::numeric_limits<RealType>::epsilon()
                           * RealType(1e4);

    auto x_guess = etf::trapezoidal_rule_prepartition(pdf,
                                                      RealType(0.0), xtail, n);

    auto dpdf = [m](RealType x) {
        return x==RealType(0.0)?
            0.0 : (m - RealType(0.5)*x)
                  *std::exp(std::log(x)*(m - RealType(1.0)) - RealType(0.5)*x);
    };
    auto p = etf::newton_partition_monotonic(
        pdf, dpdf,
        x_guess.begin(), x_guess.end(),
        rel_tol);
    
    *static_cast<Parent*>(this) = etf::make_distribution<RealType, W, N>(
            p.x.begin(), p.x.end(), p.finf.begin(), p.fsup.begin(),
            pdf, tail_dist, tail_pdf, tail_area);
}

#endif // ETF_CHI_SQUARED_HPP

