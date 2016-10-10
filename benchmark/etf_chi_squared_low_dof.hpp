
#ifndef ETF_CHI_SQUARED_LOW_DOF_HPP
#define ETF_CHI_SQUARED_LOW_DOF_HPP

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


// Outer composite distribution for chi-squared distributions with less than
// 2 degrees of freedom.
template<typename RealType, std::size_t W>
class ChiSquaredOuterDistribution
{
private:
    class LeftDist
    {
    public:
        LeftDist() = default;

        LeftDist(RealType k, RealType x0) : x0_(x0), p_(RealType(2.0)/k) {}

        // Sample a distribution with PDF ~ x^(k/2-1) in [0,x0).
        template<class G>
        RealType operator()(G& g) {
            return x0_*std::pow(etf::generate_random_real<RealType, W>(g), p_);
        }

    private:
        RealType x0_;
        RealType p_;
    };


public:
    ChiSquaredOuterDistribution() = default;

    ChiSquaredOuterDistribution(RealType k, RealType x0, RealType xtail)
        : left_dist_(k, x0), right_dist_(xtail, 1.0, 2.0, 0.0)
    {
        // Compute the sampling probability of each distribution.
        const RealType right_dist_area = RealType(2.0)
            *std::pow(xtail, RealType(0.5)*k - RealType(1.0))*std::exp(RealType(-0.5)*xtail);
        const RealType left_dist_area =
            RealType(2.0)/k * std::pow(x0, RealType(0.5)*k);
        area_ = left_dist_area + right_dist_area;
        dist_switch_ = left_dist_area/area_;
    }

    template<class G>
    RealType operator()(G& g) {
        if (etf::generate_random_real<RealType, W>(g) < dist_switch_)
            // Left distribution sampling.
            return left_dist_(g);
        else
            // Right (tail) distribution sampling.
            return right_dist_(g);
    }
   
    RealType total_non_normalized_area() {
        return area_;
    }

private:
    RealType area_;
    RealType dist_switch_;
    LeftDist left_dist_;
    etf::weibull_tail_distribution<RealType, W> right_dist_;
};


// PDF of the outer composite distribution for chi-squared distribution with
// less than 2 degrees of freedom.
template<typename RealType>
class ChiSquaredOuterPdf
{
public:
    ChiSquaredOuterPdf() = default;

    ChiSquaredOuterPdf(RealType k, RealType x0, RealType xtail)
        : m_(RealType(0.5)*k - RealType(1.0)),
          x_switch_(RealType(0.5)*(x0 + xtail)),
          right_pdf_(1.0, 2.0, 0.0, RealType(2.0)*std::pow(xtail, m_)) {}

    RealType operator()(RealType x) {
        return x < x_switch_ ? std::pow(x, m_) : right_pdf_(x);
    }


private:
    RealType m_;
    RealType x_switch_;
    etf::weibull_pdf<RealType> right_pdf_;
};


// ETF-based chi-squared distribution with less than 2 degrees of freedom.
template<typename RealType, std::size_t W, std::size_t N>
class EtfChiSquaredLowDofDistribution
    : public etf::distribution<RealType, W, N,
                               ChiSquaredPdf<RealType>,
                               ChiSquaredOuterDistribution<RealType, W>,
                               ChiSquaredOuterPdf<RealType>>
{
private:
    using Parent =
        etf::distribution<RealType, W, N,
                          ChiSquaredPdf<RealType>,
                          ChiSquaredOuterDistribution<RealType, W>,
                          ChiSquaredOuterPdf<RealType>>;


public:
    EtfChiSquaredLowDofDistribution() = default;

    EtfChiSquaredLowDofDistribution(RealType k, RealType x0, RealType xtail);
};


template<typename RealType, std::size_t W, std::size_t N>
EtfChiSquaredLowDofDistribution<RealType, W, N>::EtfChiSquaredLowDofDistribution(
    RealType k, RealType x0, RealType xtail)
{
    const std::size_t n = std::size_t(1) << N;
    const RealType m = RealType(0.5)*k - RealType(1.0);
    
    // Compute the quantiles.
    ChiSquaredPdf<RealType> pdf(k);
    const double rel_tol = std::numeric_limits<RealType>::epsilon()
                           * RealType(1e4);

    auto x_guess = etf::trapezoidal_rule_prepartition(pdf, x0, xtail, n);

    auto dpdf = [m](RealType x) {
        return (m - RealType(0.5)*x)*std::exp(std::log(x)*(m - RealType(1.0))
                                              - RealType(0.5)*x);
    };
    auto p = etf::newton_partition_monotonic(
        pdf, dpdf,
        x_guess.begin(), x_guess.end(),
        rel_tol);
    
    ChiSquaredOuterDistribution<RealType, W> outer_dist(k, x0, xtail);
    ChiSquaredOuterPdf<RealType> outer_pdf(k, x0, xtail);
    const RealType outer_area = outer_dist.total_non_normalized_area();

    //RealType upper_quadrature_area = 0.0;
    //for (std::size_t i=0; i!=n; ++i) {
    //    upper_quadrature_area +=
    //        (p.x[i+1] - p.x[i])*p.fsup[i];
    //}
    //std::cout << "Relative outer area: " << outer_area/upper_quadrature_area << std::endl;

    *static_cast<Parent*>(this) = etf::make_distribution<RealType, W, N>(
            p.x.begin(), p.x.end(), p.finf.begin(), p.fsup.begin(),
            pdf, outer_dist, outer_pdf, outer_area);
}

#endif // ETF_CHI_SQUARED_LOW_DOF_HPP

