#ifndef ETF_UTIL_HPP
#define ETF_UTIL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <iterator>
#include <utility>
#include <vector>

#include "random_digits.hpp"


/// Exclusive Top Floor namespace.
///
namespace etf {

namespace detail {

// Beware: for efficiency the diagonal and RHS terms are modified in place.
template<typename RealType>
void
solve_tridiagonal_system(const std::vector<RealType>& a,
                         std::vector<RealType>& b,
                         const std::vector<RealType>& c,
                         std::vector<RealType>& rhs,
                         std::vector<RealType>& sol) {
    using size_type = typename std::vector<RealType>::size_type;

    auto m = a.size();
    
    // Eliminate the sub-diagonal.
    for (size_type i=1; i!=m; ++i) {
        RealType pivot = a[i]/b[i-1];
        b[i] -= pivot*c[i-1];
        rhs[i] -= pivot*rhs[i-1];
    }
    
    // Solve the remaining upper bidiagonal system.
    sol[m-1] = rhs[m-1]/b[m-1];
    for (size_type i=m-2; ; --i) {
        sol[i] = (rhs[i] - c[i]*sol[i+1])/b[i];
        if (i==0) break;
    }
}        

} // namespace detail


/// Computes a partition dividing approximately evenly the area under a
/// function using the trapezoidal rule.
///
/// The trapezoidal rule is applied to function `f` over a regular grid with
/// `nb_points` grid points (including outer and inner nodes) to divide
/// interval [`x0`, `x1`] into the specified number of sub-intervals such
/// that the areas under the trapeze quandrature is the same in each
/// sub-interval.
/// The returned vector is the set of abscissae.
///
template<typename RealType, class Func>
#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
__attribute__ ((noinline))
#endif
std::vector<RealType>
trapezoidal_rule_prepartition(Func f,
                              RealType x0,
                              RealType x1,
                              std::size_t nb_intervals,
                              std::size_t nb_points) {
    using size_type = typename std::vector<RealType>::size_type;
    
    // Convenient const aliases.
    const auto& n = nb_points;
    const auto& m = nb_intervals;
    
    // Compute the curve.
    RealType dx = (x1-x0)/(n-1);
    std::vector<RealType> x(n);
    std::vector<RealType> y(n);
    for (size_type i=0; i!=(n-1); ++i) {
        x[i] = x0 + i*dx;
        y[i] = f(x[i]);
    }
    x[n-1] = x1;
    y[n-1] = f(x1);
    
    // Total area (scaled by 1/dx).
    RealType s = RealType(0.5)*(y[0] + y[n-1]);
    for (size_type i=1; i!=(n-2); ++i) {
        s += y[i];
    }
    
    // Choose abscissae that evenly split the area under the curve.
    std::vector<RealType> xp(m+1);
    xp[0] = x0;
    xp[m] = x1;   
    {
        RealType al = 0.0;
        RealType ar = RealType(0.5)*(y[0] + y[1]);
        size_type i=0;
        for (size_type j=1; j!=m; ++j) {
            RealType a = s*(static_cast<RealType>(j)/static_cast<RealType>(m));
            while (a>ar) {
                ++i;
                al = ar;
                ar += RealType(0.5)*(y[i] + y[i+1]);
            }
            xp[j] = x[i] + (x[i+1]-x[i])*((a-al)/(ar-al));
        }
    }
    
    return xp;
}


/// Computes a partition dividing approximately evenly the area under a function
/// using the trapezoidal rule.
///
/// This overload sets the number of grid points to the requested number of
/// partitions.
///
template<typename RealType, class Func>
#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
__attribute__ ((noinline))
#endif
std::vector<RealType>
trapezoidal_rule_prepartition(Func f,
                              RealType x0,
                              RealType x1,
                              std::size_t nb_intervals) {
    return trapezoidal_rule_prepartition<RealType, Func>(
        f, x0, x1, nb_intervals, nb_intervals );
}


/// A partition and the local function extrema over each sub-interval.
///
/// Partition of an interval into sub-intervals [`x[i]`, `x[i+1]`].
/// The infimum and supremum of the function over each sub-interval `i` are
/// stored as `finf[i]` and `fsup[i]`, respectively.
///
template<typename RealType>
struct partition_data
{
    std::vector<RealType> x;
    std::vector<RealType> finf;
    std::vector<RealType> fsup;
};


/// Computes an ETF partition using Newton's method.
///
/// A Newton's method (multivariate) is used to determine a partition of the
/// interval defined by the first and last point of the vector of abcissae
/// passed in argument in such a way that the rectangles making up an upper
/// Riemann sum of function `f` have equal areas.
///
/// The returned `partition_data` object includes as well the infimum and
/// supremum of the function over each sub-interval. 
///
/// For faster convergence it is recommended to provide a reasonable initial
/// estimate of the partition abcissae passed in arguments.
///
/// The derivative `df` of `f` and an ordered sequence of the inner function
/// extrema (boundary points excluded) must as well be provided.
///
/// The tolerance is the maximum relative dispersion of upper rectangle areas,
/// computed as the difference between the largest and smallest area relative
/// to the average area.
///
/// Empty tables (with size 0) are returned if the algorithm fails to
/// converge.
///
/// In order to improve convergence robustness (resp. speed), under-relaxation
/// (resp. over-relaxation) may be optionaly mandated by setting `relax` at
/// less (resp. more) than 1.
///
/// The maximum number of iterations for the Newtow method may be optionally
/// specified.
///
template<class Func, class DFunc, class InputIt1, class InputIt2>
#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
__attribute__ ((noinline))
#endif
auto
newton_partition(Func f,
                 DFunc df,
                 InputIt1 x_initial_first,
                 InputIt1 x_initial_last,
                 InputIt2 x_extremum_first,
                 InputIt2 x_extremum_last,
                 typename std::iterator_traits<InputIt1>::value_type tol,
                 typename std::iterator_traits<InputIt1>::value_type relax = 1,
                 unsigned int max_iter = 100)
-> partition_data<typename std::iterator_traits<InputIt1>::value_type> {

    using RealType = typename std::iterator_traits<InputIt1>::value_type;
    using size_type = typename std::vector<RealType>::size_type;

    // Partition object and convenient aliases.
    partition_data<RealType> p;
    auto& x    = p.x;
    auto& finf = p.finf;
    auto& fsup = p.fsup;
    
    // Initialization of the quadrature and extrema vectors.
    x.assign(x_initial_first, x_initial_last);
    const auto n = x.size() - 1;
    finf.resize(n);
    fsup.resize(n);
    
    std::vector<std::pair<RealType, RealType>> extrema;
    while (x_extremum_first!=x_extremum_last) {
        if ((*x_extremum_first-x.front())*(*x_extremum_first-x.back())<=0.0) {
            extrema.push_back(
                std::pair<RealType, RealType>(*x_extremum_first, f(*x_extremum_first)) );
        }
        ++x_extremum_first;
    }
    
    // define the main vectors and pre-compute edge values
    std::vector<RealType> y(n+1);
    std::vector<RealType> dx(n-1);
    std::vector<RealType> dy_dx(n+1);
    std::vector<RealType> dfsup_dxl(n), dfsup_dxr(n);
    std::vector<RealType> minus_s(n-1);
    std::vector<RealType> ds_dxc(n-1), ds_dxl(n-1), ds_dxr(n-1);
    
    y.front() = f(x.front());
    y.back()  = f(x.back());
    dy_dx.front() = 0.0;
    dy_dx.back()  = 0.0;
    
    unsigned int iter = 0;
    while(true)
    {
        // Compute the values at inner points.
        for (size_type i=1; i!=n; ++i) {
            y[i] = f(x[i]);
            dy_dx[i] = df(x[i]);
        }
        
        // Determine the supremum fsup of y in the range (x[i], x[i+1]),
        // the partial derivatives of fsup with respect to x[i] and x[i+1],
        // the minimum and maximum partition areas and the total area.
        auto extremum = extrema.begin();
        RealType max_area = 0.0;
        RealType min_area = std::numeric_limits<RealType>::max();
        RealType sum_area = 0.0;
        for (size_type i=0; i!=n; ++i) {
            if (y[i]>y[i+1]) {
                fsup[i] = y[i];
                dfsup_dxl[i] = dy_dx[i];
                dfsup_dxr[i] = 0.0;
            }
            else {
                fsup[i] = y[i+1];
                dfsup_dxl[i] = 0.0;
                dfsup_dxr[i] = dy_dx[i+1];
            }
            
            // Check if there are extrema within the (x[i], x[i+1]) range.
            while ( extremum!=extrema.end() &&
                    (extremum->first-x[i])*(extremum->first-x[i+1])<=0.0 )
            {
                if (extremum->second > fsup[i]) {
                    fsup[i] = extremum->second;
                    dfsup_dxl[i] = 0.0;
                    dfsup_dxr[i] = 0.0;
                }
                ++extremum;
            }
            
            RealType area = fsup[i]*std::abs(x[i+1]-x[i]);
            max_area = std::max(area, max_area);
            min_area = std::min(area, min_area);
            sum_area += area;
        }
        
        // Check convergence.
        if ((max_area-min_area)<tol*(sum_area/n)) {
            // Determine the infimum finf of y in the range (x[i], x[i+1]).
            extremum = extrema.begin();
            for (size_type i=0; i!=n; ++i) {
                if (y[i]>y[i+1]) {
                    finf[i] = y[i+1];
                }
                else {
                    finf[i] = y[i];
                }
                
                // Check if there are extrema within the (x[i], x[i+1]) range.
                while (extremum!=extrema.end() &&
                    (extremum->first-x[i])*(extremum->first-x[i+1])<=0.0 )
                {
                    if (extremum->second < finf[i]) {
                        finf[i] = extremum->second;
                    }
                    ++extremum;
                }                
            }
            break;
        }
        
        if (++iter>max_iter) {
            p.x.clear();
            p.finf.clear();
            p.fsup.clear();
            
            break;
        }
        
        // Area difference between neigboring rectangles and partial
        // derivatives of s with respect to x[i], x[i+1] and x[i+2].
        for (size_type i=0; i!=(n-1); ++i) {
            minus_s[i] = fsup[i]*(x[i+1]-x[i]) - fsup[i+1]*(x[i+2]-x[i+1]);
            
            ds_dxl[i] = fsup[i] - (x[i+1]-x[i])*dfsup_dxl[i];
            ds_dxc[i] = (x[i+2]-x[i+1])*dfsup_dxl[i+1]
                      - (x[i+1]-x[i])*dfsup_dxr[i]
                      - (fsup[i] + fsup[i+1]);
            ds_dxr[i] = fsup[i+1] + (x[i+2]-x[i+1])*dfsup_dxr[i+1];
        }
        
        // Solve the tri-diagonal system S + (dS/dX)*dX = 0 with:
        //         | ds0/dx1 ds0/dx2    0     ...                    0     |
        //         | ds1/dx1 ds1/dx2 ds1/dx3    0     ...            0     |
        // dS/dX = |    0    ds2/dx2 ds2/dx3 ds2/dx4    0     ...    0     |
        //         |                       ...                             |
        //         |    0     ...     0    ds(n-2)/dx(n-2) ds(n-2)/dx(n-2) |
        //
        //
        // and:
        //      | dx1     |         | minus_s0     |
        // dX = | ...     |    -S = | ...    |
        //      | dx(n-1) |         | minus_s(n-2) |
        detail::solve_tridiagonal_system<RealType>(
            ds_dxl, ds_dxc, ds_dxr, minus_s, dx );
        
        // For the sake of stability, updated positions are constrained within
        // the bounds set by former neighbors positions.
        {
            RealType x0 = x[0];
            for (size_type i=1; i!=n; ++i) {
                std::pair<RealType, RealType> x_range = std::minmax(x0, x[i+1]);
                x0 = x[i];
                x[i] = std::max(x[i] + relax*dx[i-1], x_range.first);
                x[i] = std::min(x[i], x_range.second);
            }
        }
    }
    
    // Voila.
    return p;
}


/// Computes an ETF partition using Newton's method.
///
/// This overload can be used if the function is monotonic over the specified
/// interval.
///
template<class Func, class DFunc, class InputIt1>
#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
__attribute__ ((noinline))
#endif
auto
newton_partition_monotonic(
    Func f,
    DFunc df,
    InputIt1 x_initial_first,
    InputIt1 x_initial_last,
    typename std::iterator_traits<InputIt1>::value_type tol,
    typename std::iterator_traits<InputIt1>::value_type relax = 1,
    unsigned int max_iter = 100)
-> partition_data<typename std::iterator_traits<InputIt1>::value_type> {
    
    typename std::iterator_traits<InputIt1>::value_type* dummy_ptr = 0;
    return newton_partition(f, df, x_initial_first, x_initial_last,
        dummy_ptr, dummy_ptr, tol, relax, max_iter);
}


/// Tail of a 3-parameter Weibull distribution generated by inversion sampling.
///
/// Generates the tail of a shifted Weibull distribution such that:
///
///  `f(x|a,b,c) = s*((x-c)/b)^(a-1)*exp[-((x-c)/b)^a]` if `x/b > x0/b`
///
/// and otherwise:
///
///  `f(x|a,b,c) = 0`
///
/// where `a` is strictly positive. The scale parameter `b` may be positive for
/// a tail extending to `+inf` and negative for a tail extending to `-inf`.
/// Parameter `c` is the so-called location parameter.
/// The (positive) normalization constant `s` need not be specified.
///
/// Template parameter `W` sets the requested precision (in bits) for the
/// generation of floating point random number.
///
template<typename RealType, std::size_t W>
class weibull_tail_distribution
{
public:
    using result_type = RealType;
    using param_type = weibull_tail_distribution<RealType, W>;
    
    
    weibull_tail_distribution(RealType x0=0.0,
                              RealType a=1.0, RealType b=1.0, RealType c=0.0)
    : inv_a_(RealType(1.0)/a), b_(b), c_(c), x0_(x0), 
      alpha_(std::pow((x0-c)/b, a))
    {}
    
    
    /// Returns a random variate using the random number generator passed as
    /// argument.
    ///
    /// This method is thread-safe as long as non-const methods are not used.
    template<class RngType>
    result_type operator()(RngType& g) const {
        RealType r = generate_random_real<RealType, W>(g);
        return c_ + b_*std::pow(alpha_ - std::log(RealType(1.0)-r), inv_a_);
    }
    
    
    /// Resets the distribution.
    ///
    /// This method is a no-op; it is defined for the sake of compatibility with
    /// distributions of the standard library.
    void reset() {}
    
    
    /// Returns an object containing the distribution parameters.
    ///
    param_type param() const {
        return *this;
    }
    
    
    /// Initializes the distribution with new distribution parameters.
    ///
    void param(const param_type& params)
    {
        *this = params;
    }
    
    /// Returns the smallest value potentially returned by `operator()`.
    result_type min() const {
        constexpr RealType MinusInf =
            std::numeric_limits<RealType>::is_iec559 ?
                 -std::numeric_limits<RealType>::infinity()
                : std::numeric_limits<RealType>::min();
        return b_<RealType(0.0) ? MinusInf : x0_;
    }
    
    
    /// Returns the greatest value potentially returned by `operator()`.
    result_type max() const {
        constexpr RealType PlusInf =
            std::numeric_limits<RealType>::is_iec559 ?
                  std::numeric_limits<RealType>::infinity()
                : std::numeric_limits<RealType>::min();
        return b_>RealType(0.0) ? PlusInf : x0_;
    }
    
    
    /// Returns distribution parameter `a`.
    result_type a() const {
        return RealType(1.0)/inv_a_;
    }
    
    
    /// Returns distribution parameter `b`.
    result_type b() const {
        return b_;
    }
    
    
    /// Returns distribution parameter `c`.
    result_type c() const {
        return c_;
    }
    
    
private:
    RealType inv_a_;
    RealType b_;
    RealType c_;
    RealType x0_;
    RealType alpha_;
};




/// 3-parameter Weibull probability density function with optional weighting.
///
/// The function is defined as:
///
///  `f(x|a,b,c) = w*a/|b|*((x-c)/b)^(a-1)*exp[-((x-c)/b)^a]` if `x/b > c/b`
///
/// and otherwise:
///   
///  `f(x|a,b,c) = 0`
///
/// where `a` is strictly positive, `b` may be negative and `w` is an optional
/// positive weighting factor. When `w` is 1 the function is the normalized
/// Weibull probability density function.
///
template<typename RealType>
class weibull_pdf
{
public:
    using result_type = RealType;
    
    weibull_pdf(RealType a=1.0, RealType b=1.0, RealType c=0.0, RealType w=1.0)
    : a_(a), inv_b_(RealType(1.0)/b), c_(c), s_(w*std::abs(a/b))
    {}
    
    
    /// Returns the value at `x`.
    ///
    result_type operator()(RealType x) const {
        RealType y = (x - c_)*inv_b_;
        if (y<RealType(0.0)) return RealType(0.0);
        RealType z = std::pow(y, a_-RealType(1.0));
        
        return s_*z*std::exp(-y*z);
    }
    
    /// Returns the total area under the function.
    ///
    /// This is equal to construction parameter `w`.
    ///
    result_type total_area() const {
        return s_/std::abs(a_*inv_b_);
    }
    
    /// Returns the area under the function from `x0` to `sign(b)*infinity`.
    ///
    result_type tail_area(RealType x0) const {
        RealType z0 = std::pow((x0 - c_)*inv_b_, a_);
        
        return s_*std::exp(-z0)/(a_*inv_b_);
    }
    
    
private:
    RealType a_;
    RealType inv_b_;
    RealType c_;
    RealType s_;
};



} // namespace etf


#endif // ETF_UTIL_HPP

