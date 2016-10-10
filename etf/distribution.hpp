#ifndef ETF_DISTRIBUTION_HPP
#define ETF_DISTRIBUTION_HPP

#include <cstddef>
#include <stdexcept>

#include "implementation.hpp"


/// Exclusive Top Floor namespace.
///
namespace etf {

/// Asymmetric ETF distribution with a rejection-sampled tail.
///
template<typename RealType, std::size_t W, std::size_t N,
         typename Func, typename OuterDist=void, typename OuterFunc=void>
class distribution : public
    detail::builder<RealType, W, N,
        detail::asymmetric<RealType, W, N,
            detail::rejection_composite<RealType, W, Func,
                                        OuterDist, OuterFunc>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::asymmetric<RealType, W, N,
                detail::rejection_composite<RealType, W,
                    Func, OuterDist, OuterFunc>>>;
    
public:
    distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    distribution(InputIt1 x_first, InputIt1 x_last,
                 InputIt2 finf_first,
                 InputIt3 fsup_first,
                 Func func,
                 OuterDist outer_dist, OuterFunc outer_func,
                 RealType outer_area)
    : Parent(func, outer_dist, outer_func) {
        this->build(0.0, x_first, x_last, finf_first, fsup_first, outer_area);
    }

private:
    using Parent::build;
};


/// Create a distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func, typename OuterDist, typename OuterFunc> // implicit
inline         
auto make_distribution(InputIt1 x_first, InputIt1 x_last,
                       InputIt2 finf_first,
                       InputIt3 fsup_first,
                       Func func, OuterDist outer_dist, OuterFunc outer_func,
                       RealType outer_area)
-> etf::distribution<RealType, W, N, Func, OuterDist, OuterFunc> {
    return etf::distribution<RealType, W, N, Func, OuterDist, OuterFunc>(
        x_first, x_last, finf_first, fsup_first, func,
        outer_dist, outer_func, outer_area);
}


/// Asymmetric ETF distribution with a user-provided tail distribution.
///
template<typename RealType, std::size_t W, std::size_t N,
         typename Func, typename OuterDist>
class distribution<RealType, W, N, Func, OuterDist, void> : public
    detail::builder<RealType, W, N,
        detail::asymmetric<RealType, W, N,
            detail::composite<RealType, W, Func, OuterDist>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::asymmetric<RealType, W, N,
                detail::composite<RealType, W, Func, OuterDist>>>;
    
public:
    distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    distribution(InputIt1 x_first, InputIt1 x_last,
                 InputIt2 finf_first,
                 InputIt3 fsup_first,
                 Func func, OuterDist outer_dist, RealType outer_area)
    : Parent(func, outer_dist) {
        this->build(0.0, x_first, x_last, finf_first, fsup_first, outer_area);
    }


private:
    using Parent::build;
};


/// Create a distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func, typename OuterDist> // implicit
inline         
auto make_distribution(InputIt1 x_first, InputIt1 x_last,
                       InputIt2 finf_first,
                       InputIt3 fsup_first,
                       Func func, OuterDist outer_dist,
                       RealType outer_area)
-> etf::distribution<RealType, W, N, Func, OuterDist> {
    return etf::distribution<RealType, W, N, Func, OuterDist>(
        x_first, x_last, finf_first, fsup_first, func, outer_dist, outer_area);
}


/// Asymmetric ETF distribution defined on a bounded interval.
///
template<typename RealType, std::size_t W, std::size_t N, typename Func>
class distribution<RealType, W, N, Func, void, void> : public
    detail::builder<RealType, W, N,
        detail::asymmetric<RealType, W, N,
            detail::bounded<RealType, W, Func>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::asymmetric<RealType, W, N,
                detail::bounded<RealType, W, Func>>>;
    
public:
    distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    distribution(InputIt1 x_first, InputIt1 x_last,
                 InputIt2 finf_first,
                 InputIt3 fsup_first,
                 Func func)
    : Parent(func) {
        this->build(0.0, x_first, x_last, finf_first, fsup_first, 0.0);
    }

private:
    using Parent::build;
};


/// Create a distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func> // implicit
inline         
auto make_distribution(InputIt1 x_first, InputIt1 x_last,
                       InputIt2 finf_first,
                       InputIt3 fsup_first,
                       Func func)
-> etf::distribution<RealType, W, N, Func> {
    return etf::distribution<RealType, W, N, Func>(
        x_first, x_last, finf_first, fsup_first, func);
}


/// Central ETF distribution with a rejection-sampled tail.
///
/// This is an efficient specialization of `symmetric_distribution` for
/// distributions that are symmetric about x=0.
///
template<typename RealType, std::size_t W, std::size_t N,
         typename Func, typename OuterDist=void, typename OuterFunc=void>
class central_distribution : public
    detail::builder<RealType, W, N,
        detail::central<RealType, W, N,
            detail::rejection_composite<RealType, W, Func,
                                        OuterDist, OuterFunc>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::central<RealType, W, N,
                detail::rejection_composite<RealType, W,
                    Func, OuterDist, OuterFunc>>>;
    
public:
    central_distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    central_distribution(InputIt1 x_first, InputIt1 x_last,
                         InputIt2 finf_first,
                         InputIt3 fsup_first,
                         Func func,
                         OuterDist outer_dist, OuterFunc outer_func,
                         RealType outer_area)
    : Parent(func, outer_dist, outer_func) {
        this->build(0.0, x_first, x_last, finf_first, fsup_first, outer_area);
    }

private:
    using Parent::build;
};


/// Create a central_distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func, typename OuterDist, typename OuterFunc> // implicit
inline         
auto make_central_distribution(InputIt1 x_first, InputIt1 x_last,
                               InputIt2 finf_first,
                               InputIt3 fsup_first,
                               Func func,
                               OuterDist outer_dist, OuterFunc outer_func,
                               RealType outer_area)
-> etf::central_distribution<RealType, W, N, Func, OuterDist, OuterFunc> {
    return etf::central_distribution<RealType, W, N, Func,
                                     OuterDist, OuterFunc>(
        x_first, x_last, finf_first, fsup_first, func,
        outer_dist, outer_func, outer_area);
}


/// Central ETF distribution with a user-provided tail distribution.
///
/// This is an efficient specialization of `symmetric_distribution` for
/// distributions that are symmetric about x=0.
///
template<typename RealType, std::size_t W, std::size_t N,
         typename Func, typename OuterDist>
class central_distribution<RealType, W, N, Func, OuterDist, void> : public
    detail::builder<RealType, W, N,
        detail::central<RealType, W, N,
            detail::composite<RealType, W, Func, OuterDist>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::central<RealType, W, N,
                detail::composite<RealType, W, Func, OuterDist>>>;
    

public:
    central_distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    central_distribution(InputIt1 x_first, InputIt1 x_last,
                         InputIt2 finf_first,
                         InputIt3 fsup_first,
                         Func func, OuterDist outer_dist,
                         RealType outer_area)
    : Parent(func, outer_dist) {
        this->build(0.0, x_first, x_last, finf_first, fsup_first, outer_area);
    }

private:
    using Parent::build;
};


/// Create a central_distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func, typename OuterDist> // implicit
inline         
auto make_central_distribution(InputIt1 x_first, InputIt1 x_last,
                               InputIt2 finf_first,
                               InputIt3 fsup_first,
                               Func func, OuterDist outer_dist, 
                               RealType outer_area)
-> etf::central_distribution<RealType, W, N, Func, OuterDist> {
    return etf::central_distribution<RealType, W, N, Func, OuterDist>(
        x_first, x_last, finf_first, fsup_first, func, outer_dist, outer_area);
}


/// Central ETF distribution defined on a bounded interval.
///
/// This is an efficient specialization of `symmetric_distribution` for
/// distributions that are symmetric about x=0.
///
template<typename RealType, std::size_t W, std::size_t N, typename Func>
class central_distribution<RealType, W, N, Func, void, void> : public
    detail::builder<RealType, W, N,
        detail::central<RealType, W, N,
            detail::bounded<RealType, W, Func>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::central<RealType, W, N,
                detail::bounded<RealType, W, Func>>>;

public:
    central_distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    central_distribution(InputIt1 x_first, InputIt1 x_last,
                         InputIt2 finf_first,
                         InputIt3 fsup_first,
                         Func func)
    : Parent(func) {
        this->build(0.0, x_first, x_last, finf_first, fsup_first, 0.0);
    }

private:
    using Parent::build;
};


/// Create a central_distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func> // implicit
inline         
auto make_central_distribution(InputIt1 x_first, InputIt1 x_last,
                               InputIt2 finf_first,
                               InputIt3 fsup_first,
                               Func func)
-> etf::central_distribution<RealType, W, N, Func> {
    return etf::central_distribution<RealType, W, N, Func>(
        x_first, x_last, finf_first, fsup_first, func);
}


/// Symmetric ETF distribution with a rejection-sampled tail.
///
template<typename RealType, std::size_t W, std::size_t N,
         typename Func, typename OuterDist=void, typename OuterFunc=void>
class symmetric_distribution : public
    detail::builder<RealType, W, N,
        detail::symmetric<RealType, W, N,
            detail::rejection_composite<RealType, W, Func,
                                        OuterDist, OuterFunc>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::symmetric<RealType, W, N,
                detail::rejection_composite<RealType, W,
                    Func, OuterDist, OuterFunc>>>;
    
public:
    symmetric_distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    symmetric_distribution(RealType x0,
                           InputIt1 x_first, InputIt1 x_last,
                           InputIt2 finf_first,
                           InputIt3 fsup_first,
                           Func func,
                           OuterDist outer_dist, OuterFunc outer_func,
                           RealType outer_area)
    : Parent(x0, func, outer_dist, outer_func) {
        this->build(x0, x_first, x_last, finf_first, fsup_first, outer_area);
    }

private:
    using Parent::build;
};


/// Create a symmetric_distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func, typename OuterDist, typename OuterFunc> // implicit
auto make_symmetric_distribution(RealType x0,
                                 InputIt1 x_first, InputIt1 x_last,
                                 InputIt2 finf_first,
                                 InputIt3 fsup_first,
                                 Func func,
                                 OuterDist outer_dist, OuterFunc outer_func,
                                 RealType outer_area)
-> etf::symmetric_distribution<RealType, W, N, Func, OuterDist, OuterFunc> {
    return etf::symmetric_distribution<RealType, W, N,
        Func, OuterDist, OuterFunc>(
            x0, x_first, x_last, finf_first, fsup_first,
            func, outer_dist, outer_func, outer_area);
}


/// Symmetric ETF distribution with a user-provided tail distribution.
///
template<typename RealType, std::size_t W, std::size_t N,
         typename Func, typename OuterDist>
class symmetric_distribution<RealType, W, N, Func, OuterDist, void> : public
    detail::builder<RealType, W, N,
        detail::symmetric<RealType, W, N,
            detail::composite<RealType, W, Func, OuterDist>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::symmetric<RealType, W, N,
                detail::composite<RealType, W, Func, OuterDist>>>;
    

public:
    symmetric_distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    symmetric_distribution(RealType x0,
                           InputIt1 x_first, InputIt1 x_last,
                           InputIt2 finf_first,
                           InputIt3 fsup_first,
                           Func func, OuterDist outer_dist,
                           RealType outer_area)
    : Parent(x0, func, outer_dist) {
        this->build(x0, x_first, x_last, finf_first, fsup_first, outer_area);
    }

private:
    using Parent::build;
};


/// Create a symmetric_distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func, typename OuterDist> // implicit
auto make_symmetric_distribution(RealType x0,
                                 InputIt1 x_first, InputIt1 x_last,
                                 InputIt2 finf_first,
                                 InputIt3 fsup_first,
                                 Func func, OuterDist outer_dist,
                                 RealType outer_area)
-> etf::symmetric_distribution<RealType, W, N, Func, OuterDist> {
    return etf::symmetric_distribution<RealType, W, N, Func, OuterDist>(
        x0, x_first, x_last, finf_first, fsup_first, func,
        outer_dist, outer_area);
}


/// Symmetric ETF distribution defined on a bounded interval.
///
template<typename RealType, std::size_t W, std::size_t N, typename Func>
class symmetric_distribution<RealType, W, N, Func, void, void> : public
    detail::builder<RealType, W, N,
        detail::symmetric<RealType, W, N,
            detail::bounded<RealType, W, Func>>>
{
private:
    using Parent =
        detail::builder<RealType, W, N,
            detail::symmetric<RealType, W, N,
                detail::bounded<RealType, W, Func>>>;

public:
    symmetric_distribution() = default;

    template<typename InputIt1, typename InputIt2, typename InputIt3>
    symmetric_distribution(RealType x0,
                           InputIt1 x_first, InputIt1 x_last,
                           InputIt2 finf_first,
                           InputIt3 fsup_first,
                           Func func)
    : Parent(x0, func) {
        this->build(x0, x_first, x_last, finf_first, fsup_first, 0.0);
    }

private:
    using Parent::build;
};


/// Create a symmetric_distribution object, deducing trailing types.
///
template<typename RealType, std::size_t W, std::size_t N, // explicit
         typename InputIt1, typename InputIt2, typename InputIt3, // implicit
         typename Func> // implicit
auto make_symmetric_distribution(RealType x0,
                                 InputIt1 x_first, InputIt1 x_last,
                                 InputIt2 finf_first,
                                 InputIt3 fsup_first,
                                 Func func)
-> etf::symmetric_distribution<RealType, W, N, Func> {
    return etf::symmetric_distribution<RealType, W, N, Func>(
        x0, x_first, x_last, finf_first, fsup_first, func);
}


} // namespace etf

#endif // ETF_DISTRIBUTION_HPP

