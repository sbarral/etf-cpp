#ifndef ETF_IMPLEMENTATION_HPP
#define ETF_IMPLEMENTATION_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "exceptions.hpp"
#include "random_digits.hpp"


namespace etf {

namespace detail {

template<typename RealType, std::size_t W>
class data {
public:
    using result_type = RealType;

protected:
    using UIntType = typename etf::integer_traits<W>::uint_least_t;

private:
    struct Datum
    {
        UIntType scaled_fratio;
        RealType scaled_fsup;
        RealType scaled_dx;
    };

protected:
    std::vector<RealType> x_;
    std::vector<Datum> data_;
};


template<typename RealType, std::size_t W, typename Func>
class bounded : public data<RealType, W>
{
private:
    using Parent = data<RealType, W>;

protected:
    using UIntType = typename Parent::UIntType;

public:
    template<typename=void>
    void reset() {}

protected:
    bounded() = default;

    bounded(Func func) : func_(func) {}

    void define_outer_switch(UIntType) {}; // never called

    UIntType outer_switch() { return 0; } // never called

    template<class RngType>
    bool sample_outer(RngType& g, RealType& x) { return false; } // never called

    template<typename=void>
    RealType outer_min() const { return 0.0; } // never called
    
    template<typename=void>
    RealType outer_max() const { return 0.0; } // never called
    
protected:
    Func func_;
    static constexpr bool HasOuter = false;
    static constexpr bool HasRejection = false;
};



template<typename RealType, std::size_t W, typename Func, typename OuterDist>
class composite : public bounded<RealType, W, Func>
{
private:
    using Parent = bounded<RealType, W, Func>;

protected:
    using UIntType = typename Parent::UIntType;

public:
    template<typename=void>
    void reset() {
        outer_dist_.reset();
    }

protected:
    composite() = default;

    composite(Func func, OuterDist outer_dist)
    : Parent(func), outer_dist_(outer_dist) {}

    void define_outer_switch(UIntType outer_switch) {
        outer_switch_ = outer_switch;
    }

    UIntType outer_switch() {
        return outer_switch_;
    }

    template<class RngType>
    bool sample_outer(RngType& g, RealType& x)
    {
        x = outer_dist_(g);
        return true;
    }

    template<typename=void>
    RealType outer_min() const {
        return outer_dist_.min();
    }
    
    template<typename=void>
    RealType outer_max() const {
        return outer_dist_.max();
    }
    
protected:
    UIntType outer_switch_;
    OuterDist outer_dist_;
    static constexpr bool HasOuter = true;
    static constexpr bool HasRejection = false;
};


template< typename RealType, std::size_t W,
          typename Func, typename OuterDist, typename OuterFunc>
class rejection_composite : public composite<RealType, W, Func, OuterDist>
{
private:
    using Parent = composite<RealType, W, Func, OuterDist>;

protected:
    using UIntType = typename Parent::UIntType;

protected:
    rejection_composite() = default;

    rejection_composite(Func func, OuterDist outer_dist, OuterFunc outer_func)
    : Parent(func, outer_dist), outer_func_(outer_func)
    {}

    template<class RngType>
    bool sample_outer(RngType& g, RealType& x)
    {
        RealType r = generate_random_real<RealType, W>(g);
        x = this->outer_dist_(g);
        return r*outer_func_(x) <= this->func_(x);
    }

protected:
    OuterFunc outer_func_;
    static constexpr bool HasRejection = true;
};


template<typename RealType, std::size_t W, std::size_t N, class Category>
class asymmetric : public Category
{
protected:
    using UIntType = typename Category::UIntType;

public:
    /// Returns a random number.
    ///
    template<class RngType>
    RealType operator()(RngType& g) {
        while (true)
        {
            // Generate a table index and a positive value from a single random
            // number.
            auto r = generate_random_integer<UIntType, W>(g);
            // Mantissa is made of bits 0:(W-N-1).
            constexpr UIntType m_mask = (UIntType(1) << (W - N)) - 1;
            UIntType u = r & m_mask;
            // The table index is made of bits (W-N):(W-1).
            auto i = std::size_t(r >> (W - N));
            
            const auto& d = this->data_[i];
            // Note that the following test will also fail if 'u' is greater or
            // equal to the outer switch value since all 'fratio' values are
            // lower than the switch value.
            if (u<d.scaled_fratio)
                return this->x_[i] + d.scaled_dx*u;
            
            // Should the outer distribution be sampled?
            if (Category::HasOuter && u>=this->outer_switch()) {
                RealType x;
                bool acceptance = this->sample_outer(g, x);
                if (!Category::HasRejection || acceptance)
                    return x;
            }

            // Otherwise it is a wedge, test y<f(x) for rejection sampling.
            RealType v = generate_random_real<RealType, W>(g); // v in [0,1)
            RealType x = this->x_[i] + v*(this->x_[i+1] - this->x_[i]);
            if ((u*d.scaled_fsup) < this->func_(x))
                return x;
        }
    }

    template<typename=void>
    RealType min() const {
        RealType m = std::min(this->x_.front(), this->x_.back());
        return Category::HasOuter? std::min(m, this->outer_min()) : m;
    }
    
    template<typename=void>
    RealType max() const {
        RealType m = std::max(this->x_.front(), this->x_.back());
        return Category::HasOuter? std::max(m, this->outer_max()) : m;
    }
    
protected:
    using Category::Category;

protected:
    static constexpr bool IsSymmetric = false;
};


template<typename RealType, std::size_t W, std::size_t N, class Category>
class central : public Category
{
protected:
    using UIntType = typename Category::UIntType;

public:
    template<class RngType>
    RealType operator()(RngType& g) {
        while (true)
        {
            // Generate a table index, a sign bit and a positive value from a
            // single random number.
            auto r = generate_random_integer<UIntType, W>(g);
            // Mantissa is made of bits 0:(W-N-2).
            constexpr UIntType m_mask = (UIntType(1) << (W - N - 1)) - 1;
            UIntType u = r & m_mask;
            // The table index is made of bits (W-N-1):(W-2).
            constexpr std::size_t i_mask = (std::size_t(1) << N) - 1;
            auto i = std::size_t(r >> (W - N - 1)) & i_mask;
            // Sign is bit (W-1).
            int s = r >> (W - 1) ? 1 : -1;
            
            const auto& d = this->data_[i];
            // Note that the following test will also fail if 'u' is greater or
            // equal to the outer switch value since all 'fratio' values are
            // lower than the switch value.
            if (u<d.scaled_fratio)
                return s*(this->x_[i] + d.scaled_dx*u);
            
            // Should the outer distribution be sampled?
            if (Category::HasOuter && u>=this->outer_switch()) {
                RealType x;
                bool acceptance = this->sample_outer(g, x);
                if (!Category::HasRejection || acceptance)
                    return s*x;
            }

            // Otherwise it is a wedge, test y<f(x) for rejection sampling.
            RealType v = generate_random_real<RealType, W>(g); // v in [0,1)
            RealType x = this->x_[i] + v*(this->x_[i+1] - this->x_[i]);
            if ((u*d.scaled_fsup) < this->func_(x))
                return s*x;
        }
    }

    template<typename=void>
    RealType min() const {
        auto mm = std::minmax(this->x_.front(), this->x_.back());
        RealType m = std::min(mm.first, -mm.second);
        if (Category::HasOuter) {
            RealType t = std::min(this->outer_min(), -this->outer_max());
            return std::min(m, t);
        }
        else
            return m;
    }
    
    template<typename=void>
    RealType max() const {
        auto mm = std::minmax(this->x_.front(), this->x_.back());
        RealType m = std::max(-mm.first, mm.second);
        if (Category::HasOuter) {
            RealType t = std::max(this->outer_max(), -this->outer_min());
            return std::max(m, t);
        }
        else
            return m;
    }
    
protected:
    using Category::Category;

protected:
    static constexpr bool IsSymmetric = true;
};


template<typename RealType, std::size_t W, std::size_t N, class Category>
class symmetric : public Category
{
protected:
    using UIntType = typename Category::UIntType;

public:
    template<class RngType>
    RealType operator()(RngType& g) {
        while (true)
        {
            // Generate a table index, a sign bit and a positive value from a
            // single random number.
            auto r = generate_random_integer<UIntType, W>(g);
            // Mantissa is made of bits 0:(W-N-2).
            constexpr UIntType m_mask = (UIntType(1) << (W - N - 1)) - 1;
            UIntType u = r & m_mask;
            // The table index is made of bits (W-N-1):(W-2).
            constexpr std::size_t i_mask = (std::size_t(1) << N) - 1;
            auto i = std::size_t(r >> (W - N - 1)) & i_mask;
            // Sign is bit (W-1).
            int s = r >> (W - 1) ? 1 : -1;
            
            const auto& d = this->data_[i];
            // Note that the following test will also fail if 'u' is greater or
            // equal to the outer switch value since all 'fratio' values are
            // lower than the switch value.
            if (u<d.scaled_fratio)
                return x_origin_ + s*(this->x_[i] + d.scaled_dx*u);
            
            // Should the outer distribution be sampled?
            if (Category::HasOuter && u>=this->outer_switch()) {
                RealType x;
                bool acceptance = this->sample_outer(g, x);
                if (!Category::HasRejection || acceptance)
                    return x_origin_ + s*(x - x_origin_);
            }

            // Otherwise it is a wedge, test y<f(x) for rejection sampling.
            RealType v = generate_random_real<RealType, W>(g); // v in [0,1)
            RealType x = this->x_[i] + v*(this->x_[i+1] - this->x_[i]);
            if ((u*d.scaled_fsup) < this->func_(x + x_origin_))
                return x_origin_ + s*x;
        }
    }

    template<typename=void>
    RealType min() const {
        auto mm = std::minmax(this->x_.front(), this->x_.back());
        RealType m = std::min(x_origin_ + mm.first, x_origin_ - mm.second);
        if (Category::HasOuter) {
            RealType t = std::min(this->outer_min(),
                                  2*x_origin_ - this->outer_max());
            return std::min(m, t);
        }
        else
            return m;
    }
    
    template<typename=void>
    RealType max() const {
        auto mm = std::minmax(this->x_.front(), this->x_.back());
        RealType m = std::max(x_origin_ - mm.first, x_origin_ + mm.second);
        if (Category::HasOuter) {
            RealType t = std::max(this->outer_max(),
                                  2*x_origin_ - this->outer_min());
            return std::max(m, t);
        }
        else
            return m;
    }

protected:
    symmetric() = default;

    template<typename... Args>
    symmetric(RealType x0, Args... args): Category(args...), x_origin_(x0) {}

protected:
    RealType x_origin_;
    static constexpr bool IsSymmetric = true;
};


template<typename RealType, std::size_t W, std::size_t N, class Shape>
struct builder : public Shape {
private:
    using typename Shape::UIntType;

protected:
    using Shape::Shape;

    template<class InputIt1, class InputIt2, class InputIt3>
    void build(RealType x_origin,
               InputIt1 x_first, InputIt1 x_last,
               InputIt2 finf_first, InputIt3 fsup_first,
               RealType outer_area);

private:
    using Shape::define_outer_switch;
    using Shape::outer_switch;
    using Shape::sample_outer;
    using Shape::outer_min;
    using Shape::outer_max;

private:
    using Shape::x_;
    using Shape::data_;

};


template<typename RealType, std::size_t W, std::size_t N, class Shape>
template<class InputIt1, class InputIt2, class InputIt3>
void builder<RealType, W, N, Shape>::build(RealType x_origin,
                                           InputIt1 x_first, InputIt1 x_last,
                                           InputIt2 finf_first,
                                           InputIt3 fsup_first,
                                           RealType outer_area) {
    using UIntType = typename Shape::UIntType;

    // Bit-width of the sign, if any.
    constexpr std::size_t S = builder::IsSymmetric ? 1 : 0;

    // Resize tables and assign x coordinates.
    const std::size_t n = std::size_t(1) << N;
    this->data_.resize(n);
    this->x_.assign(x_first, x_last);
    if (this->x_.size()!=(n + 1))
        throw invalid_table_size();
    
    // Translate x to make it relative to the origin if necessary.
    if (x_origin!=RealType(0.0)) {
        for (auto& x: this->x_)
            x -= x_origin;
    }

    // Assign fsup but do not perform any scaling for now.
    for (std::size_t i=0; i!=n; ++i)
       this->data_[i].scaled_fsup = *fsup_first++; 

    // Compute the outer switch, i.e. an integer threshold such that when
    // drawing a random integer r, the probability:
    //  P(r>=switch)
    // expresses the probability to sample the outer distribution.
    UIntType outer_switch;
    if (builder::HasOuter) {
        RealType upper_quadrature_area = 0.0;
        for (std::size_t i=0; i!=n; ++i) {
            upper_quadrature_area +=
                (this->x_[i+1] - this->x_[i])*this->data_[i].scaled_fsup;
        }
        outer_switch = static_cast<UIntType>(
            std::round(RealType(UIntType(1) << (W - N - S)) *
            (upper_quadrature_area/(outer_area + upper_quadrature_area))));
        this->define_outer_switch(outer_switch);
    } else {
        outer_switch = UIntType(1) << (W - N - S);
    }

    // Compute the tables.
    for (std::size_t i=0; i!=n; ++i) {
        RealType fratio = (*finf_first++)/this->data_[i].scaled_fsup;
        if (fratio>=RealType(0.5)) // will we loose at most 1 bit of accuracy?
            this->data_[i].scaled_fratio =
                static_cast<UIntType>(fratio*outer_switch);
        else // otherwise, force wedge sampling to ensure high quality samples
            this->data_[i].scaled_fratio = 0;
        this->data_[i].scaled_fsup /= outer_switch;
        this->data_[i].scaled_dx = (this->x_[i+1] - this->x_[i])/
                                   this->data_[i].scaled_fratio;
    }
}

} // namespace detail


} // namespace etf

#endif // ETF_IMPLEMENTATION_HPP

