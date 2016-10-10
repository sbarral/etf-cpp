#ifndef ETF_RANDOM_DIGITS_HPP
#define ETF_RANDOM_DIGITS_HPP

#include <cstdint>
#include <limits>



/// Exclusive Top Floor namespace.
///
namespace etf {

namespace detail {

// Compile-time computation of the bit width.
template<typename UIntType>
constexpr std::size_t bit_width(UIntType value, std::size_t shift = 0) {
    return value==UIntType(0) ? shift : bit_width(value >> 1, shift + 1);
}


// Number of significant digits of a random number generator with range 2^N.
template<typename RngType>
constexpr std::size_t rng_digits() {
    return bit_width(RngType::max() - RngType::min());
}


// Check that a number can be written as 2^N-1, i.e. all lower bits are set.
template<typename UIntType>
constexpr bool is_power_of_2_less_1(UIntType value) {
    return value==1 ? true
         : ((value & 1)==0 ? false
         : is_power_of_2_less_1(value >> 1));
}



// Check the range of the random number generator.
//
// For efficiency reason, many operations assume that the range of the random
// generator is [0, 2^N-1] so this is checked at compile-time.
// As a small, pragmatic exception to this rule, random number generators which
// never produce 0 and/or 2^N-1 are tolerated.
template<typename RngType>
constexpr std::size_t check_rng_range() {
    using RngResultType = decltype(RngType::min());
    static_assert( RngType::min()==0 || RngType::min()==1,
                  "random number generator min value is not 0 or 1");
    static_assert(   is_power_of_2_less_1(RngType::max())
                  || is_power_of_2_less_1(RngType::max() | RngResultType(1)),
                  "random number generator max value is not 2^N-1 or 2^N-2");
    return 0; // meaningless, C++11 does not accept void return type
}


// Generate P random bits.
template<typename UIntType, std::size_t P>
struct random_bits
{
    template<typename RngType>
    static UIntType generate(RngType& rng) {
        constexpr std::size_t N = rng_digits<RngType>();
        constexpr std::size_t R = P>N ? P-N : 0;
        constexpr std::size_t S = N>P ? N-P : 0;
        return random_bits<UIntType, R>::generate_next(rng,
            (static_cast<UIntType>(rng() >> S)));
    }

    template<typename RngType>
    static UIntType generate_next(RngType& rng, UIntType u) {
        constexpr std::size_t N = rng_digits<RngType>();
        constexpr std::size_t R = P>N ? P-N : 0;
        constexpr std::size_t S = N>P ? N-P : 0;
        constexpr std::size_t T = N<P ? N : P;
        return random_bits<UIntType, R>::generate_next(rng,
            (u << T) | (static_cast<UIntType>(rng()) >> S));
    }
};

template<typename UIntType>
struct random_bits<UIntType, 0>
{
    template<typename RngType>
    static UIntType generate_next(RngType& rng, UIntType u) {
        return u;
    }
};


struct BITFIELD_TOO_WIDE_ERROR;

template<std::size_t N>
struct bitfield_category
{
    using int_fast_t = BITFIELD_TOO_WIDE_ERROR;
    using int_least_t = BITFIELD_TOO_WIDE_ERROR;
    using uint_fast_t = BITFIELD_TOO_WIDE_ERROR;
    using uint_least_t = BITFIELD_TOO_WIDE_ERROR;
};

template<>
struct bitfield_category<1>
{
    using int_fast_t = std::int_fast16_t;
    using int_least_t = std::int_least16_t;
    using uint_fast_t = std::uint_fast16_t;
    using uint_least_t = std::uint_least16_t;
};

template<>
struct bitfield_category<2>
{
    using int_fast_t = std::int_fast32_t;
    using int_least_t = std::int_least32_t;
    using uint_fast_t = std::uint_fast32_t;
    using uint_least_t = std::uint_least32_t;
};

template<>
struct bitfield_category<3>
{
    using int_fast_t = std::int_fast64_t;
    using int_least_t = std::int_least64_t;
    using uint_fast_t = std::uint_fast64_t;
    using uint_least_t = std::uint_least64_t;
};

template<>
struct bitfield_category<4>
{
    using int_fast_t = long long;
    using int_least_t = long long;
    using uint_fast_t = unsigned long long;
    using uint_least_t = unsigned long long;
};

template<>
struct bitfield_category<5>
{
    using int_fast_t = std::intmax_t;
    using int_least_t = std::intmax_t;
    using uint_fast_t = std::uintmax_t;
    using uint_least_t = std::uintmax_t;
};

} // end namespace detail


/// Traits specifying convenient integer types able to hold W bits.
///
template<std::size_t W>
struct integer_traits
{
private:
    static const int value =
           W<=16 ? 1
        : (W<=32 ? 2
        : (W<=64 ? 3
        : ((  W<=std::numeric_limits<long long>::digits
           && W<=std::numeric_limits<unsigned long long>::digits) ? 4
        : ((  W<=std::numeric_limits<std::intmax_t>::digits
           && W<=std::numeric_limits<std::uintmax_t>::digits) ? 5
        : 0)))); // invalid bit width, category value 0 generates an error

public:
    using int_fast_t =
        typename detail::bitfield_category<value>::int_fast_t;
    using int_least_t =
        typename detail::bitfield_category<value>::int_least_t;
    using uint_fast_t =
        typename detail::bitfield_category<value>::uint_fast_t;
    using uint_least_t =
        typename detail::bitfield_category<value>::uint_least_t;
};


/// Generate a W-bit precision random integer equidistributed in [0, 2^W-1].
///
/// As many random numbers are generated as necessary to fill W random bits.
///
template<typename UIntType, std::size_t W, typename RngType>
inline
UIntType generate_random_integer(RngType& rng) {
    detail::check_rng_range<RngType>();
    return detail::random_bits<UIntType, W>::generate(rng);
}


/// Generate a W-bit precision floating point value in [0,1).
///
/// This a drop-in relacement for std::generate_canonical with the following
/// caveats:
///  * it is faster,
///  * it is hopefully less buggy, meaning it should never return 1.0
///    and should have a better bin distribution quality,
///  * it is less flexible: the RNG ideally needs to generate values between 0
///    and a power of 2 less 1, though a minimum value of 1 and/or a maximum
///    value equal to a power of 2 less 2 are tolerated,
///  * it is less generic: the platform must have an unsigned type able to
///    store W bits or the number of significant bits of the floating point
///    type, whichever is less.
/// 
/// Just like std::generate_canonical, it may actually generate less bits than
/// the W bits requested if the floating point type has less that W significant
/// bits. In any event, it will generate as many random numbers as necessary to
/// create the floating point number.
///
template<typename RealType, std::size_t W, typename RngType>
inline
RealType generate_random_real(RngType& rng) {
    constexpr std::size_t N = std::numeric_limits<RealType>::digits;
    constexpr std::size_t M = N<W ? N : W;
    using UIntType = typename integer_traits<M>::uint_fast_t;
    constexpr RealType S = RealType(1)/(RealType(UIntType(1) << M/2)
                                       *RealType(UIntType(1) << (M - M/2)));
    return S*static_cast<RealType>(generate_random_integer<UIntType, M>(rng));
}


} // namespace etf

#endif // ETF_RANDOM_DIGITS_HPP

