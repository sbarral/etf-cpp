#ifndef ETF_EXCEPTIONS_HPP
#define ETF_EXCEPTIONS_HPP

#include <stdexcept>


/// Exclusive Top Floor namespace.
///
namespace etf {


/// Exception thrown when an ETF distribution with a N-bit table index is
/// constructed using an `x` table with a size different from 2^N+1.
///
class invalid_table_size : public std::invalid_argument {
public:
    invalid_table_size() : std::invalid_argument("Invalid ETF table size") {}
};

} // namespace etf

#endif // ETF_EXCEPTIONS_HPP
