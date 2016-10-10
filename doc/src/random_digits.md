
# <etf/random_digits.hpp>

The `<etf/random_digits.hpp>` header contains fixed-precision random number
generation functions and fixed-precision integer traits.

These facilities are used internally in the ETF algorithm implementation but may
also be useful for the implementation of user-supplied outer distributions. In
particular, the `generate_real()` function provides an efficient and sound(er)
replacement for the standard library's `generate_canonical()` function.
