# *\<etf\>*: a C++ implementation of the Exclusive Top Floor algorithm

The [ETF algorithm](https://sbarral.github.io/etf) is a blazingly fast
algorithm to sample arbitrary univariate continuous probability distributions
solely from their probability density function and using marginally more than
one random number per sample.

The *\<etf\>* library is a ready-to-use implementation of the ETF algorithm with
optimized specializations for asymmetric, symmetric and central distributions.


## Documentation

The *\<etf\>* library API is documented [here](https://sbarral.github.io/etf-cpp).
You may read more about the algorithm itself in the
[online overview](https://sbarral.github.io/etf) of the method.


## Installation

This is a header-only library: once installed locally, it is enough to make
sure that the *etf* directory is in the include path.


## Benchmarking

The [benchmarking directory](benchmark) contains a
[timing benchmark](benchmark/timing.cpp) as well as a
[statistical collision test benchmark](benchmark/collision.cpp).
These are discussed in details in the
[ETF overview](https://sbarral.github.io/etf).

Note that the timing benchmark requires the
[nonius](https://nonius.io) microbenchmarking library which in turn needs a
couple of the Boost libraries.

## Examples

A short tutorial will be coming soon. In the meantime, the
[benchmarking directory](benchmark) contains example applications of the
*\<etf\>* library for the fast generation of normal and chi-squared variates.


## License

This software is licensed under the
[Apache License version 2.0](LICENSE-APACHE), the
[MIT license](LICENSE-MIT) or the
[Boost license version 1.0](LICENSE-BOOST), at your option.

Copyright © 2016 Serge Barral.

