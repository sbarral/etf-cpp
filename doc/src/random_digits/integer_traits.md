
## Integer traits

Unlike the standard C++ library where random variates precision is largely
implementation-defined, the `<etf>` library recognizes the need to control
the trade-off between speed and quality and requires thus the specification of
the number of random digits *W* to be used for random variates generation.

This additional flexibility brings the need for a generic way to select
integers able to hold the generated *W* random bits, a need that is addressed
by the `integer_trait` class:


```c++
template<std::size_t W>
struct integer_traits;
```


### Template parameters

 Parameter   | Description
-------------|-----------------------------------------------------------------
 `W`         | minimum number of digits (in bits) of the member types; if no signed integer or no unsigned integer can be found that can hold *W* bits, a static assertion is triggered (note that the availability of such integers for *W*>64 is platform-dependent)


### Member types

 Type           | Description
----------------|--------------------------------------------------------------
 `int_fast_t`   | A fast signed integer type with at least *W* bits
 `int_least_t`  | A small signed integer type with at least *W* bits
 `uint_fast_t`  | A fast unsigned integer type with at least *W* bits
 `uint_least_t` | A small unsigned integer type with at least *W* bits

