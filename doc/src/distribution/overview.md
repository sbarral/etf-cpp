
## Distributions overview

Three distribution type families are provided:

* `distribution<...>` for arbitrarily shaped distributions,

* `central_distribution<...>` for distributions that are symmetric about *x*=0,

* `symmetric_distribution<...>` for distributions that are symmetric about some
  arbitrary abcissa *x*.


### Template parameters

All distribution families accept 4 to 6 explicit template parameters:

* the 4-parameters distribution types are for distributions defined over a
  bounded interval \[*x₀*, *x₁*\],

* the 5-parameters distribution types are composite distributions: the ETF
  algorithm is used to sample the distribution over a bounded interval while a
  user-supplied outer distribution is used to sample over an arbitrary outer
  interval (e.g. an infinite tail),

* the 6-parameters distributions are composite distributions which use
  rejection sampling over the outer interval; the user must provide an
  appropriately scaled probability function that majorizes the reference
  probability density function over the outer interval, together with the
  corresponding proposal outer distribution.

The template parameters are provided in the below table:

 Parameter   | Description
-------------|-----------------------------------------------------------------
 `RealType`  | Floating point type
 `W`         | Number of digits (in bits) of random numbers used for the generation of distribution variates; if the number of digits does not match that of the RNG passed in argument when requesting a variate, the elementary random numbers produced by the RNG will be appropriately concatenated or truncated; for this reason, it is most efficient to set W equal to the number of digits produced by the RNG passed to `operator()`, or to a multiple thereof
 `N`         | Number of digits (in bits) of the lookup table index; lookup tables have 2*ᴺ* elements
 `Func`      | Probability density function type; instances of `Func` need *not* be normalized and are only required to be proportional to the normalized probability density function
 `OuterDist` | Distribution type used when sampling over the outer interval [*only for composite distributions*]
 `OuterFunc` | Type of the non-normalized majorizing probability distribution in case the outer distribution is used in conjunction with rejection sampling [*only for composite distributions with rejection sampling over the outer interval*]




### Class declarations

Note that only template parameters to be explicitly provided are enlisted.
Even though the template definition of each distribution family counts 6
template parameters, the 5-th and 6-th parameters default to `void` and the
proper template specialization is selected based on the non-void parameters.

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func>
class distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist>
class distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist, typename OuterFunc>
class distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func>
class central_distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist>
class central_distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist, typename OuterFunc>
class central_distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func>
class symmetric_distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist>
class symmetric_distribution;
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist, typename OuterFunc>
class symmetric_distribution;
```
