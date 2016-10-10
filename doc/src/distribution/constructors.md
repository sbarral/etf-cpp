
## Constructors

Distributions are generally constructed from pre-computed tables containing
the partition abcissae as well as the infinum and supremum of the probability
density function over the sub-interval.
Composite distributions also require an outer distribution and, possibly,
the probability density function used for rejection sampling on the outer
interval.

All distributions also declare default, copy and move constructors, which may
or may not be deleted depending on the types of associated functions and outer
distributions.

Note that using a default-constructed distribution through the `min()`, `max()`
or `operator()` members functions is undefined behavior.


### Mandatory arguments

The following arguments are mandated to construct a distribution from partition
data:

 Argument                   | Description
----------------------------|--------------------------------------------------
 `x_first`, `x_last`        | Input iterator range of a sequence of abcissae defining a partition of a bounded interval for which the rectangles making up an upper Riemann sum of `func` have equal areas; the number of elements of the sequence must be equal to 2*ᴺ*+1 where *N* is the third class template parameter
 `finf_first`, `fsup_first` | Input iterators pointing to the beginning of sequences of the infimum and supremum of `func` over each sub-interval of the partition; the number of elements of both sequences must be equal to 2*ᴺ*
 `func`                     | Function or functor proportional to the probability density function of the distribution to be sampled; `func` does not need to be normalized

### Additional argument for composite distributions

The generic 5-parameter distribution types (composite distributions)
additionaly require:

 Argument     | Description
--------------|----------------------------------------------------------------
 `outer_dist` | User-provided distribution for outer interval sampling
 `outer_area` | Total area under `func` over the outer interval

whereas the generic 6-parameter distribution types (rejection-sampled composite
distributions) require:

 Argument     | Description
--------------|----------------------------------------------------------------
 `outer_dist` | User-provided proposal distribution for rejection sampling within the outer interval
 `outer_func` | Majorizing, non-normalized probability distribution function corresponding to the provided `outer_dist`; its scaling must be consistent with that of `func`, i.e. `outer_func` may not be less than `func` anywhere within the outer interval
 `outer_area` | Total area under `outer_func` over the outer interval

Finally, the `symmetric_distribution<...>` family requires as first argument:

 Argument | Description
----------|--------------------------------------------------------------------
`x0`      | Location of the symmetry axis.

### Exceptions

The following exception is thrown if a distribution with a *N*-bit table index
is constructed with a sequence `x_first`, `x_last` which length differs from
2*ᴺ*+1:

```c++
class invalid_table_size : public std::invalid_argument {...};
```


### Constructor declarations for the *distribution<...>* family

(default, copy and move constructors omitted)

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func>
template<typename InputIt1, typename InputIt2, typename InputIt3>
distribution<RealType, W, N, Func>::distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func);
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist>
template<typename InputIt1, typename InputIt2, typename InputIt3>
distribution<RealType, W, N, Func, OuterDist>::distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist,
    RealType outer_area);
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist, typename OuterFunc>
template<typename InputIt1, typename InputIt2, typename InputIt3>
distribution<RealType, W, N, Func, OuterDist, OuterFunc>::distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist, OuterFunc outer_func,
    RealType outer_area);
```



### Constructor declarations for the *central_distribution<...>* family

(default, copy and move constructors omitted)

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func>
template<typename InputIt1, typename InputIt2, typename InputIt3>
central_distribution<RealType, W, N, Func>::central_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func);
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist>
template<typename InputIt1, typename InputIt2, typename InputIt3>
central_distribution<RealType, W, N, Func, OuterDist>::central_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist,
    RealType outer_area);
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist, typename OuterFunc>
template<typename InputIt1, typename InputIt2, typename InputIt3>
central_distribution<RealType, W, N, Func,
                     OuterDist, OuterFunc>::central_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist, OuterFunc outer_func,
    RealType outer_area);
```

### Constructor declarations for the *symmetric_distribution<...>* family

(default, copy and move constructors omitted)

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func>
template<typename InputIt1, typename InputIt2, typename InputIt3>
symmetric_distribution<RealType, W, N, Func>::symmetric_distribution(
    RealType x0,
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func);
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist>
template<typename InputIt1, typename InputIt2, typename InputIt3>
symmetric_distribution<RealType, W, N, Func,
                       OuterDist>::symmetric_distribution(
    RealType x0,
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist,
    RealType outer_area);
```

```c++
template<typename RealType, std::size_t W, std::size_t N, typename Func,
         typename OuterDist, typename OuterFunc>
template<typename InputIt1, typename InputIt2, typename InputIt3>
symmetric_distribution<RealType, W, N, Func,
                       OuterDist, OuterFunc>::symmetric_distribution(
    RealType x0,
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist, OuterFunc outer_func,
    RealType outer_area);
```

