
## Partitioning

The ETF algorithm requires a pre-processing step whereby a majorizing function
for the probability density function is sought in the form of equal-area
contiguous vertical rectangles sitting on the *x* axis.
The support of these rectangles over the *x* axis defines the ETF partition.

At the moment, the only partition solvers implemented are based on a
multi-variate Newton method:

```c++
template<class Func, class DFunc, class InputIt1, class InputIt2>
partition_data<typename std::iterator_traits<InputIt1>::value_type>
newton_partition(
    Func f,
    DFunc df,
    InputIt1 x_initial_first,
    InputIt1 x_initial_last,
    InputIt2 x_extremum_first,
    InputIt2 x_extremum_last,
    typename std::iterator_traits<InputIt1>::value_type eps,
    typename std::iterator_traits<InputIt1>::value_type relax = 1,
    unsigned int max_iter = 100);
```

```c++
template<class Func, class DFunc, class InputIt1>
partition_data<typename std::iterator_traits<InputIt1>::value_type>
newton_partition_monotonic(
    Func f,
    DFunc df,
    InputIt1 x_initial_first,
    InputIt1 x_initial_last,
    typename std::iterator_traits<InputIt1>::value_type eps,
    typename std::iterator_traits<InputIt1>::value_type relax = 1,
    unsigned int max_iter = 100);
```

The second solver is a specialization for the case of monotonic probability
density functions.

These solvers compute an ETF partition of the interval bounded by the first and
last points of the vector of abcissae passed in argument. The returned object
include the partition as well the infimum and supremum of the function over
each sub-interval.

For faster convergence it is recommended to provide a reasonable initial
estimate of the partition abcissae passed in arguments.

The derivative `df` of `f` must be provided. For non-monotonic probability
density functions, an ordered sequence of the *x* abcissae of inner function
extrema (boundary points excluded) must as well be provided.

The tolerance is the maximum relative dispersion of upper rectangle areas,
computed as the difference between the largest and smallest area relative
to the average area.

An object with empty tables (size 0) is returned if the algorithm fails to
converge.

In order to improve convergence robustness (respectively convergence rate),
successive under-relaxation (respectively over-relaxation) may be optionally
requested by setting the relaxation factor to less (respectively more)
than 1.

The maximum number of iterations may be optionally specified but the default
value is rather conservative and should prove adequate for most cases unless
a very small relaxation factor is used.

The complexity of these solvers is Ο(*N*) where *N* is the number of
sub-intervals. For sufficiently well-behaved functions the convergence rate
is quadratic with respect to tolerance `eps`.


### Arguments

 Argument                              | Description
---------------------------------------|---------------------------------------
 `f`                                   | Function proportional to the probability density function
 `df`                                  | Derivative of `f`
 `x_initial_first`, `x_initial_last`   | Input iterator range of an ordered sequence of abcissae defining a partition of a bounded interval for which the rectangles making up an upper Riemann sum of `func` have, preferably, approximately equal areas
 `x_extremum_first`, `x_extremum_last` | Input iterator range of an ordered sequence of the abcissae *x* at which `f` has a local extremum, excluding boundary points
 `eps`                                 | Tolerance, defined as the maximum dispersion of upper rectangle areas relatively to the average rectangle area
 `relax`                               | Relaxation factor for the iterative Newton solver
 `max_iter`                            | Maximum number of iteration of the Newton method before giving up


### Return value

The return value is a struct with 3 public member variables:

```c++
template<typename RealType>
struct partition_data
{
    std::vector<RealType> x;
    std::vector<RealType> finf;
    std::vector<RealType> fsup;
};
```

Upon successful convergence of the solver, the members of the return value
satisfy the following post-conditions:

* `size(x)==std::distance(x_initial_first, x_initial_last)`

* `size(finf)==size(x)-1`

* `size(fsup)==size(x)-1`

If the solver fails to converge after `max_iter` iterations, however, the
post-conditions become:

* `size(x)==0`

* `size(finf)==0`

* `size(fsup)==0`

 Member variable | Description
-----------------|-------------------------------------------------------------
 `x`             | Sequence of abcissae defining the ETF partition
 `finf`          | Sequence of infima for each of the partition sub-interval
 `fsup`          | Sequence of suprema for each of the partition sub-interval

