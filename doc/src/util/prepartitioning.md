
## Pre-partitioning

The ETF algorithm requires a pre-processing step whereby a majorizing function
for the probability density function is sought in the form of equal-area
contiguous vertical rectangles sitting on the *x* axis.
The support of these rectangles over the *x* axis defines the ETF partition.

Pre-partitioning consists in finding a first approximation of the ETF partition
of the *x* axis to be used by the final partition solver.

At the moment, a trivial (but fast) pre-partitioning algorithm based on
the trapezoidal rule is available:

```c++
template<typename RealType, class Func>
std::vector<RealType>
trapezoidal_rule_prepartition(Func f,
                              RealType x0,
                              RealType x1,
                              std::size_t nb_intervals,
                              std::size_t nb_points);
```

```c++
template<typename RealType, class Func>
std::vector<RealType>
trapezoidal_rule_prepartition(Func f,
                              RealType x0,
                              RealType x1,
                              std::size_t nb_intervals);
```

where the second overload is equivalent to setting `nb_points` to
`nb_partitions`.

These functions apply the trapezoidal rule to the supplied function `f` over a
regular grid with `nb_points` grid points (including boundary nodes) to divide
interval [`x0`, `x1`] into the specified number of sub-intervals such that the
areas under the trapeze quadrature are the same in each sub-interval.  The
returned vector is the set of abscissae.


### Arguments

 Argument        | Description
-----------------|-----------------------------------------------------------------
 `f`             | Function proportional to the probability density function
 `x0`, `x1`      | Boundaries of the interval
 `nb_intervals`  | Number of sub-intervals in the generated partition
 `nb_points`     | Number of points (including boundary nodes) to be used for the trapezoidal rule


### Return value

The set of abcissae defining the partition.
