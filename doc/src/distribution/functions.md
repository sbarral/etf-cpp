
## Non-member functions

These convenience functions create distribution objects, deducing trailing
argument types to select the distribution type.

The required explicit template parameters are `RealType`, `W` and `N` (see
the [distributions overview](distribution/overview.html) for template
parameters and function arguments descriptions).


### Non-member function declarations for the *distribution<...>* family

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func>
etf::distribution<RealType, W, N, Func>
make_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func);
```

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func, typename OuterDist>
etf::distribution<RealType, W, N, Func, OuterDist>
make_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist,
    RealType outer_area);
```

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func, typename OuterDist, typename OuterFunc>
etf::distribution<RealType, W, N, Func, OuterDist, OuterFunc>
make_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist, OuterFunc outer_func,
    RealType outer_area);
```


### Non-member function declarations for the *central_distribution<...>* family

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func>
etf::central_distribution<RealType, W, N, Func>
make_central_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func);
```

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func, typename OuterDist>
etf::central_distribution<RealType, W, N, Func, OuterDist>
make_central_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist,
    RealType outer_area);
```

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func, typename OuterDist, typename OuterFunc>
etf::central_distribution<RealType, W, N, Func, OuterDist, OuterFunc>
make_central_distribution(
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist, OuterFunc outer_func,
    RealType outer_area);
```

### Non-member function declarations for the *symmetric_distribution<...>* family

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func>
etf::symmetric_distribution<RealType, W, N, Func>
make_symmetric_distribution(
    RealType x0,
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func);
```

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func, typename OuterDist>
etf::symmetric_distribution<RealType, W, N, Func, OuterDist>
make_symmetric_distribution(
    RealType x0,
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist,
    RealType outer_area);
```

```c++
template<typename RealType, std::size_t W, std::size_t N,
         typename InputIt1, typename InputIt2, typename InputIt3,
         typename Func, typename OuterDist, typename OuterFunc>
etf::symmetric_distribution<RealType, W, N, Func, OuterDist, OuterFunc>
make_symmetric_distribution(
    RealType x0,
    InputIt1 x_first, InputIt1 x_last,
    InputIt2 finf_first,
    InputIt3 fsup_first,
    Func func,
    OuterDist outer_dist, OuterFunc outer_func,
    RealType outer_area);
```

