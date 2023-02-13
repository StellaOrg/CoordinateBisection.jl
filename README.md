# CoordinateBisection

Fast Recursive Axis-Aligned Coordinate Bisection for Load Balancing.


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://anicusan.github.io/CoordinateBisection.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anicusan.github.io/CoordinateBisection.jl/dev/)
[![Build Status](https://github.com/anicusan/CoordinateBisection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anicusan/CoordinateBisection.jl/actions/workflows/CI.yml?query=branch%3Amain)


Split a given set of points with optional weights into N balanced axis-aligned domains.
Balanced means either:
- If the weights are `nothing`, the number of points in each domain is equalised.
- If the weights are a `Vector` for each point, the sum of weights on each domain is equalised.
- If the weights are a `Matrix` for each dimension of each 3D point, the sum of coordinate-wise
  weights is equalised.

The resulting `RCB` tree stores each domain's neighbours, each domain's span, and the spatial limits of the points
used to construct it.

The input `points` should have shape `(3, NumPoints)`. When computing the neighbours, an optional
"skin" distance may be specified, such that if any domains are less than "skin" apart (even if they
are not physically touching), they are considered neighbours.

Once the initial `RCB` tree is constructed, it can be iteratively improved in `optimize`
iterations. This is a mathematically terrible problem, so no absolute guarantees can be offered at
the moment as to whether they will _always_ improve the splitting; however, in all our tests - even
with skewed points distributions - it behaves well.


## Examples

Split 10 points into 3 domains:

```julia
using CoordinateBisection
using Random

Random.seed!(123)
points = rand(3, 10)

rcb = RCB{3, Float64}(points)
display(rcb)

#output
RCB{3, Float64}
  x, y, z::NTuple{3, Vector{(NeighborIndex, SplitCoordinate)}}
  xspan, yspan, zspan::NTuple{3, (MinCoordinate, MaxCoordinate)}
  xlimit, ylimit, zlimit::Tuple{MinCoordinate, MaxCoordinate}
```

You can plot results if `Makie` is loaded:

```julia
using GLMakie

# Possible rcbplot arguments:
rcbplot(rcb)
rcbplot(rcb, points)
rcbplot(rcb, points, weights)
```

Interactive plots depicting each optimisation iteration can be shown with `GLMakie`:

```julia
using GLMakie

interactive_rcbplot(points, weights)
```


## References

While this implementation is original - including original weighted point formulation and domain
optimisations - and does not follow any other code, we build upon the ideas of:

> [1] Berger, M.J. and Bokhari, S.H., 1987. A partitioning strategy for nonuniform problems on multiprocessors. IEEE Transactions on Computers, 36(05), pp.570-580.


## License

`CoordinateBisection.jl` is GPL v3.0 licensed.