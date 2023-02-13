

using Random
using Distributions
using CoordinateBisection


Random.seed!(123)


# Generate 8 points in a box
points = [0 0 0
          1 0 0
          0 1 0
          0 0 1
          1 1 0
          1 0 1
          0 1 1
          1 1 1]
points = Matrix{Float64}(points')


n = 5
np = 16n

# points = [
#     rand(Distributions.Chi(1), np)'
#     rand(Distributions.Normal(), np)'
#     rand(Distributions.Beta(), np)'
# ]

points = rand(3, np)

# weights = nothing
# weights = rand(size(points, 2))
weights = rand(size(points)...)


rcb = RCB{n, Float64}(points, weights, 0.1, 10)
display(rcb)
println(rcb)

println("Domain weights:")
display(domain_weights(rcb, points, weights))


using GLMakie
# rcbplot(rcb, points, weights)
interactive_rcbplot(points, weights, Val(n), 5n; markersize=20)

