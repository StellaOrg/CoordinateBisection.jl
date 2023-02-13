using CoordinateBisection
using Test

using Random


@testset "RCB - Basic" begin
    Random.seed!(123)

    points = rand(3, 10)
    rcb = RCB{3, Float64}(points, nothing, 0., 0)

    @test ndims(rcb) == 3
    @test eltype(rcb) == Float64
    @test length(rcb.x) == 3
    @test length(rcb.xspan) == 3
    @test length(rcb.xlimit) == 2
end


# @testset "domain_weights"
# @testset "weighted_quantile"
# @testset "interpolate"
