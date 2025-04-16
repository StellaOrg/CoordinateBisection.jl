# Optional, loaded only if GLMakie is available
using .Makie


# Optional exports
export interactive_rcbplot


function box_quads(lower, upper)
    corners = [
        lower[1] lower[2] lower[3]
        upper[1] lower[2] lower[3]
        upper[1] upper[2] lower[3]
        lower[1] upper[2] lower[3]

        lower[1] lower[2] upper[3]
        upper[1] lower[2] upper[3]
        upper[1] upper[2] upper[3]
        lower[1] upper[2] upper[3]
    ]

    quad_faces = [
        1 2 3 4
        1 2 6 5
        2 3 7 6
        3 4 8 7
        4 1 5 8
        5 6 7 8
    ]
    
    corners, quad_faces
end


function side_colors(rp, i)

    rcb = rp[1][]
    colors = Vector{Any}(undef, 6)

    colors[1] = isinf(rcb.zspan[i][1]) ? rp.outer[] : rp.inner[]
    colors[2] = isinf(rcb.yspan[i][1]) ? rp.outer[] : rp.inner[]
    colors[3] = isinf(rcb.xspan[i][2]) ? rp.outer[] : rp.inner[]
    colors[4] = isinf(rcb.yspan[i][2]) ? rp.outer[] : rp.inner[]
    colors[5] = isinf(rcb.xspan[i][1]) ? rp.outer[] : rp.inner[]
    colors[6] = isinf(rcb.zspan[i][2]) ? rp.outer[] : rp.inner[]

    colors
end


function domain_box(rp, i)

    rcb = rp[1][]

    lower = [
        isinf(rcb.xspan[i][1]) ? rcb.xlimit[1] : rcb.xspan[i][1]
        isinf(rcb.yspan[i][1]) ? rcb.ylimit[1] : rcb.yspan[i][1]
        isinf(rcb.zspan[i][1]) ? rcb.zlimit[1] : rcb.zspan[i][1]
    ]
    upper = [
        isinf(rcb.xspan[i][2]) ? rcb.xlimit[2] : rcb.xspan[i][2]
        isinf(rcb.yspan[i][2]) ? rcb.ylimit[2] : rcb.yspan[i][2]
        isinf(rcb.zspan[i][2]) ? rcb.zlimit[2] : rcb.zspan[i][2]
    ]

    corners, faces = box_quads(lower, upper)
    colors = side_colors(rp, i)

    corners, faces, colors
end


# Create Makie plotting recipe
@recipe(RCBPlot) do scene
    Attributes(
        inner = :blue,
        outer = :grey,
        alpha = 0.1,
        marker = theme(scene, :marker),
        markersize = theme(scene, :markersize),
    )
end


# Different signatures for rcbplot, with 1 / 2 / 3 optional arguments
function Makie.plot!(rp::RCBPlot{<:Tuple{RCB}})
    plot(rp, rp[1], nothing, nothing)
end


function Makie.plot!(rp::RCBPlot{<:Tuple{RCB, Union{Nothing, AbstractArray}}})
    plot(rp, rp[1], rp[2], nothing)
end


function Makie.plot!(rp::RCBPlot{<:Tuple{RCB,
                                         Union{Nothing, AbstractArray},
                                         Union{Nothing, AbstractArray}}})
    plot(rp, rp[1], rp[2], rp[3])
end


function plot(rp, rcb, points, weights)
    # Plot each domain as a box
    N = ndims(rcb[])
    for i in 1:N
        corners, faces, colors = domain_box(rp, i)

        # Plot individual faces
        for j in 1:6
            mesh!(rp,
                  @view(corners[@view(faces[j, :]), :]),
                  [1 2 3; 3 4 1],
                  color = (colors[j], rp.alpha[]),
                  transparency = true)
        end
    end

    # Plot points if given
    if !isnothing(points)
        points = points[]
        isnothing(weights) || (weights = weights[])

        @assert ndims(points) == 2
        @assert size(points, 1) == 3

        # Make sizes proportional to weights if given
        if isnothing(weights)
            markersize = rp.markersize
        elseif ndims(weights) == 1
            @assert size(weights, 1) == size(points, 2)
            markersize = rp.markersize[] * weights / maximum(weights)
        elseif ndims(weights) == 2
            @assert size(weights, 2) == size(points, 2)
            merged = weights[1, :] + weights[2, :] + weights[3, :]
            markersize = rp.markersize[] * merged / maximum(merged)
        end

        scatter!(rp, points', markersize=markersize, marker=rp.marker)
    end

    rp
end


function interactive_rcbplot(
    points,
    weights,
    ::Val{N},
    maxopt=10N;
    verbose=true,
    kwargs...
) where N
    # Extract type of elements in the points array
    T = eltype(points)

    # Create first RCB for 0 optimisations
    rcb0 = RCB{N, T}(points, weights, zero(T), 0)
    if verbose
        println("Domain weights for 0 optimisations:")
        display(domain_weights(rcb0, points, weights))
    end

    # Create figure and first plot of RCB
    fig = Figure()
    ax, plot_obj = rcbplot(fig[1, 1], rcb0, points, weights; kwargs...)

    # Create slider grid for number of optimisation passes to do
    sgrid = SliderGrid(
        fig[2, 1],
        (label="Optimisations", range=0:1:maxopt, startvalue=0),
    )

    # Recompute RCB and replace plot when changing the slider grid values
    on(sgrid.sliders[1].value) do nopt
        delete!(ax, plot_obj)

        rcb = RCB{N, T}(points, weights, zero(T), nopt)

        if verbose
            println("Domain weights for $nopt optimisations:")
            display(domain_weights(rcb, points, weights))
        end

        plot_obj = rcbplot!(ax, rcb, points, weights; kwargs...)
    end

    fig
end
