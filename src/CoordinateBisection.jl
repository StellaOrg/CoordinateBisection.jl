module CoordinateBisection

# Exports
export RCB, domain_weights


# Internal imports
using StaticArrays
using Requires
using DocStringExtensions


# Optional code
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("plotting.jl")
end


"""
    $(TYPEDEF)

Split a given set of points with optional weights into N balanced axis-aligned domains.
Balanced means:
- If the weights are `nothing`, the number of points in each domain is equalised.
- If the weights are a `Vector` for each point, the sum of weights on each domain is equalised.
- If the weights are a `Matrix` for each dimension of each 3D point, the sum of coordinate-wise
  weights is equalised.

The resulting `RCB` tree stores each domain's neighbours in the `x`, `y` and `z` members, each
domain's span in the `xspan`, `yspan` and `zspan` members, and the spatial limits of the points
used to construct it in the `xlimit`, `ylimit` and `zlimit` members.

The input `points` should have shape `(3, NumPoints)`. When computing the neighbours, an optional
"skin" distance may be specified, such that if any domains are less than "skin" apart (even if they
are not physically touching), they are considered neighbours.

Once the initial `RCB` tree is constructed, it can be iteratively improved in `optimize`
iterations. This is a mathematically terrible problem, so no absolute guarantees can be offered at
the moment as to whether they will _always_ improve the splitting; however, in all our tests - even
with skewed points distributions - it behaves well.

# Fields
    $(TYPEDFIELDS)


# Methods
    RCB{N, T}(points, weights=nothing, skin=T(0.), optimize=2N) where {N, T <: Real}
    ndims(::RCB{N, T}) where {N, T} = N
    eltype(::RCB{N, T}) where {N, T} = T


# Examples

Split 10 points into 3 domains:

```julia
using CoordinateBisection
using Random

Random.seed!(123)
points = rand(3, 10)

rcb = RCB{3, Float64}(points)
```

You can plot results if `Makie` is loaded:

```julia
using GLMakie

# Possible rcbplot arguments:
rcbplot(rcb)
rcbplot(rcb, points)
rcbplot(rcb, points, weights)
```

Interactive plots showing optimisation iterations can be shown with `GLMakie`:

```julia
using GLMakie

interactive_rcbplot(points, weights)
```

"""
struct RCB{N, T <: Real}
    "Each domain's neighbours in the X/Y/Z dimension, as `Tuple{DomainIndex, SplitCoordinate}`"
    x::NTuple{N, Vector{Tuple{Int64, T}}}
    y::NTuple{N, Vector{Tuple{Int64, T}}}
    z::NTuple{N, Vector{Tuple{Int64, T}}}

    "Each domain's span in the X/Y/Z dimension, as `Tuple{MinCoordinate, MaxCoordinate}`"
    xspan::NTuple{N, Tuple{T, T}}
    yspan::NTuple{N, Tuple{T, T}}
    zspan::NTuple{N, Tuple{T, T}}

    "Limits of points over which the RCB was constructed, as `Tuple{MinCoordinate, MaxCoordinate}`"
    xlimit::Tuple{T, T}
    ylimit::Tuple{T, T}
    zlimit::Tuple{T, T}
end


Base.ndims(::RCB{N, T}) where {N, T} = N
Base.eltype(::RCB{N, T}) where {N, T} = T


function Base.show(io::IO, ::MIME"text/plain", rcb::RCB{N, T}) where {N, T}
    # Multiline display
    print(io, """
              RCB{$N, $T}
                x, y, z::NTuple{$N, Vector{(NeighborIndex, SplitCoordinate)}}
                xspan, yspan, zspan::NTuple{$N, (MinCoordinate, MaxCoordinate)}
                xlimit, ylimit, zlimit::Tuple{MinCoordinate, MaxCoordinate}
              """)
end


function Base.show(io::IO, rcb::RCB{N, T}) where {N, T}
    # Short single-line display
    print(io, "RCB{$N, $T}(x, y, z, xspan, yspan, zspan, xlimit, ylimit, zlimit)")
end


function RCB{N, T}(points, weights=nothing, skin=T(0.), optimize=5N) where {N, T <: Real}
    # Error checking
    N >= 1 || throw(ArgumentError("N (the number of domains) must be larger or equal to 1"))
    @assert ndims(points) == 2
    @assert size(points, 1) == 3
    if !isnothing(weights)
        @assert 1 <= ndims(weights) <= 2
        ndims(weights) == 1 && @assert size(weights, 1) == size(points, 2)
        ndims(weights) == 2 && @assert size(weights, 2) == size(points, 2)
    end

    # TODO: no domain can have a size smaller than skin
    # TODO: any smarter way to ensure correct splitting for very few points?
    @assert size(points, 2) >= 2N "Must have at least 2N points for N domains"

    # Find limits of points
    xlimit, ylimit, zlimit = points_limits(points)

    # Trivial case for a single node: no neighbours
    N == 1 && return RCB{N, T}((Tuple{Int64, T}[],),
                               (Tuple{Int64, T}[],),
                               (Tuple{Int64, T}[],),
                               ((typemin(T), typemax(T)),),
                               ((typemin(T), typemax(T)),),
                               ((typemin(T), typemax(T)),),
                               xlimit, ylimit, zlimit)

    # Indices for each dimension
    xdim = 1
    ydim = 2
    zdim = 3

    # Encode binary tree as a parent vector: each element in the vector contains its parent's index
    num_leaves = N
    num_nodes = 2N - 1

    parents = MVector{num_nodes, Int64}(undef)                      # Indices of parents
    leaves = MVector{num_leaves, Int64}(undef)                      # Indices of leaves
    spans = MVector{num_nodes, NTuple{3, Tuple{T, T}}}(undef)       # Spatial span of each node

    # Root node encompasses all Cartesian space - (-Inf, +Inf)
    parents[1] = 0
    leaves[1] = 1
    spans[1] = ((typemin(T), typemax(T)),
                (typemin(T), typemax(T)),
                (typemin(T), typemax(T)))

    # Need to do N - 1 splits
    for i in 2:N

        # Go over each leaf to find the one with maximum total weight
        bestleaf = 1
        besttotal = typemin(T)

        for j in 1:i - 1
            xtotal, ytotal, ztotal = leaf_weight(points, weights, spans[leaves[j]])

            if xtotal > besttotal
                bestleaf = j
                besttotal = xtotal
            end
            if ytotal > besttotal
                bestleaf = j
                besttotal = ytotal
            end
            if ztotal > besttotal
                bestleaf = j
                besttotal = ztotal
            end
        end

        # Found leaf to split; now find dimension with largest weighted 5-95% distribution
        isplit = leaves[bestleaf]

        xspread = leaf_spread(points, weights, spans[isplit], xdim)
        yspread = leaf_spread(points, weights, spans[isplit], ydim)
        zspread = leaf_spread(points, weights, spans[isplit], zdim)

        bestdim = xdim
        bestspread = xspread

        if yspread > bestspread
            bestspread = yspread
            bestdim = ydim
        end
        if zspread > bestspread
            bestspread = zspread
            bestdim = zdim
        end

        # Now split `leaves[bestleaf]` along `bestdim` by its weighted median
        if bestdim == xdim
            xsplit = leaf_median(points, weights, spans[isplit], xdim)
            spans[2i - 2] = ((spans[isplit][xdim][1], xsplit),
                             (spans[isplit][ydim][1], spans[isplit][ydim][2]),
                             (spans[isplit][zdim][1], spans[isplit][zdim][2]))

            spans[2i - 1] = ((xsplit, spans[isplit][xdim][2]),
                             (spans[isplit][ydim][1], spans[isplit][ydim][2]),
                             (spans[isplit][zdim][1], spans[isplit][zdim][2]))
        end

        if bestdim == ydim
            ysplit = leaf_median(points, weights, spans[isplit], ydim)
            spans[2i - 2] = ((spans[isplit][xdim][1], spans[isplit][xdim][2]),
                             (spans[isplit][ydim][1], ysplit),
                             (spans[isplit][zdim][1], spans[isplit][zdim][2]))

            spans[2i - 1] = ((spans[isplit][xdim][1], spans[isplit][xdim][2]),
                             (ysplit, spans[isplit][ydim][2]),
                             (spans[isplit][zdim][1], spans[isplit][zdim][2]))
        end

        if bestdim == zdim
            zsplit = leaf_median(points, weights, spans[isplit], zdim)
            spans[2i - 2] = ((spans[isplit][xdim][1], spans[isplit][xdim][2]),
                             (spans[isplit][ydim][1], spans[isplit][ydim][2]),
                             (spans[isplit][zdim][1], zsplit))

            spans[2i - 1] = ((spans[isplit][xdim][1], spans[isplit][xdim][2]),
                             (spans[isplit][ydim][1], spans[isplit][ydim][2]),
                             (zsplit, spans[isplit][zdim][2]))
        end

        # Sprout two leaves: set their parent's index
        parents[2i - 2] = isplit
        parents[2i - 1] = isplit

        # Indices of current active leaves: replace current domain we split and add a new one
        leaves[bestleaf] = 2i - 2
        leaves[i] = 2i - 1
    end

    # Extract leaf spans into dimension-wise spans, e.g. [(xmin, xmax), ...]
    xspan = MVector{N, Tuple{T, T}}(undef)
    yspan = MVector{N, Tuple{T, T}}(undef)
    zspan = MVector{N, Tuple{T, T}}(undef)

    for i in 1:N
        xspan[i] = spans[leaves[i]][1]
        yspan[i] = spans[leaves[i]][2]
        zspan[i] = spans[leaves[i]][3]
    end

    # Pre-allocate vectors of neighbours
    x = Vector{Vector{Tuple{Int64, T}}}(undef, N)
    y = Vector{Vector{Tuple{Int64, T}}}(undef, N)
    z = Vector{Vector{Tuple{Int64, T}}}(undef, N)

    # Iteratively improve domain sizes
    # TODO: optimisation with weights doesn't do anything?

    if N >= 3 && optimize > 0
        # Ideal weight per domain is equally-split
        ideal_totals = compute_ideal_totals(points, weights, N, T)

        dimspans = (xspan, yspan, zspan)
        totals = MVector{N, Tuple{T, T, T}}(undef)      # Dimension-wise weights sum inside leaves

        # Keep track of number of improvements done on each domain; reuse `leaves` MVector
        num_adjusts = leaves
        num_adjusts .= 0

        for _ in 1:optimize

            # Find each domain's neighbours and compute sum of weights inside them
            find_neighbors!(x, y, z, xspan, yspan, zspan)
            evaluate_leaves!(totals, points, weights, xspan, yspan, zspan)

            # Adjust new domains with worst weights compared to ideals, one per dimension
            xleaf, yleaf, zleaf = worst_totals(totals, ideal_totals, num_adjusts)

            num_adjusts[xleaf] += 1
            num_adjusts[yleaf] += 1
            num_adjusts[zleaf] += 1

            adjust_leaf!(points, weights, dimspans, x[xleaf], totals, xdim, xleaf)
            adjust_leaf!(points, weights, dimspans, y[yleaf], totals, ydim, yleaf)
            adjust_leaf!(points, weights, dimspans, z[zleaf], totals, zdim, zleaf)
        end
    end

    # Find each domain's neighbours, including a skin distance / tolerance between them
    find_neighbors!(x, y, z, xspan, yspan, zspan, skin)

    RCB{N, T}(
        x |> NTuple{N, Vector{Tuple{Int64, T}}},
        y |> NTuple{N, Vector{Tuple{Int64, T}}},
        z |> NTuple{N, Vector{Tuple{Int64, T}}},
        xspan |> NTuple{N, Tuple{T, T}},
        yspan |> NTuple{N, Tuple{T, T}},
        zspan |> NTuple{N, Tuple{T, T}},
        xlimit, ylimit, zlimit,
    )
end


"""
    domain_weights(
        rcb::RCB{N, T},
        points,
        weights=nothing,
    ) where {N, T} -> MVector{N, Tuple{T, T, T}}

Compute the total points' weight in each domain of `rcb`.
"""
function domain_weights(rcb::RCB{N, T}, points, weights=nothing) where {N, T}
    @assert ndims(points) == 2
    @assert size(points, 1) == 3
    if !isnothing(weights)
        @assert 1 <= ndims(weights) <= 2
        ndims(weights) == 1 && @assert size(weights, 1) == size(points, 2)
        ndims(weights) == 2 && @assert size(weights, 2) == size(points, 2)
    end

    # Dimension-wise weights sum inside leaves
    totals = MVector{N, Tuple{T, T, T}}(undef)
    evaluate_leaves!(totals, points, weights, rcb.xspan, rcb.yspan, rcb.zspan)

    totals
end


function points_limits(points)
    xmin = ymin = zmin = typemax(eltype(points))
    xmax = ymax = zmax = typemin(eltype(points))

    for i in 1:size(points, 2)
        point = @view points[:, i]

        point[1] < xmin && (xmin = point[1])
        point[1] > xmax && (xmax = point[1])

        point[2] < ymin && (ymin = point[2])
        point[2] > ymax && (ymax = point[2])

        point[3] < zmin && (zmin = point[3])
        point[3] > zmax && (zmax = point[3])
    end

    (xmin, xmax), (ymin, ymax), (zmin, zmax)
end


function leaf_weight(points, weights, span::NTuple{3, Tuple{T, T}}) where T
    # points: (3, NP)
    # weights: nothing / (NP,) / (3, NP)
    # span: ((xmin, xmax), (ymin, ymax), (zmin, zmax))

    xtotal = zero(T)
    ytotal = zero(T)
    ztotal = zero(T)

    for i in 1:size(points, 2)
        point = @view points[:, i]

        # Select only points within this leaf's box
        contained = (span[1][1] <= point[1] < span[1][2] &&
                     span[2][1] <= point[2] < span[2][2] &&
                     span[3][1] <= point[3] < span[3][2])
        
        if contained
            # If we don't have any weights, just count number of elements in each leaf
            if isnothing(weights)
                xtotal += 1
                ytotal += 1
                ztotal += 1
            # Single weight per point
            elseif ndims(weights) == 1
                xtotal += weights[i]
                ytotal += weights[i]
                ztotal += weights[i]
            # Dimension-wise weights
            else
                xtotal += weights[1, i]
                ytotal += weights[2, i]
                ztotal += weights[3, i]
            end
        end
    end

    xtotal, ytotal, ztotal
end


function leaf_spread(points, weights, span, dim)
    # points: (3, NP)
    # weights: nothing / (NP,) / (3, NP)
    # span: ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    # dim: 1 / 2 / 3
    T = eltype(points)

    # Select only points within this leaf's box
    select = (
        (span[1][1] .< @view(points[1, :])) .& (@view(points[1, :]) .< span[1][2]) .&
        (span[2][1] .< @view(points[2, :])) .& (@view(points[2, :]) .< span[2][2]) .&
        (span[3][1] .< @view(points[3, :])) .& (@view(points[3, :]) .< span[3][2])
    )

    select_coords = points[dim, select]

    # If we don't have weights, simply sort coordinates
    if isnothing(weights)
        sort!(select_coords)
        select_weights = nothing

    # Otherwise sort both coordinates and their corresponding weights
    else
        select_weights = ndims(weights) == 1 ? weights[select] : weights[dim, select]

        permutations = sortperm(select_coords)
        select_coords = select_coords[permutations]
        select_weights = select_weights[permutations]

        # Compute normalised cumulative sum of weights
        cumsum!(select_weights, select_weights)
        select_weights ./= select_weights[end]
    end

    # Spread is between 90% - 10% quantiles
    q05 = weighted_quantile(select_coords, select_weights, T(0.1))
    q95 = weighted_quantile(select_coords, select_weights, T(0.9))

    q95 - q05
end


function leaf_median(points, weights, span, dim)
    # points: (3, NP)
    # weights: nothing / (NP,) / (3, NP)
    # span: ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    # dim: 1 / 2 / 3
    T = eltype(points)

    # Select only points within this leaf's box
    select = (
        (span[1][1] .< @view(points[1, :])) .& (@view(points[1, :]) .< span[1][2]) .&
        (span[2][1] .< @view(points[2, :])) .& (@view(points[2, :]) .< span[2][2]) .&
        (span[3][1] .< @view(points[3, :])) .& (@view(points[3, :]) .< span[3][2])
    )

    select_coords = points[dim, select]

    # If we don't have weights, simply sort coordinates
    if isnothing(weights)
        sort!(select_coords)
        select_weights = nothing

    # Otherwise sort both coordinates and their corresponding weights
    else
        select_weights = ndims(weights) == 1 ? weights[select] : weights[dim, select]

        permutations = sortperm(select_coords)
        select_coords = select_coords[permutations]
        select_weights = select_weights[permutations]

        # Compute normalised cumulative sum of weights
        cumsum!(select_weights, select_weights)
        select_weights ./= select_weights[end]
    end

    weighted_quantile(select_coords, select_weights, T(0.5))
end


@inline function weighted_quantile(coords, weights, q)
    # Weights must have been accumulated and normalised!

    # No weights - interpolate between coordinates at indices adjacent to q
    if isnothing(weights)
        qmid = 1 + q * (length(coords) - 1)
        qindex = qmid |> ceil |> Int64
        if qindex == 1
            return coords[1]
        end

        return interpolate(qmid, qindex - 1, qindex, coords[qindex - 1], coords[qindex])
    end

    if q <= 0.5
        qcoord = coords[begin]
        for i in firstindex(weights):lastindex(weights)
            if weights[i] > q
                if i > firstindex(weights)
                    qcoord = interpolate(q, weights[i - 1], weights[i], coords[i - 1], coords[i])
                end
                break
            end
        end
    else
        qcoord = coords[end]
        for i in lastindex(weights):firstindex(weights)
            if weights[i] < q
                if i < lastindex(weights)
                    qcoord = interpolate(q, weights[i], weights[i + 1], coords[i], coords[i + 1])
                end
                break
            end
        end
    end

    qcoord
end


function interpolate(xmid, x0, x1, y0, y1)
    y0 + (y1 - y0) / (x1 - x0) * (xmid - x0)
end


function find_neighbors!(
    x, y, z,
    xspan::MVector{N, Tuple{T, T}},
    yspan::MVector{N, Tuple{T, T}},
    zspan::MVector{N, Tuple{T, T}},
    skin=zero(T),
) where {N, T}

    # Initialise vectors of neighbours
    for i in 1:N
        x[i] = Tuple{Int64, T}[]
        y[i] = Tuple{Int64, T}[]
        z[i] = Tuple{Int64, T}[]
    end

    # Find each leaf's neighbours
    for i in 1:N - 1
        for j in i + 1:N
            # Check if they share a boundary in one dimension and both other dimensions overlap
            if xspan[i][1] == xspan[j][2] && both_overlap(yspan, zspan, i, j, skin)
                push!(x[i], (j, xspan[i][1]))
                push!(x[j], (i, xspan[j][2]))
            end
            if xspan[i][2] == xspan[j][1] && both_overlap(yspan, zspan, i, j, skin)
                push!(x[i], (j, xspan[i][2]))
                push!(x[j], (i, xspan[j][1]))
            end

            if yspan[i][1] == yspan[j][2] && both_overlap(xspan, zspan, i, j, skin)
                push!(y[i], (j, yspan[i][1]))
                push!(y[j], (i, yspan[j][2]))
            end
            if yspan[i][2] == yspan[j][1] && both_overlap(xspan, zspan, i, j, skin)
                push!(y[i], (j, yspan[i][2]))
                push!(y[j], (i, yspan[j][1]))
            end

            if zspan[i][1] == zspan[j][2] && both_overlap(xspan, yspan, i, j, skin)
                push!(z[i], (j, zspan[i][1]))
                push!(z[j], (i, zspan[j][2]))
            end
            if zspan[i][2] == zspan[j][1] && both_overlap(xspan, yspan, i, j, skin)
                push!(z[i], (j, zspan[i][2]))
                push!(z[j], (i, zspan[j][1]))
            end
        end
    end

    nothing
end


function both_overlap(
    spans1::AbstractVector{Tuple{T, T}},
    spans2::AbstractVector{Tuple{T, T}},
    i, j,
    skin=zero(T),
) where T
    overlap(spans1[i], spans1[j], skin) && overlap(spans2[i], spans2[j], skin)
end


function overlap(a::Tuple{T, T}, b::Tuple{T, T}, skin=zero(T)) where T
    a[2] + skin > b[1] && a[1] < b[2] + skin
end


function evaluate_leaves!(
    totals::MVector{N, Tuple{T, T, T}},
    points, weights,
    xspan, yspan, zspan,
) where {N, T}

    # Initialise totals
    for i in 1:N
        totals[i] = (zero(T), zero(T), zero(T))
    end

    # Find each point's containing leaf and add weights to the corresponding total
    for i in 1:size(points, 2)
        point = @view points[:, i]
        for j in 1:N
            contained = (
                xspan[j][1] <= point[1] < xspan[j][2] &&
                yspan[j][1] <= point[2] < yspan[j][2] &&
                zspan[j][1] <= point[3] < zspan[j][2]
            )

            if contained
                if isnothing(weights)
                    totals[j] = (totals[j][1] + 1,
                                 totals[j][2] + 1,
                                 totals[j][3] + 1)
                elseif ndims(weights) == 1
                    totals[j] = (totals[j][1] + weights[i],
                                 totals[j][2] + weights[i],
                                 totals[j][3] + weights[i])
                elseif ndims(weights) == 2
                    totals[j] = (totals[j][1] + weights[1, i],
                                 totals[j][2] + weights[2, i],
                                 totals[j][3] + weights[3, i])
                end
            end
        end
    end

    nothing
end


function compute_ideal_totals(points, weights, N, T)
    if isnothing(weights)
        num_points = T(size(points, 2) / N)
        return (num_points, num_points, num_points)
    elseif ndims(weights) == 1
        total = T(sum(weights) / N)
        return (total, total, total)
    else
        xtotal, ytotal, ztotal = sum(weights, dims=2)
        return (T(xtotal / N), T(ytotal / N), T(ztotal / N))
    end
end


function worst_totals(
    totals::MVector{N, Tuple{T, T, T}},
    ideal_totals::Tuple{T, T, T},
    num_adjusts::MVector{N, Int64},
) where {N, T}

    # Find leaves with smallest number of previous adjustments and most different total weights
    # to the ideal value, one per dimension
    # TODO: only prioritise new leaves if there isn't a leaf with less than 0.5 ideal
    xleaf = 1
    yleaf = 1
    zleaf = 1

    xadjusts = typemax(N)
    yadjusts = typemax(N)
    zadjusts = typemax(N)

    xdiff = typemin(T)
    ydiff = typemin(T)
    zdiff = typemin(T)

    # Find new leaf with largest total difference in the X dimension
    for leaf_index in 1:N
        leaf_xtotal, leaf_ytotal, leaf_ztotal = totals[leaf_index]
        leaf_xdiff = abs(leaf_xtotal - ideal_totals[1])

        # If any leaf is in danger of being reduced to nothing, adjust that
        if leaf_xtotal < T(0.5) * ideal_totals[1]
            xleaf = leaf_index
            xadjusts = num_adjusts[leaf_index]
            xdiff = leaf_xdiff
            break
        end

        # Prioritise leaves with least adjustments, or equal adjustments but worse weight
        if num_adjusts[leaf_index] < xadjusts ||
                (num_adjusts[leaf_index] == xadjusts && leaf_xdiff > xdiff)
            xleaf = leaf_index
            xadjusts = num_adjusts[leaf_index]
            xdiff = leaf_xdiff
        end
    end

    # Same, but find different leaf to `xleaf`
    for leaf_index in 1:N
        (leaf_index == xleaf) && continue

        leaf_xtotal, leaf_ytotal, leaf_ztotal = totals[leaf_index]
        leaf_ydiff = abs(leaf_ytotal - ideal_totals[2])

        if leaf_ytotal < T(0.5) * ideal_totals[2]
            yleaf = leaf_index
            yadjusts = num_adjusts[leaf_index]
            ydiff = leaf_ydiff
            break
        end

        if num_adjusts[leaf_index] < yadjusts || 
                (num_adjusts[leaf_index] == yadjusts && leaf_ydiff > ydiff)
            yleaf = leaf_index
            yadjusts = num_adjusts[leaf_index]
            ydiff = leaf_ydiff
        end
    end

    # Same, but find different leaf to `xleaf` and `yleaf`
    for leaf_index in 1:N
        (leaf_index == xleaf || leaf_index == yleaf) && continue

        leaf_xtotal, leaf_ytotal, leaf_ztotal = totals[leaf_index]
        leaf_zdiff = abs(leaf_ztotal - ideal_totals[3])

        if leaf_ztotal < T(0.5) * ideal_totals[3]
            zleaf = leaf_index
            zadjusts = num_adjusts[leaf_index]
            zdiff = leaf_zdiff
            break
        end

        if num_adjusts[leaf_index] < zadjusts ||
                (num_adjusts[leaf_index] == zadjusts && leaf_zdiff > zdiff)
            zleaf = leaf_index
            zadjusts = num_adjusts[leaf_index]
            zdiff = leaf_zdiff
        end
    end

    xleaf, yleaf, zleaf
end


function adjust_leaf!(
    points, weights,
    dimspans::NTuple{3, MVector{N, Tuple{T, T}}},
    neighbors::Vector{Tuple{Int64, T}},
    totals::MVector{N, Tuple{T, T, T}},
    dim, leaf_index,
) where {N, T}

    # Adjust a domain `leaf_index` along dimension `dim` against neighbour with largest summed
    # weight difference

    # If there are no neighbours in this dimension, skip
    if length(neighbors) == 0
        return nothing
    end

    # Find index and boundary coordinate of neighbour with largest summed weight difference
    neigh_index, split_coord = largest_difference_neighbor(totals, neighbors, dim, leaf_index)

    # Find new boundary coordinate that splits weights evenly
    leaf_spans = (dimspans[1][leaf_index], dimspans[2][leaf_index], dimspans[3][leaf_index])
    neigh_spans = (dimspans[1][neigh_index], dimspans[2][neigh_index], dimspans[3][neigh_index])

    # Update span of all leaves with an upper / lower boundary at the old split coordinate
    new_split = leaves_median(leaf_spans, neigh_spans, points, weights, dim)
    span = dimspans[dim]

    if span[leaf_index][1] == split_coord                       # Split at lower boundary
        update_spans!(span, span[leaf_index][1], new_split)
    elseif span[leaf_index][2] == split_coord                   # Split at upper boundary
        update_spans!(span, span[leaf_index][2], new_split)
    else
        # TODO: extensive tests, then remove this exception
        throw(ArgumentError("Split was not the same!"))
    end

    nothing
end


function largest_difference_neighbor(
    totals::MVector{N, Tuple{T, T, T}},
    neighbors::Vector{Tuple{Int64, T}},
    dim, leaf_index,
) where {N, T}

    # Find domain (leaf_index) neighbour with largest weight difference
    @assert length(neighbors) > 0

    leaf_total = totals[leaf_index][dim]

    neigh_index = 1
    neigh_split = T(NaN)
    neigh_diff = typemin(T)

    for (current_index, split_coord) in neighbors
        current_diff = totals[current_index][dim] - leaf_total

        if current_diff > neigh_diff
            neigh_index = current_index
            neigh_split = split_coord
            neigh_diff = current_diff
        end
    end

    neigh_index, neigh_split
end


function leaves_median(span1, span2, points, weights, dim)
    # Find median from points found in either of domain span1 or span2
    # span: ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    # points: (3, NP)
    # weights: nothing / (NP,) / (3, NP)
    # dim: 1 / 2 / 3

    # Select only points within this leaf's box
    select1 = (
        (span1[1][1] .< @view(points[1, :])) .& (@view(points[1, :]) .< span1[1][2]) .&
        (span1[2][1] .< @view(points[2, :])) .& (@view(points[2, :]) .< span1[2][2]) .&
        (span1[3][1] .< @view(points[3, :])) .& (@view(points[3, :]) .< span1[3][2])
    )
    select2 = (
        (span2[1][1] .< @view(points[1, :])) .& (@view(points[1, :]) .< span2[1][2]) .&
        (span2[2][1] .< @view(points[2, :])) .& (@view(points[2, :]) .< span2[2][2]) .&
        (span2[3][1] .< @view(points[3, :])) .& (@view(points[3, :]) .< span2[3][2])
    )

    select = select1 .| select2
    select_coords = points[dim, select]

    # If we don't have weights, simply sort coordinates
    if isnothing(weights)
        sort!(select_coords)
        select_weights = nothing

    # Otherwise sort both coordinates and their corresponding weights
    else
        select_weights = ndims(weights) == 1 ? weights[select] : weights[dim, select]

        permutations = sortperm(select_coords)
        select_coords = select_coords[permutations]
        select_weights = select_weights[permutations]

        # Compute normalised cumulative sum of weights
        cumsum!(select_weights, select_weights)
        select_weights ./= select_weights[end]
    end

    weighted_quantile(select_coords, select_weights, 0.5)
end


function update_spans!(span::MVector{N, Tuple{T, T}}, old_split, new_split) where {N, T}
    # Update all spans that share one boundary at the old split
    for i in 1:N
        if span[i][1] == old_split
            span[i] = (new_split, span[i][2])
        elseif span[i][2] == old_split
            span[i] = (span[i][1], new_split)
        end
    end

    nothing
end










# function sorted_dimensions(points, weights::AbstractArray{N, T}) where {N, T}
# 
#     if N < 1 || N > 2 || (N == 2 && size(weights, 1) != 3)
#         throw(DimensionMismatch("weights must be a 1D vector or 2D matrix with 3 columns"))
#     end
# 
#     permutations = similar(points, Int64, size(points, 2))
# 
#     sortperm!(permutations, @view points[1, :])
#     xsort = points[1, permutations]
#     xweights = N == 1 ? weights[permutations] : weights[1, permutations]
# 
#     sortperm!(permutations, @view points[2, :])
#     ysort = points[2, permutations]
#     yweights = N == 1 ? weights[permutations] : weights[2, permutations]
# 
#     sortperm!(permutations, @view points[3, :])
#     zsort = points[3, permutations]
#     zweights = N == 1 ? weights[permutations] : weights[3, permutations]
# 
#     xsort, ysort, zsort, xweights, yweights, zweights
# end
# 
# 
# 
# 
# 
# function weighted_quantile(span, coords, weights, q)
#     select = (coords .> span[1]) .& (coords .< span[2])
# 
#     select_coords = coords[select]
# 
#     if isnothing(weights)
#         select_weights = nothing
#     else
#         select_weights = weights[select]
#         cumsum!(select_weights, select_weights)
#         select_weights ./= select_weights[end]
#     end
# 
#     weighted_quantile(select_coords, select_weights, q)
# end
# 
# 
# 
# 
# 
# 
# 
# 
# function expand_leaf!(span::MVector{N, Tuple{T, T}},
#                       neighbors::Vector{Vector{Tuple{Int64, T}}},
#                       totals::MVector{N, Tuple{T, T, T}},
#                       extents::MVector{N, Tuple{T, T, T}},
#                       dim, leaf_index) where {N, T}
# 
#     # Expand a domain `leaf_index` along dimension `dim`. The expansion amount Δ is calculated
#     # relative to the largest neighbour and the domain's span
#     leaf_neighbors = neighbors[leaf_index]
# 
#     # If there are no neighbours in this dimension, skip
#     if length(leaf_neighbors) == 0
#         return nothing
#     end
# 
#     # Find index of neighbour with largest summed weight
#     neigh_index = largest_neighbor(totals, dim, leaf_neighbors)
# 
#     # Compute relative difference in weight totals
#     leaf_total = totals[leaf_index][dim]
#     neigh_total = totals[neigh_index][dim]
# 
#     reldiff = T(1) - leaf_total / neigh_total
#     if reldiff < T(0.01)
#         return nothing
#     end
# 
#     # @show leaf_index neigh_index leaf_total neigh_total
#     # println()
# 
#     # Expansion amount as a fraction of current leaf span
#     Δ = T(0.2) * extents[leaf_index][dim] * reldiff
# 
#     # Decrease lower bounds and increase upper bounds
#     update_spans!(span, span[leaf_index][1], span[leaf_index][1] - Δ)
#     update_spans!(span, span[leaf_index][2], span[leaf_index][2] + Δ)
# 
#     nothing
# end
# 
# 
# function contract_leaf!(span::MVector{N, Tuple{T, T}},
#                         neighbors::Vector{Vector{Tuple{Int64, T}}},
#                         totals::MVector{N, Tuple{T, T, T}},
#                         extents::MVector{N, Tuple{T, T, T}},
#                         dim, leaf_index) where {N, T}
# 
#     # Contract a domain `leaf_index` along dimension `dim`. The contraction amount Δ is calculated
#     # relative to the smallest neighbour and the domain's span
#     leaf_neighbors = neighbors[leaf_index]
# 
#     # If there are no neighbours in this dimension, skip
#     if length(leaf_neighbors) == 0
#         return nothing
#     end
# 
#     # Find index of neighbour with smallest summed weight
#     neigh_index = smallest_neighbor(totals, dim, leaf_neighbors)
# 
#     # Compute relative difference in weight totals
#     leaf_total = totals[leaf_index][dim]
#     neigh_total = totals[neigh_index][dim]
# 
#     # If the relative difference is negative or too small, skip
#     reldiff = T(1) - neigh_total / leaf_total
#     if reldiff < T(0.01)
#         return nothing
#     end
# 
#     # @show leaf_index neigh_index leaf_total neigh_total
#     # println()
# 
#     # Contraction amount as a fraction of current leaf span
#     Δ = T(0.2) * extents[leaf_index][dim] * reldiff
# 
#     # Increase lower bounds and decrease upper bounds
#     update_spans!(span, span[leaf_index][1], span[leaf_index][1] + Δ)
#     update_spans!(span, span[leaf_index][2], span[leaf_index][2] - Δ)
# 
#     nothing
# end
# 
# 
# 
# 
# function smallest_total(totals::MVector{N, Tuple{T, T, T}}) where {N, T}
#     # Find leaf with smallest total *in any dimension*
#     leaf_index = 1
#     leaf_total = typemax(T)
# 
#     for ileaf in 1:N
#         xtotal, ytotal, ztotal = totals[ileaf]
#         if xtotal < leaf_total
#             leaf_index = ileaf
#             leaf_total = xtotal
#         end
#         if ytotal < leaf_total
#             leaf_index = ileaf
#             leaf_total = ytotal
#         end
#         if ztotal < leaf_total
#             leaf_index = ileaf
#             leaf_total = ztotal
#         end
#     end
# 
#     leaf_index
# end
# 
# 
# function largest_total(totals::MVector{N, Tuple{T, T, T}}) where {N, T}
#     # Find leaf with largest total *in any dimension*
#     leaf_index = 1
#     leaf_total = typemin(T)
# 
#     for ileaf in 1:N
#         xtotal, ytotal, ztotal = totals[ileaf]
#         if xtotal > leaf_total
#             leaf_index = ileaf
#             leaf_total = xtotal
#         end
#         if ytotal > leaf_total
#             leaf_index = ileaf
#             leaf_total = ytotal
#         end
#         if ztotal > leaf_total
#             leaf_index = ileaf
#             leaf_total = ztotal
#         end
#     end
# 
#     leaf_index
# end
# 
# 
# function smallest_neighbor(totals::MVector{N, Tuple{T, T, T}}, dim,
#                            leaf_neighbors::Vector{Tuple{Int64, T}}) where {N, T}
#     # Find neighbor with smallest total in a given dimension
#     neigh_index = 1
#     neigh_total = typemax(T)
# 
#     for (ineigh, _) in leaf_neighbors
#         if totals[ineigh][dim] < neigh_total
#             neigh_index = ineigh
#             neigh_total = totals[ineigh][dim]
#         end
#     end
# 
#     neigh_index
# end
# 
# 
# function largest_neighbor(totals::MVector{N, Tuple{T, T, T}}, dim,
#                           leaf_neighbors::Vector{Tuple{Int64, T}}) where {N, T}
#     # Find neighbor with maximum total in a given dimension 
#     neigh_index = 1
#     neigh_total = typemin(T)
# 
#     for (ineigh, _) in leaf_neighbors
#         if totals[ineigh][dim] > neigh_total
#             neigh_index = ineigh
#             neigh_total = totals[ineigh][dim]
#         end
#     end
# 
#     neigh_index
# end
# 
# 
# 
# 
# function increase!(span, split_coord, ileaf, Δ)
#     if split_coord == span[ileaf][1]
#         # Enlarge current `ileaf` lower boundary
#         update_spans!(span, span[ileaf][1], span[ileaf][1] - Δ)
#     elseif split_coord == span[ileaf][2]
#         # Enlarge current `ileaf` upper boundary
#         update_spans!(span, span[ileaf][2], span[ileaf][2] + Δ)
#     else
#         throw(AssertionError("split_coord was not the same! [1]"))
#     end
# 
#     nothing
# end
# 
# 
# function decrease!(span, split_coord, ileaf, Δ)
#     if split_coord == span[ileaf][1]
#         # Contract current `ileaf` lower boundary
#         update_spans!(span, span[ileaf][1], span[ileaf][1] + Δ)
#     elseif split_coord == span[ileaf][2]
#         # Contract current `ileaf` upper boundary
#         update_spans!(span, span[ileaf][2], span[ileaf][2] - Δ)
#     else
#         throw(AssertionError("split_coord was not the same! [2]"))
#     end
# 
#     nothing
# end





end     # module CoordinateBisection
