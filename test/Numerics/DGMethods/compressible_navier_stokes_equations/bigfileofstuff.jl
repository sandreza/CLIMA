using ClimateMachine
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Elements
import ClimateMachine.Mesh.Elements: baryweights
using ClimateMachine.Mesh.Grids: polynomialorders
using GaussQuadrature
using Base.Threads


# Depending on CliMa version 
# old, should return a tuple of polynomial orders
# polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = Tuple([N for i in 1:dim])
# new, should return a tuple of polynomial orders
# polynomialorders(::DiscontinuousSpectralElementGrid{T, dim, N}) where {T, dim, N} = N

# utils.jl
"""
function cellaverage(Q; M = nothing)
# Description
Compute the cell-average of Q given the mass matrix M.
Assumes that Q and M are the same size
# Arguments
`Q`: MPIStateArrays (array)
# Keyword Arguments
`M`: Mass Matrix (array)
# Return
The cell-average of Q
"""
function cellaverage(Q; M = nothing)
    if M == nothing
        return nothing
    end
    return (sum(M .* Q, dims = 1) ./ sum(M, dims = 1))[:]
end

"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
# Description
Gets the (x,y,z) coordinates corresponding to the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- `x, y, z`: views of x, y, z coordinates
"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end

"""
function cellcenters(Q; M = nothing)
# Description
Get the cell-centers of every element in the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- Tuple of cell-centers
"""
function cellcenters(grid::DiscontinuousSpectralElementGrid)
    x, y, z = coordinates(grid)
    M = view(grid.vgeo, :, grid.Mid, :)  # mass matrix
    xC = cellaverage(x, M = M)
    yC = cellaverage(y, M = M)
    zC = cellaverage(z, M = M)
    return xC[:], yC[:], zC[:]
end

# find_element.jl
# 3D version
function findelement(xC, yC, zC, location, p, lin)
    ex, ey, ez = size(lin)
    # i 
    currentmin = ones(1)
    minind = ones(Int64, 1)
    currentmin[1] = abs.(xC[p[lin[1, 1, 1]]] .- location[1])
    for i in 2:ex
        current = abs.(xC[p[lin[i, 1, 1]]] .- location[1])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    i = minind[1]
    # j 
    currentmin[1] = abs.(yC[p[lin[1, 1, 1]]] .- location[2])
    minind[1] = 1
    for i in 2:ey
        current = abs.(yC[p[lin[1, i, 1]]] .- location[2])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    j = minind[1]
    # k 
    currentmin[1] = abs.(zC[p[lin[1, 1, 1]]] .- location[3])
    minind[1] = 1
    for i in 2:ez
        current = abs.(zC[p[lin[1, 1, i]]] .- location[3])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    k = minind[1]
    return p[lin[i, j, k]]
end

# 2D version
function findelement(xC, yC, location, p, lin)
    ex, ey = size(lin)
    # i 
    currentmin = ones(1)
    minind = ones(Int64, 1)
    currentmin[1] = abs.(xC[p[lin[1, 1]]] .- location[1])
    for i in 2:ex
        current = abs.(xC[p[lin[i, 1]]] .- location[1])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    i = minind[1]
    # j 
    currentmin[1] = abs.(yC[p[lin[1, 1]]] .- location[2])
    minind[1] = 1
    for i in 2:ey
        current = abs.(yC[p[lin[1, i]]] .- location[2])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    j = minind[1]
    return p[lin[i, j]]
end

# gridhelper.jl

struct InterpolationHelper{S, T}
    points::S
    quadrature::S
    interpolation::S
    cartesianindex::T
end

function InterpolationHelper(g::DiscontinuousSpectralElementGrid)
    porders = polynomialorders(g)
    if length(porders) == 3
        npx, npy, npz = porders
        rx, wx = GaussQuadrature.legendre(npx + 1, both)
        ωx = baryweights(rx)
        ry, wy = GaussQuadrature.legendre(npy + 1, both)
        ωy = baryweights(ry)
        rz, wz = GaussQuadrature.legendre(npz + 1, both)
        ωz = baryweights(rz)
        linlocal = reshape(
            collect(1:((npx + 1) * (npy + 1) * (npz + 1))),
            (npx + 1, npy + 1, npz + 1),
        )
        return InterpolationHelper(
            (rx, ry, rz),
            (wx, wy, wz),
            (ωx, ωy, ωz),
            linlocal,
        )
    elseif length(porders) == 2
        npx, npy = porders
        rx, wx = GaussQuadrature.legendre(npx + 1, both)
        ωx = baryweights(rx)
        ry, wy = GaussQuadrature.legendre(npy + 1, both)
        ωy = baryweights(ry)
        linlocal =
            reshape(collect(1:((npx + 1) * (npy + 1))), (npx + 1, npy + 1))
        return InterpolationHelper((rx, ry), (wx, wy), (ωx, ωy), linlocal)
    else
        println("Not supported")
        return nothing
    end
    return nothing
end

struct ElementHelper{S, T, U, Q, V, W}
    cellcenters::S
    coordinates::T
    cartesiansizes::U
    polynomialorders::Q
    permutation::V
    cartesianindex::W
end

addup(xC, tol) = sum(abs.(xC[1] .- xC) .≤ tol)

# only valid for cartesian domains
function ElementHelper(g::DiscontinuousSpectralElementGrid)
    porders = polynomialorders(g)
    x, y, z = coordinates(g)
    xC, yC, zC = cellcenters(g)
    ne = size(x)[2]
    ex = round(Int64, ne / addup(xC, 10^4 * eps(maximum(abs.(x)))))
    ey = round(Int64, ne / addup(yC, 10^4 * eps(maximum(abs.(y)))))
    ez = round(Int64, ne / addup(zC, 10^4 * eps(maximum(abs.(z)))))

    check = ne == ex * ey * ez
    check ? true : error("improper counting")
    p = getperm(xC, yC, zC, ex, ey, ez)
    # should use dispatch ...
    if length(porders) == 3
        npx, npy, npz = porders
        lin = reshape(collect(1:length(xC)), (ex, ey, ez))
        return ElementHelper(
            (xC, yC, zC),
            (x, y, z),
            (ex, ey, ez),
            porders,
            p,
            lin,
        )
    elseif length(porders) == 2
        npx, npy = porders
        lin = reshape(collect(1:length(xC)), (ex, ey))
        check = ne == ex * ey
        check ? true : error("improper counting")
        return ElementHelper((xC, yC), (x, y), (ex, ey), porders, p, lin)
    else
        println("no constructor for polynomial order = ", porders)
        return nothing
    end
    return nothing
end

struct GridHelper{S, T, V}
    interpolation::S
    element::T
    grid::V
end

function GridHelper(g::DiscontinuousSpectralElementGrid)
    return GridHelper(InterpolationHelper(g), ElementHelper(g), g)
end

function getvalue(f, location, gridhelper::GridHelper)
    ih = gridhelper.interpolation
    eh = gridhelper.element
    porders = gridhelper.element.polynomialorders
    if length(porders) == 3
        npx, npy, npz = gridhelper.element.polynomialorders
        fl = reshape(f, (npx + 1, npy + 1, npz + 1, prod(eh.cartesiansizes)))
        ip = getvalue(
            fl,
            eh.cellcenters...,
            location,
            eh.permutation,
            eh.cartesianindex,
            ih.cartesianindex,
            eh.coordinates...,
            ih.points...,
            ih.interpolation...,
        )
        return ip
    elseif length(porders) == 2
        npx, npy = gridhelper.element.polynomialorders
        fl = reshape(f, (npx + 1, npy + 1, prod(eh.cartesiansizes)))
        ip = getvalue(
            fl,
            eh.cellcenters...,
            location,
            eh.permutation,
            eh.cartesianindex,
            ih.cartesianindex,
            eh.coordinates...,
            ih.points...,
            ih.interpolation...,
        )
        return ip
    end
    return nothing
end

# lagrange_interpolation.jl
function checkgl(x, rx)
    for i in eachindex(rx)
        if abs(x - rx[i]) ≤ eps(rx[i])
            return i
        end
    end
    return 0
end

function lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
    icheck = checkgl(newx, rx)
    jcheck = checkgl(newy, ry)
    kcheck = checkgl(newz, rz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
        if kcheck == 0
            Δz = (newz .- rz[k])
            polez = ωz[k] ./ Δz
            kk = k
        else
            polez = 1.0
            k = eachindex(rz)[end]
            kk = kcheck
        end
        for j in eachindex(ry)
            if jcheck == 0
                Δy = (newy .- ry[j])
                poley = ωy[j] ./ Δy
                jj = j
            else
                poley = 1.0
                j = eachindex(ry)[end]
                jj = jcheck
            end
            for i in eachindex(rx)
                if icheck == 0
                    Δx = (newx .- rx[i])
                    polex = ωx[i] ./ Δx
                    ii = i
                else
                    polex = 1.0
                    i = eachindex(rx)[end]
                    ii = icheck
                end
                numerator[1] += f[ii, jj, kk] * polex * poley * polez
                denominator[1] += polex * poley * polez
            end
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval(f, newx, newy, rx, ry, ωx, ωy)
    icheck = checkgl(newx, rx)
    jcheck = checkgl(newy, ry)
    numerator = zeros(1)
    denominator = zeros(1)
    for j in eachindex(ry)
        if jcheck == 0
            Δy = (newy .- ry[j])
            poley = ωy[j] ./ Δy
            jj = j
        else
            poley = 1.0
            j = eachindex(ry)[end]
            jj = jcheck
        end
        for i in eachindex(rx)
            if icheck == 0
                Δx = (newx .- rx[i])
                polex = ωx[i] ./ Δx
                ii = i
            else
                polex = 1.0
                i = eachindex(rx)[end]
                ii = icheck
            end
            numerator[1] += f[ii, jj] * polex * poley
            denominator[1] += polex * poley
        end
    end
    return numerator[1] / denominator[1]
end


function lagrange_eval_nocheck(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
        Δz = (newz .- rz[k])
        polez = ωz[k] ./ Δz
        for j in eachindex(ry)
            Δy = (newy .- ry[j])
            poley = ωy[j] ./ Δy
            for i in eachindex(rx)
                Δx = (newx .- rx[i])
                polex = ωx[i] ./ Δx
                numerator[1] += f[i, j, k] * polex * poley * polez
                denominator[1] += polex * poley * polez
            end
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval_nocheck(f, newx, newy, rx, ry, ωx, ωy)
    numerator = zeros(1)
    denominator = zeros(1)
    for j in eachindex(ry)
        Δy = (newy .- ry[j])
        poley = ωy[j] ./ Δy
        for i in eachindex(rx)
            Δx = (newx .- rx[i])
            polex = ωx[i] ./ Δx
            numerator[1] += f[i, j] * polex * poley
            denominator[1] += polex * poley
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval_nocheck(f, newx, rx, ωx)
    numerator = zeros(1)
    denominator = zeros(1)
    for i in eachindex(rx)
        Δx = (newx .- rx[i])
        polex = ωx[i] ./ Δx
        numerator[1] += f[i] * polex * poley
        denominator[1] += polex * poley
    end
    return numerator[1] / denominator[1]
end

# 3D, only valid for rectangles
function getvalue(
    fl,
    xC,
    yC,
    zC,
    location,
    p,
    lin,
    linlocal,
    x,
    y,
    z,
    rx,
    ry,
    rz,
    ωx,
    ωy,
    ωz,
)
    e = findelement(xC, yC, zC, location, p, lin)
    # need bounds to rescale, only value for cartesian

    xmax = x[linlocal[length(rx), 1, 1], e]
    xmin = x[linlocal[1, 1, 1], e]
    ymax = y[linlocal[1, length(ry), 1], e]
    ymin = y[linlocal[1, 1, 1], e]
    zmax = z[linlocal[1, 1, length(rz)], e]
    zmin = z[linlocal[1, 1, 1], e]

    # rescale new point to [-1,1]³
    newx = 2 * (location[1] - xmin) / (xmax - xmin) - 1
    newy = 2 * (location[2] - ymin) / (ymax - ymin) - 1
    newz = 2 * (location[3] - zmin) / (zmax - zmin) - 1

    return lagrange_eval(
        view(fl, :, :, :, e),
        newx,
        newy,
        newz,
        rx,
        ry,
        rz,
        ωx,
        ωy,
        ωz,
    )
end

# 2D
function getvalue(fl, xC, yC, location, p, lin, linlocal, x, y, rx, ry, ωx, ωy)
    e = findelement(xC, yC, location, p, lin)
    # need bounds to rescale
    xmax = x[linlocal[length(rx), 1, 1], e]
    xmin = x[linlocal[1, 1, 1], e]
    ymax = y[linlocal[1, length(ry), 1], e]
    ymin = y[linlocal[1, 1, 1], e]

    # rescale new point to [-1,1]²
    newx = 2 * (location[1] - xmin) / (xmax - xmin) - 1
    newy = 2 * (location[2] - ymin) / (ymax - ymin) - 1

    return lagrange_eval(view(fl, :, :, e), newx, newy, rx, ry, ωx, ωy)
end

# permutations.jl
function getperm(xC, yC, zC, ex, ey, ez)
    pz = sortperm(zC)
    tmpY = reshape(yC[pz], (ex * ey, ez))
    tmp_py = [sortperm(tmpY[:, i]) for i in 1:ez]
    py = zeros(Int64, length(pz))
    for i in eachindex(tmp_py)
        n = length(tmp_py[i])
        ii = (i - 1) * n + 1
        py[ii:(ii + n - 1)] .= tmp_py[i] .+ ii .- 1
    end
    tmpX = reshape(xC[pz][py], (ex, ey * ez))
    tmp_px = [sortperm(tmpX[:, i]) for i in 1:(ey * ez)]
    px = zeros(Int64, length(pz))
    for i in eachindex(tmp_px)
        n = length(tmp_px[i])
        ii = (i - 1) * n + 1
        px[ii:(ii + n - 1)] .= tmp_px[i] .+ ii .- 1
    end
    p = [pz[py[px[i]]] for i in eachindex(px)]
    return p
end

# ScalarFields

using Base.Threads, LinearAlgebra
import Base: getindex, materialize!, broadcasted

abstract type AbstractField end
struct ScalarField{S, T} <: AbstractField
    data::S
    grid::T
end

function (ϕ::ScalarField)(x::Tuple)
    return getvalue(ϕ.data, x, ϕ.grid)
end

function (ϕ::ScalarField)(x::Number, y::Number, z::Number)
    return getvalue(ϕ.data, (x, y, z), ϕ.grid)
end

function (ϕ::ScalarField)(x::Number, y::Number)
    return getvalue(ϕ.data, (x, y), ϕ.grid)
end

getindex(ϕ::ScalarField, i::Int) = ϕ.data[i]

materialize!(ϕ::ScalarField, f::Base.Broadcast.Broadcasted) =
    materialize!(ϕ.data, f)
broadcasted(identity, ϕ::ScalarField) = broadcasted(Base.identity, ϕ.data)

function (ϕ::ScalarField)(
    xlist::StepRangeLen,
    ylist::StepRangeLen,
    zlist::StepRangeLen;
    threads = false,
)
    newfield = zeros(length(xlist), length(ylist), length(zlist))
    if threads
        @threads for k in eachindex(zlist)
            for j in eachindex(ylist)
                for i in eachindex(xlist)
                    newfield[i, j, k] =
                        getvalue(ϕ.data, (xlist[i], ylist[j], zlist[k]), ϕ.grid)
                end
            end
        end
    else
        for k in eachindex(zlist)
            for j in eachindex(ylist)
                for i in eachindex(xlist)
                    newfield[i, j, k] =
                        getvalue(ϕ.data, (xlist[i], ylist[j], zlist[k]), ϕ.grid)
                end
            end
        end
    end
    return newfield
end

function (ϕ::ScalarField)(
    xlist::StepRangeLen,
    ylist::StepRangeLen;
    threads = false,
)
    newfield = zeros(length(xlist), length(ylist))
    if threads
        @threads for j in eachindex(ylist)
            for i in eachindex(xlist)
                newfield[i, j] = getvalue(ϕ.data, (xlist[i], ylist[j]), ϕ.grid)
            end
        end
    else
        for j in eachindex(ylist)
            for i in eachindex(xlist)
                newfield[i, j] = getvalue(ϕ.data, (xlist[i], ylist[j]), ϕ.grid)
            end
        end
    end
    return newfield
end

function uniform_grid(Ω::AbstractDomain; resolution = (32, 32, 32))
    dims = ndims(Ω)
    resolution = resolution[1:dims]
    uniform = []
    for i in 1:dims
        push!(uniform, range(Ω[i].min, Ω[i].max, length = resolution[i]))
    end
    return Tuple(uniform)
end

# Vizinanigans

using GLMakie, Statistics, Printf

"""
visualize(states::AbstractArray; statenames = string.(1:length(states)), quantiles = (0.1, 0.99), aspect = (1,1,1), resolution = (1920, 1080), statistics = false, title = "Field = ")
# Description 
Visualize 3D states 
# Arguments
- `states`: Array{Array{Float64,3},1}. An array of arrays containing different fields
# Keyword Arguments
- `statenames`: Array{String,1}. An array of stringnames
- `aspect`: Tuple{Int64,Int64,Float64}. Determines aspect ratio of box for volumes
- `resolution`: Resolution of preliminary makie window
- `statistics`: boolean. toggle for displaying statistics 
# Return
- `scene`: Scene. A preliminary scene object for manipulation
"""
function visualize(
    states::AbstractArray;
    statenames = string.(1:length(states)),
    units = ["" for i in eachindex(states)],
    aspect = (1, 1, 1),
    resolution = (1920, 1080),
    statistics = false,
    title = "Field = ",
    bins = 300,
)
    # Create scene
    scene, layout = layoutscene(resolution = resolution)
    lscene = layout[2:4, 2:4] = LScene(scene)
    width = round(Int, resolution[1] / 4) # make menu 1/4 of preliminary resolution

    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    if statistics
        llscene =
            layout[4, 1] = Axis(
                scene,
                xlabel = @lift(statenames[$statenode] * units[$statenode]),
                xlabelcolor = :black,
                ylabel = "pdf",
                ylabelcolor = :black,
                xlabelsize = 40,
                ylabelsize = 40,
                xticklabelsize = 25,
                yticklabelsize = 25,
                xtickcolor = :black,
                ytickcolor = :black,
                xticklabelcolor = :black,
                yticklabelcolor = :black,
            )
        layout[3, 1] = Label(scene, "Statistics", width = width, textsize = 50)
    end

    # x,y,z are for determining the aspect ratio of the box
    if (typeof(aspect) <: Tuple) & (length(aspect) == 3)
        x, y, z = aspect
    else
        x, y, z = size(states[1])
    end

    # Clim sliders
    upperclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node = upperclim_slider.value
    lowerclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node = lowerclim_slider.value

    # Lift Nodes
    state = @lift(states[$statenode])
    statename = @lift(statenames[$statenode])
    clims = @lift((
        quantile($state[:], $lowerclim_node),
        quantile($state[:], $upperclim_node),
    ))
    cmap_rgb = @lift(to_colormap($colornode))
    titlename = @lift(title * $statename) # use padding and appropriate centering

    # Statistics
    if statistics
        histogram_node = @lift(histogram($state, bins = bins))
        xs = @lift($histogram_node[1])
        ys = @lift($histogram_node[2])
        pdf = GLMakie.AbstractPlotting.barplot!(
            llscene,
            xs,
            ys,
            color = :red,
            strokecolor = :red,
            strokewidth = 1,
        )
        @lift(GLMakie.AbstractPlotting.xlims!(llscene, extrema($state)))
        @lift(GLMakie.AbstractPlotting.ylims!(
            llscene,
            extrema($histogram_node[2]),
        ))
        vlines!(
            llscene,
            @lift($clims[1]),
            color = :black,
            linewidth = width / 100,
        )
        vlines!(
            llscene,
            @lift($clims[2]),
            color = :black,
            linewidth = width / 100,
        )
    end

    # Volume Plot 
    volume!(
        lscene,
        0..x,
        0..y,
        0..z,
        state,
        camera = cam3d!,
        colormap = cmap_rgb,
        colorrange = clims,
    )
    # Camera
    cam = cameracontrols(scene.children[1])
    eyeposition = Float32[2, 2, 1.3]
    lookat = Float32[0.82, 0.82, 0.1]
    # Title
    supertitle =
        layout[1, 2:4] = Label(scene, titlename, textsize = 50, color = :black)


    # Menus
    statemenu = Menu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = Menu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end
    lowerclim_string = @lift(
        "lower clim quantile = " *
        @sprintf("%0.2f", $lowerclim_node) *
        ", value = " *
        @sprintf("%0.1e", $clims[1])
    )
    upperclim_string = @lift(
        "upper clim quantile = " *
        @sprintf("%0.2f", $upperclim_node) *
        ", value = " *
        @sprintf("%0.1e", $clims[2])
    )
    # depends on makie version, vbox for old, vgrid for new
    layout[2, 1] = vgrid!(
        Label(scene, "State", width = nothing),
        statemenu,
        Label(scene, "Color", width = nothing),
        colormenu,
        Label(scene, lowerclim_string, width = nothing),
        lowerclim_slider,
        Label(scene, upperclim_string, width = nothing),
        upperclim_slider,
    )
    layout[1, 1] = Label(scene, "Menu", width = width, textsize = 50)

    # Modify Axis

    # axis = scene.children[1][OldAxis]
    # axis[:names][:axisnames] = ("↓ Zonal [m] ", "Meriodonal [m]↓ ", "Depth [m]↓ ")
    # axis[:names][:axisnames] = ("↓", "↓ ", "↓ ")
    # axis[:names][:align] = ((:left, :center), (:right, :center), (:right, :center))
    # need to adjust size of ticks first and then size of axis names
    # axis[:names][:textsize] = (50.0, 50.0, 50.0)
    # axis[:ticks][:textsize] = (00.0, 00.0, 00.0)
    # axis[:ticks][:ranges_labels].val # current axis labels
    #=
    xticks = collect(range(-0, aspect[1], length = 2))
    yticks = collect(range(-0, aspect[2], length = 6))
    zticks = collect(range(-0, aspect[3], length = 2))
    ticks = (xticks, yticks, zticks)
    axis[:ticks][:ranges] = ticks
    xtickslabels = [@sprintf("%0.1f", (xtick)) for xtick in xticks]
    xtickslabels[end] = "1e6"
    ytickslabels = ["", "south", "", "", "north", ""]
    ztickslabels = [@sprintf("%0.1f", (xtick)) for xtick in xticks]
    labels = (xtickslabels, ytickslabels, ztickslabels)
    axis[:ticks][:labels] = labels
    =#

    display(scene)
    # Change the default camera position after the fact
    # note that these change dynamically as the plot is manipulated
    return scene
end


"""
histogram(array; bins = 100)
# Description
return arrays for plotting histogram
"""
function histogram(
    array;
    bins = minimum([100, length(array)]),
    normalize = true,
)
    tmp = zeros(bins)
    down, up = extrema(array)
    down, up = down == up ? (down - 1, up + 1) : (down, up) # edge case
    bucket = collect(range(down, up, length = bins + 1))
    normalization = normalize ? length(array) : 1
    for i in eachindex(array)
        # normalize then multiply by bins
        val = (array[i] - down) / (up - down) * bins
        ind = ceil(Int, val)
        # handle edge cases
        ind = maximum([ind, 1])
        ind = minimum([ind, bins])
        tmp[ind] += 1 / normalization
    end
    return (bucket[2:end] + bucket[1:(end - 1)]) .* 0.5, tmp
end

# 2D visualization
function visualize(
    states::Array{Array{S, 2}, 1};
    statenames = string.(1:length(states)),
    units = ["" for i in eachindex(states)],
    aspect = (1, 1, 1),
    resolution = (2412, 1158),
    title = "2D  ",
    xlims = (0, 1),
    ylims = (0, 1),
    bins = 300,
) where {S}
    # Create scene
    scene, layout = layoutscene(resolution = resolution)
    lscene =
        layout[2:4, 2:4] = Axis(
            scene,
            xlabel = "x ",
            xlabelcolor = :black,
            ylabel = "y ",
            ylabelcolor = :black,
            xlabelsize = 40,
            ylabelsize = 40,
            xticklabelsize = 25,
            yticklabelsize = 25,
            xtickcolor = :black,
            ytickcolor = :black,
            xticklabelcolor = :black,
            yticklabelcolor = :black,
            titlesize = 50,
        )
    width = round(Int, resolution[1] / 4) # make menu 1/4 of preliminary resolution

    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    interpolationlabels = ["contour", "heatmap"]
    interpolationchoices = [true, false]
    interpolationnode = Node(interpolationchoices[1])

    # Statistics
    llscene =
        layout[4, 1] = Axis(
            scene,
            xlabel = @lift(statenames[$statenode] * " " * units[$statenode]),
            xlabelcolor = :black,
            ylabel = "pdf",
            ylabelcolor = :black,
            xlabelsize = 40,
            ylabelsize = 40,
            xticklabelsize = 25,
            yticklabelsize = 25,
            xtickcolor = :black,
            ytickcolor = :black,
            xticklabelcolor = :black,
            yticklabelcolor = :black,
        )
    layout[3, 1] = Label(scene, "Statistics", width = width, textsize = 50)

    # Clim sliders
    upperclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node = upperclim_slider.value
    lowerclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node = lowerclim_slider.value

    #ylims = @lift(range($lowerval, $upperval, length = $))
    # Lift Nodes
    state = @lift(states[$statenode])
    statename = @lift(statenames[$statenode])
    unit = @lift(units[$statenode])
    oclims = @lift((
        quantile($state[:], $lowerclim_node),
        quantile($state[:], $upperclim_node),
    ))
    cmap_rgb = colornode
    clims = @lift(
        $oclims[1] != $oclims[2] ? (minimum($oclims), maximum($oclims)) :
        (minimum($oclims) - 1, maximum($oclims) + 1)
    )
    xlims = Array(range(xlims[1], xlims[2], length = 4)) #collect(range(xlims[1], xlims[2], length = size(states[1])[1]))
    ylims = Array(range(ylims[1], ylims[2], length = 4)) #@lift(collect(range($lowerval], $upperval, length = size($state)[2])))
    # newrange = @lift(range($lowerval, $upperval, length = 4))
    # lscene.yticks = @lift(Array($newrange))
    titlename = @lift(title * $statename) # use padding and appropriate centering
    layout[1, 2:4] = Label(scene, titlename, textsize = 50)
    # heatmap 
    heatmap1 = heatmap!(
        lscene,
        xlims,
        ylims,
        state,
        interpolate = interpolationnode,
        colormap = cmap_rgb,
        colorrange = clims,
    )


    # statistics
    histogram_node = @lift(histogram($state, bins = bins))
    xs = @lift($histogram_node[1])
    ys = @lift($histogram_node[2])
    pdf = GLMakie.AbstractPlotting.barplot!(
        llscene,
        xs,
        ys,
        color = :red,
        strokecolor = :red,
        strokewidth = 1,
    )
    @lift(GLMakie.AbstractPlotting.xlims!(llscene, extrema($state)))
    @lift(GLMakie.AbstractPlotting.ylims!(llscene, extrema($histogram_node[2])))
    vlines!(llscene, @lift($clims[1]), color = :black, linewidth = width / 100)
    vlines!(llscene, @lift($clims[2]), color = :black, linewidth = width / 100)

    # Menus
    statemenu = Menu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = Menu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end

    interpolationmenu =
        Menu(scene, options = zip(interpolationlabels, interpolationchoices))
    on(interpolationmenu.selection) do s
        interpolationnode[] = s
        heatmap1 = heatmap!(
            lscene,
            xlims,
            ylims,
            state,
            interpolate = s,
            colormap = cmap_rgb,
            colorrange = clims,
        )
    end

    newlabel = @lift($statename * " " * $unit)
    cbar = Colorbar(scene, heatmap1, label = newlabel)
    cbar.width = Relative(1 / 3)
    cbar.height = Relative(5 / 6)
    cbar.halign = :center
    # cbar.flipaxisposition = true
    # cbar.labelpadding = -350
    cbar.labelsize = 50

    lowerclim_string = @lift(
        "clim quantile = " *
        @sprintf("%0.2f", $lowerclim_node) *
        ", value = " *
        @sprintf("%0.1e", $clims[1])
    )
    upperclim_string = @lift(
        "clim quantile = " *
        @sprintf("%0.2f", $upperclim_node) *
        ", value = " *
        @sprintf("%0.1e", $clims[2])
    )

    # depends on makie version, vbox for old, vgrid for new
    layout[2, 1] = vgrid!(
        Label(scene, "State", width = nothing),
        statemenu,
        Label(
            scene,
            "plotting options",
            width = width,
            textsize = 30,
            padding = (0, 0, 10, 0),
        ),
        interpolationmenu,
        Label(scene, "Color", width = nothing),
        colormenu,
        Label(scene, lowerclim_string, width = nothing),
        lowerclim_slider,
        Label(scene, upperclim_string, width = nothing),
        upperclim_slider,
    )

    layout[2:4, 5] = vgrid!(
        Label(
            scene,
            "Color Bar",
            width = width / 2,
            textsize = 50,
            padding = (25, 0, 0, 00),
        ),
        cbar,
    )
    layout[1, 1] = Label(scene, "Menu", width = width, textsize = 50)
    display(scene)
    return scene
end

function volumeslice(
    states::AbstractArray;
    statenames = string.(1:length(states)),
    units = ["" for i in eachindex(states)],
    aspect = (1, 1, 32 / 192),
    resolution = (2678, 1030),
    statistics = false,
    title = "Volume plot of ",
    bins = 300,
    statlabelsize = (20, 20),
)
    scene, layout = layoutscene(resolution = resolution)
    volumescene = layout[2:4, 2:4] = LScene(scene)
    menuwidth = round(Int, 350)
    layout[1, 1] = Label(scene, "Menu", width = menuwidth, textsize = 50)

    slice_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.0)
    slice_node = slice_slider.value

    directionindex = [1, 2, 3]
    directionnames = ["x-slice", "y-slice", "z-slice"]
    directionnode = Node(directionindex[1])

    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    layout[1, 2:4] =
        Label(scene, @lift(title * statenames[$statenode]), textsize = 50)

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    state = @lift(states[$statenode])
    statename = @lift(statenames[$statenode])
    unit = @lift(units[$statenode])
    nx = @lift(size($state)[1])
    ny = @lift(size($state)[2])
    nz = @lift(size($state)[3])
    nr = @lift([$nx, $ny, $nz])

    nslider = 100
    xrange = range(0.00, aspect[1], length = nslider)
    yrange = range(0.00, aspect[2], length = nslider)
    zrange = range(0.00, aspect[3], length = nslider)
    constx = collect(reshape(xrange, (nslider, 1, 1)))
    consty = collect(reshape(yrange, (1, nslider, 1)))
    constz = collect(reshape(zrange, (1, 1, nslider)))
    matx = zeros(nslider, nslider, nslider)
    maty = zeros(nslider, nslider, nslider)
    matz = zeros(nslider, nslider, nslider)
    matx .= constx
    maty .= consty
    matz .= constz
    sliceconst = [matx, maty, matz]
    planeslice = @lift(sliceconst[$directionnode])

    upperclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node = upperclim_slider.value
    lowerclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node = lowerclim_slider.value

    clims = @lift((
        quantile($state[:], $lowerclim_node),
        quantile($state[:], $upperclim_node),
    ))

    volume!(
        volumescene,
        0..aspect[1],
        0..aspect[2],
        0..aspect[3],
        state,
        overdraw = false,
        colorrange = clims,
        colormap = @lift(to_colormap($colornode)),
    )


    alpha_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.5)
    alphanode = alpha_slider.value

    slicecolormap = @lift(cgrad(:viridis, alpha = $alphanode))
    v = volume!(
        volumescene,
        0..aspect[1],
        0..aspect[2],
        0..aspect[3],
        planeslice,
        algorithm = :iso,
        isorange = 0.005,
        isovalue = @lift($slice_node * aspect[$directionnode]),
        transparency = true,
        overdraw = false,
        visible = true,
        colormap = slicecolormap,
        colorrange = [-1, 0],
    )

    # Volume histogram

    layout[3, 1] = Label(scene, "Statistics", textsize = 50)
    hscene =
        layout[4, 1] = Axis(
            scene,
            xlabel = @lift(statenames[$statenode] * " " * units[$statenode]),
            xlabelcolor = :black,
            ylabel = "pdf",
            ylabelcolor = :black,
            xlabelsize = 40,
            ylabelsize = 40,
            xticklabelsize = statlabelsize[1],
            yticklabelsize = statlabelsize[2],
            xtickcolor = :black,
            ytickcolor = :black,
            xticklabelcolor = :black,
            yticklabelcolor = :black,
        )

    histogram_node = @lift(histogram($state, bins = bins))
    vxs = @lift($histogram_node[1])
    vys = @lift($histogram_node[2])
    pdf = GLMakie.AbstractPlotting.barplot!(
        hscene,
        vxs,
        vys,
        color = :red,
        strokecolor = :red,
        strokewidth = 1,
    )

    @lift(GLMakie.AbstractPlotting.xlims!(hscene, extrema($vxs)))
    @lift(GLMakie.AbstractPlotting.ylims!(hscene, extrema($vys)))
    vlines!(
        hscene,
        @lift($clims[1]),
        color = :black,
        linewidth = menuwidth / 100,
    )
    vlines!(
        hscene,
        @lift($clims[2]),
        color = :black,
        linewidth = menuwidth / 100,
    )


    # Slice
    sliceupperclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    sliceupperclim_node = sliceupperclim_slider.value
    slicelowerclim_slider =
        Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    slicelowerclim_node = slicelowerclim_slider.value


    slicexaxislabel = @lift(["y", "x", "x"][$directionnode])
    sliceyaxislabel = @lift(["z", "z", "y"][$directionnode])

    slicexaxis = @lift([[1, $ny], [1, $nx], [1, $nx]][$directionnode])
    sliceyaxis = @lift([[1, $nz], [1, $nz], [1, $ny]][$directionnode])

    slicescene =
        layout[2:4, 5:6] =
            Axis(scene, xlabel = slicexaxislabel, ylabel = sliceyaxislabel)

    sliced_state1 = @lift( $state[
        round(Int, 1 + $slice_node * (size($state)[1] - 1)),
        1:size($state)[2],
        1:size($state)[3],
    ])
    sliced_state2 = @lift( $state[
        1:size($state)[1],
        round(Int, 1 + $slice_node * (size($state)[2] - 1)),
        1:size($state)[3],
    ])
    sliced_state3 = @lift( $state[
        1:size($state)[1],
        1:size($state)[2],
        round(Int, 1 + $slice_node * (size($state)[3] - 1)),
    ])
    sliced_states = @lift([$sliced_state1, $sliced_state2, $sliced_state3])
    sliced_state = @lift($sliced_states[$directionnode])

    oclims = @lift((
        quantile($sliced_state[:], $slicelowerclim_node),
        quantile($sliced_state[:], $sliceupperclim_node),
    ))
    slicecolormapnode = @lift($oclims[1] < $oclims[2] ? $colornode : $colornode)
    sliceclims = @lift(
        $oclims[1] != $oclims[2] ? (minimum($oclims), maximum($oclims)) :
        (minimum($oclims) - 1, maximum($oclims) + 1)
    )

    heatmap1 = heatmap!(
        slicescene,
        slicexaxis,
        sliceyaxis,
        sliced_state,
        interpolate = true,
        colormap = slicecolormapnode,
        colorrange = sliceclims,
    )

    # Colorbar
    newlabel = @lift($statename * " " * $unit)
    cbar = Colorbar(scene, heatmap1, label = newlabel)
    cbar.width = Relative(1 / 3)
    # cbar.height = Relative(5/6)
    cbar.halign = :left
    # cbar.flipaxisposition = true
    # cbar.labelpadding = -250
    cbar.labelsize = 50

    @lift(GLMakie.AbstractPlotting.xlims!(slicescene, extrema($slicexaxis)))
    @lift(GLMakie.AbstractPlotting.ylims!(slicescene, extrema($sliceyaxis)))

    sliceindex = @lift([
        round(Int, 1 + $slice_node * ($nx - 1)),
        round(Int, 1 + $slice_node * ($ny - 1)),
        round(Int, 1 + $slice_node * ($nz - 1)),
    ][$directionnode])
    slicestring =
        @lift(directionnames[$directionnode] * " of " * statenames[$statenode])
    layout[1, 5:6] = Label(scene, slicestring, textsize = 50)


    axis = scene.children[1][OldAxis]
    # axis[:names][:axisnames] = ("↓", "↓ ", "↓ ")
    axis[:names][:align] =
        ((:left, :center), (:right, :center), (:right, :center))
    axis[:names][:textsize] = (50.0, 50.0, 50.0)
    axis[:ticks][:textsize] = (00.0, 00.0, 00.0)


    # Menus
    statemenu = Menu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = Menu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end


    # Slice Statistics
    layout[1, 7] = Label(scene, "Slice Menu", width = menuwidth, textsize = 50)
    layout[3, 7] = Label(scene, "Slice Statistics", textsize = 50)
    hslicescene =
        layout[4, 7] = Axis(
            scene,
            xlabel = @lift(statenames[$statenode] * " " * units[$statenode]),
            xlabelcolor = :black,
            ylabel = "pdf",
            ylabelcolor = :black,
            xlabelsize = 40,
            ylabelsize = 40,
            xticklabelsize = statlabelsize[1],
            yticklabelsize = statlabelsize[2],
            xtickcolor = :black,
            ytickcolor = :black,
            xticklabelcolor = :black,
            yticklabelcolor = :black,
        )

    slicehistogram_node = @lift(histogram($sliced_state, bins = bins))
    xs = @lift($slicehistogram_node[1])
    ys = @lift($slicehistogram_node[2])
    pdf = GLMakie.AbstractPlotting.barplot!(
        hslicescene,
        xs,
        ys,
        color = :blue,
        strokecolor = :blue,
        strokewidth = 1,
    )

    @lift(GLMakie.AbstractPlotting.xlims!(hslicescene, extrema($xs)))
    @lift(GLMakie.AbstractPlotting.ylims!(hslicescene, extrema($ys)))
    vlines!(
        hslicescene,
        @lift($sliceclims[1]),
        color = :black,
        linewidth = menuwidth / 100,
    )
    vlines!(
        hslicescene,
        @lift($sliceclims[2]),
        color = :black,
        linewidth = menuwidth / 100,
    )

    interpolationnames = ["contour", "heatmap"]
    interpolationchoices = [true, false]
    interpolationnode = Node(interpolationchoices[1])
    interpolationmenu =
        Menu(scene, options = zip(interpolationnames, interpolationchoices))

    on(interpolationmenu.selection) do s
        interpolationnode[] = s
        # hack
        heatmap!(
            slicescene,
            slicexaxis,
            sliceyaxis,
            sliced_state,
            interpolate = s,
            colormap = slicecolormapnode,
            colorrange = sliceclims,
        )
    end

    directionmenu = Menu(scene, options = zip(directionnames, directionindex))

    on(directionmenu.selection) do s
        directionnode[] = s
    end

    slicemenustring = @lift(
        directionnames[$directionnode] *
        " at index " *
        string(round(Int, 1 + $slice_node * ($nr[$directionnode] - 1)))
    )
    lowerclim_string = @lift(
        "quantile = " *
        @sprintf("%0.2f", $lowerclim_node) *
        ", value = " *
        @sprintf("%0.1e", $clims[1])
    )
    upperclim_string = @lift(
        "quantile = " *
        @sprintf("%0.2f", $upperclim_node) *
        ", value = " *
        @sprintf("%0.1e", $clims[2])
    )
    alphastring = @lift("Slice alpha = " * @sprintf("%0.2f", $alphanode))
    layout[2, 1] = vgrid!(
        Label(scene, "State", width = nothing),
        statemenu,
        Label(scene, "Color", width = nothing),
        colormenu,
        Label(scene, "Slice Direction", width = nothing),
        directionmenu,
        Label(scene, alphastring, width = nothing),
        alpha_slider,
        Label(scene, slicemenustring, width = nothing),
        slice_slider,
        Label(scene, lowerclim_string, width = nothing),
        lowerclim_slider,
        Label(scene, upperclim_string, width = nothing),
        upperclim_slider,
    )

    slicelowerclim_string = @lift(
        "quantile = " *
        @sprintf("%0.2f", $slicelowerclim_node) *
        ", value = " *
        @sprintf("%0.1e", $sliceclims[1])
    )
    sliceupperclim_string = @lift(
        "quantile = " *
        @sprintf("%0.2f", $sliceupperclim_node) *
        ", value = " *
        @sprintf("%0.1e", $sliceclims[2])
    )

    layout[2, 7] = vgrid!(
        Label(scene, "Contour Plot Type", width = nothing),
        interpolationmenu,
        Label(scene, slicelowerclim_string, width = nothing),
        slicelowerclim_slider,
        Label(scene, sliceupperclim_string, width = nothing),
        sliceupperclim_slider,
        cbar,
    )

    display(scene)
    return scene
end

# quickviz
function visualize(simulation::Simulation; statenames =[string(i) for i in 1:size(simulation.state)[2]], resolution = (32,32,32) )
    a_, statesize, b_ = size(simulation.state)
    mpistate = simulation.state
    grid = simulation.model.grid
    grid_helper = GridHelper(grid)
    r  = coordinates(grid)
    states = []
    ϕ = ScalarField(copy(r[1]), grid_helper)
    r = uniform_grid(Ω, resolution = resolution)
    # statesymbol = vars(Q).names[i] # doesn't work for vectors
    for i in 1:statesize
        ϕ .= mpistate[:, i, :]
        ϕnew = ϕ(r...)
        push!(states, ϕnew)
    end
    visualize([states...], statenames = statenames)
end
