using ClimateMachine, MPI
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Topologies
using ClimateMachine.MPIStateArrays

import ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
import ClimateMachine.Mesh.Topologies: StackedBrickTopology, BrickTopology
include("domains.jl")
include("vizinanigans.jl")
include("bigfileofstuff.jl")
# some convenience functions
function convention(a::NamedTuple{(:vertical, :horizontal), T}, ::Val{3}) where T
    return (a.horizontal, a.horizontal, a.vertical)
end

function convention(a::Number, ::Val{3})
    return (a, a, a)
end

function convention(a::NamedTuple{(:vertical, :horizontal), T}, ::Val{2}) where T
    return (a.horizontal, a.vertical)
end

function convention(a::Number, ::Val{2})
    return (a, a)
end

function convention(a::Tuple, b)
    return a
end

# brick range brickbuilder
function uniformbrickbuilder(Ω, elements; FT = Float64)
    dimension = ndims(Ω)
    tuple_ranges = []
    for i in 1:dimension
        push!(tuple_ranges, range(FT(Ω[i].a); length = elements[i] + 1, stop = FT(Ω[i].b)))
    end
    brickrange = Tuple(tuple_ranges)
    return brickrange
end

# Grid Constructor
"""
function DiscontinuousSpectralElementGrid(Ω::ProductDomain; elements = nothing, polynomialorder = nothing)
# Description 
Computes a DiscontinuousSpectralElementGrid as specified by a product domain
# Arguments
-`Ω`: A product domain object
# Keyword Arguments 
-`elements`: A tuple of integers ordered by (Nx, Ny, Nz) for number of elements
-`polynomialorder`: A tupe of integers ordered by (npx, npy, npz) for polynomial order
-`FT`: floattype, assumed Float64 unless otherwise specified
-`topology`: default = StackedBrickTopology
-`mpicomm`: default = MPI.COMM_WORLD
-`array`: default = ClimateMachine.array_type()
-`brickbuilder`: default = uniformbrickbuilder, 
  brickrange=uniformbrickbuilder(Ω, elements)
# Return 
A DiscontinuousSpectralElementGrid object
"""
function DiscontinuousSpectralElementGrid(
    Ω::ProductDomain; 
    elements = nothing, 
    polynomialorder = nothing, 
    FT = Float64,         
    mpicomm = MPI.COMM_WORLD, 
    array = ClimateMachine.array_type(),
    topology = StackedBrickTopology,
    brickbuilder = uniformbrickbuilder
    )

    if elements==nothing
        error_message = "Please specify the number of elements as a tuple whose size is commensurate with the domain,"
        error_message *= " e.g., a 3 dimensional domain would need a specification like elements = (10,10,10)."
        error_message *= " or elements = (vertical = 8, horizontal = 5)"
        @error(error_message)
        return nothing
    end

    if polynomialorder==nothing
        error_message = "Please specify the polynomial order as a tuple whose size is commensurate with the domain,"
        error_message *= "e.g., a 3 dimensional domain would need a specification like polynomialorder = (3,3,3)."
        error_message *= " or polynomialorder = (vertical = 8, horizontal = 5)"
        @error(error_message)
        return nothing
    end

    dimension = ndims(Ω)

    if (dimension < 2) || (dimension > 3)
        error_message = "SpectralElementGrid only works with dimensions 2 or 3. "
        error_message *= "The current dimension is " * string(ndims(Ω))
        println("The domain is ", Ω)
        @error(error_message)
        return nothing
    end

    elements = convention(elements, Val(dimension))
    if ndims(Ω) != length(elements)
        @error("Incorrectly specified elements for the dimension of the domain")
        return nothing
    end

    polynomialorder = convention(polynomialorder, Val(dimension))
    if ndims(Ω) != length(polynomialorder)
        @error("Incorrectly specified polynomialorders for the dimension of the domain")
        return nothing
    end

    brickrange = brickbuilder(Ω, elements, FT = FT)

    if dimension == 2
        boundary = ((1,2), (3,4))
    else
        boundary = ((1,2), (3,4), (5,6))
    end

    periodicity = periodicityof(Ω)

    topl = topology(
        mpicomm,
        brickrange;
        periodicity = periodicity,
        boundary = boundary
    )

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = array,
        polynomialorder = polynomialorder,
    )

    return grid
end

## 
#=
# perhaps return wrapper to dg_grid instead
ClimateMachine.init()
Ω = Periodic(0,1) × Interval(0,1) × Periodic(0,1)
dggrid = DiscontinuousSpectralElementGrid(
    Ω, 
    elements = (3,4,5),
    polynomialorder = (1, 2, 3), 
    topology = BrickTopology,
)
#
Ω = Periodic(0,1) × Interval(0,1) × Periodic(0,1)
dggrid = DiscontinuousSpectralElementGrid(
    Ω, 
    elements = (vertical = 1, horizontal = 2),
    polynomialorder = 3, 
)
#
Ω = Periodic(0,1) × Interval(0,1) × Periodic(0,1)
dggrid = DiscontinuousSpectralElementGrid(
    Ω, 
    elements = (vertical = 1, horizontal = 4),
    polynomialorder = (vertical = 2, horizontal = 3), 
    topology = BrickTopology,
)
#
Ω = Periodic(0,1) × Periodic(0,1)
dggrid = DiscontinuousSpectralElementGrid(
    Ω, 
    elements = (vertical = 1, horizontal = 2),
    polynomialorder = 1, 
)
=#


"""
function visualize(g::DiscontinuousSpectralElementGrid)
# Description
- Use Makie to visualize a 3D grid
# Arguments
- `g`: A DiscontinuousSpectralElementGrid
"""
function visualize(g::DiscontinuousSpectralElementGrid)
    x, y, z = coordinates(g)

    e = collect(1:size(x)[2])
    nfaces = 6
    faces = collect(1:nfaces)
    opacities = [0.0, 0.1, 1.0]

    scene, layout = layoutscene()

    element = Node{Any}(e[1])
    face = Node{Any}(faces[6])
    opacity = Node{Any}(opacities[1])

    tmpx = @lift(x[:, $element])
    tmpy = @lift(y[:, $element])
    tmpz = @lift(z[:, $element])

    tmpxf = @lift(x[g.vmap⁻[:, $face, $element]])
    tmpyf = @lift(y[g.vmap⁻[:, $face, $element]])
    tmpzf = @lift(z[g.vmap⁻[:, $face, $element]])

    total_color = @lift(RGBAf0(0,0,0, $opacity))

    lscene = layout[1:4, 2:4] = LScene(scene)
    GLMakie.scatter!(lscene, x[:], y[:], z[:], color = total_color, markersize = 100.0, strokewidth = 0)
    GLMakie.scatter!(lscene, tmpx, tmpy, tmpz, color = RGBAf0(0,0,0,0.5), markersize = 100.0, strokewidth = 0, camera = cam3d!)
    GLMakie.scatter!(lscene, tmpxf, tmpyf, tmpzf, color = RGBAf0(1,0,0,1.0), markersize = 100.0, strokewidth = 0, camera = cam3d!)
    supertitle = layout[1,2] = Label(scene, " "^10 * " Gauss-Lobatto Points " * " "^10, textsize = 50, color = :black)

    menu  = Menu(scene, options = zip(e,e))
    menu2 = Menu(scene, options = zip(faces,faces))
    menu3 = Menu(scene, options = zip(opacities,opacities))
    layout[1, 1] = vgrid!(
        Label(scene, "Element", width = nothing),
        menu,
        Label(scene, "Face", width = nothing),
        menu2,
        Label(scene, "Opacity", width = nothing),
        menu3,
    )
    on(menu.selection) do s
        element[] = s
    end
    on(menu2.selection) do s
        face[] = s
    end
    on(menu3.selection) do s
        opacity[] = s
    end
    display(scene)
    return scene
end