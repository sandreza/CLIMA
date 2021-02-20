abstract type AbstractDissipation end

struct Laplacian{S} <: AbstractDissipation 
    model::S
end