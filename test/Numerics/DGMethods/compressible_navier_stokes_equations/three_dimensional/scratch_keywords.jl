

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

function convention(a::Tuple)
    return a
end




elements = (vertical = 4, horizontal = 8)
convention(elements, Val(3))
convention(elements, Val(2))
