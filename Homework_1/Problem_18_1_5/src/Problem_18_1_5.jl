module Problem_18_1_5

export SimTridiag

"""
    SimTridiag(gd, sd)

Data type for a symmetric tridiagonal matrix with the given main diagonal `gd` and sub-diagonal `sd`.
"""
struct SimTridiag
    gd  # main diagonal
    sd  # sub-diagonal (both lower and upper, since it's symmetric)
end

import Base: getindex, setindex!, firstindex, lastindex, *


"""
    getindex(T::SimTridiag, i::Int, j::Int)

Get the element at position (i, j) of the symmetric tridiagonal matrix `T`.
If `i` and `j` are out of bounds, a `BoundsError` is thrown.
"""
function getindex(T::SimTridiag, i::Int, j::Int)
    n = length(T.gd)
    if i < 1 || i > n || j < 1 || j > n
        throw(BoundsError(T, (i, j)))
    end

    if i == j
        return T.gd[i]
    elseif i == j + 1 || i + 1 == j
        return T.sd[min(i, j)]
    else
        return 0.0  # outside the tridiagonal band
    end
end


"""
    setindex!(T::SimTridiag, value, i::Int, j::Int)
Set the element at position (i, j) of the symmetric tridiagonal matrix `T` to `value`.
If `i` and `j` are out of bounds, a `BoundsError` is thrown.
"""
function setindex!(T::SimTridiag, value, i::Int, j::Int)
    n = length(T.gd)
    if i < 1 || i > n || j < 1 || j > n
        throw(BoundsError(T, (i, j)))
    end

    if i == j
        T.gd[i] = value
    elseif i == j + 1 || i + 1 == j
        T.sd[min(i, j)] = value
    else
        throw(BoundsError("Cannot set value outside the tridiagonal band. Only main and sub-diagonal elements can be set."))
    end
end

firstindex(T::SimTridiag, d::Int) = 1
lastindex(T::SimTridiag, d::Int) = length(T.gd)

function *(T::SimTridiag, x::Vector)
    n = length(T.gd)
    if length(x) != n
        throw(DimensionMismatch("Matrix and vector dimensions do not match."))
    end

    result = zeros(n)
    for i in 1:n
        result[i] = T.gd[i] * x[i]
        if i > 1
            result[i] += T.sd[i - 1] * x[i - 1]
        end
        if i < n
            result[i] += T.sd[i] * x[i + 1]
        end
    end
    return result
end

function *(T::SimTridiag, x::Matrix)
    n = length(T.gd)
    m = size(x, 2)
    if size(x, 1) != n
        throw(DimensionMismatch("Matrix dimensions do not match."))
    end

    result = zeros(n, m)
    for j in 1:m
        for i in 1:n
            result[i, j] = T.gd[i] * x[i, j]
            if i > 1
                result[i, j] += T.sd[i - 1] * x[i - 1, j]
            end
            if i < n
                result[i, j] += T.sd[i] * x[i + 1, j]
            end
        end
    end
    return result
end

end # module Problem_18_1_5