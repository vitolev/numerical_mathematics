module Homework_1

export SimTridiag, ZgornjeTridiag

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

"""
    *(T::SimTridiag, x::Vector)

Multiply the symmetric tridiagonal matrix `T` by a vector `x`.
Returns a new vector containing the result of the multiplication.
"""
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

"""
    *(T::SimTridiag, x::Matrix)

Multiply the symmetric tridiagonal matrix `T` by a matrix `x`.
Returns a new matrix containing the result of the multiplication.
"""
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

export ZgornjeTridiag, Givens

"""
    ZgornjeTridiag(diag, superdiag)

Data type for an upper triangular matrix with the given diagonal `diag` and superdiagonal `superdiag`.
"""
struct ZgornjeTridiag
    diag::Vector{Float64}
    superdiag::Vector{Float64}
end

"""
    Givens(rotations)

Data type for a Givens rotation matrix, represented by a list of rotations. 
The order of rotations is such that the first rotation (first element in the list) is applied first, then second, and so on.
Each rotation is a tuple (c, s, i, j) where c is the cosine, s is the sine, and (i,j) are the rows / columns being rotated.
If passed i > j it will store the rotation as (c, s, j, i) instead.
"""
struct Givens
    rotations::Vector{Tuple{Float64, Float64, Int, Int}}

    function Givens(rotations)
        new_rotations = Vector{Tuple{Float64, Float64, Int, Int}}(undef, length(rotations))
        for (k, rotation) in enumerate(rotations)
            c, s, i, j = rotation
            if !isapprox(c^2 + s^2, 1.0)
                throw(ArgumentError("Cosine and sine values must satisfy c^2 + s^2 = 1 (got c=$c, s=$s)"))
            end
            if i == j
                throw(ArgumentError("Cannot rotate the same row/column."))
            end
            if i > j
                i, j = j, i
            end
            new_rotations[k] = (c, s, i, j)
        end
        new(new_rotations)
    end
end

"""
    *(R::ZgornjeTridiag, x::Vector)

Multiply the upper triangular matrix `R` by a vector `x`.
Returns a new vector containing the result of the multiplication.
"""
function Base.:*(R::ZgornjeTridiag, x::Vector)
    n = length(R.diag)
    if length(x) != n
        throw(DimensionMismatch("Matrix and vector dimensions do not match."))
    end
    y = zeros(n)
    for i in 1:n
        y[i] = R.diag[i] * x[i]
        if i < n
            y[i] += R.superdiag[i] * x[i+1]
        end
    end
    return y
end

"""
    *(R::ZgornjeTridiag, A::Matrix)

Multiply the upper triangular matrix `R` by a matrix `A`.
Returns a new matrix containing the result of the multiplication.
"""
function Base.:*(R::ZgornjeTridiag, A::Matrix)
    n = length(R.diag)
    if size(A, 1) != n
        throw(DimensionMismatch("Matrix dimensions do not match."))
    end
    m = size(A, 2)
    B = zeros(n, m)
    for j in 1:m
        for i in 1:n
            B[i, j] = R.diag[i] * A[i, j]
            if i < n
                B[i, j] += R.superdiag[i] * A[i+1, j]
            end
        end
    end
    return B
end

"""
    *(Q::Givens, x::Vector)
    
Multiply the Givens rotation matrix `Q` by a vector `x`.
Returns a new vector containing the result of the multiplication.
"""
function Base.:*(Q::Givens, x::Vector)
    y = copy(x)
    for (c, s, i, j) in Q.rotations
        yi = y[i]
        yj = y[j]
        y[i] = c * yi - s * yj
        y[j] = s * yi + c * yj
    end
    return y
end

"""
    *(Q::Givens, A::Matrix)

Multiply the Givens rotation matrix `Q` by a matrix `A`.
Returns a new matrix containing the result of the multiplication.
"""
function Base.:*(Q::Givens, A::Matrix)
    B = copy(A)
    for (c, s, i, j) in Q.rotations
        for col in 1:size(B, 2)
            a = B[i, col]
            b = B[j, col]
            B[i, col] = c * a - s * b
            B[j, col] = s * a + c * b
        end
    end
    return B
end

end # module Homework_1
