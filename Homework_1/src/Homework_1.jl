module Homework_1

export SimTridiag, Tridiag, ZgornjeTridiag, Givens, qr, eigen

import Base: getindex, setindex!, firstindex, lastindex, *, isapprox, show, convert

########################################
############## SimTridiag ##############    
########################################
"""
    SimTridiag(gd, sd)

Data type for a symmetric tridiagonal matrix with the given main diagonal `gd` and sub-diagonal `sd`.
"""
struct SimTridiag
    gd  # main diagonal
    sd  # sub-diagonal (both lower and upper, since it's symmetric)
end

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

##########################################
############## Tridiag ###################
##########################################
"""
    Tridiag(sd, d, zd)

Data type for a tridiagonal matrix with the given sub-diagonal `sd`, main diagonal `gd`, and super-diagonal `zd`.
"""
struct Tridiag
    sd
    gd
    zd

    function Tridiag(sd, gd, zd)
        n = length(gd)
        if length(sd) != n - 1 || length(zd) != n - 1
            throw(DimensionMismatch("Diagonal and off-diagonal lengths do not match."))
        end
        new(sd, gd, zd)
    end
end

"""
    getindex(T::Tridiag, i::Int, j::Int)

Get the element at position (i, j) of the tridiagonal matrix `T`.
"""
function getindex(T::Tridiag, i::Int, j::Int)
    n = length(T.d)
    if i < 1 || i > n || j < 1 || j > n
        throw(BoundsError(T, (i, j)))
    end

    if i == j
        return T.d[i]
    elseif i == j + 1
        return T.sd[j]
    elseif i == j - 1
        return T.zd[i]
    else
        return 0.0
    end
end

"""
    setindex!(T::Tridiag, value, i::Int, j::Int)

Set the element at position (i, j) of the tridiagonal matrix `T` to `value`.
If `i` and `j` are out of bounds, a `BoundsError` is thrown.
"""
function setindex!(T::Tridiag, v, i::Int, j::Int)
    n = length(T.d)
    if i < 1 || i > n || j < 1 || j > n
        throw(BoundsError(T, (i, j)))
    end

    if i == j
        T.d[i] = v
    elseif i == j + 1
        T.sd[j] = v
    elseif i == j - 1
        T.zd[i] = v
    else
        throw(ArgumentError("Only diagonal and off-diagonal elements can be set."))
    end
end

"""
    *(T::Tridiag, x::Vector)

Multiply the tridiagonal matrix `T` by a vector `x`.
Returns a new vector containing the result of the multiplication.
"""
function *(T::Tridiag, v::Vector)
    n = length(T.d)
    if length(v) != n
        throw(DimensionMismatch("Vector length must be equal to the number of columns in the matrix."))
    end

    b = zero(v)
    b[1] = T[1,1] * v[1] + T[1,2] * v[2]
    for i = 2:n-1
        b[i] = T[i,i-1] * v[i-1] + T[i,i] * v[i] + T[i,i+1] * v[i+1]
    end
    b[n] = T[n,n-1] * v[n-1] + T[n,n] * v[n]
    return b
end

"""
    convert(::Type{Tridiag}, T::SimTridiag)

Convert a `SimTridiag` matrix to a `Tridiag` matrix.
"""
function Base.convert(::Type{Tridiag}, T::SimTridiag)
    n = length(T.gd)
    Tridiag(copy(T.sd), copy(T.gd), copy(T.sd))  # sd is used as both sub- and super-diagonal
end

###############################################
############## ZgornjeTridiag #################
###############################################
"""
    ZgornjeTridiag(gd, sd, sd2)

Data type for an upper triangular matrix with the given diagonal `gd` first superdiagonal `sd` and second superdiagonal `sd2`.
"""
struct ZgornjeTridiag
    gd
    sd
    sd2

    function ZgornjeTridiag(gd, sd, sd2)
        n = length(gd)
        if length(sd) != n - 1 || length(sd2) != n - 2
            throw(DimensionMismatch("Diagonal and superdiagonal lengths do not match."))
        end
        new(gd, sd, sd2)
    end
end

"""
    getindex(R::ZgornjeTridiag, i::Int, j::Int)

Get the element at position (i, j) of the upper triangular matrix `R`.
If `i` and `j` are out of bounds, a `BoundsError` is thrown.
"""
function getindex(R::ZgornjeTridiag, i::Int, j::Int)
    n = length(R.gd)
    if i < 1 || i > n || j < 1 || j > n
        throw(BoundsError(R, (i, j)))
    end

    if i == j
        return R.gd[i]
    elseif i == j - 1
        return R.sd[i]
    elseif i == j - 2
        return R.sd2[i]
    else
        return 0.0
    end
end

"""
    setindex!(R::ZgornjeTridiag, value, i::Int, j::Int)

Set the element at position (i, j) of the upper triangular matrix `R` to `value`.
If `i` and `j` are out of bounds, a `BoundsError` is thrown.
"""
function setindex!(R::ZgornjeTridiag, v, i::Int, j::Int)
    n = length(R.gd)
    if i < 1 || i > n || j < 1 || j > n
        throw(BoundsError(R, (i, j)))
    end

    if i == j
        R.gd[i] = v
    elseif i == j - 1
        R.sd[i] = v
    elseif i == j - 2
        R.sd2[i] = v
    else
        throw(ArgumentError("Only diagonal and superdiagonal elements can be set."))
    end
end

"""
    *(R::ZgornjeTridiag, x::Vector)

Multiply the upper triangular matrix `R` by a vector `x`.
Returns a new vector containing the result of the multiplication.
"""
function Base.:*(R::ZgornjeTridiag, x::Vector)
    n = length(R.gd)
    if length(x) != n
        throw(DimensionMismatch("Matrix and vector dimensions do not match."))
    end
    y = zeros(n)
    for i in 1:n
        y[i] = R.gd[i] * x[i]
        if i < n
            y[i] += R.sd[i] * x[i+1]
        end
        if i < n - 1
            y[i] += R.sd2[i] * x[i+2]
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
    n = length(R.gd)
    if size(A, 1) != n
        throw(DimensionMismatch("Matrix dimensions do not match."))
    end
    m = size(A, 2)
    B = zeros(n, m)
    for j in 1:m
        for i in 1:n
            B[i, j] = R.gd[i] * A[i, j]
            if i < n    # First superdiagonal contribution
                B[i, j] += R.sd[i] * A[i+1, j]
            end
            if i < n - 1 # Second superdiagonal contribution
                B[i, j] += R.sd2[i] * A[i+2, j]
            end
        end
    end
    return B
end

"""
    isapprox(R1::ZgornjeTridiag, R2::ZgornjeTridiag)

Check if two upper triangular matrices `R1` and `R2` are approximately equal.
Returns `true` if they are approximately equal, `false` otherwise.
"""
function isapprox(R1::ZgornjeTridiag, R2::ZgornjeTridiag)
    if length(R1.gd) != length(R2.gd) || length(R1.sd) != length(R2.sd) || length(R1.sd2) != length(R2.sd2)
        return false
    end
    for i in 1:length(R1.gd)
        if !isapprox(R1.gd[i], R2.gd[i])
            return false
        end
    end
    for i in 1:length(R1.sd)
        if !isapprox(R1.sd[i], R2.sd[i])
            return false
        end
    end
    for i in 1:length(R1.sd2)
        if !isapprox(R1.sd2[i], R2.sd2[i])
            return false
        end
    end
    return true
end

################################################
################# Givens #######################
################################################
"""
    Givens(rotations)

Data type for a Givens rotation matrix, represented by a list of rotations. 
The order of rotations is such that if we express the Givens rotation matrix as a product of rotation matrices,
    then the order of the rotations in the list is the order from left to right in the product.
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
    *(Q::Givens, x::Vector)
    
Multiply the Givens rotation matrix `Q` by a vector `x`.
Returns a new vector containing the result of the multiplication.
"""
function Base.:*(Q::Givens, x::Vector)
    y = copy(x)
    for (c, s, i, j) in reverse(Q.rotations)
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
    for (c, s, i, j) in reverse(Q.rotations)
        for col in 1:size(B, 2)
            a = B[i, col]
            b = B[j, col]
            B[i, col] = c * a - s * b
            B[j, col] = s * a + c * b
        end
    end
    return B
end

"""
    *(A::Matrix, Q::Givens)

Multiply a matrix `A` by the Givens rotation matrix `Q`.
Returns a new matrix containing the result of the multiplication.
"""
function Base.:*(A::Matrix, Q::Givens)
    B = copy(A)
    for (c, s, i, j) in Q.rotations
        for col in 1:size(B, 1)
            a = B[col, i]
            b = B[col, j]
            B[col, i] = c * a + s * b
            B[col, j] = -s * a + c * b
        end
    end
    return B
end

"""
    isapprox(G1::Givens, G2::Givens)

Check if two Givens rotation matrices `G1` and `G2` are approximately equal.
Returns `true` if they are approximately equal, `false` otherwise.
"""
function isapprox(G1::Givens, G2::Givens; atol=1e-10, rtol=1e-6)
    if length(G1.rotations) != length(G2.rotations)
        return false
    end
    for (r1, r2) in zip(G1.rotations, G2.rotations)
        if !isapprox(r1[1], r2[1]; atol=atol, rtol=rtol) || !isapprox(r1[2], r2[2]; atol=atol, rtol=rtol)
            return false
        end
        # Indices i and j in Givens rotations can be swapped
        if r1[3] == r2[3] && r1[4] == r2[4]
            continue
        elseif r1[3] == r2[4] && r1[4] == r2[3]
            continue
        else
            return false
        end
    end
    return true
end

####################################################################

"""
    *(Q::Givens, R::ZgornjeTridiag)

Multiply the Givens rotation matrix `Q` by the upper triangular matrix `R`.
Returns a Tridiag matrix if product of Givens and ZgornjeTridiag results in a tridiagonal matrix,
    otherwise returns a full matrix. 
"""
function Base.:*(Q::Givens, R::ZgornjeTridiag)
    n = length(R.gd)
    B = [R[i, j] for i in 1:n, j in 1:n]
    A = Q * B
    # Check if the result is tridiagonal
    eps = 1e-12
    is_tridiagonal = true
    for i in 1:n
        for j in 1:n
            if abs(i - j) > 1 && abs(A[i, j]) > eps
                is_tridiagonal = false
                break
            end
        end
        if !is_tridiagonal
            break
        end
    end
    if is_tridiagonal
        return Tridiag([A[i+1, i] for i in 1:n-1], [A[i, i] for i in 1:n], [A[i, i+1] for i in 1:n-1])
    else
        return A
    end
end

"""
    *(R::ZgornjeTridiag, Q::Givens)

Multiply the upper triangular matrix `R` by the Givens rotation matrix `Q`.
Returns a Tridiag matrix if product of Givens and ZgornjeTridiag results in a tridiagonal matrix,
    otherwise returns a full matrix. 
"""
function Base.:*(R::ZgornjeTridiag, Q::Givens)
    n = length(R.gd)
    B = [R[i, j] for i in 1:n, j in 1:n]
    A = B * Q
    # Check if the result is tridiagonal
    eps = 1e-12
    is_tridiagonal = true
    for i in 1:n
        for j in 1:n
            if abs(i - j) > 1 && abs(A[i, j]) > eps
                is_tridiagonal = false
                break
            end
        end
        if !is_tridiagonal
            break
        end
    end
    if is_tridiagonal
        return Tridiag([A[i+1, i] for i in 1:n-1], [A[i, i] for i in 1:n], [A[i, i+1] for i in 1:n-1])
    else
        return A
    end
end

"""
    qr(T::Tridiag)

Perform QR decomposition of the tridiagonal matrix `T` using Givens rotations.
Returns a tuple (Q, R) where `Q` is an orthogonal matrix of type Givens and `R` is a matrix of type ZgornjeTridiag.
"""
function qr(T::Tridiag)
    n = length(T.gd)
    gd = Float64.(T.gd)   # Convert to Float64
    sd = Float64.(T.zd)   # Convert to Float64
    sd2 = zeros(Float64, n - 2)  # Allocate zeros as Float64
    R = ZgornjeTridiag(gd, sd, sd2)
    rotations = Vector{Tuple{Float64, Float64, Int, Int}}()

    for i in 1:n-1
        a = R.gd[i]
        b = T.sd[i]
        if iszero(b)
            continue
        end

        # Compute Givens rotation
        r = sqrt(a^2 + b^2)
        c, s = (a/r, b/r)
        push!(rotations, (c, s, i, i+1))

        R[i, i] = c * a + s * b
        old_R_i_i1 = R[i, i+1]
        old_R_i1_i1 = R[i+1, i+1]
        R[i, i+1] = c * old_R_i_i1 + s * old_R_i1_i1
        R[i+1, i+1] = -s * old_R_i_i1 + c * old_R_i1_i1

        if i < n - 1
            old_R_i1_i2 = R[i+1, i+2]
            R[i+1, i+2] = c * old_R_i1_i2
            R[i, i+2] = s * old_R_i1_i2
        end
    end

    Q = Givens(rotations)
    return Q, R
end

"""
    qr(T::SimTridiag)

Perform QR decomposition of the symmetric tridiagonal matrix `T` using Givens rotations, by converting it to a `Tridiag` matrix
and then applying the QR decomposition.
Returns a tuple (Q, R) where `Q` is an orthogonal matrix of type Givens and `R` is a matrix of type ZgornjeTridiag.
"""
function qr(T::SimTridiag)
    qr(convert(Tridiag, T))
end

"""
    eigen(T::SimTridiag; max_iter=1000)

Compute the eigenvalues and eigenvectors of the symmetric tridiagonal matrix `T` using the QR algorithm.
Returns a tuple (eigenvalues, eigenvectors) where `eigenvalues` is a vector of eigenvalues and `eigenvectors` is a list of vectors.
"""
function eigen(T::SimTridiag; max_iter=1000)
    n = length(T.gd)
    V = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]  # Initialize eigenvector matrix as identity

    # Because T is symmetric tridiagonal, the product R * Q is also symmetric tridiagonal.
    Q, R = qr(T)
    for _ in 1:max_iter
        T = R * Q
        V = V * Q
        Q, R = qr(T)
    end
    T = R * Q
    V = V * Q
    eigenvalues = T.gd

    eigenvectors_list = [V[:, i] for i in 1:n]

    return eigenvalues, eigenvectors_list
end

end # module Homework_1
