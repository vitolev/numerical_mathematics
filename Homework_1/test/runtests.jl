using Test
using Homework_1
import LinearAlgebra: norm

@testset "Data type 'SimTriag' is defined" begin
    gd = [1, 2, 3]
    sd = [4, 5]
    T = SimTridiag(gd, sd)
    @test T.gd == gd
    @test T.sd == sd
end

@testset "getindex function works correctly" begin
    gd = [1, -7, 3, 1]
    sd = [4, 5, 9]
    T = SimTridiag(gd, sd)
    
    @test T[1, 1] == 1
    @test T[2, 2] == -7
    @test T[3, 2] == 5
    @test T[2, 3] == 5
    # BoundsError test
    @test_throws BoundsError T[10, 1]
    @test_throws BoundsError T[1, -3]
end

@testset "setindex! function works correctly" begin
    gd = [3, -8, -2, 5]
    sd = [6, 6, -1]
    T = SimTridiag(gd, sd)
    
    T[1, 1] = 10
    @test T[1, 1] == 10
    
    T[1, 1] = -5
    @test T[1, 1] == -5

    T[4, 3] = 7
    @test T[4, 3] == 7
    @test T[3, 4] == 7
    
    # BoundsError test
    @test_throws BoundsError T[10, 1] = 5
    @test_throws BoundsError T[1, -1] = -1
end

@testset "firstindex function works correctly" begin
    gd = [1, 2, 3]
    sd = [4, 5]
    T = SimTridiag(gd, sd)
    
    @test firstindex(T, 1) == 1
    @test firstindex(T, 2) == 1
end

@testset "lastindex function works correctly" begin
    gd = [1, 2, 3]
    sd = [4, 5]
    T = SimTridiag(gd, sd)
    
    @test lastindex(T, 1) == 3
    @test lastindex(T, 2) == 3
end

@testset "Multiplication from right with vector" begin
    gd = [3,2,6,-5]
    sd = [1,0,1]
    T = SimTridiag(gd, sd)
    x1 = [1,2,3,4]
    x2 = [2,0,-1,-1]
    @test T * x1 == [5, 5, 22, -17]
    @test T * x2 == [6, 2, -7, 4]

    # DimensionMismatch test
    @test_throws DimensionMismatch T * [1, 2, 3]
end

@testset "Multiplication from right with matrix" begin
    gd = [3,2,6,-5]
    sd = [1,0,1]
    T = SimTridiag(gd, sd)
    x1 = [1 2; 3 4; 5 6; 7 8]
    x2 = [2 0; -1 -1; 0 1; 1 2]
    @test T * x1 == [6 10; 7 10; 37 44; -30 -34]
    @test T * x2 == [5 -1; 0 -2; 1 8; -5 -9]

    # DimensionMismatch test
    @test_throws DimensionMismatch T * [1 2; 2 0; 1 3]
end

@testset "Data type 'ZgornjeTridiag' is defined" begin
    diag = [1, 2, 3]
    superdiag = [4, 5]
    superdiag2 = [6]
    T = ZgornjeTridiag(diag, superdiag, superdiag2)
    @test T.gd == diag
    @test T.sd == superdiag
    @test T.sd2 == superdiag2
    @test T[1, 1] == 1
    @test T[2, 3] == 5
    @test T[1, 3] == 6
    T[1, 2] = 10
    @test T[1, 2] == 10

    superdiag = [4, 5, 6]
    @test_throws DimensionMismatch ZgornjeTridiag(diag, superdiag, superdiag2)
    superdiag = [6, 1]
    superdiag2 = [4, 5]
    @test_throws DimensionMismatch ZgornjeTridiag(diag, superdiag, superdiag2)


end

@testset "Data type 'Givens' is defined" begin
    c1 = 1
    s1 = 0
    c2 = sqrt(2)/2
    s2 = sqrt(2)/2
    rotations = [(c1, s1, 3, 1), (c2, s2, 4, 1)]
    G = Givens(rotations)
    @test length(G.rotations) == 2

    rotations = [(c1, s1, 1, 1), (c2, s2, 3, 4)]
    @test_throws ArgumentError Givens(rotations)

    c1 = 2
    s1 = 1
    rotations = [(c1, s1, 2, 1)]
    @test_throws ArgumentError Givens(rotations)

    c1 = 0.5
    s1 = 0.5
    rotations = [(c1, s1, 3, 2)]
    @test_throws ArgumentError Givens(rotations)
end

@testset "ZgornjeTridiag multiplication works correctly" begin
    diag = [1, 2, 3, 4]
    superdiag = [5, 6, 7]
    superdiag2 = [8, 9]
    T = ZgornjeTridiag(diag, superdiag, superdiag2)
    x1 = [1, 2, 3, 4]
    x2 = [2, 0, -1, 0]
    @test T * x1 == [35, 58, 37, 16]
    @test T * x2 == [-6, -6, -3, 0]

    # DimensionMismatch test
    @test_throws DimensionMismatch T * [1, 2]

    A = [1 2; 0 3; 2 1; 1 0]
    @test T * A == [17 25; 21 12; 13 3; 4 0]

    # DimensionMismatch test
    @test_throws DimensionMismatch T * [1 2; 3 4; 5 6]
end

@testset "Givens multiplication works correctly" begin
    c1 = 0
    s1 = 1
    G1 = Givens([(c1, s1, 2, 1)])
    x1 = [1, 2, 3]
    @test G1 * x1 == [-2, 1, 3]

    c2 = sqrt(3)/2
    s2 = 1/2
    G2 = Givens([(c1, s1, 2, 1), (c2, s2, 3, 2)])
    x2 = [1.0, 0.0, 1.0]
    @test isapprox(G2 * x2, [0.5, 1, 0.5*sqrt(3)])

    A = [1.0 0.0; 0 1; 2 1]
    At = [1.0 0 2; 0 1 1]
    @test isapprox(G2 * A, [1 0.5*(1-sqrt(3)); 1 0; sqrt(3) 0.5*(1+sqrt(3))])
    @test isapprox(At * G2, [0 0.5*(2-sqrt(3)) 0.5+sqrt(3); 1 0.5 0.5*sqrt(3)])

    c1 = sqrt(2)/2
    s1 = sqrt(2)/2
    G3 = Givens([(c1, s1, 2, 1)])
    A = [1.0 5.0 8.0 0.0; 0 0 0 -4; 0 0 3 7; 0 2 6 9]
    @test isapprox(G3 * A, [sqrt(2)/2 5*sqrt(2)/2 8*sqrt(2)/2 4*sqrt(2)/2; sqrt(2)/2 5*sqrt(2)/2 8*sqrt(2)/2 -4*sqrt(2)/2; 0 0 3 7; 0 2 6 9])
    @test isapprox(A * G3, [3*sqrt(2) 2*sqrt(2) 8 0; 0 0 0 -4; 0 0 3 7; sqrt(2) sqrt(2) 6 9])
end

@testset "Givens multiplication with ZgornjeTridiag works correctly" begin
    diag = [1, 2, 3, 4]
    superdiag = [5, 6, 7]
    superdiag2 = [8, 9]
    R = ZgornjeTridiag(diag, superdiag, superdiag2)
    c1 = 0
    s1 = 1
    G1 = Givens([(c1, s1, 2, 4)])
    @test G1 * R == [1 5 8 0; 0 0 0 -4; 0 0 3 7; 0 2 6 9]
    c2 = sqrt(2)/2
    s2 = sqrt(2)/2
    G2 = Givens([(c2, s2, 1, 2), (c1, s1, 2, 4)])
    @test G2 * R == [sqrt(2)/2 5*sqrt(2)/2 8*sqrt(2)/2 4*sqrt(2)/2; sqrt(2)/2 5*sqrt(2)/2 8*sqrt(2)/2 -4*sqrt(2)/2; 0 0 3 7; 0 2 6 9]
end

@testset "QR decomposition of (Sim)Tridiag works correctly" begin
    gd = [1,1,1]
    zd = [0,1]
    sd = [1,0]
    T = Tridiag(sd, gd, zd)
    Q, R = qr(T)
    @test isapprox(Q, Givens([(sqrt(2)/2, sqrt(2)/2, 1, 2)]))
    @test isapprox(R, ZgornjeTridiag([sqrt(2), sqrt(2)/2, 1], [sqrt(2)/2, sqrt(2)/2], [sqrt(2)/2]))

    gd = [1,2,2,1,1]
    zd = [1,1,1,2]
    sd = [1,1,1,2]
    T = Tridiag(sd, gd, zd)
    Q, R = qr(T)
    @test isapprox(Q, Givens([(1/sqrt(2),1/sqrt(2),1,2),(sqrt(3)/3,sqrt(6)/3,2,3),(0.5,sqrt(3)/2,3,4),(0,1,4,5)]))
    @test isapprox(R, ZgornjeTridiag([sqrt(2), sqrt(3/2), 2/sqrt(3), 2, -1], [3/sqrt(2), 5/sqrt(6), 2/sqrt(3), 1], [1/sqrt(2), sqrt(2/3), sqrt(3)]))

    ST = SimTridiag(gd, sd)
    Q, R = qr(ST)
    @test isapprox(Q, Givens([(1/sqrt(2),1/sqrt(2),1,2),(sqrt(3)/3,sqrt(6)/3,2,3),(0.5,sqrt(3)/2,3,4),(0,1,4,5)]))
    @test isapprox(R, ZgornjeTridiag([sqrt(2), sqrt(3/2), 2/sqrt(3), 2, -1], [3/sqrt(2), 5/sqrt(6), 2/sqrt(3), 1], [1/sqrt(2), sqrt(2/3), sqrt(3)]))
end

@testset "QR iterations work correctly" begin
    gd = [3,2,1]
    sd = [1,1]
    T = SimTridiag(gd, sd)

    eigenvalues, eigenvectors = eigen(T)
    @test isapprox(eigenvalues, [2+sqrt(3), 2, 2-sqrt(3)])
    v1 = [2+sqrt(3), 1+sqrt(3), 1]
    v2 = [-1, 1, 1]
    v3 = [2-sqrt(3), 1-sqrt(3), 1]
    v1 = v1 / norm(v1)
    v2 = v2 / norm(v2)
    v3 = v3 / norm(v3)
    @test isapprox(eigenvectors[1], v1)
    @test isapprox(eigenvectors[2], v2)
    @test isapprox(eigenvectors[3], v3)
end