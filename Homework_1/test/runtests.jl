using Test
using Homework_1

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

@testset "ZgornjeTridiag is defined" begin
    diag = [1, 2, 3]
    superdiag = [4, 5]
    T = ZgornjeTridiag(diag, superdiag)
    @test T.diag == diag
    @test T.superdiag == superdiag
end

@testset "Givens is defined" begin
    c1 = 1
    s1 = 0
    c2 = sqrt(2)/2
    s2 = sqrt(2)/2
    rotations = [(c1, s1, 3, 1), (c2, s2, 4, 1)]
    G = Givens(rotations)
    @test length(G.rotations) == 2

    rotations = [(c1, s1, 1, 2), (c2, s2, 3, 4)]
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