using Test
using Problem_18_1_5

@testset "Data type 'SimTriag' is defined" begin
    gd = [1, 2, 3]
    sd = [4, 5]
    T = SimTridiag(gd, sd)
    @test T.gd == gd
    @test T.sd == sd
end