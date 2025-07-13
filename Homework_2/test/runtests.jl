using Test
using Homework_2

# ----------------------------
# Helper tests
# ----------------------------

@testset "Simpson Weight Tests" begin
    @test simpson_weight(0, 10) == 1
    @test simpson_weight(10, 10) == 1
    @test simpson_weight(4, 10) == 2
    @test simpson_weight(3, 10) == 4

    @test simpson_weight6D((0, 1, 2, 3, 4, 5), 10) ==
        simpson_weight(0, 10) *
        simpson_weight(1, 10) *
        simpson_weight(2, 10) *
        simpson_weight(3, 10) *
        simpson_weight(4, 10) *
        simpson_weight(5, 10)
end

@testset "Force Integrands" begin
    # Check that force_integrand6D and force_integrand3D give correct values in simple cases
    x1, y1, z1 = 0.0, 0.0, 0.0
    x2, y2, z2 = 1.0, 0.0, 0.0
    @test isapprox(force_integrand6D(x1, y1, z1, x2, y2, z2), -1.0)

    x, y, z = -1.0, 0.0, 0.0
    @test isapprox(force_integrand3D(x, y, z), -1.0)
end

@testset "Volume Overlap" begin
    @test isapprox(volume_overlap(-2.0, 0.0, 0.0), 1.0) # Cubes fully overlap
    @test isapprox(volume_overlap(-3.0, 0.0, 0.0), 0.0) # No overlap
end

@testset "Custom Integrand" begin
    f_val = custom_integrand(-2.0, 0.0)
    @test f_val < 0  # x is negative, force is negative
end

# ----------------------------
# Integration method tests
# ----------------------------

@testset "Integration Consistency" begin
    # We expect these to be roughly equal, not exactly the same
    n_simpson = 4  # low resolution for speed
    n_mc = 10_000  # sufficient for Monte Carlo convergence

    sim6 = simpson6D(n_simpson)
    sim3 = simpson3D(n_simpson)
    sim2 = simpson2D(n_simpson)
    mc   = monte_carlo_integration(n_mc)

    # All methods should give force in the same direction (negative x)
    @test sim6 < 0
    @test sim3 < 0
    @test sim2 < 0
    @test mc < 0

    tol = 0.1   # Tolerance for relative error, just approximate to ensure methods are consistent
    ref = sim2
    @test isapprox(sim3, ref; rtol=tol)
    @test isapprox(sim6, ref; rtol=tol)
    @test isapprox(mc, ref; rtol=tol)
end