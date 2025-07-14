using Homework_3
using Test

@testset "distances to Earth and Moon" begin
    X = [-0.1, 0.0, 0.0, 0, 0, 0]
    mu = 0.012150585609624  # Gravitational parameter for the Earth-Moon system

    @test isapprox(distance_to_M(X, mu), 0.1 - mu)
    @test isapprox(distance_to_m(X, mu), 1.1 - mu)
end

@testset "derivative vector calculation" begin
    X = [-0.1, 0.0, 0.0, 0.01, 0.02, 0.03]
    mu = 0.012150585609624
    R = distance_to_M(X, mu)
    r = distance_to_m(X, mu)

    der = derivative_vector(X, mu, R, r)
    @test length(der) == 6
end

@testset "Runge-Kutta 4th order step" begin
    X = [-0.1, 0.0, 0.0, 0.01, 0.02, 0.03]
    mu = 0.012150585609624
    dt = 0.0001

    new_X = rk4_step(X, mu, dt)
    @test length(new_X) == 6
    @test new_X != X  # The state vector has to change, beacause with such initial state vector the derivative is not 0.
end

@testset "Runge-Kutta 4th order method" begin
    X0 = [-0.02, 0.0, 0.0, 0.0, 10.0, 0.00]
    mu = 0.012150585609624
    dt = 0.0001
    n_steps = 100

    trajectory, min_dist_m, steps = rk4_method(X0, mu, dt, n_steps)
    @test all(length(X) == 6 for X in trajectory)  # Each state vector should have 6 components
end