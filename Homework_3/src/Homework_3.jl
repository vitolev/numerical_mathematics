module Homework_3

export rk4_method, distance_to_M, distance_to_m, derivative_vector, rk4_step

"""
    distance_to_M(X, mu)

Calculate the distance from the spacecraft to the Earth.
Arguments:
- `X`: State vector of the spacecraft.
- `mu`: Gravitational parameter of the Earth-Moon system.
Returns the distance in Earth-Moon normalized units.
"""
function distance_to_M(X, mu)
    x, y, z, _, _, _ = X
    return sqrt((x+mu)^2 + y^2 + z^2)
end

"""
    distance_to_m(X, mu)

Calculate the distance from the spacecraft to the Moon.
Arguments:
- `X`: State vector of the spacecraft.
- `mu`: Gravitational parameter of the Earth-Moon system.
Returns the distance in Earth-Moon normalized units.
"""
function distance_to_m(X, mu)
    x, y, z, _, _, _ = X
    return sqrt((x-(1-mu))^2 + y^2 + z^2)
end

"""
    derivative_vector(X, mu, R, r)

Calculate the derivatives of the state vector.
Arguments:
- `X`: State vector of the spacecraft.
- `mu`: Gravitational parameter of the Earth-Moon system.
- `R`: Distance to the Earth.
- `r`: Distance to the Moon.
Returns the derivatives of the state vector.
"""
function derivative_vector(X, mu, R, r)
    x, y, z, alpha, beta, gamma = X
    der_alpha = x + 2*beta - (1-mu)*(x+mu)/R^3 - mu*(x-(1-mu))/r^3
    der_beta = y - 2*alpha - (1-mu)*y/R^3 - mu*y/r^3
    der_gamma = -(1-mu)*z/R^3 - mu*z/r^3
    der_x = alpha
    der_y = beta
    der_z = gamma
    return [der_x, der_y, der_z, der_alpha, der_beta, der_gamma]
end

"""
    rk4_step(X, mu, dt)

Perform one Runge-Kutta 4th order step.
Arguments:
- `X`: Current state vector of the spacecraft.
- `mu`: Gravitational parameter of the Earth-Moon system.
- `dt`: Time step for the RK4 method.
Returns the updated state vector after one RK4 step.
"""
function rk4_step(X, mu, dt)
    # Helper to compute derivatives
    function dX(X)
        R = distance_to_M(X, mu)
        r = distance_to_m(X, mu)
        return derivative_vector(X, mu, R, r)
    end

    k1 = dX(X)
    k2 = dX(X .+ 0.5dt .* k1)
    k3 = dX(X .+ 0.5dt .* k2)
    k4 = dX(X .+ dt .* k3)

    return X .+ dt/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
end

"""
    rk4_method(X0, mu, dt, n_steps; check_collision=false)

Run the RK4 method for `n_steps` with initial state `X0`.
Arguments:
- `X0`: Initial state vector of the spacecraft.
- `mu`: Gravitational parameter of the Earth-Moon system.
- `dt`: Time step for the RK4 method.
- `n_steps`: Number of RK4 steps to perform.
- `check_collision`: If true, checks for collisions with Earth or Moon and stops if a collision occurs.
Returns the trajectory, minimum distance to the Moon, and number of steps taken.
"""
function rk4_method(X0, mu, dt, n_steps; check_collision=false)
    X = X0
    trajectory = [X]
    min_dist_m = 10e6  # Large initial value
    steps = n_steps
    for i in 1:n_steps
        X = rk4_step(X, mu, dt)
        R = distance_to_M(X, mu)
        r = distance_to_m(X, mu)
        if r < min_dist_m
            min_dist_m = r
        end
        if check_collision
            if R < 6400/384400
                println("Collision detected with the Earth at step $i")
                steps = i
                break
            end
            if r < 1700/384400
                println("Collision detected with the Moon at step $i")
                steps = i
                break
            end
        end
        push!(trajectory, X)
    end
    return trajectory, min_dist_m * 384400, steps
end

end # module Homework_3
