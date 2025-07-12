module Homework_3

export euler_method, rk4_method

function distance_to_M(X, mu)
    x, y, z, _, _, _ = X
    return sqrt((x+mu)^2 + y^2 + z^2)
end

function distance_to_m(X, mu)
    x, y, z, _, _, _ = X
    return sqrt((x-(1-mu))^2 + y^2 + z^2)
end

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

function euler_step(X, mu, dt)
    R = distance_to_M(X, mu)
    r = distance_to_m(X, mu)
    derivatives = derivative_vector(X, mu, R, r)
    return X .+ dt .* derivatives
end

function euler_method(X0, mu, dt, n_steps)
    X = X0
    trajectory = [X]
    for _ in 1:n_steps
        X = euler_step(X, mu, dt)
        push!(trajectory, X)
    end
    return trajectory
end

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

function rk4_method(X0, mu, dt, n_steps; check_collision=false)
    X = X0
    trajectory = [X]
    min_dist_m = 10e6  # Large initial value
    for i in 1:n_steps
        X = rk4_step(X, mu, dt)
        if check_collision
            R = distance_to_M(X, mu)
            r = distance_to_m(X, mu)
            if R < 6400/384400 || r < 1700/384400
                println("Collision detected at step $i: R = $R, r = $r")
                break
            end
            if r < min_dist_m
                min_dist_m = r
            end
        end
        push!(trajectory, X)
    end
    return trajectory, min_dist_m * 384400
end


end # module Homework_3
