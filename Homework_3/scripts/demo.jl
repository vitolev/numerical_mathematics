using Homework_3
using Plots

# Function to run one simulation with given parameters
function run_simulation(fi, v, dt, t_max; check_collision=false)
    mu = 0.012150585609624  # Gravitational parameter for the Earth-Moon system
    r = 6600 / 384400  # Initial orbiting radius of the spacecraft
    X0 = [-mu + r*cos(fi), r*sin(fi), 0.0, v*-sin(fi), v*cos(fi), 0.0]  # Initial state vector

    n_steps = floor(Int, t_max / dt)
    trajectory, min_dist_m, steps = rk4_method(X0, mu, dt, n_steps, check_collision=check_collision)
    println("Minimum distance to the Moon: $min_dist_m km")
    time = steps * dt * 375600 / 3600 / 24  # Convert to days
    println("Time of flight: $time days")
    return trajectory, min_dist_m, steps
end

# Helper function to plot the trajectory
function plot_trajectory(trajectory; label="Trajectory")
    xs = [X[1] for X in trajectory]
    ys = [X[2] for X in trajectory]

    # Plot the trajectory
    plot!(xs, ys, aspect_ratio=:equal, label=label)
end

# Helper function to plot filled circles for Earth and Moon
function plot_filled_circle(x, y, r; color=:blue)
    θ = range(0, 2π, length=100)
    x_circle = x .+ r * cos.(θ)
    y_circle = y .+ r * sin.(θ)
    plot!([Shape(x_circle, y_circle)], color=color, linewidth=0, label="")
end

# Function to compute free return orbits shown in report
function compute_free_return_orbits()
    fi = [227.3, 210, -5, 30]
    v = [10.6459, -10.6823, 10.64555, -10.68257]
    mu = 0.012150585609624
    for i in 1:length(fi)
        traj, min_dist_m, steps = run_simulation(fi[i] * pi / 180, v[i], 0.000001, 8; check_collision=true)

        plot(legend=false, aspect_ratio=:equal, size=(480, 320))
        plot_trajectory(traj)
        plot_filled_circle(-mu, 0.0, 6400/384400; color=:blue)     # Earth (blue)
        plot_filled_circle(1.0 - mu, 0.0, 1700/384400; color=:gray) # Moon (gray)
        
        savefig("free_return_orbit_$i.png")
    end
end

# Function to test the effect of time step on the trajectory as presented in the report
function test_time_step_effects()
    fi = 227.3 * pi / 180
    v = 10.6459
    mu = 0.012150585609624
    dt_values = 10 .^ range(log10(0.00001), log10(0.0000001), length=100)
    x = []
    y = []
    for dt in dt_values
        traj, min_dist_m, steps = run_simulation(fi, v, dt, 4; check_collision=true)
        # Extract the coordinates of last point
        last_point = traj[end]
        push!(x, last_point[1])
        push!(y, last_point[2])
    end
    # Plot the results
    plot(dt_values, x, label="x coordinate", xlabel="h", ylabel="x", legend=false, size=(480, 320), xscale=:log10)
    savefig("x_coordinate_vs_dt.png")
    plot(dt_values, y, label="y coordinate", xlabel="h", ylabel="y", legend=false, size=(480, 320), xscale=:log10)
    savefig("y_coordinate_vs_dt.png")
end

# Function to test how initial conditions affect the trajectory as presented in the report
function test_stability()
    fi = 227.3 * pi / 180
    v = 10.6459
    mu = 0.012150585609624
    dt = 0.000001
    dfi = 0.01 * pi / 180  # Small change in initial angle
    dv = 0.0001  # Small change in initial velocity

    min_dist_matrix = zeros(11, 11)

    for i in -4:4
        for j in -4:4
            traj, min_dist_m, steps = run_simulation(fi + i * dfi, v + j * dv, dt, 4; check_collision=true)
            row = i + 5
            col = j + 5
            if min_dist_m < 1700
                min_dist_matrix[row, col] = 0  # Collision
            elseif steps == 4 / dt
                min_dist_matrix[row, col] = 10e8  # Flyby
            else
                min_dist_matrix[row, col] = min_dist_m  # Normal
            end
        end
    end

    println("Minimum distance matrix:")
    println(min_dist_matrix)
end

compute_free_return_orbits()
#test_time_step_effects()
#test_stability()

