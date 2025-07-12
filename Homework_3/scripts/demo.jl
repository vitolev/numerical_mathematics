using Homework_3
using Plots

function run_simulation(fi, v, dt, t_max; check_collision=false, output_file="trajectory.png")
    mu = 0.012150585609624  # Gravitational parameter for the Earth-Moon system
    r = 6600 / 384400  # Initial orbiting radius of the spacecraft
    X0 = [-mu + r*cos(fi), r*sin(fi), 0.0, v*-sin(fi), v*cos(fi), 0.0]  # Initial state vector

    n_steps = floor(Int, t_max / dt)
    trajectory, min_dist_m = rk4_method(X0, mu, dt, n_steps, check_collision=check_collision)
    println("Minimum distance to the Moon: $min_dist_m km")

    xs = [X[1] for X in trajectory]
    ys = [X[2] for X in trajectory]

    # Plot the trajectory
    plot(xs, ys, legend=false, color=:black, aspect_ratio=:equal)  # Trajectory in black

    function plot_filled_circle(x, y, r; color=:blue)
        θ = range(0, 2π, length=100)
        x_circle = x .+ r * cos.(θ)
        y_circle = y .+ r * sin.(θ)
        plot!([Shape(x_circle, y_circle)], color=color, linewidth=0)
    end

    # Plot Earth and Moon as filled shapes
    plot_filled_circle(-mu, 0.0, 6400/384400; color=:blue)     # Earth (blue)
    plot_filled_circle(1.0 - mu, 0.0, 1700/384400; color=:gray) # Moon (gray)

    savefig(output_file)
end

fi = [227.3, 210, -5, 30]
v = [10.6459, -10.6823, 10.64555, -10.68257]
for i in 1:length(fi)
    run_simulation(fi[i] * pi / 180, v[i], 0.00001, 8; check_collision=true, output_file="trajectory_$i.png")
end