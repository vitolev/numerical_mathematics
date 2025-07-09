using Homework_3
using Plots

# Initial conditions
mu = 0.012150585609624
fi = 235 * pi / 180  # Initial angle in radians
v = 10.6285   # Initial velocity
r = 6600 / 384400  # Initial orbiting radius of the spacecraft
X0 = [-mu + r*cos(fi), r*sin(fi), 0.0, v*-sin(fi), v*cos(fi), 0.0]  # Initial state vector

dt = 0.00001  # Time step
t_max = 2.0
n_steps = floor(Int, t_max / dt)
trajectory = rk4_method(X0, mu, dt, n_steps, check_collision=true)

xs = [X[1] for X in trajectory]
ys = [X[2] for X in trajectory]

# Plot the trajectory
plot(xs, ys, legend=false, color=:black)  # Trajectory in black

function plot_filled_circle(x, y, r; color=:blue)
    θ = range(0, 2π, length=100)
    x_circle = x .+ r * cos.(θ)
    y_circle = y .+ r * sin.(θ)
    # Use `shape` to draw a filled circle
    plot!([Shape(x_circle, y_circle)], color=color, linewidth=0)
end

# Plot Earth and Moon as filled shapes
plot_filled_circle(-mu, 0.0, 6400/384400; color=:blue)     # Earth (blue)
plot_filled_circle(1.0 - mu, 0.0, 1700/384400; color=:gray) # Moon (gray)

savefig("trajectory.png")
