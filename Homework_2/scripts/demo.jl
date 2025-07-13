using Homework_2
using Plots

"""
    compute_force(n)
Computes the force vector from cube 2 to cube 1 using Simpson's 2D method.
Arguments:
- `n`: The number of subdivisions in each dimension for the Simpson's 2D method. Must be even. If not +1 is added to make it even.
Prints the force vector and elapsed time.
"""
function compute_force(n)
    println("Computing the force vector from cube 2 to cube 1 using Simpson's 2D method...")
    t = @elapsed F = simpson2D(n)
    println("Force: ", F)
    println("Elapsed time: ", t, " seconds")
end

"""
    compare_convergence(n_values)
Compares the convergence of Simpson's 2D, 3D, and 6D methods and Monte Carlo integration.
Arguments:
- `n_values`: A range of values for `n` to test the convergence. Must contain only even numbers.
Plots the results and saves the figure as "force_convergence.png".
"""
function compare_convergence(n_values)
    simp2D = [simpson2D(n) for n in n_values]
    println("Simpson's 2D results: ", simp2D)
    simp3D = [simpson3D(n) for n in n_values]
    println("Simpson's 3D results: ", simp3D)
    simp6D = [simpson6D(n) for n in n_values]
    println("Simpson's 6D results: ", simp6D)
    monteCarlo = [monte_carlo_integration(n) for n in n_values]

    # Plot the results
    plot(n_values, simp2D, label="Simpson 2D", xlabel="n", ylabel="Sila")
    plot!(n_values, simp3D, label="Simpson 3D")
    plot!(n_values, simp6D, label="Simpson 6D")
    plot!(n_values, monteCarlo, label="Monte Carlo", linestyle=:dash)
    savefig("force_convergence.png")
end

"""
    compare_time(n_max2D, n_max3D, n_max6D, n_maxMC)
Compares the time taken by Simpson's 2D, 3D, and 6D methods and Monte Carlo integration.
Arguments:
- `n_max2D`: Maximum value of `n` for the 2D method.
- `n_max3D`: Maximum value of `n` for the 3D method.
- `n_max6D`: Maximum value of `n` for the 6D method
- `n_maxMC`: Maximum value of `n` for the Monte Carlo method.
Plots the times and saves the figure as "force_time.png".
"""
function compare_time(n_max2D, n_max3D, n_max6D, n_maxMC)
    println("Running 2D method")
    times2D = [@elapsed simpson2D(n) for n in 2:2:n_max2D]
    println("Running 3D method")
    times3D = [@elapsed simpson3D(n) for n in 2:2:n_max3D]
    println("Running 6D method")
    times6D = [@elapsed simpson6D(n) for n in 2:2:n_max6D]
    println("Running Monte Carlo method")
    timesMC = [@elapsed monte_carlo_integration(n) for n in 2:2:n_maxMC]

    # Plot the times
    plot(2:2:n_max2D, times2D, label="Simpson 2D", xlabel="n", ylabel="Čas [s]", yscale=:log10, xscale=:log10, legend=:topleft)
    plot!(2:2:n_max3D, times3D, label="Simpson 3D")
    plot!(2:2:n_max6D, times6D, label="Simpson 6D")
    plot!(2:2:n_maxMC, timesMC, label="Monte Carlo", linestyle=:dash)
    savefig("force_time.png")
end

#compute_force(200000)
compare_convergence(2:2:20)
compare_time(1000, 300, 20, 1000)