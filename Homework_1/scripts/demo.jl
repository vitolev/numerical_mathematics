# Script to compare the performace of our QR decomposition with the one from LinearAlgebra

using LinearAlgebra
using Statistics
using Homework_1
using Plots
using Random

possible_sizes = [10, 50, 100, 250, 500, 750, 1000]

# Create a random symmetric tridiagonal matrix of size n*n
function create_tridiagonal(n; seed=42)
    Random.seed!(seed)  # Set seed for reproducibility
    
    # Generate random diagonal and sub-diagonal elements
    diag = randn(n)  # Main diagonal
    sub_diag = randn(n - 1)  # Sub-diagonal (one less than main diagonal)
    
    # Create the tridiagonal matrix using the SimTridiag type and the full dense matrix representation for comparison
    tridiagonal_matrix = zeros(n, n)  # Initialize a dense matrix
    for i in 1:n
        tridiagonal_matrix[i, i] = diag[i]  # Main diagonal
        if i > 1
            tridiagonal_matrix[i, i - 1] = sub_diag[i - 1]  # Sub-diagonal
            tridiagonal_matrix[i - 1, i] = sub_diag[i - 1]  # Super-diagonal (same as sub-diagonal)
        end
    end
    return SimTridiag(diag, sub_diag), tridiagonal_matrix
end

# Function to measure the time taken for QR decomposition
function measure_time(n)
    # Create a random symmetric tridiagonal matrix of size n x n
    A, full_matrix = create_tridiagonal(n)
    
    # Measure the time taken for QR decomposition using the custom implementation
    custom_time = @elapsed Q_custom, R_custom = Homework_1.qr(A)

    # Measure the time taken for QR decomposition using LinearAlgebra
    linear_time = @elapsed Q_linear, R_linear = LinearAlgebra.qr(full_matrix)
    
    return custom_time, linear_time
end

"""
Function to compare the performance of custom QR decomposition with LinearAlgebra's QR decomposition
"""
function compare_qr_decompositions()
    means_custom = Float64[]
    uncertainties_custom = Float64[]
    means_linear = Float64[]
    uncertainties_linear = Float64[]

    for n in possible_sizes
        custom_times = Float64[]
        linear_times = Float64[]
        for _ in 1:100 # Repeat the measurement 100 times for each size
            custom_time, linear_time = measure_time(n)
            push!(custom_times, custom_time)
            push!(linear_times, linear_time)
        end
        avg_custom_time = mean(custom_times)
        avg_linear_time = mean(linear_times)
        std_custom_time = std(custom_times)
        std_linear_time = std(linear_times)

        custom_time_un = std_custom_time / sqrt(length(custom_times))
        linear_time_un = std_linear_time / sqrt(length(linear_times))
        
        push!(means_custom, avg_custom_time)
        push!(uncertainties_custom, custom_time_un)
        push!(means_linear, avg_linear_time)
        push!(uncertainties_linear, linear_time_un)

        println("Size: $n")
        println("Custom QR time: $avg_custom_time ± $custom_time_un seconds")
        println("LinearAlgebra QR time: $avg_linear_time ± $linear_time_un seconds")
    end

    # Plotting the results

    gr(size=(600, 350))
    # Plot custom QR times with error bars
    plot(possible_sizes, means_custom; 
        yerror=uncertainties_custom, 
        label="Prilagojen QR", 
        marker=:circle, 
        linewidth=2, 
        linestyle=:solid,
        legend=:topleft)

    # Plot LinearAlgebra QR times with error bars
    plot!(possible_sizes, means_linear; 
        yerror=uncertainties_linear, 
        label="LinearAlgebra QR", 
        marker=:square, 
        linewidth=2, 
        linestyle=:solid)

    xlabel!("Velikost matrike")
    ylabel!("Čas (sekunde)")
    savefig("qr_performance_comparison.png")
end

"""
Function to perform qr iteration on example and plots the convergence
"""
function qr_iteration()
    T = SimTridiag([4.0, 2.0, 1.0], [1.0, 2.0])

    # The same function as in Homework_1.qr, but here we store intermediate results for plotting
    eigenvalue_1 = [T.gd[1]]
    eigenvalue_2 = [T.gd[2]]
    eigenvalue_3 = [T.gd[3]]
    max_iter = 19
    n = length(T.gd)
    V = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]  # Initialize eigenvector matrix as identity
    Q, R = Homework_1.qr(T)
    for _ in 1:max_iter
        T = R * Q
        eigenvalue_1 = push!(eigenvalue_1, T.gd[1])
        eigenvalue_2 = push!(eigenvalue_2, T.gd[2])
        eigenvalue_3 = push!(eigenvalue_3, T.gd[3])
        V = V * Q
        Q, R = Homework_1.qr(T)
    end
    T = R * Q
    V = V * Q
    eigenvalue_1 = push!(eigenvalue_1, T.gd[1])
    eigenvalue_2 = push!(eigenvalue_2, T.gd[2])
    eigenvalue_3 = push!(eigenvalue_3, T.gd[3])

    # Print final eigenvalues
    println("Final eigenvalues after QR iteration:")
    println("Eigenvalue 1: ", T.gd[1])
    println("Eigenvalue 2: ", T.gd[2])
    println("Eigenvalue 3: ", T.gd[3])

    # Plot the convergence of eigenvalues
    gr(size=(600, 350))
    plot(0:max_iter+1, eigenvalue_1, label="Eigenvalue 1", marker=:circle, linewidth=2, legend=:right)
    plot!(0:max_iter+1, eigenvalue_2, label="Eigenvalue 2", marker=:circle, linewidth=2)
    plot!(0:max_iter+1, eigenvalue_3, label="Eigenvalue 3", marker=:circle, linewidth=2)
    xlabel!("Iteration")
    ylabel!("Eigenvalue")
    savefig("qr_iteration_convergence.png")
end

# Run the comparison function
compare_qr_decompositions()
qr_iteration()