# Script to compare the performace of our QR decomposition with the one from LinearAlgebra

using LinearAlgebra
using Statistics
using Homework_1
using Plots

possible_sizes = [10, 50, 100, 250, 500, 750, 1000]

# Create a random symmetric tridiagonal matrix of size n*n
function create_tridiagonal(n)
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

# Run the comparison function
compare_qr_decompositions()