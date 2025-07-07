module Homework_2

using Random

export force_integrand6D, simpson_weight, simpson_weight6D, simpson6D, monte_carlo_integration, force_integrand3D, simpson_weights_1D, volume_overlap, simpson3D, simpson2D, custom_integrand

############################################################
############# Naive 6D Integration of Force ################
############################################################

# This section contains functions for calculating gravitational force between two cubes in 3D space using both Simpson's rule and Monte Carlo integration.
# It computes the 6D integral, which is quickly turned out to be infeasible for large n, so we introduce better approaches later.
# The section is kept for educational purposes ann to display the progress of the homework and for comparison with the later approaches.

"""
    force_integrand6D(x1, y1, z1, x2, y2, z2)

Calculate the integrand for the force between two points in 3D space.
Arguments:
- `x1`, `y1`, `z1`: Coordinates of the first point.
- `x2`, `y2`, `z2`: Coordinates of the second point.
Returns:
- An x-component of the integrand, as the resulting force will be only in the x-direction.
"""
function force_integrand6D(x1, y1, z1, x2, y2, z2)
    dist_sq = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
    return (x1 - x2) / dist_sq^(3/2)
end

"""
    simpson_weight(i, n)

Calculate the Simpson's rule weight in 1D for the i-th index in an n-point integration.
Returns:
- 1 for the endpoints (i = 0 or n)
- 2 for even indices (i = 2, 4, ...)
- 4 for odd indices (i = 1, 3, ...)
"""
function simpson_weight(i, n)
    if i == 0 || i == n
        return 1    # Endpoints have weight 1
    elseif iseven(i)
        return 2    # Even indices have weight 2
    else
        return 4    # Odd indices have weight 4
    end
end

"""
    simpson_weight6D(indices, n)

Calculate the weight for a point in 6D integral using Simpson's rule.
The weight is the product of the weights for each dimension.
Arguments:
- `indices`: A vector of indices for each dimension (length 6).
- `n`: The number of points in each dimension (assumed to be the same for all dimensions).
Returns:
- The weight for a point in the 6D integral.
"""
function simpson_weight6D(indices, n)
    # Weight is w_1 * w_2 * ... * w_6, where w_i is the weight for the i-th dimension
    w = 1
    for i in indices
        w *= simpson_weight(i, n)
    end
    return w
end

"""
    simpson6D(n)

Calculate the 6D integral of the force between two cubes in 3D space using Simpson's rule.
Arguments:
- `n`: The number of subdivisions in each dimension for the Simpson's rule integration. Must be even. If not, it will be incremented by 1 to ensure evenness.
Returns:
- An x-component of the total force vector from cube 2 to cube 1, as the resulting force will be only in the x-direction.
"""
function simpson6D(n)
    # Ensure n is even. This is necessary for Simpson's rule as the weights are 1, 4, 2, 4, ..., 2, 4, 1.
    # Hence, we need odd number of function evaluations, therefore we need even number of subdivisions, so even n.
    if n % 2 != 0
        n += 1
    end
    m = n / 2  # midpoint needed for later

    # Cube 1 bounds
    # X:
    x_1_min = 0.0
    x_1_max = 1.0
    # Y:
    y_1_min = 0.0
    y_1_max = 1.0
    # Z:
    z_1_min = 0.0
    z_1_max = 1.0
    # Cube 2 bounds
    # X:
    x_2_min = 2.0
    x_2_max = 3.0
    # Y:
    y_2_min = 0.0
    y_2_max = 1.0
    # Z:
    z_2_min = 0.0
    z_2_max = 1.0

    h = 1.0 / n
    factor = (h / 3)^6
    F_total = 0.0

    # We iterate through entire first cube and 1/8 of the second (specifically along x-axis and triangle in y-z plane that is formed by center, vertex and midpoint of an edge).
    # This is because the force is symmetric, so we can calculate only 1/8 of the integral and multiply by 8.
    for i1 in 0:n, i2 in 0:n, i3 in 0:n, i4 in 0:n, i5 in 0:m, i6 in 0:i5
        x1 = x_1_min + i1 * h
        y1 = y_1_min + i2 * h
        z1 = z_1_min + i3 * h
        x2 = x_2_min + i4 * h
        y2 = y_2_min + i5 * h
        z2 = z_2_min + i6 * h

        if x1 == x2 && y1 == y2 && z1 == z2
            continue  # Skip if points are the same to avoid division by zero
        end

        W = simpson_weight6D((i1, i2, i3, i4, i5, i6), n)
        if i5 == i6 || (i5 == m)
            if i6 == m
                W = W / 8   # If we are at the exact middle in y-z plane there are 8 triangles in that point, so we divide by 8
            else
                W = W / 2   # If not we are either at diagonal edge or at the vertical edge, so we divide by 2, as there are 2 triangles that share this edge
            end
        end

        f = force_integrand6D(x1, y1, z1, x2, y2, z2)
        F_total += W * f
    end

    return 8 * factor * F_total # Multiply by 8 to account for the symmetry of the force vector
end

"""
    monte_carlo_integration(n; seed=42)

Perform Monte Carlo integration to estimate the integral of the force between two cubes in 3D space.
Arguments:
- `n`: The number of random points to sample in the integration.
- `seed`: (Optional) Seed for the random number generator. Defaults to 42.
Returns:
- An x-component of the total force vector from cube 2 to cube 1, as the resulting force will be only in the x-direction.
"""
function monte_carlo_integration(n; seed=42)
    rng = MersenneTwister(seed)  # Use fixed seed for reproducibility
    F_total = 0.0

    for _ in 1:n
        x1 = rand(rng)  # Random x1 in [0, 1]
        y1 = rand(rng)  # Random y1 in [0, 1]
        z1 = rand(rng)  # Random z1 in [0, 1]
        x2 = rand(rng) + 2.0  # Random x2 in [2, 3]
        y2 = rand(rng)  # Random y2 in [0, 1]
        z2 = rand(rng)  # Random z2 in [0, 1]

        f = force_integrand6D(x1, y1, z1, x2, y2, z2)
        F_total += f
    end

    return F_total / n
end

#############################################################
############ Smart integration over 3D space ################
#############################################################

# This section computes the force between two cubes in 3D space using a more efficient approach. Instead of integrating over the entire 6D space (3D for each cube),
# we integrate just once over the 3D space of relative positions of the cubes. This is significant improvement over the naive 6D integration.

"""
    force_integrand3D(x, y, z)

Calculate the integrand for the force between two points in 3D space.
Arguments: 
- `x`, `y`, `z`: Coordinates of the relative position vector from cube 2 to cube 1.
Returns:
- An x-component of the integrand, as the resulting force will be only in the x-direction.
"""
function force_integrand3D(x, y, z)
    dist_sq = x^2 + y^2 + z^2
    return x / dist_sq^(3/2)
end

"""
    simpson_weights_1D(n)

Generate a vector of Simpson's rule weights for 1D integration.
Arguments:
- `n`: The number of subdivisions in the 1D integration.
Returns:
- A vector of weights for each point in the 1D integration (1, 4, 2, 4, ..., 2, 4, 1).
"""
function simpson_weights_1D(n)
    weights = fill(2.0, n + 1)
    weights[1] = 1.0
    weights[end] = 1.0
    for i in 2:2:n
        weights[i] = 4.0
    end
    return weights
end

"""
    volume_overlap(x, y, z)

Calculate the volume of overlap between two cubes in 3D space based on the relative position vector.
Arguments:
- `x`, `y`, `z`: Coordinates of the relative position vector from cube 2 to cube 1. x is expected to be in [-3, -1] and y, z to be in [-1, 1].
Returns:
- The volume of overlap between the two cubes, which is a function of the relative position vector.
"""
function volume_overlap(x, y, z)
    return (min(1.0, 3.0+x) - max(0.0, 2.0+x))*(1-abs(y))*(1-abs(z))
end

"""
    simpson3D(n)

Calculate the 3D integral of the force between two cubes in 3D space using Simpson's rule, based on the relative position vector.
Arguments:
- `n`: The number of subdivisions in each dimension for the Simpson's rule integration. Must be even. If not, it will be incremented by 1 to ensure evenness.
Returns:
- An x-component of the total force vector from cube 2 to cube 1, as the resulting force will be only in the x-direction.
"""
function simpson3D(n)
    # Ensure n is even. This is necessary for Simpson's rule as the weights are 1, 4, 2, 4, ..., 2, 4, 1.
    # Hence, we need odd number of function evaluations, therefore we need even number of subdivisions, so even n.
    if n % 2 != 0
        n += 1
    end

    h = 2.0 / n # x range is [-3, -1], y and z ranges are [-1, 1], so difference is 2 for all dimensions
    factor = (h / 3)^3
    F_total = 0.0
    w = simpson_weights_1D(n)

    # Due to symmetry in y and z, we can integrate only over [-1, 0] and then multiply by 2 for y and z.
    for i1 in 0:n, i2 in 0:(n÷2), i3 in 0:(n÷2)
        x = -3.0 + i1 * h
        y = -1.0 + i2 * h
        z = -1.0 + i3 * h

        W = w[i1 + 1] * w[i2 + 1] * w[i3 + 1]
        f = force_integrand3D(x, y, z)
        V = volume_overlap(x, y, z)

        symmetric_factor = 1.0
        symmetric_factor *= (i2 == n÷2) ? 1.0 : 2.0 
        symmetric_factor *= (i3 == n÷2) ? 1.0 : 2.0
        F_total += W * f * V * symmetric_factor
    end

    return factor * F_total # Multiply by the factor to account for the volume and the step size
end

##################################################################
###### Smart integration with partial analytical solution ########
##################################################################

# This section computes the force between two cubes in 3D space using the same approach as the previous section, 
# but first we analytical integrate over z dimenstion, to reduce the integral to 2D integral over x and y dimensions.

"""
    simpson2D(n)

Calculate the 2D integral of the force between two cubes in 3D space using Simpson's rule, based on the relative position vector and
accounting for the analytical integration over the z dimension.
Arguments:
- `n`: The number of subdivisions in each dimension for the Simpson's rule integration. Must be even. If not, it will be incremented by 1 to ensure evenness.
Returns:
- An x-component of the total force vector from cube 2 to cube 1, as the resulting force will be only in the x-direction.
"""
function simpson2D(n)
    # Ensure n is even. This is necessary for Simpson's rule as the weights are 1, 4, 2, 4, ..., 2, 4, 1.
    # Hence, we need odd number of function evaluations, therefore we need even number of subdivisions, so even n.
    if n % 2 != 0
        n += 1
    end

    h = 2.0 / n # x range is [-3, -1] and y is [-1, 1], so difference is 2 for both dimensions
    factor = (h / 3)^2
    F_total = 0.0
    w = simpson_weights_1D(n)   # Precompute Simpson's weights for 1D integration

    # Due to symmetry in y, we can integrate only over [-1, 0] and then multiply by 2 for y.
    for i1 in 0:n, i2 in 0:(n÷2)
        x = -3.0 + i1 * h
        y = -1.0 + i2 * h

        W = w[i1 + 1] * w[i2 + 1]   # Access the precomputed weights
        f = custom_integrand(x, y)  # Use the custom integrand that accounts for analytical integration over z

        symmetric_factor = (i2 == n÷2) ? 1.0 : 2.0  # If we are at the midpoint in y, we only count it once, otherwise we count it twice
        F_total += W * f * symmetric_factor
    end
    return factor * F_total # Multiply by the factor to account for the volume and the step size
end

"""
    custom_integrand(x, y)

Calculate the integrand for the force between two points in 3D space, accounting for the analytical integration over z dimension.
Arguments:
- `x`, `y`: Coordinates of the relative position vector from cube 2 to cube 1.
Returns:
- An x-component of the integrand, as the resulting force will be only in the x-direction.
"""
function custom_integrand(x, y)
    sq_sum = x^2 + y^2
    volume = (min(1.0, 3.0+x) - max(0.0, 2.0+x)) * (1 + y)
    factor = 1 / (sq_sum*(sq_sum + 1)^(1/2)) + (sq_sum + 1)^(-1/2) - sq_sum^(-1/2)
    return 2*x * volume * factor
end


end # module Homework_2
