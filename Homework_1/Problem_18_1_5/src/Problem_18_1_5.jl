module Problem_18_1_5

export SimTridiag

"""
    SimTridiag(gd, sd)

Data type for a symmetric tridiagonal matrix with the given main diagonal `gd` and sub-diagonal `sd`.
"""
struct SimTridiag
    gd  # main diagonal
    sd  # sub-diagonal (both lower and upper, since it's symmetric)
end

end # module Problem_18_1_5