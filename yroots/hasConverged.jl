module HasConverged

function hasConverged(coeff::AbstractArray{T}, coeff2::AbstractArray{T}, tol::T) where T
    """Determine whether the high-degree coefficients of a Chebyshev approximation have converged
    to machine epsilon.

    Parameters
    ----------
    coeff : AbstractArray{T}
        Absolute values of chebyshev coefficients of degree n approximation.
    coeff2 : AbstractArray{T}
        Absolute values of chebyshev coefficients of degree 2n+1 approximation.
    tol : T
        Tolerance (distance from zero) used to determine whether the coefficients have converged.
    
    Returns
    -------
    hasConverged : Bool
        True if all the values of coeff and coeff2 are within tol of each other; False otherwise
    """
    coeff3 = copy(coeff2)
    # Subtract off coeff from coeff2 elementwise and ensure all elements are then less than tol
    #TODO: for optimiaztion, it may be faster to iterate through each individually and break early if it's not within the tolerance level rather than calling maximum 
	coeff3[CartesianIndices(coeff)] .-= coeff 
    return maximum(abs.(coeff3)) < tol
end

end