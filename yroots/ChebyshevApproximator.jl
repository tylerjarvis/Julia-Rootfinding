function transformpoints(x,a,b)
    """Transforms points from the interval [-1, 1] to the interval [a, b].

    Parameters
    ----------
    x : array
        The points to be tranformed.
    a : array
        The lower bounds on the interval.
    b : array
        The upper bounds on the interval.

    Returns
    -------
    transformed_pts : array
        The transformed points.
    """

    return ((b-a).*x .+(b+a))/2
end

# function hasconverged(coeff,coeff2,tol)
#     """Determine whether the high-degree coefficients of a Chebyshev approximation have converged
#     to machine epsilon.

#     Parameters
#     ----------
#     coeff : array
#         Absolute values of chebyshev coefficients of degree n approximation.
#     coeff2 : array
#         Absolute values of chebyshev coefficients of degree 2n+1 approximation.
#     tol : float
#         Tolerance (distance from zero) used to determine wheher the coefficients have converged.

#     Returns
#     -------
#     hasConverged : bool
#         True if all the values of coeff and coeff2 are within tol of each other; False otherwise
#     """
# end

function getfinal_degree(coeff,tol)
    """Finalize the degree of Chebyshev approximation to use along one particular dimension.

    This function is called after the coefficients have started converging at degree n. A degree
    2n+1 approximation is passed in. Assuming that the coefficients have fully converged by degree 
    3n/2, the cutoff epsVal is calculated as twice the max coefficient of degree at least 3n/2.
    The final degree is then set as the largest coefficient with magnitude greater than epsVal.
    
    The rate of convergence is calculated assuming that the coefficients converge geometrically
    starting from the largest coefficient until machine epsilon is reached. This is a lower bound, as
    in practice, the coefficients usually slowly decrease at first but drop off fast at the end.

    Parameters
    ----------
    coeff : array
        Absolute values of chebyshev coefficients.
    
    Returns
    -------
    degree : int
        The numerical degree of the approximation
    epsVal : float
        The epsilon value to which the coefficients have converged
    rho : float
        The geometric rate of convergence of the coefficients
    """
    # # Set the final degree to the position of the last coefficient greater than convergence value
    # convergedDeg = int(3 * (len(coeff) - 1) / 4) # Assume convergence at degree 3n/2.
    # epsVal = 2*np.max(coeff[convergedDeg:]) # Set epsVal to 2x the largest coefficient past degree 3n/2
    # nonZeroCoeffs = np.where(coeff > epsVal)[0]
    # degree = 1 if len(nonZeroCoeffs) == 0 else max(1, nonZeroCoeffs[-1])

    # # Set degree to 0 for constant functions (all coefficients but first are less than tol)
    # if np.all(coeff[1:] < tol):
    #     degree = 0
    
    # # Calculate the rate of convergence
    # maxSpot = np.argmax(coeff)
    # if epsVal == 0: #Avoid divide by 0. epsVal shouldn't be able to shrink by more than 1e-24 cause floating point.
    #      epsVal = coeff[maxSpot] * 1e-24
    # rho = (coeff[maxSpot]/epsVal)**(1/((degree - maxSpot) + 1)) 
    # return degree, epsVal, rho
end

function startedconverging(coefflist,tol)
    """Determine whether the high-degree coefficients of a given Chebyshev approximation are near 0.

    Parameters
    ----------
    coeffList : array
        Absolute values of chebyshev coefficients.
    tol : float
        Tolerance (distance from zero) used to determine whether coeffList has started converging.
    
    Returns
    -------
    startedConverging : bool
        True if the last 5 coefficients of coeffList are less than tol; False otherwise
    """
    return all(x -> x < tol, coefflist[end-4:end])
end