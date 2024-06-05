#import Pkg; Pkg.add("IterTools")
using IterTools

function getApproxError(degs, epsilons, rhos)
    """
    Computes an upper bound for the error of the Chebyshev approximation.

    Using the epsilon values and rates of geometric convergence calculated in getChebyshevDegrees,
    calculates the infinite sum of the coefficients past those used in the approximation.
    
    Parameters
    ----------
    degs :
        The degrees used in each dimension of the approximation.
    epsilons : 
        The values to which the approximation converged in each dimension.
    rhos : 
        The calculated rate of convergence in each dimension.
    
    Returns
    -------
    approxError : Float64
        An upper bound on the approximation error
    """

    approxError = 0.0
    
    # Create a partition of coefficients where idxs[i]=1 represents coefficients being greater than
    # degs[i] in dimension i and idxs[i]=0 represents coefficients being less than [i] in dimension i.
    for idxs in Iterators.product(Iterators.repeated((0,1), length(degs))...)
        # Skip the set of all 0's, corresponding to the terms actually included in the approximation.
        if sum(idxs) == 0
            continue
        end
        
        s = 1.0
        thisEps = 0.0
        
        for (i, used) in enumerate(idxs)
            if Bool(used)
                # multiply by infinite sum of coeffs past the degree at which the approx stops in dim i
                #1/rhos[i] is the rate, so this is (1/rhos[i]) / (1 - 1/rhos[i]) = 1/(rhos[i]-1)
                s /= (rhos[i] - 1)
                # The points in this section are going to be < max(epsilons[i] that contribute to it)
                thisEps = max(thisEps, epsilons[i])
            else
                # multiply by the number of coefficients in the approximation along dim i
                s *= (degs[i] + 1)
            end
        end
        
        # Append to the error
        approxError += s * thisEps
    end
    
    return approxError
end

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

    # Set the final degree to the position of the last coefficient greater than convergence value
    converged_deg = Int64(div((3 * (length(coeff) - 1) / 4),1)) + 1 # Assume convergence at degree 3n/2.
    epsval = 2*maximum(coeff[converged_deg:end]) # Set epsVal to 2x the largest coefficient past degree 3n/2
    nonzero_coeffs_index = [i for i in 1:length(coeff) if coeff[i]>epsval]
    if isempty(nonzero_coeffs_index) 
        degree = 1
    else
        degree = max(1,nonzero_coeffs_index[end]-1)
    end

    # Set degree to 0 for constant functions (all coefficients but first are less than tol)
    if all(x -> x < tol, coeff[2:end])
        degree = 0
    end
    
    # Calculate the rate of convergence
    maxspot = argmax(coeff)
    if epsval == 0 #Avoid divide by 0. epsVal shouldn't be able to shrink by more than 1e-24 cause floating point.
         epsval = coeff[maxspot] * 1e-24
    end

    rho = (coeff[maxspot]/epsval)^(1/(degree - (maxspot[2]) + 2))
    return degree, epsval, rho
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

function check_constant_in_dimension(f,a,b,currdim,tol)
    """Check to see if the output of f is not dependent on the input coordinate of a dimension.
    
    Uses predetermined random numbers to find a point x in the interval where f(x) != 0 and checks
    whether the value of f changes as the dimension currDim coordinate of x changes. Repeats twice.

    Parameters
    ----------
    f : function
        The function being evaluated.
    a : array
        The lower bound on the interval.
    b : array
        The upper bound on the interval.
    currDim : int
        The dimension being examined.
    
    Returns
    -------
    is_constant : bool
        Whether the dimension is constant in dimension currDim. Returns False if the test is
        indeterminate or f is seen to vary with different values of x[dim]. Returns True otherwise.
    """
    dim = length(a)
    currdim = currdim + 1
    # First test point x1
    x1 = transformpoints([0.8984743990614998^(val) for val in 1:dim]',a,b)
    eval1 = f(x1...)
    if isapprox(eval1,0,rtol=tol)
        return false
    end
    # Test how changing x_1[dim] changes the value of f for several values         
    for val in transformpoints([-0.7996847717584993 0.18546110255464776 -0.13975937255055182 0. 1. -1.]',a[currdim],b[currdim])
        x1[currdim] = val
        eval2 = f(x1...)
        if !isapprox(eval1,eval2,rtol=tol) # Corresponding points gave different values for f(x)
            return false
        end
    end

    # Second test point x_2
    x2 = transformpoints([(-0.2598647169391334*(val)/(dim))^2 for val in 1:dim]',a,b)
    eval1 = f(x2...)
    if isapprox(eval1,0,rtol=tol) # Make sure f(x_2) != 0 (unlikely)
        return false
    end

    for val in transformpoints([-0.17223860129797386,0.10828286380141305,-0.5333148248321931,0.46471703497219596]',a[currdim],b[currdim])
        x2[currdim] = val
        eval2 = f(x2...)
        if !isapprox(eval1,eval2,rtol=tol)
            return false # Corresponding points gave different values for f(x)
        end
    end
    # Both test points had not zeros of f and had no variance along dimension currDim.
    return true
end

function hasConverged(coeff, coeff2, tol)
    """Determine whether the high-degree coefficients of a Chebyshev approximation have converged
    to machine epsilon.

    Parameters
    ----------
    coeff : array
        Absolute values of chebyshev coefficients of degree n approximation.
    coeff2 : array
        Absolute values of chebyshev coefficients of degree 2n+1 approximation.
    tol : float
        Tolerance (distance from zero) used to determine whether the coefficients have converged.
    
    Returns
    -------
    hasConverged : Bool
        True if all the values of coeff and coeff2 are within tol of each other; False otherwise
    """
    coeff3 = copy(coeff2)
	coeff3[CartesianIndices(coeff)] .-= coeff 
    return maximum(abs.(coeff3)) < tol
end

function create_meshgrid(point_arrays...)
    num_arrays = length(point_arrays)
    matrix_lengths = [length(point_array) for point_array in point_arrays]
    outputs = []

    if num_arrays == 1
        return point_arrays[1]
    end

    for i in 1:num_arrays
        arr = []
        point_array = point_arrays[i]
        if i == 1
            repeat = prod(matrix_lengths[2:end])
            for item in point_array
                for j in 1:repeat
                push!(arr,item)
                end
            end
            push!(outputs,reshape(arr,Tuple(matrix_lengths)))
        elseif i == num_arrays
            for j in 1:prod(matrix_lengths[1:i-1])
                for item in point_array
                    push!(arr,item)
                end
            end
            push!(outputs,reshape(arr,Tuple(matrix_lengths)))
        else
            repeat = prod(matrix_lengths[i+1:end])
            for j in 1:product(matrix_lengths[1:i-1])
                for item in point_array
                    for k in 1:repeat
                    push!(arr,item)
                    end
                end
            end
        end
        push!(outputs,reshape(arr,Tuple(matrix_lengths)))
    end
    return outputs
end

function create_meshgrid2(arrays...)
    """ Takes arguments x1,x2,...,xn, each xi being a row vector. Does meshgrid. 
    Output is in the format [X,Y,Z,...], where each element in the list is a matrix.
    Example: create_meshgrid2([1 2],[3 4]) -> [[1;1;;2;2],[3;4;;3;4]].
        Note: This would output as [[1 2;1 2],[3 3;4 4]], which looks wrong, but for our 
        purposes it will be easier to think of the matricies in the first format instead of the printing one.
        This way, accessing elements and slices is exactly like in Python but reversed.
    """
    dims = []
    for array in arrays
        push!(dims,length(array))
    end
    finals = []
    for iter in range(1,length(arrays))
        # Get product of array sizes before and after the current array
        reps = 1
        full_reps=1
        try 
            reps = prod(dims[iter+1:end])
        catch end
        try 
            full_reps = prod(dims[1:iter-1])
        catch end
        # Repeat the current array into a matrix with height and width determined by the sizes before and after
        newArray = arrays[iter]
        endArray = repeat((arrays[iter])',full_reps,reps)
        # Reshape it to the meshgrid array and push it onto our final list
        push!(finals,reshape(reduce(vcat,endArray';init=[0])[2:end],Tuple(reverse(dims)))) 
    end
    return finals
end


#f = (x,y,z) -> x+y+z
#A = [1 2 3;4 5 6;1 1 1]
#X = mapslices(x -> f(x...),A;dims=2)
#println(X)
#X = f(eachrow(A))
#println(X)
#X = [f(A[i,:]...) for i in range(1,2)]
#println(X)
#v = [0 1 2; 3 4 5;;; 6 7 8; 9 10 11;;; 12 13 14; 15 16 17;;; 18 19 20; 21 22 23]
#w =  [0;3;;1;4;;2;5;;;6;7;;8;9;;10;11;;;;0;3;;1;4;;2;5;;;6;7;;8;9;;10;11]
#println(w)
#println(w[1,:,2,:])
#println(size(w))
#println(v[5,3,2])
#x = [1;2;;3;4]
#display(v)
#println(reverse([1 2 3]))
x = [1 2]
y = [2 3 4]
z = [3 1 6 7]
#println(create_meshgrid2(x,y,z))

