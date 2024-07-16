# TODO: import from a library like this one instead of crowding our sourcecode with pre-written code https://github.com/JeffreySarnoff/ErrorfreeArithmetic.jl/blob/main/src/sum.jl
function twoSum(a,b)
    """Returns x,y such that a+b=x+y exactly, and a+b=x in floating point."""
    x = a+b
    z = x-a
    y = (a-(x-z)) + (b-z)
    return x, y
end

# TODO: import from a library instead
function split(a)
    """Returns x,y such that a = x+y exactly and a = x in floating point."""
    c = (2^27 + 1) * a
    x = c-(c-a)
    y = a-x
    return x,y
end

# TODO: import from a library instead
function twoProd(a,b)
    """Returns x,y such that a*b=x+y exactly and a*b=x in floating point."""
    x = a*b
    a1,a2 = split(a)
    b1,b2 = split(b)
    y=a2*b2-(((x-a1*b1)-a2*b1)-a1*b2)
    return x,y
end

# TODO: import from a library instead
function TwoProdWithSplit(a,b,a1,a2)
    """Returns x,y such that a*b = x+y exactly and a*b = x in floating point but with a already split."""
    x = a*b
    b1,b2 = split(b)
    y=a2*b2-(((x-a1*b1)-a2*b1)-a1*b2)
    return x,y
end 

function getLinearTerms(M)
    """Gets the linear terms of the Chebyshev coefficient tensor M.

    Uses the fact that the linear terms are located at
    M[(0,0, ... ,0,1)]
    M[(0,0, ... ,1,0)]
    ...
    M[(0,1, ... ,0,0)]
    M[(1,0, ... ,0,0)]
    which are indexes
    1, M.shape[-1], M.shape[-1]*M.shape[-2], ... when looking at M.ravel().

    Parameters
    ----------
    M : array
        The coefficient array to get the linear terms from

    Returns
    -------
    A: array
        An array with the linear terms of M
    """
    A = []
    spot = 1
    
    for i in size(M)
        push!(A, (i) == 1 ? 0 : reshape(M,(1,length(M)))[spot+1])
        spot *= (i)
    end

    return reverse(A) # Return linear terms in dimension order.
end

function linearCheck1(totalErrs,A,consts)
    """Takes A, the linear terms of each function approximation, and makes any possible reduction 
        in the interval based on the totalErrs.


    Parameters
    ----------
    totalErrs : array
        gives bounds for the function using error in our approximation and coefficients
    A : array 
        each row represents a function with the linear coefficients of each dimension as the columns
    consts : array
        constant terms for each function

    Returns
    -------
    a : array
        lower bound
    b : array
        lower bound
        
    """
    dim = length(A[1,:])
    a = -ones(dim) * Inf
    b = ones(dim) * Inf
    for row in 1:dim
        for col in 1:dim
            if A[col,row] != 0 #Don't bother running the check if the linear term is too small.
                v1 = totalErrs[row] / abs(A[col,row]) - 1
                v2 = 2 * consts[row] / A[col,row]
                if v2 >= 0
                    c = -v1
                    d = v1-v2
                else
                    c = -v2-v1
                    d = v1
                end
                a[col] = max(a[col], c)
                b[col] = min(b[col], d)
            end
        end
    end
    return a, b
end

function reduceSolvedDim(Ms, errors, trackedInterval, dim)
    # dim is python dim (starting from 0)

    val = (trackedInterval.interval[1,dim+1] + trackedInterval.interval[2,dim+1]) / 2
    # Get the matrix of the linear terms
    A = []
    for M in Ms
        if isempty(A)
            A = getLinearTerms(M)
        else 
            A = hcat(A,getLinearTerms(M))
        end
    end
    # Figure out which linear approximation is flattest in the dimension we want to reduce
    dot_prods = (A./sqrt.(sum(A.^2, dims = 1)))[dim+1,:]
    func_idx = argmax(dot_prods)
    # Remove that polynomial from Ms and errors
    deleteat!(Ms,func_idx)
    new_errors = Base.copy(errors)
    deleteat!(new_errors,func_idx)

    # Evaluate other polynomials on the solved dimension
    # Ms are already scaled, so we just want to evaluate the T_i(x_dim)'s at x_dim = 0
    final_Ms = []
    for M in Ms
        # Get the dimensions of each M
        degs = reverse(size(M))
        total_dim = length(degs)
        x = zeros(degs[dim+1])
        x[1:2:end] = [(-1)^i for i in collect(0:(degs[dim+1]+1)//2-1)]
        new_M = reshape(mapslices(v->sum(x.*v),M,dims=total_dim-dim),reverse(degs[1:end .!= dim+1]))
        push!(final_Ms,new_M)
    end

    # Remove the point dimension from the tracked interval
    trackedInterval.ndim -= 1
    trackedInterval.interval = trackedInterval.interval[:,1:end .!=dim+1]
    # trackedInterval.interval = np.delete(trackedInterval.interval,dim,axis=0)
    push!(trackedInterval.reducedDims,dim)
    push!(trackedInterval.solvedVals,val)

    return final_Ms,new_errors,trackedInterval
end

function transformChebInPlace1D1D(coeffs,alpha,beta)
    """
    Does the 1d coeff array case for transformChebInPlace1D
    """
    coeffs_shape = size(coeffs)
    last_dim_length = coeffs_shape[end]
    transformedCoeffs = zeros(coeffs_shape)
    # Initialize three arrays to represent subsequent columns of the transformation matrix.
    arr1 = zeros(last_dim_length)
    arr2 = zeros(last_dim_length)
    arr3 = zeros(last_dim_length)

    #The first column of the transformation matrix C. Since T_0(alpha*x + beta) = T_0(x) = 1 has 1 in the top entry and 0's elsewhere.
    arr1[1] = 1.

    transformedCoeffs[1] = coeffs[1] # arr1[0] * coeffs[0] (matrix multiplication step)
    #The second column of C. Note that T_1(alpha*x + beta) = alpha*T_1(x) + beta*T_0(x).
    arr2[1] = beta
    arr2[2] = alpha
    transformedCoeffs[1] += beta * coeffs[2] # arr2[0] * coeffs[1] (matrix muliplication)
    transformedCoeffs[2] += alpha * coeffs[2] # arr2[1] * coeffs[1] (matrix multiplication)
    maxRow = 2
    for col in 2:last_dim_length-1 # For each column, calculate each entry and do matrix mult
        thisCoeff = coeffs[col+1] # the row of coeffs corresponding to the column col of C (for matrix mult)
        # The first entry
        arr3[1] = -arr1[1] + alpha.*arr2[2] + 2*beta*arr2[1]
        transformedCoeffs[1] += thisCoeff * arr3[1]
        # The second entry
        if maxRow > 2
            arr3[2] = -arr1[2] + alpha*(2*arr2[1] + arr2[3]) + 2*beta*arr2[2]
            transformedCoeffs[2] += thisCoeff * arr3[2]
        end

        # All middle entries
        for i in 2:maxRow-2
            arr3[i+1] = -arr1[i+1] + alpha*(arr2[i] + arr2[i+2]) + 2*beta.*arr2[i+1]
            transformedCoeffs[i+1] += thisCoeff * arr3[i+1]
        end

        # The second to last entry
        i = maxRow -1
        arr3[i+1] = -arr1[i+1] + (i == 1 ? 2 : 1)*alpha*(arr2[i]) + 2*beta*arr2[i+1]
        transformedCoeffs[i+1] += thisCoeff * arr3[i+1]
        #The last entry
        finalVal = alpha*arr2[i+1]
        # This final entry is typically very small. If it is essentially machine epsilon,
        # zero it out to save calculations.
        if abs(finalVal) > 1e-16 #TODO: Justify this val!
            arr3[maxRow+1] = finalVal
            transformedCoeffs[maxRow+1] += thisCoeff * finalVal
            maxRow += 1 # Next column will have one more entry than the current column.
        end

        # Save the values of arr2 and arr3 to arr1 and arr2 to get ready for calculating the next column.
        arr = arr1
        arr1 = arr2
        arr2 = arr3
        arr3 = arr
    end
    return transformedCoeffs[1:maxRow]
end

function transformChebInPlace1D(coeffs,alpha,beta)
    """Applies the transformation alpha*x + beta to one dimension of a Chebyshev approximation.

    Recursively finds each column of the transformation matrix C from the previous two columns
    and then performs entrywise matrix multiplication for each entry of the column, thus enabling
    the transformation to occur while only retaining three columns of C in memory at a time.

    Parameters
    ----------
    coeffs : array
        The coefficient array
    alpha : double
        The scaler of the transformation
    beta : double
        The shifting of the transformation

    Returns
    -------
    transformedCoeffs : array
        The new coefficient array following the transformation
    """
    coeffs_shape = size(coeffs)
    if length(coeffs_shape) == 1
        return transformChebInPlace1D1D(coeffs,alpha,beta)
    end
    last_dim_length = coeffs_shape[end]
    transformedCoeffs = zeros(coeffs_shape)
    # Initialize three arrays to represent subsequent columns of the transformation matrix.
    arr1 = zeros(last_dim_length)
    arr2 = zeros(last_dim_length)
    arr3 = zeros(last_dim_length)

    #The first column of the transformation matrix C. Since T_0(alpha*x + beta) = T_0(x) = 1 has 1 in the top entry and 0's elsewhere.
    arr1[1] = 1.

    # Get the correct number of colons for indexing transformedCoeffs
    idxs = []
    dims = length(coeffs_shape)
    for i in 1:dims-1
        push!(idxs,:)
    end

    transformedCoeffs[idxs...,1] = coeffs[idxs...,1] # arr1[0] * coeffs[0] (matrix multiplication step)
    #The second column of C. Note that T_1(alpha*x + beta) = alpha*T_1(x) + beta*T_0(x).
    arr2[1] = beta
    arr2[2] = alpha
    transformedCoeffs[idxs...,1] .+= beta .* coeffs[idxs...,2] # arr2[0] * coeffs[1] (matrix muliplication)
    transformedCoeffs[idxs...,2] .+= alpha .* coeffs[idxs...,2] # arr2[1] * coeffs[1] (matrix multiplication)
    maxRow = 2
    for col in 2:last_dim_length-1 # For each column, calculate each entry and do matrix mult
        thisCoeff = coeffs[idxs...,col+1] # the row of coeffs corresponding to the column col of C (for matrix mult)
        # The first entry
        arr3[1] = -arr1[1] + alpha.*arr2[2] + 2*beta.*arr2[1]
        transformedCoeffs[idxs...,1] += thisCoeff .* arr3[1]
        # The second entry
        if maxRow > 2
            arr3[2] = -arr1[2] + alpha*(2*arr2[1] + arr2[3]) + 2*beta.*arr2[2]
            transformedCoeffs[idxs...,2] += thisCoeff .* arr3[2]
        end

        # All middle entries
        for i in 2:maxRow-2
            arr3[i+1] = -arr1[i+1] + alpha.*(arr2[i] + arr2[i+2]) + 2*beta.*arr2[i+1]
            transformedCoeffs[idxs...,i+1] += thisCoeff .* arr3[i+1]
        end

        # The second to last entry
        i = maxRow -1
        arr3[i+1] = -arr1[i+1] + (i == 1 ? 2 : 1)*alpha.*(arr2[i]) + 2*beta.*arr2[i+1]
        transformedCoeffs[idxs...,i+1] += thisCoeff .* arr3[i+1]
        #The last entry
        finalVal = alpha*arr2[i+1]
        # This final entry is typically very small. If it is essentially machine epsilon,
        # zero it out to save calculations.
        if abs(finalVal) > 1e-16 #TODO: Justify this val!
            arr3[maxRow+1] = finalVal
            transformedCoeffs[idxs...,maxRow+1] += thisCoeff * finalVal
            maxRow += 1 # Next column will have one more entry than the current column.
        end

        # Save the values of arr2 and arr3 to arr1 and arr2 to get ready for calculating the next column.
        arr = arr1
        arr1 = arr2
        arr2 = arr3
        arr3 = arr
    end
    return transformedCoeffs[idxs...,1:maxRow]
end

function TransformChebInPlaceND(coeffs, dim, alpha, beta, exact)
    """Transforms a single dimension of a Chebyshev approximation for a polynomial.

    Parameters
    ----------
    coeffs : array
        The coefficient tensor to transform
    dim : int
        The index of the dimension to transform (numpy index for now)
    alpha: double
        The scaler of the transformation
    beta: double
        The shifting of the transformation
    exact: bool
        Whether to perform the transformation with higher precision to minimize error (currently unimplemented)

    Returns
    -------
    transformedCoeffs : array
        The new coefficient array following the transformation
    """

    #TODO: Could we calculate the allowed error beforehand and pass it in here?
    #TODO: Make this work for the power basis polynomials
    if (((alpha == 1.0) && (beta == 0.0)) || (size(coeffs)[end-dim] == 1))
        return coeffs # No need to transform if the degree of dim is 0 or transformation is the identity.
    end

    transformFunc = transformChebInPlace1D

    if dim == 0
        return transformFunc(coeffs, alpha, beta)
    else # Need to transpose the matrix to line up the multiplication for the current dim
        ndim = length(size(coeffs))
        # Move the current dimension to the dim 0 spot in the np array.
        python_order = vcat([dim+1], collect(1:dim), collect(dim+2:ndim))
        order = (ndim+1) .- reverse(python_order)
        # Then transpose with the inverted order after the transformation occurs.
        python_backOrder = zeros(Int,ndim)
        python_backOrder[python_order] = collect(1:length(size(coeffs)))
        backOrder = (ndim+1) .- reverse(python_backOrder)
        return permutedims(transformFunc(permutedims(coeffs,order), alpha, beta),backOrder)
    end
end

function getTransformationError(M,dim)
    """Returns an upper bound on the error of transforming the Chebyshev approximation M

    In the transformation of dimension dim in M, the matrix multiplication of M by the transformation
    matrix C has each element of M involved in n element multiplications, where n is the number of rows
    in C, which is equal to the degree of approximation of M in dimension dim, or M.shape[dim].

    Parameters
    ----------
    M : array
        The Chebyshev approximation coefficient tensor being transformed
    dim : int
        The dimension of M being transformed

    Returns
    -------
    error : float
        The upper bound for the error associated with the transformation of dimension dim in M
    """

    machEps = 2^-52
    error = reverse(size(M))[dim+1] * machEps * sum(abs.(M))
    return error #TODO: Figure out a more rigurous bound!
end

function transformCheb(M,alphas,betas,error,exact)
    """Transforms an entire Chebyshev coefficient matrix using the transformation xHat = alpha*x + beta.

    Parameters
    ----------
    M : array
        The chebyshev coefficient matrix
    alphas : iterable
        The scalers in each dimension of the transformation.
    betas : iterable
        The offset in each dimension of the transformation.
    error : float
        A bound on the error of the chebyshev approximation
    exact : bool
        Whether to perform the transformation with higher precision to minimize error

    Returns
    -------
    M : numpy array
        The coefficient matrix transformed to the new interval
    error : float
        An upper bound on the error of the transformation
    """
    #This just does the matrix multiplication on each dimension. Except it's by a tensor.
    ndim = length(size(M))
    for (dim,n,alpha,beta) in zip(0:ndim-1,size(M),alphas,betas)
        error += getTransformationError(M, dim)
        M = TransformChebInPlaceND(M,dim,alpha,beta,exact)
    end
    return M, error
end

function transformChebToInterval(Ms, alphas, betas, errors, exact)
    """Transforms an entire list of Chebyshev approximations to a new interval xHat = alpha*x + beta.

    Parameters
    ----------
    Ms : list of arrays
        The chebyshev coefficient matrices
    alphas : iterable
        The scalers of the transformation we are doing.
    betas : iterable
        The offsets of the transformation we are doing.
    errors : array
        A bound on the error of each Chebyshev approximation
    exact : bool
        Whether to perform the transformation with higher precision to minimize error

    Returns
    -------
    newMs : list of arrays
        The coefficient matrices transformed to the new interval
    newErrors : array
        The new errors associated with the transformed coefficient matrices
    """
    #Transform the chebyshev polynomials
    newMs = []
    newErrors = []
    for (M,e) in zip(Ms, errors)
        newM, newE = transformCheb(M, alphas, betas, e, exact)
        push!(newMs,newM)
        push!(newErrors,(newE))
    end
    return newMs, newErrors
end

function findVertices()
    return 0
end