include("QuadraticCheck.jl")
include("StructsWithTheirFunctions/SolverOptions.jl")
include("StructsWithTheirFunctions/TrackedInterval.jl")
using LinearAlgebra
using Logging
using RecursiveArrayTools

# TODO: import from a library like this one instead of crowding our sourcecode with pre-written code https://github.com/JeffreySarnoff/ErrorfreeArithmetic.jl/blob/main/src/sum.jl
function twoSum(a,b)
    """Returns x,y such that a+b=x+y exactly, and a+b=x in floating point."""
    x = a .+ b
    z = x .- a
    y = (a .- (x .- z)) .+ (b .- z)
    return x, y
end

# TODO: import from a library instead
function homebrewSplit(a)
    """Returns x,y such that a = x+y exactly and a = x in floating point."""
    c = (2^27 + 1) * a
    x = c-(c-a)
    y = a-x
    return x,y
end

# TODO: import from a library instead
function twoProd(a,b)
    """Returns x,y such that a*b=x+y exactly and a*b=x in floating point."""
    x = a .* b
    a1,a2 = homebrewSplit(a)
    b1,b2 = homebrewSplit(b)
    y=a2 .*b2-(((x-a1 .*b1)-a2 .*b1)-a1 .*b2)
    return x,y
end

# TODO: import from a library instead
function TwoProdWithSplit(a,b,a1,a2)
    """Returns x,y such that a*b = x+y exactly and a*b = x in floating point but with a already split."""
    x = a*b
    b1,b2 = homebrewSplit(b)
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
    newM = reshape(M,(1,length(M)))
    
    for i in size(M)
        push!(A, (i) == 1 ? 0 : newM[spot+1])
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

function transformChebInPlace1D1D(coeffs,alpha,beta,precision)
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
    for col in 3:last_dim_length # For each column, calculate each entry and do matrix mult
        thisCoeff = coeffs[col] # the row of coeffs corresponding to the column col of C (for matrix mult)
        # The first entry
        arr3[1] = -arr1[1] + alpha*arr2[2] + 2*beta*arr2[1]
        transformedCoeffs[1] += thisCoeff * arr3[1]
        # The second entry
        if maxRow > 2
            arr3[2] = -arr1[2] + alpha*(2*arr2[1] + arr2[3]) + 2*beta*arr2[2]
            transformedCoeffs[2] += thisCoeff * arr3[2]
        end

        # All middle entries
        for i in 3:maxRow-1
            arr3[i] = -arr1[i] + alpha*(arr2[i-1] + arr2[i+1]) + 2*beta*arr2[i]
            transformedCoeffs[i] += thisCoeff * arr3[i]
        end

        # The second to last entry
        i = maxRow
        arr3[i] = -arr1[i] + (i == 2 ? 2 : 1)*alpha*(arr2[i-1]) + 2*beta*arr2[i]
        transformedCoeffs[i] += thisCoeff * arr3[i]
        #The last entry
        finalVal = alpha*arr2[i]
        # This final entry is typically very small. If it is essentially machine epsilon,
        # zero it out to save calculations.
        if abs(finalVal) > 2.0^(-precision) #TODO: Justify this val!
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

function transformChebInPlace1D(coeffs,alpha,beta,precision)
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
        return transformChebInPlace1D1D(coeffs,alpha,beta,precision)
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
    # idxs = []
    dims = length(coeffs_shape)
    # for i in 1:dims-1
    #     push!(idxs,:)
    # end

    # oldSlices = mapslices(x->[x],coeffs,dims = 1:dims-1)
    oldSlices = collect(eachslice(coeffs,dims=dims))
    # newSlices = mapslices(x->[x],transformedCoeffs,dims = 1:dims-1)
    newSlices = collect(Array.(eachslice(transformedCoeffs,dims=dims)))
    newSlices[1] = oldSlices[1] # arr1[0] * coeffs[0] (matrix multiplication step)
    #The second column of C. Note that T_1(alpha*x + beta) = alpha*T_1(x) + beta*T_0(x).
    arr2[1] = beta
    arr2[2] = alpha
    newSlices[1] += beta * oldSlices[2] # arr2[0] * coeffs[1] (matrix muliplication)
    newSlices[2] += alpha * oldSlices[2] # arr2[1] * coeffs[1] (matrix multiplication)
    maxRow = 2
    for col in 2:last_dim_length-1 # For each column, calculate each entry and do matrix mult
        thisCoeff = oldSlices[col+1] # the row of coeffs corresponding to the column col of C (for matrix mult)
        # The first entry
        arr3[1] = -arr1[1] + alpha*arr2[2] + 2*beta*arr2[1]
        newSlices[1] += thisCoeff * arr3[1]
        # The second entry
        if maxRow > 2
            arr3[2] = -arr1[2] + alpha*(2*arr2[1] + arr2[3]) + 2*beta*arr2[2]
            newSlices[2] += thisCoeff * arr3[2]
        end

        # All middle entries
        for i in 3:maxRow-1
            arr3[i] = -arr1[i] + alpha*(arr2[i-1] + arr2[i+1]) + 2*beta*arr2[i]
            # newSlices[i] += thisCoeff .* arr3[i]
            newSlices[i] += thisCoeff * arr3[i]
        end

        # The second to last entry
        i = maxRow
        arr3[i] = -arr1[i] + (i == 2 ? 2 : 1)*alpha*(arr2[i-1]) + 2*beta*arr2[i]
        newSlices[i] += thisCoeff * arr3[i]
        #The last entry
        finalVal = alpha*arr2[i]
        # This final entry is typically very small. If it is essentially machine epsilon,
        # zero it out to save calculations.
        if abs(finalVal) > 2.0^(-precision) #TODO: Justify this val!
            arr3[maxRow+1] = finalVal
            newSlices[maxRow+1] += thisCoeff * finalVal
            maxRow += 1 # Next column will have one more entry than the current column.
        end

        # Save the values of arr2 and arr3 to arr1 and arr2 to get ready for calculating the next column.
        arr = arr1
        arr1 = arr2
        arr2 = arr3
        arr3 = arr
    end
    VA = VectorOfArray(newSlices[1:maxRow])
    return convert(Array,VA)
    # return cat(newSlices[1:maxRow]...,dims = dims)
end

function TransformChebInPlaceND(coeffs, dim, alpha, beta, exact, precision)
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
    if (((alpha == 1.0) && (beta == 0.0)) || (size(coeffs)[end-(dim-1)] == 1))
        return coeffs # No need to transform if the degree of dim is 0 or transformation is the identity.
    end

    transformFunc = transformChebInPlace1D

    if dim == 1
        return transformFunc(coeffs, alpha, beta,precision)
    else # Need to transpose the matrix to line up the multiplication for the current dim
        ndim = length(size(coeffs))
        # Move the current dimension to the dim 0 spot in the np array.
        python_order = vcat([dim], collect(1:dim-1), collect(dim+1:ndim))
        order = (ndim+1) .- reverse(python_order)
        # Then transpose with the inverted order after the transformation occurs.
        python_backOrder = zeros(Int,ndim)
        python_backOrder[python_order] = collect(1:length(size(coeffs)))
        backOrder = (ndim+1) .- reverse(python_backOrder)
        return permutedims(transformFunc(permutedims(coeffs,order), alpha, beta, precision),backOrder)
    end
end

function getTransformationError(M,dim,precision)
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

    machEps = 2.0^(-(precision-1.))
    error = reverse(size(M))[dim] * machEps * sum(abs.(M))
    return error #TODO: Figure out a more rigurous bound!
end

function transformCheb(M,alphas,betas,error,exact,precision)
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
    for (dim,alpha,beta) in zip(1:ndim,alphas,betas)
        error += getTransformationError(M, dim,precision)
        M = TransformChebInPlaceND(M,dim,alpha,beta,exact,precision)
    end
    return M, error
end

function transformChebToInterval(Ms, alphas, betas, errors, exact, precision)
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
        newM, newE = transformCheb(M, alphas, betas, e, exact,precision)
        push!(newMs,newM)
        push!(newErrors,(newE))
    end
    return newMs, newErrors
end

function findVertices(A,b,errors)
    dim = length(b)
    V = reshape([(-1)^div(2^j*(i-1),2^(dim)) for i=1:2^dim for j=1:dim],(dim,2^dim))
    W = V .* e .+ z
    A = hcat(A',W)
    B = rref(A)
    verts = B[:,dim+1:end]
    bs = maximum(verts,dims=2)
    as = minimum(verts,dims=2)
    return as,bs
end

function boundingIntervalLinearSystem(Ms, errors, finalStep)
    """Finds a smaller region in which any root must be.

    Parameters
    ----------
    Ms : list of numpy arrays
        Each numpy array is the coefficient tensor of a chebyshev polynomials
    errors : iterable of floats
        The maximum error of chebyshev approximations
    finalStep : bool
        Whether we are in the final step of the algorithm

    Returns
    -------
    newInterval : numpy array
        The smaller interval where any root must be
    changed : bool
        Whether the interval has shrunk at all
    should_stop : bool
        Whether we should stop subdividing
    throwout :
        Whether we should throw out the interval entirely
    """
    if finalStep
        errors = zeros(size(errors))
    end
    s=1
    dim = length(Ms)
    #Some constants we use here
    widthToAdd = 1e-10 #Add this width to the new intervals we find to avoid rounding error throwing out roots
    minZoomForChange = 0.99 #If the volume doesn't shrink by this amount say that it hasn't changed
    minZoomForBaseCaseEnd = 0.4^dim #If the volume doesn't change by at least this amount when running with no error, stop
    #Get the matrix of the linear terms
    A = Matrix{Float64}(reduce(hcat,[getLinearTerms(M) for M in Ms]))
    #Get the Vector of the constant terms
    consts = [M[1] for M in Ms]'
    #Get the Error of everything else combined.
    totalErrs = [sum(abs.(Ms[i])) + errors[i] for i = 1:dim]'
    linear_sums = sum(abs.(A),dims=1)
    err = totalErrs - abs.(consts) - linear_sums

    #Scale all the polynomials relative to one another
    errors = Base.copy(errors')
    errors_0 = Base.copy(errors)
    for i in 1:dim
        scaleVal = maximum(abs.(A[:,i]))
        if scaleVal > 0
            s = 2. ^Int(floor(log2(abs(scaleVal))))
            A[:,i] /= s
            consts[i] /= s
            totalErrs[i] /= s
            linear_sums[i] /= s
            err[i] /= s
            errors[i] /= s
        end
    end
    #Precondition the columns. (AP)X = B -> A(PX) = B. So scale columns, solve, then scale the solution.
    colScaler = ones(dim)
    for i in 1:dim
        scaleVal = maximum(abs.(A[i,:]))
        if scaleVal > 0
            s = 2^(-floor(log2(abs(scaleVal))))
            colScaler[i] = s
            totalErrs += abs.(A[i,:])' * (s - 1)
            A[i,:] *= s
        end
    end

    #Run linear algorithm for shrinking or deciding whether to subdivide.
    #This loop will only execute the second time if the interval was not changed on the first iteration and it needs to run again with tighter errors
    #Use the other interval shrinking method
    a0, b0 = linearCheck1(totalErrs, A, consts)
    a_orig = a0
    b_orig = b0
    U,S,Vh = svd(A')
    Ainv = ((1 ./ S).*Vh')' * (U')
    center = -Ainv*consts'
    wellConditioned = S[1] > 0 && S[end]/S[1] > 1e-10
    for i = 0:1
        #Now do the linear solve check
        #We use the matrix inverse to find the width, so might as well use it both spots. Should be fine as dim is small.
        if wellConditioned #Make sure conditioning is ok.
            width = abs.(Ainv)*err'
            a1 = center-width
            b1 = center + width
            a = mapslices(x->maximum(x),hcat(a0,a1),dims = 2)
            b = mapslices(x->minimum(x),hcat(b0,b1),dims = 2)
        else
            a = a0
            b = b0
        end
        #Undo the column preconditioning
        a .*= colScaler
        b .*= colScaler
        #Add error and bound
        a .-= widthToAdd
        b .+= widthToAdd
        throwOut = any(a .> b) || any(a .> 1) || any(b .< -1)
        a[a .< -1] .= -1
        b[b .< -1] .= -1
        a[a .> 1] .= 1
        b[b .> 1] .= 1

        forceShouldStop = finalStep && !wellConditioned
        # Calculate the "changed" variable
        newRatio = prod(b - a) ./ 2^dim
        changed = false
        if throwOut
            changed = true
        elseif i == 0
            changed = newRatio < minZoomForChange
        else
            changed = newRatio < minZoomForBaseCaseEnd
        end

        if i == 0 && changed
            #If it is the first time through the loop and there was a change, return the interval it shrunk down to and set "is_done" to false
            return hcat(a,b)', changed, forceShouldStop, throwOut
        elseif i == 0 && !changed
            #If it is the first time through the loop and there was not a change, save the a and b as the original values to return,
            #and then try running through the loop again with a tighter error to see if we shrink then
            a_orig = a
            b_orig = b
            err = errors
        elseif changed
            #If it is the second time through the loop and it did change, it means we didn't change on the first time,
            #but that the interval did shrink with tighter errors. So return the original interval with changed = False and is_done = False
            return hcat(a_orig, b_orig)', false, forceShouldStop, false
        else
            #If it is the second time through the loop and it did NOT change, it means we will not shrink the interval even if we subdivide,
            #so return the original interval with changed = False and is_done = wellConditioned
            return hcat(a_orig,b_orig)', false, wellConditioned || forceShouldStop, false
        end
    end
end

function zoomInOnIntervalIter(Ms, errors, trackedInterval, exact, precision)
    """One iteration of shrinking an interval that may contain roots.

    Calls BoundingIntervaLinearSystem which determines a smaller interval in which any roots are
    bound to lie. Then calls transformChebToInterval to transform the current coefficient
    approximations to the new interval.

    Parameters
    ----------
    Ms : list of arrays
        The Chebyshev coefficient tensors of each approximation
    errors : array
        An upper bound on the error of each Chebyshev approximation
    trackedInterval : TrackedInterval
        The current interval for which the Chebyshev approximations are valid
    exact : bool
        Whether the transformation should be done with higher precision to minimize error

    Returns
    -------
    Ms : list of arrays
        The chebyshev coefficient matrices transformed to the new interval
    errors : array
        The new errors associated with the transformed coefficient matrices
    trackedInterval : TrackedInterval
        The new interval that the transformed coefficient matrices are valid for
    changed : bool
        Whether or not the interval shrunk significantly during the iteration
    should_stop : bool
        Whether or not to continue subdiviing after the iteration of shrinking is completed
    """

    dim = length(Ms)
    #Zoom in on the current interval
    interval, changed, should_stop, throwOut = boundingIntervalLinearSystem(Ms, errors, trackedInterval.finalStep)
    #Don't zoom in if we're already at a point
    for curr_dim in 1:dim
        if trackedInterval.interval[1,curr_dim] == trackedInterval.interval[2,curr_dim]
            interval[1,curr_dim] = -1.
            interval[2,curr_dim] = 1.
        end
    end
    #We can't throw out on the final step
    if throwOut &&  ! canThrowOut(trackedInterval)
        throwOut = false
        should_stop = true
        changed = true
    end
    #Check if we can throw out the whole thing
    if throwOut
        trackedInterval.empty = true
        return Ms, errors, trackedInterval, true, true
    end
    #Check if we are done iterating
    if !changed
        return Ms, errors, trackedInterval, changed, should_stop
    end
    #Transform the chebyshev polynomials
    addTransform(trackedInterval,interval)
    Ms, errors = transformChebToInterval(Ms, getLastTransform(trackedInterval)[:,1],getLastTransform(trackedInterval)[:,2], errors, exact,precision)
    #We should stop in the final step once the interval has become a point
    if trackedInterval.finalStep && isPoint(trackedInterval)
        should_stop = true
        changed = false
    end

    return Ms, errors, trackedInterval, changed, should_stop
end

function getSubdivisionDims(Ms,trackedInterval,level)
    """Decides which dimensions to subdivide in and in what order.
    
    Parameters
    ----------
    Ms : list of arrays
        The chebyshev coefficient matrices
    trackedInterval : trackedInterval
        The interval to be subdivided
    level : int
        The current depth of subdivision from the original interval

    Returns
    -------
    allDims : numpy array
        The ith row gives the dimensions in which Ms[i] should be subdivided, in order.
    """
    dim = length(Ms)
    dims_to_consider = collect(1:dim)
    new_dims = []
    for i in 1:dim
        if !isapprox(trackedInterval.interval[1,i], trackedInterval.interval[2,i],rtol=1e-5,atol=1e-8) || (i == dim && length(new_dims) == 0)
            push!(new_dims,dims_to_consider[i])
        end
    end
    dims_to_consider = new_dims
    if level > 5
        idxs_by_dim = [reverse(dims_to_consider[sortperm(reverse(collect(size(M)))[(dims_to_consider)])]) for M in Ms]
        return idxs_by_dim
    else
        dim_lengths = dimSize(trackedInterval)
        max_length = maximum(dim_lengths[dims_to_consider])
        dims_to_consider = filter(x -> dim_lengths[(x)] > max_length/5,dims_to_consider)
        if length(dims_to_consider) > 1
            shapes = reverse(reduce(hcat,[collect(size(M)) for M in Ms]))
            degree_sums = sum(shapes,dims=2)
            total_sum = sum(shapes)
            for i in Base.copy(dims_to_consider)
                if length(dims_to_consider) > 1 && (degree_sums[i] < floor(total_sum/(dim+1)))
                    dims_to_consider = deleteat!(dims_to_consider, findall(x->x==i,dims_to_consider))
                end
            end
        end
        if length(dims_to_consider) == 0
            dims_to_consider = [1]
        end
        idxs_by_dim = [reverse(dims_to_consider[sortperm(reverse(collect(size(M)))[(dims_to_consider)])]) for M in Ms]
        return idxs_by_dim
    end
end

function getSubdivisionDim(shapes,dim) # largest degree subdivide
    return Int.(((argmax(shapes)-1)%dim) .* ones(dim)')
end

# function getSubdivisionDim(interval) # largest size subdivide
#     (argmax(interval[1,:] - interval[2,:]) - 1) .* ones(Int(length(interval)/2))'
# end

function getInverseOrder(order)
    """Gets a particular order of matrices needed in getSubdivisionIntervals (helper function).

    Takes the order of dimensions in which a Chebyshev coefficient tensor M was subdivided and gets
    the order of the indexes that will arrange the list of resulting transformed matrices as if the
    dimensions had bee subdivided in standard index order. For example, if dimensions 0, 3, 1 were
    subdivided in that order, this function returns the order [0,2,1,3,4,6,5,7] corresponding to the
    indices of currMs such that when arranged in this order, it appears as if the dimensions were
    subdivided in order 0, 1, 3.

    Parameters
    ----------
    order : array
        The order of dimensions along which a coefficient tensor was subdivided

    Returns
    -------
    invOrder : array
        The order of indices of currMs (in the function getSubdivisionIntervals) that arranges the
        matrices resulting from the subdivision as if the original matrix had been subdivided in
        numerical order
    """
    t = ones(length(order))
    t[sortperm(order)] = collect(0:length(t)-1)
    order = t
    len = length(order)-1
    order = Int.([2^(len-x) for x in order])
    combinations = Iterators.product(fill([0,1],len+1)...)
    newOrder_matrix = [collect(reverse(i))'*order for i in combinations]
    newOrder = reshape(newOrder_matrix,(1,length(newOrder_matrix)))
    invOrder = ones(length(newOrder))
    invOrder[newOrder .+ 1] = collect(1:length(newOrder))
    return Tuple(Int.(invOrder))
end

function getSubdivisionIntervals(Ms,errors,trackedInterval,exact,level,precision;oneDimension=false)
    """Gets the matrices, error bounds, and intervals for the next iteration of subdivision.

    Parameters
    ----------
    Ms : list of arrays
        The chebyshev coefficient matrices
    errors : array
        An upper bound on the error of each Chebyshev approximation
    trackedInterval : trackedInterval
        The interval to be subdivided
    exact : bool
        Whether transformations should be completed with higher precision to minimize error
    level : int
        The current depth of subdivision from the original interval

    Returns
    -------
    allMs : list of arrays
        The transformed coefficient matrices associated with each new interval
    allErrors : array
        A list of upper bounds for the errors associated with each transformed coefficient matrix
    allIntervals : list of TrackedIntervals
        The intervals from the subdivision (corresponding one to one with the matrices in allMs)
    """
    
    if oneDimension
        subdivisionDims = getSubdivisionDim(reduce(vcat,[[size(M)...] for M in Ms]),length(Ms))
    else
        subdivisionDims = getSubdivisionDims(Ms,trackedInterval,level)
    end
    dimSet = sort(subdivisionDims[1])
    allMs = []
    allErrors = []
    idx = 0
    dim = length(Ms)
    for (M,error,order_num) in zip(Ms, errors, collect(1:dim))
        order = subdivisionDims[order_num]
        idx += 1
        #Iterate through the dimensions, highest degree first.
        currMs, currErrs = [M],[error]
        for thisDim in order
            newMidpoint = trackedInterval.nextTransformPoints[thisDim]
            alpha, beta = (newMidpoint+1)/2, (newMidpoint-1)/2
            tempMs = []
            tempErrs = []
            for (T,E) in zip(currMs, currErrs)
                #Transform the polys
                P1, P2 = TransformChebInPlaceND(T, thisDim, alpha, beta, thisDim,precision), TransformChebInPlaceND(T, thisDim, -beta, alpha, exact,precision)
                E1 = getTransformationError(T, thisDim,precision)
                push!(tempMs,P1)
                push!(tempMs,P2)
                push!(tempErrs,E1+E)
                push!(tempErrs,E1+E)
            end
            currMs = tempMs
            currErrs = tempErrs
        end
        if ndims(M) == 1
            push!(allMs,currMs) #Already ordered because there's only 1.
            push!(allErrors,currErrs) #Already ordered because there's only 1.
        else
            #Order the polynomials so they match the intervals in subdivideInterval
            invOrder = getInverseOrder(order)
            push!(allMs,[currMs[i] for i in invOrder])
            push!(allErrors,([currErrs[i] for i in invOrder]))
        end
    end
    allMs = [[allMs[i][j] for i in eachindex(allMs)] for j in eachindex(allMs[1])]
    allErrors = [[allErrors[i][j] for i in eachindex(allErrors)] for j in eachindex(allErrors[1])]
    #Get the intervals
    allIntervals = [trackedInterval]
    
    for thisDim in dimSet
        newMidpoint = trackedInterval.nextTransformPoints[thisDim]
        newSubinterval = ones(size(trackedInterval.interval)) #TODO: Make this outside for loop
        newSubinterval[1,:] .= -1.
        newIntervals = []
        for oldInterval in allIntervals
            newInterval1 = copyInterval(oldInterval)
            newInterval2 = copyInterval(oldInterval)
            newSubinterval[:,thisDim] = [-1.; newMidpoint]
            addTransform(newInterval1,newSubinterval)
            newSubinterval[:,thisDim] = [newMidpoint; 1.]
            addTransform(newInterval2,newSubinterval)
            newInterval1.nextTransformPoints[thisDim] = 0
            newInterval2.nextTransformPoints[thisDim] = 0
            push!(newIntervals,newInterval1)
            push!(newIntervals,newInterval2)
        end
        allIntervals = newIntervals
    end
    return allMs, allErrors, allIntervals
end

function isExteriorInterval(originalInterval,trackedInterval)
    """Determines if the current interval is exterior to its original interval."""
    return any(getIntervalForCombining(trackedInterval) .== getIntervalForCombining(originalInterval))
end

function trimMs(Ms, errors, relApproxTol=1e-3, absApproxTol=0)
    """Reduces the degree of each chebyshev approximation M when doing so has negligible error.

    The coefficient matrices are trimmed in place. This function iteratively looks at the highest
    degree coefficient row of each M along each dimension and trims it as long as the error introduced
    is less than the allowed error increase for that dimension.

    Parameters
    ----------
    Ms : list of arrays
        The chebyshev approximations of the functions
    errors : array
        The max error of the chebyshev approximation from the function on the interval
    relApproxTol : double
        The relative error increase allowed
    absApproxTol : double
        The absolute error increase allowed
    """
    dim = ndims(Ms[1])
    for polyNum in 1:dim #Loop through the polynomials
        allowedErrorIncrease = absApproxTol + errors[polyNum] * relApproxTol
        #Use slicing to look at a slice of the highest degree in the dimension we want to trim
        slices = []
        for i in 1:dim
            push!(slices,:)
        end
        # [: for i in 1:dim] # equivalent to selecting everything
        for currDim in dim:-1:1
            slices[currDim] = size(Ms[polyNum])[currDim] # Now look at just the last row of the current dimension's approximation
            lastSum = sum(abs.(Ms[polyNum][slices...]))
            # Iteratively eliminate the highest degree row of the current dimension if
            # the sum of its approximation coefficients is of low error, but keep deg at least 2
            while (lastSum < allowedErrorIncrease) && (size(Ms[polyNum])[currDim] > 3)
                # Trim the polynomial
                slices[currDim] = 1:(slices[currDim]-1)
                Ms[polyNum] = Ms[polyNum][slices...]
                # Update the remaining error increase allowed an the error of the approximation.
                allowedErrorIncrease -= lastSum
                errors[polyNum] += lastSum
                # Reset for the next iteration with the next highest degree of the current dimension.
                slices[currDim] = size(Ms[polyNum])[currDim]
                lastSum = sum(abs.(Ms[polyNum][slices...]))
            end
            # Reset to select all of the current dimension when looking at the next dimension.
            slices[currDim] = 1:size(Ms[polyNum])[currDim]
        end
    end
end

function solvePolyRecursive(Ms,trackedInterval,errors,solverOptions,precision)
    """Recursively shrinks and subdivides the given interval to find the locations of all roots.

    Parameters
    ----------
    Ms : list of arrays
        The chebyshev approximations of the functions
    trackedInterval : TrackedInterval
        The information about the interval we are solving on.
    errors : array
        An upper bound for the error of the Chebyshev approximation of the function on the interval
    solverOptions : SolverOptions
        Desired settings for running interval checks, transformations, and subdivision.

    Returns
    -------
    boundingBoxesInterior : list of arrays (optional)
        Each element of the list is an interval in which there may be a root. The interval is on the interior of the current
        interval
    boundingBoxesExterior : list of arrays (optional)
        Each element of the list is an interval in which there may be a root. The interval is on the exterior of the current
        interval
    """

    #TODO: Check if trackedInterval.interval has width 0 in some dimension, in which case we should get rid of that dimension.
    #If the interval is a point, return it
    if isPoint(trackedInterval)
        return [], [trackedInterval]
    end

    #If we ever change the options in this function, we will need to do a copy here.
    #Should be cheap, but as we never change them for now just avoid the copy
    solverOptions = copySO(solverOptions)
    solverOptions.level += 1

    #Constant term check, runs at the beginning of the solve and before each subdivision
    #If the absolute value of the constant term for any of the chebyshev polynomials is greater than the sum of the
    #absoulte values of any of the other terms, it will return that there are no zeros on that interval
    if solverOptions.constant_check
        consts = [M[1] for M in Ms]
        err = [sum(abs.(M))-abs(c)+e for (M,e,c) in zip(Ms,errors,consts)]
        if any(abs.(consts) .> err)
            return [], []
        end
    end

    #Runs quadratic check after constant check, only for dimensions 2 and 3 by default
    #More expensive than constant term check, but testing show it saves time in lower dimensions
    if (solverOptions.low_dim_quadratic_check && ndims(Ms[1]) <= 3) || solverOptions.all_dim_quadratic_check
        for i in eachindex(Ms)
            if quadraticCheck(Ms[i], errors[i])
                return [], []
            end
        end
    end

    #Trim
    Ms = copy(Ms)
    originalMs = copy(Ms)
    trackedInterval = copyInterval(trackedInterval)
    errors = copy(errors)
    trimMs(Ms, errors)

    #Solve
    dim = ndims(Ms[1])
    changed = true
    zoomCount = 0
    originalInterval = copyInterval(trackedInterval)
    originalIntervalSize = sizeOfInterval(trackedInterval)
    #Zoom in while we can
    global should_stop = false
    lastSizes = dimSize(trackedInterval)
    while changed && zoomCount <= solverOptions.maxZoomCount
        #Zoom in until we stop changing or we hit machine epsilon
        Ms, errors, trackedInterval, changed, should_stop = zoomInOnIntervalIter(Ms, errors, trackedInterval, solverOptions.exact,precision)

        if trackedInterval.empty #Throw out the interval
            return [], []
        end
        #Only count in towards the max is we don't cut the interval in half
        newSizes = dimSize(trackedInterval)
        if all(newSizes .>= (lastSizes ./ 2)) #Check all dims and use >= to account for a dimension being 0.
            zoomCount += 1
        end
        lastSizes = newSizes
    end

    if should_stop
        #Start the final step if the is in the options and we aren't already in it.
        if trackedInterval.finalStep || !solverOptions.useFinalStep
            if solverOptions.verbose
                print("*")
            end
            if isExteriorInterval(originalInterval, trackedInterval)
                #input("val:")
                return [], [trackedInterval]
            else
                #input("val:")
                return [trackedInterval], []
            end
        else
            startFinalStep(trackedInterval)
            return solvePolyRecursive(Ms, trackedInterval, errors, solverOptions,precision)
        end

    elseif trackedInterval.finalStep
        trackedInterval.canThrowOutFinalStep = true
        allMs, allErrors, allIntervals = getSubdivisionIntervals(Ms, errors, trackedInterval, solverOptions.exact, solverOptions.level,precision;oneDimension=false)
        resultsAll = []
        for (newMs, newErrs, newInt) in zip(allMs, allErrors, allIntervals)
            newInterior, newExterior = solvePolyRecursive(newMs, newInt, newErrs, solverOptions,precision)
            append!(resultsAll, newInterior)
            append!(resultsAll,newExterior)
        end
        if length(resultsAll) == 0
            #Can't throw out final step! This might not actually be a root though!
            trackedInterval.possibleExtraRoot = true
            if isExteriorInterval(originalInterval, trackedInterval)
                return [], [trackedInterval]
            else
                return [trackedInterval], []
            end
        else
            #Combine all roots that converged to the same point.
            allFoundRoots = Set([])
            tempResults = []
            for result in resultsAll
                point = Tuple(result.interval[1,:])
                if point in allFoundRoots
                    continue
                end
                push!(allFoundRoots,point)
                push!(tempResults,result)
            end
            for result in tempResults
                if length(result.possibleDuplicateRoots) > 0
                    append!(trackedInterval.possibleDuplicateRoots,result.possibleDuplicateRoots)
                else
                    push!(trackedInterval.possibleDuplicateRoots,getFinalPoint(result))
                end
            end
            if isExteriorInterval(originalInterval, trackedInterval)
                return [], [trackedInterval]
            else
                return [trackedInterval], []
            end
            #TODO: Don't subdivide in the final step in dimensions that are already points!
        end
    else 
        #Otherwise, Subdivide
        if solverOptions.level == 15
            @warn "HIGH SUBDIVISION DEPTH!\nSubdivision on the search interval has now reached recursion depth 15. Runtime may be long."
        elseif solverOptions.level == 25
            @warn "HIGH SUBDIVISION DEPTH!\nExtreme subdivision depth!\nSubdivision on the search interval has now reached" *
                        " at least depth 25, which is unusual. The solver may not finish running." *
                        "Ensure the input functions meet the requirements of being continuous, smooth," *
                        "and having only finitely many simple roots on the search interval."
        end
        resultInterior, resultExterior = [], []
        #Get the new intervals and polynomials
        allMs, allErrors, allIntervals = getSubdivisionIntervals(Ms, errors, trackedInterval, solverOptions.exact, solverOptions.level,precision;oneDimension=false)
        #Run each interval
        for (newMs, newErrs, newInt) in zip(allMs, allErrors, allIntervals)
            newInterior, newExterior = solvePolyRecursive(newMs, newInt, newErrs, solverOptions,precision)
            append!(resultInterior, newInterior)
            append!(resultExterior, newExterior)
        end
        #Rerun the touching intervals
        idx1 = 1
        idx2 = 2
        #Combine any touching intervals and throw them at the end. Flip a bool saying rerun them
        #If changing this code, test it by defaulting the nextTransformationsInterals to 0, so roots lie on the boundary more.
        #TODO: Make the combining intervals it's own function!!!
        for tempInterval in resultExterior
            tempInterval.reRun = false
        end
        while idx1 <= length(resultExterior)
            while idx2 <= length(resultExterior)
                if overlapsWith(resultExterior[idx1],resultExterior[idx2])
                    #Combine, throw at the back. Set reRun to true.
                    combinedInterval = copyInterval(originalInterval)
                    if combinedInterval.finalStep
                        combinedInterval.interval = copyInterval(combinedInterval.preFinalInterval)
                        combinedInterval.transforms = copyInterval(combinedInterval.preFinalTransforms)
                    end
                    newAs = minimum(reduce(hcat,[getIntervalForCombining(resultExterior[idx1])[1,:], getIntervalForCombining(resultExterior[idx2])[1,:]]),dims=2)
                    newBs = maximum(reduce(hcat,[getIntervalForCombining(resultExterior[idx1])[2,:], getIntervalForCombining(resultExterior[idx2])[2,:]]),dims=2)
                    final1 = getFinalInterval(resultExterior[idx1])
                    final2 = getFinalInterval(resultExterior[idx2])
                    newAsFinal = minimum(reduce(hcat,[final1[1,:], final2[1,:]]),dims=2)
                    newBsFinal = maximum(reduce(hcat,[final1[2,:], final2[2,:]]),dims=2)
                    oldAs = originalInterval.interval[1,:]
                    oldBs = originalInterval.interval[2,:]
                    oldAsFinal, oldBsFinal = getFinalInterval(originalInterval)[1,:],getFinalInterval(originalInterval)[2,:]
                    #Find the final A and B values exactly. Then do the currSubinterval calculation exactly.
                    #Look at what was done on the example that's failing and see why.
                    equalMask = oldBsFinal .== oldAsFinal
                    oldBsFinal[equalMask] = oldBsFinal[equalMask] .+ 1 #Avoid a divide by zero on the next line
                    currSubinterval = ((2 .*reduce(hcat,[newAsFinal, newBsFinal]) .- oldAsFinal .- oldBsFinal)./(oldBsFinal .- oldAsFinal))'
                    #If the interval is exactly -1 or 1, make sure that shows up as exact.
                    currSubinterval[1,equalMask] .= -1
                    currSubinterval[2,equalMask] .= 1
                    currSubinterval[1,:][reshape(oldAs.==newAs,length(oldAs))] .= -1
                    currSubinterval[2,:][reshape(oldBs.==newBs,length(newBs))] .= 1
                    #Update the current subinterval. Use the best transform we can get here, but use the exact combined
                    #interval for tracking
                    addTransform(combinedInterval,currSubinterval)
                    combinedInterval.interval = reduce(hcat,[newAs, newBs])'
                    combinedInterval.reRun = true
                    deleteat!(resultExterior,idx2)
                    deleteat!(resultExterior,idx1)
                    push!(resultExterior,combinedInterval)
                    idx2 = idx1 + 1
                else
                    idx2 += 1
                end
            end
            idx1 += 1
            idx2 = idx1 + 1
        end
        #Rerun, check if still on exterior
        newResultExterior = []
        for tempInterval in resultExterior
            if tempInterval.reRun
                if isapprox(tempInterval.interval,originalInterval.interval)
                    push!(newResultExterior,tempInterval)
                else
                    #Project the MS onto the interval, then recall the function.
                    #TODO: Instead of using the originalMs, use Ms, and then don't use the original interval, use the one
                    #we started subdivision with.
                    lastTransform = getLastTransform(tempInterval)
                    tempMs, tempErrors = transformChebToInterval(originalMs, lastTransform[:,1],lastTransform[:,2], errors, solverOptions.exact,precision)
                    tempResultsInterior, tempResultsExterior = solvePolyRecursive(tempMs, tempInterval, tempErrors, solverOptions,precision)
                    #We can assume that nothing in these has to be recombined
                    append!(resultInterior,tempResultsInterior)
                    append!(newResultExterior,tempResultsExterior)
                end
            elseif isExteriorInterval(originalInterval, tempInterval)
                push!(newResultExterior,tempInterval)
            else
                push!(resultInterior,tempInterval)
            end
        end
        return resultInterior, newResultExterior
    end
end

function solveChebyshevSubdivision(Ms, errors; verbose = false, returnBoundingBoxes = false, exact = false, constant_check = true, low_dim_quadratic_check = true, all_dim_quadratic_check = false,precision=53.)
    """
    Initiates shrinking and subdivision recursion and returns the roots and bounding boxes.

    Parameters
    ----------
    Ms : Vector{Array}
    The Chebyshev approximations of the functions on the interval given to CombinedSolver
    errors : Vector{Float64}
    The max error of the Chebyshev approximation from the function on the interval
    verbose : Bool
    Defaults to false. Whether or not to output progress of solving to the terminal.
    returnBoundingBoxes : Bool (Optional)
    Defaults to false. If true, returns the bounding boxes around each root as well as the roots.
    exact : Bool
    Whether transformations should be done with higher precision to minimize error.
    constant_check : Bool
    Defaults to true. Whether or not to run constant term check after each subdivision.
    low_dim_quadratic_check : Bool
    Defaults to true. Whether or not to run quadratic check in dim 2, 3.
    all_dim_quadratic_check : Bool
    Defaults to false. Whether or not to run quadratic check in dim ≥ 4.

    Returns
    -------
    roots : Vector
    The roots of the system of functions on the interval given to Combined Solver
    boundingBoxes : Vector{Array} (optional)
    List of intervals for each root in which the root is bound to lie.
    """

    # Assert that we have n nD polys
    if any(ndims(M) != length(Ms) for M in Ms)
        throw(ArgumentError("Solver takes in N polynomials of dimension N!"))
    end
    if length(Ms) != length(errors)
        throw(ArgumentError("Ms and errors must be same length!"))
    end

    # Solve
    ndim = length(Ms)
    #set precision of TrackedInterval
    if precision <= 11
        polys = Vector{Array{Float16,dim}}(undef, dim)
        errs = Float16.(errs)
        a = Float16.(a)
        b = Float16.(b)
    elseif precision <= 24
        polys = Vector{Array{Float32,dim}}(undef, dim)
        errs = Float32.(errs)
        a = Float32.(a)
        b = Float32.(b)
    elseif precision <= 53
        originalInterval = TrackedInterval(hcat(fill(-1.,ndim), fill(1.,ndim))')
    else
        setprecision(precision)
        originalInterval = TrackedInterval(hcat(BigFloat(fill(-1.,ndim)), BigFloat(fill(1.,ndim)))')
    end
    solverOptions = SolverOptions()
    solverOptions.verbose = verbose
    solverOptions.exact = exact
    solverOptions.constant_check = constant_check
    solverOptions.low_dim_quadratic_check = low_dim_quadratic_check
    solverOptions.all_dim_quadratic_check = all_dim_quadratic_check
    solverOptions.useFinalStep = true

    if verbose
        println("Finding roots...")
    end

    b1, b2 = solvePolyRecursive(Ms, originalInterval, errors, solverOptions,precision)

    boundingIntervals = append!(b1, b2)
    roots = []
    hasDupRoots = false
    hasExtraRoots = false

    for interval in boundingIntervals
        #TODO: Figure out the best way to return the bounding intervals.
        #Right now interval.finalInterval is the interval where we say the root is.
        getFinalInterval(interval)
        if interval.possibleExtraRoot
            hasExtraRoots = true
        end
        if !isempty(interval.possibleDuplicateRoots)
            append!(roots, interval.possibleDuplicateRoots)
            hasDupRoots = true
        else
            push!(roots, getFinalPoint(interval))
        end
    end

    # Warn if extra or duplicate roots
    if hasExtraRoots
        @warn "Might have extra roots! See bounding boxes for details!"
    end
    if hasDupRoots
        @warn "Might have duplicate roots! See bounding boxes for details!"
    end

    # Return
    roots = collect(roots)
    if verbose
        println("\nFound $(length(roots)) root(s)\n")
    end

    if returnBoundingBoxes
        return roots, boundingIntervals
    else
        return roots
    end
end