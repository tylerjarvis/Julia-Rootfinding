#using LinearAlgebra


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

function BoundingIntervalLinearSystem(Ms, errors, finalStep)
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
    for i = 0:1
        println(err)
        #Now do the linear solve check
        U,S,Vh = svd(A')
        wellConditioned = S[1] > 0 && S[end]/S[1] > 1e-10
        #We use the matrix inverse to find the width, so might as well use it both spots. Should be fine as dim is small.
        if wellConditioned #Make sure conditioning is ok.
            Ainv = ((1 ./ S).*Vh')' * (U')

            center = -Ainv*consts'
            width = abs.(Ainv)*err'
            a1 = center-width
            b1 = center + width
            a = mapslices(x->maximum(x),hcat(a0,a1),dims = 2)
            b = mapslices(x->minimum(x),hcat(b0,b1),dims = 2)
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

        println(a)
        println(b)
        println(throwOut)

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
            println("1")
            return hcat(a,b), changed, forceShouldStop, throwOut
        elseif i == 0 && !changed
            #If it is the first time through the loop and there was not a change, save the a and b as the original values to return,
            #and then try running through the loop again with a tighter error to see if we shrink then
            a_orig = a
            b_orig = b
            err = errors
            println("2")
        elseif changed
            #If it is the second time through the loop and it did change, it means we didn't change on the first time,
            #but that the interval did shrink with tighter errors. So return the original interval with changed = False and is_done = False
            #print("subdivide")
            println("3")
            return hcat(a_orig, b_orig), false, forceShouldStop, false
        else
            #If it is the second time through the loop and it did NOT change, it means we will not shrink the interval even if we subdivide,
            #so return the original interval with changed = False and is_done = wellConditioned
            #print("done")
            #print("throwout:",throwOut)
            println("4")
            return hcat(a_orig,b_orig), false, wellConditioned || forceShouldStop, false
        end
    end
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
    dims_to_consider = collect(0:dim-1)
    for i in 0:dim-1
        if isapprox(trackedInterval.interval[1,i+1], trackedInterval.interval[2,i+1],rtol=1e-5,atol=1e-8)
            if length(dims_to_consider) != 1
                dims_to_consider = deleteat!(dims_to_consider, findall(x->x==i,dims_to_consider))
            end
        end
    end
    if level > 5
        idxs_by_dim = [reverse(dims_to_consider[sortperm(reverse(collect(size(M))[(1 .+ dims_to_consider)]))]) for M in Ms]
        return reshape([item for sublist in idxs_by_dim for item in sublist],(length(dims_to_consider),dim))
    else
        dim_lengths = dimSize(trackedInterval)
        max_length = maximum([dim_lengths[i] for i in (1 .+ dims_to_consider)])
        dims_to_consider = filter(x -> dim_lengths[(1 + x)] > max_length/5,dims_to_consider)
        if length(dims_to_consider) > 1
            shapes = reverse(reduce(hcat,[collect(size(M)) for M in Ms]))
            # shapes_list = [reverse(collect(size(M))) for M in Ms]
            # shapes = reshape([item for sublist in shapes_list for item in sublist],(dim,dim))
            degree_sums = sum(shapes,dims=2)
            total_sum = sum(shapes)
            for i in Base.copy(dims_to_consider)
                if length(dims_to_consider) > 1 && (degree_sums[i+1] < floor(total_sum/(dim+1)))
                    dims_to_consider = deleteat!(dims_to_consider, findall(x->x==i,dims_to_consider))
                end
            end
        end
        if length(dims_to_consider) == 0
            dims_to_consider = [0]
        end
        idxs_by_dim = [reverse(dims_to_consider[sortperm(reverse(collect(size(M)))[(1 .+ dims_to_consider)])]) for M in Ms]
        for M in Ms
        end
        return reshape([item for sublist in idxs_by_dim for item in sublist],(length(dims_to_consider),dim))
    end
end
