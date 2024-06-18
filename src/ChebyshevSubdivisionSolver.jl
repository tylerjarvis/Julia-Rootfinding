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

    println(size(M))
    
    for i in size(M)
        println(i-1)
        push!(A, (i-1) == 0 ? 1 : reshape(M,(1,length(M)))[spot+1])
        spot *= (i)
    end

    return reverse(A) # Return linear terms in dimension order.
end