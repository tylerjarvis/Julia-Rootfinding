struct MultiPower
    """Contains the coeffs array for a MultiPower object"""
    coeff
    dim

    function MultiPower(coeff)
        dim = ndims(coeff)
        new(coeff,dim)
    end
end

struct MultiCheb
    """Contains the coeffs array for a MultiCheb object"""
    coeff
end

function eval_MultiPower(multiPower,points)
    function polyval(x, cc)
        cc = collect(eachslice(cc,dims=ndims(cc)))
        c0 = cc[end]
        for i in 1:(length(cc) - 1)
            c0 = c0.*x .+ cc[end-i]
        end
        return c0
    end
    
    if ndims(points) == 1
        if multiPower.dim > 1
            points = reshape(points,(size(points)[end],1))
        else
            points = reshape(points,(1,size(points)[end]))
        end
    end

    if size(points)[end-1] != multiPower.dim
        1/0
    end

    c = multiPower.coeff
    n = ndims(c)
    cc = reshape(c,reverse((reverse(size(c))...,ntuple(i->1, ndims(points))...)))
    c = polyval(points[1,:],cc)
    for i in 2:n
        c = polyval(points[i,:],c)
    end
    if length(c) == 1
        return c[1]
    else
        return c
    end

end