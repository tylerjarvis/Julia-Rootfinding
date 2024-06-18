mutable struct TrackedInterval 
    # This struct is implemented by passing in one argument "interval"
    # eg: TrackedInterval([-1;-3.4;0])
    topInterval # = interval (by default) 
    interval # = interval (by default) 
    transforms # = [] (by default)
    ndim # = length(interval) (by default)
    empty # = false (by default)
    finalStep # = false (by default)
    canThrowOutFinalStep # = false (by default)
    possibleDuplicateRoots # = [] (by default)
    possibleExtraRoot # = false (by default)
    nextTransformPoints #Random Point near 0
    function TrackedInterval(interval)
        ndim = length(interval)
        new(interval,interval,[],ndim,false,false,false,[],false,[0.0394555475981047]*ndim)
    end
end