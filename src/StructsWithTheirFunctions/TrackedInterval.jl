mutable struct TrackedInterval 
    """Tracks the properties of and changes to each interval as it passes through the solver.

    Parameters
    ----------
    topInterval: array
        The original interval before any changes
    interval: array
        The current interval (lower bound and upper bound for each dimension in order)
    transforms: array
        List of the alpha and beta values for all the transformations the interval has undergone
    ndim: Int
        The number of dimensions of which the interval consists
    empty: bool
        Whether the interval is known to contain no roots
    finalStep: bool
        Whether the interval is in the final step (zooming in on the bounding box to a point at the end)
    canThrowOutFinalStep: bool
        Defaults to false. Whether or not the interval should be thrown out if empty in the final step
        of solving. Changed to true if subdivision occurs in the final step.
    possibleDuplicateRoots: array
        Any multiple roots found through subdivision in the final step that would have been
        returned as just one root before the final step
    possibleExtraRoot: bool
        Defaults to false. Whether or not the interval would have been thrown out during the final step.
    nextTransformPoints: array
        Where the midpoint of the next subdivision should be for each dimension
    """

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

