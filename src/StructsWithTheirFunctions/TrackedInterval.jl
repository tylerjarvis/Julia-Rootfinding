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
    preFinalInterval # = [] (by default)
    preFinalTransforms # = [] (by default)
    reducedDims # = [] (by default)
    solvedVals # = [] (by default)
    finalInterval # = [] (by default)
    function TrackedInterval(interval)
        ndim = Int(length(interval)/2)
        new(interval,interval,[],ndim,false,false,false,[],false,fill(0.0394555475981047,ndim),[],[],[],[],[])
    end
end

"""==============================FUNCTIONS FOR TRACKED INTERVAL=============================="""

function canThrowOut(trackedInterval::TrackedInterval)
    """Ensures that an interval that has not subdivided cannot be thrown out on the final step."""
    return !trackedInterval.finalStep || trackedInterval.canThrowOutFinalStep
end

function addTransform(trackedInterval::TrackedInterval, subInterval)
    """Adds the next alpha and beta values to the list transforms and updates the current interval.

    Parameters:
    -----------
    subInterval : array
        The subinterval to which the current interval is being reduced
    """
    #Ensure the interval has non zero size; mark it empty if it doesn't
    if any(subInterval[1,:] > subInterval[2,:]) && canThrowOut(trackedInterval)
        trackedInterval.empty = true
        return
    elseif any(subInterval[1,:] > subInterval[2,:])
        #If we can't throw the interval out, it should be bounded by [-1,1].
        subInterval[1,:] = min.(subInterval[1,:], ones(length(subInterval[1,:])))
        subInterval[1,:] = max.(subInterval[1,:], -ones(length(subInterval[1,:])))
        subInterval[2,:] = min.(subInterval[2,:], ones(length(subInterval[1,:])))
        subInterval[2,:] = max.(subInterval[2,:], subInterval[1,:])
    end
    # Get the alpha and beta associated with the transformation in each dimension
    a1 = subInterval[1,:]
    b1 = subInterval[2,:] # all the lower bounds and upper bounds of the new interval, respectively
    a2 = trackedInterval.interval[1,:]
    b2 = trackedInterval.interval[2,:] # all the lower bounds and upper bounds of the original interval
    # println(a1,b1,a2,b2)
    alpha1, beta1 = (b1-a1)/2, (b1+a1)/2
    alpha2, beta2 = (b2-a2)/2, (b2+a2)/2
    # println(alpha1,beta1,alpha2,beta2) 
    push!(trackedInterval.transforms,hcat(alpha1, beta1))
    #Update the lower and upper bounds of the current interval
    for dim in 0:trackedInterval.ndim-1
        for i in 0:1
            x = subInterval[i+1,dim+1]
            #Be exact if x = +-1
            if x == -1.0
                trackedInterval.interval[i+1,dim+1] = trackedInterval.interval[1,dim+1]
            elseif x == 1.0
                trackedInterval.interval[i+1,dim+1] = trackedInterval.interval[2,dim+1]
            else
                trackedInterval.interval[i+1,dim+1] = alpha2[dim+1]*x+beta2[dim+1]
            end
        end
    end
end

function getLastTransform(trackedInterval::TrackedInterval)
    return trackedInterval.transforms[end]
end

# def getFinalInterval(self):
#     """Finds the interval that should be reported as containing a root.

#     The final interval is calculated by applying all of the recorded transformations that
#     occurred before the final step to topInterval, the original interval.

#     Returns
#     -------
#     finalInterval: numpy array
#         The final interval to be reported as containing a root
#     """
#     # TODO: Make this a seperate function so it can use njit.
#     # Make these _NoNumba calls use floats so they call call the numba functions without a seperate compile
#     finalInterval = self.topInterval.T
#     finalIntervalError = np.zeros_like(finalInterval)
#     transformsToUse = self.transforms if not self.finalStep else self.preFinalTransforms
#     for alpha,beta in transformsToUse[::-1]: # Iteratively apply each saved transform
#         finalInterval, temp = TwoProd_NoNumba(finalInterval, alpha)
#         finalIntervalError = alpha * finalIntervalError + temp
#         finalInterval, temp = TwoSum_NoNumba(finalInterval,beta)
#         finalIntervalError += temp

#     finalInterval = finalInterval.T
#     finalIntervalError = finalIntervalError.T
#     self.finalInterval = finalInterval + finalIntervalError # Add the error and save the result.
#     self.finalAlpha, alphaError = TwoSum_NoNumba(-finalInterval[:,0]/2,finalInterval[:,1]/2)
#     self.finalAlpha += alphaError + (finalIntervalError[:,1] - finalIntervalError[:,0])/2
#     self.finalBeta, betaError = TwoSum_NoNumba(finalInterval[:,0]/2,finalInterval[:,1]/2)
#     self.finalBeta += betaError + (finalIntervalError[:,1] + finalIntervalError[:,0])/2
#     return self.finalInterval

# def getFinalPoint(self):
#     """Finds the point that should be reported as the root (midpoint of the final step interval).

#     Returns
#     -------
#     root: numpy array
#         The final point to be reported as the root of the interval
#     """
#     #TODO: Make this a seperate function so it can use njit.
#     #Make these _NoNumba calls use floats so they call call the numba functions without a seperate compile
#     if not self.finalStep: #If no final step, use the midpoint of the calculated final interval.
#         self.root = (self.finalInterval[:,0] + self.finalInterval[:,1]) / 2
#     else: #If using the final step, recalculate the final interval using post-final transforms.
#         finalInterval = self.topInterval.T
#         finalIntervalError = np.zeros_like(finalInterval)
#         transformsToUse = self.transforms
#         for alpha,beta in transformsToUse[::-1]:
#             finalInterval, temp = TwoProd_NoNumba(finalInterval, alpha)
#             finalIntervalError = alpha * finalIntervalError + temp
#             finalInterval, temp = TwoSum_NoNumba(finalInterval,beta)
#             finalIntervalError += temp
#         finalInterval = finalInterval.T + finalIntervalError.T
#         self.root = (finalInterval[:,0] + finalInterval[:,1]) / 2 # Return the midpoint
#     return self.root

# def size(self):
#     """Gets the volume of the current interval."""
#     return np.product(self.interval[:,1] - self.interval[:,0])

function dimSize(trackedInterval)
    """Gets the lengths along each dimension of the current interval."""
    return trackedInterval.interval[2,:] - trackedInterval.interval[1,:]
end

function finalDimSize(trackedInterval)
    """Gets the lengths along each dimension of the current interval."""
    return trackedInterval.finalInterval[2,:] - trackedInterval.finalInterval[1,:]
end

# def finalDimSize(self):
#     """Gets the lengths along each dimension of the final interval."""
#     return self.finalInterval[:,1] - self.finalInterval[:,0]

# def copy(self):
#     """Returns a deep copy of the current interval with all changes and properties preserved."""
#     newone = TrackedInterval(self.topInterval)
#     newone.interval = self.interval.copy()
#     newone.transforms = self.transforms.copy()
#     newone.empty = self.empty
#     newone.nextTransformPoints = self.nextTransformPoints.copy()
#     if self.finalStep:
#         newone.finalStep = True
#         newone.canThrowOutFinalStep = self.canThrowOutFinalStep
#         newone.possibleDuplicateRoots = self.possibleDuplicateRoots.copy()
#         newone.possibleExtraRoot = self.possibleExtraRoot
#         newone.preFinalInterval = self.preFinalInterval.copy()
#         newone.preFinalTransforms = self.preFinalTransforms.copy()
#     return newone


function intervalCopy(trackedInterval::TrackedInterval)
    """Returns a deep copy of the current interval with all changes and properties preserved."""
    newone = TrackedInterval(trackedInterval.topInterval)
    newone.interval = copy(trackedInterval.interval)
    newone.transforms = copy(trackedInterval.transforms)
    newone.empty = trackedInterval.empty
    newone.nextTransformPoints = copy(trackedInterval.nextTransformPoints)
    if trackedInterval.finalStep
        newone.finalStep = true
        newone.canThrowOutFinalStep = trackedInterval.canThrowOutFinalStep
        newone.possibleDuplicateRoots = copy(trackedInterval.possibleDuplicateRoots)
        newone.possibleExtraRoot = trackedInterval.possibleExtraRoot
        newone.preFinalInterval = copy(trackedInterval.preFinalInterval)
        newone.preFinalTransforms = copy(trackedInterval.preFinalTransforms)
    end
    return newone
end

# def __contains__(self, point):
#     """Determines if point is contained in the current interval."""
#     return np.all(point >= self.interval[:,0]) and np.all(point <= self.interval[:,1])

"""DEFINITELY INCORRECT"""
function contains(trackedInterval::TrackedInterval, point)
    """Determines if point is contained in the current interval."""
    return all(x -> x >= point,tracked.interval[:,0]) && all(x -> x <= point, trackedInterval.interval[:,1])
end

# def overlapsWith(self, otherInterval):
#     """Determines if the otherInterval overlaps with the current interval.

#     Returns True if the lower bound of one interval is less than the upper bound of the other
#         in EVERY dimension; returns False otherwise."""
#     for (a1,b1),(a2,b2) in zip(self.getIntervalForCombining(), otherInterval.getIntervalForCombining()):
#         if a1 > b2 or a2 > b1:
#             return False
#     return True

function overlapsWith(trackedInterval::TrackedInterval, otherInterval::TrackedInterval)
    """Determines if the otherInterval overlaps with the current interval.

    Returns True if the lower bound of one interval is less than the upper bound of the other
        in EVERY dimension; returns False otherwise."""
    for (a1,b1,a2,b2) in zip(getIntervalForCombining(trackedInterval), getIntervalForCombining(otherInterval))
        if (a1 > b2) || (a2 > b1)
            return false
        end
    end
    return true
end


function isPoint(trackedInterval::TrackedInterval, macheps = 2^-52)
    """Determines if the current interval has essentially length 0 in each dimension."""
    return all(abs.(trackedInterval.interval[1,:] - trackedInterval.interval[2,:]) .< macheps)
end



# def startFinalStep(self):
#     """Prepares for the final step by saving the current interval and its transform list."""
#     self.finalStep = True
#     self.preFinalInterval = self.interval.copy()
#     self.preFinalTransforms = self.transforms.copy()

function startFinalStep(trackedInterval::TrackedInterval)
    """Prepares for the final step by saving the current interval and its transform list."""
    trackedInterval.finalStep = true
    trackedInterval.preFinalInterval = copy(trackedInterval.interval)
    trackedInterval.preFinalTransforms = copy(trackedInterval.transforms)
end

function getIntervalForCombining(trackedInterval::TrackedInterval)
    """Returns the interval to be used in combining intervals to report at the end."""
    return trackedInterval.finalStep ? trackedInterval.preFinalInterval : trackedInterval.interval
end

function toStr(trackedInterval::TrackedInterval)
    return string(trackedInterval.interval)
end