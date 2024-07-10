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
    reducedDims # = []
    solvedVals # = []
    function TrackedInterval(interval)
        ndim = length(interval)
        new(interval,interval,[],ndim,false,false,false,[],false,[0.0394555475981047]*ndim,[],[],[],[])
    end
end

"""==============================FUNCTIONS FOR TRACKED INTERVAL=============================="""

# def canThrowOut(self):
# """Ensures that an interval that has not subdivided cannot be thrown out on the final step."""
# return not self.finalStep or self.canThrowOutFinalStep


function canThrowOut(trackedInterval::TrackedInterval)
    """Ensures that an interval that has not subdivided cannot be thrown out on the final step."""
    return !trackedInterval.finalStep || trackedInterval.canThrowOutFinalStep
end

# def addTransform(self, subInterval):
#     """Adds the next alpha and beta values to the list transforms and updates the current interval.

#     Parameters:
#     -----------
#     subInterval : numpy array
#         The subinterval to which the current interval is being reduced
#     """
#     #Ensure the interval has non zero size; mark it empty if it doesn't
#     if np.any(subInterval[:,0] > subInterval[:,1]) and self.canThrowOut():
#         self.empty = True
#         return
#     elif np.any(subInterval[:,0] > subInterval[:,1]):
#         #If we can't throw the interval out, it should be bounded by [-1,1].
#         subInterval[:,0] = np.minimum(subInterval[:,0], np.ones_like(subInterval[:,0]))
#         subInterval[:,0] = np.maximum(subInterval[:,0], -np.ones_like(subInterval[:,0]))
#         subInterval[:,1] = np.minimum(subInterval[:,1], np.ones_like(subInterval[:,0]))
#         subInterval[:,1] = np.maximum(subInterval[:,1], subInterval[:,0])
#     # Get the alpha and beta associated with the transformation in each dimension
#     a1,b1 = subInterval.T # all the lower bounds and upper bounds of the new interval, respectively
#     a2,b2 = self.interval.T # all the lower bounds and upper bounds of the original interval
#     alpha1, beta1 = (b1-a1)/2, (b1+a1)/2
#     alpha2, beta2 = (b2-a2)/2, (b2+a2)/2
#     self.transforms.append(np.array([alpha1, beta1]))
#     #Update the lower and upper bounds of the current interval
#     for dim in range(self.ndim):
#         for i in range(2):
#             x = subInterval[dim][i]
#             #Be exact if x = +-1
#             if x == -1.0:
#                 self.interval[dim][i] = self.interval[dim][0]
#             elif x == 1.0:
#                 self.interval[dim][i] = self.interval[dim][1]
#             else:
#                 self.interval[dim][i] = alpha2[dim]*x+beta2[dim]

# def getLastTransform(self):
#     """Gets the alpha and beta values of the last transformation the interval underwent."""
#     return self.transforms[-1]

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

# def dimSize(self):
#     """Gets the lengths along each dimension of the current interval."""
#     return self.interval[:,1] - self.interval[:,0]

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


function copy(trackedInterval::TrackedInterval)
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

# def isPoint(self):
#     """Determines if the current interval has essentially length 0 in each dimension."""
#     return np.all(np.abs(self.interval[:,0] - self.interval[:,1]) < 1e-32)

# possibly helpful:
# result = all(((x, y) -> x + y > 5), zip(arr1, arr2)) something like this

"""DEFINITELY INCORRECT"""
function isPoint(trackedInterval::TrackedInterval)
    """Determines if the current interval has essentially length 0 in each dimension."""
    return all(abs(trackedInterval.interval[:,0] - trackedInterval.interval[:,1]) < 1e-32)
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

# def getIntervalForCombining(self):
#     """Returns the interval to be used in combining intervals to report at the end."""
#     return self.preFinalInterval if self.finalStep else self.interval

function getIntervalForCombining(trackedInterval::TrackedInterval)
    """Returns the interval to be used in combining intervals to report at the end."""
    return (trackedInterval.finalStep ? trackedInterval.preFinalInterval : trackedInterval.interval)
end

# def __repr__(self):
#     return str(self)

# def __str__(self):
#   return str(self.interval)

function toStr(trackedInterval::TrackedInterval)
    return string(trackedInterval.interval)
end