mutable struct TrackedInterval
    self.topInterval = interval
    self.interval = interval
    self.transforms = []
    self.ndim = len(self.interval)
    self.empty = False
    self.finalStep = False
    self.canThrowOutFinalStep = False
    self.possibleDuplicateRoots = []
    self.possibleExtraRoot = False
    self.nextTransformPoints = np.array([0.0394555475981047]*self.ndim) #Random Point near 0
end