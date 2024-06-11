module YRoots #This is the main module file, for versioning and test suite compatibility

# Include other source files
include("ChebyshevApproximator.jl")

# Optionally, re-export symbols from included files
using .ChebyshevApproximator

end # module YRoots
