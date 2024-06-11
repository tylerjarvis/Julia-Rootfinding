# testScript.jl
using Pkg
using Test

# Activate the project environment
Pkg.activate(".") #MAKE SURE YOU ARE IN THE /Julia-Rootfinding folder when you run this script in the terminal`
Pkg.instantiate()

include("tests/ChebyshevApproximatorTest.jl")  # Adjust the path as needed

