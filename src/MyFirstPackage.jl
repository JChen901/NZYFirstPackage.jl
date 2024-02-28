module MyFirstPackage
# import packages
using LinearAlgebra

# export interfaces
export Lorenz, integrate_step
export Point, Point2D, Point3D
export RungeKutta, Euclidean

export D2Q9, directions, flip_direction_index, Cell, density,momentum

# `include` other source files into this module
include("lorenz.jl")
include("fluid.jl")

end