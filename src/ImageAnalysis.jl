module ImageAnalysis

using FileIO, ColorTypes, Statistics

include("IO.jl")
include("utils.jl")
include("histogram.jl")

greet() = print("Hello World!")

end # module
