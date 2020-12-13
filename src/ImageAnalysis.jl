module ImageAnalysis

using FileIO, ColorTypes, FFTW

include("IO.jl")
include("utils.jl")
include("stats.jl")
include("histogram.jl")
include("clusters.jl")
include("kernels.jl")

greet() = print("Hello World!")

end # module
