module ImageAnalysis

using FileIO, ColorTypes, FFTW

include("IO.jl")
include("utils.jl")
include("stats.jl")
include("histogram.jl")
include("kernels.jl")
include("watershed.jl")
include("thresholding.jl")
include("morphology.jl")
include("CCL.jl")
include("convexHullSIMP.jl")
include("integralImages.jl")
include("entropy.jl")
include("moments.jl")

greet() = print("Hello World!")

end # module
