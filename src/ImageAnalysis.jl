module ImageAnalysis

using FileIO, ColorTypes, FFTW

include("IO.jl")
include("utils.jl")
include("stats.jl")
include("histogram.jl")
include("integralImages.jl")
include("FFTconvolutions.jl")
include("kernels.jl")
include("thresholding.jl")
include("morphology.jl")
include("CCL.jl")
include("convexHullSIMP.jl")
include("entropy.jl")
include("moments.jl")
include("integral_array_utils.jl")

greet() = print("Hello World!")

end # module
