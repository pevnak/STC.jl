module STC
using ProgressMeter
using Printf
import Base: setindex!,getindex,copy!
export additivestc, extract

function readlsbq end
function readlsb end
function trypmone end
function updatedistortion! end
abstract type AbstractImageDistortion end

ind2sub(image::AbstractImageDistortion, i::Int) = Tuple(CartesianIndices((image.height,image.width))[i])
Base.getindex(image::AbstractImageDistortion,i,j) = image.image[i,j]
Base.getindex(image::AbstractImageDistortion,i) = getindex(image,ind2sub(image, i)...)

include("lsbop.jl")
include("additive.jl")
include("common.jl")
include("randomizedtrellis.jl")
include("dynamic.jl")
include("minchanges.jl")

include("SUniward/SUniward.jl")

end