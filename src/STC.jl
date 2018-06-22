module STC
using ProgressMeter
import Base: setindex!,getindex,copy!
export additivestc, extract

function readlsbq end
function readlsb end
function trypmone end
function updatedistortion! end
abstract type AbstractImageDistortion end

Base.getindex(image::T,i,j) where {T<:AbstractImageDistortion} = image.image[i,j]
Base.getindex(image::T,i) where {T<:AbstractImageDistortion} = getindex(image,ind2sub((image.height,image.width),i)...)

include("lsbop.jl")
include("additive.jl")
include("common.jl")
include("randomizedtrellis.jl")
include("dynamic.jl")
include("minchanges.jl")

include("SUniward/SUniward.jl")

end