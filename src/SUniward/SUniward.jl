module SUniward
import Base: setindex!,getindex,copy!
using Images
using OffsetArrays
import ..STC: matchlsb!,readlsb,readlsbq,LSBcostfun, AbstractImageDistortion, trypmone

include("incrementalsuniward.jl")
include("suniwardadd.jl")
include("convolutions.jl")

const hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768];
const lpdf = (-1).^(0:length(hpdf)-1).*hpdf[end:-1:1];
const flts=[lpdf*hpdf',hpdf*lpdf',hpdf*hpdf'];
const rotflts=rot180.(flts)
const wetCost = 10^8;

function suniwardcosts(cover::Array{T},sigma::T=T(1)) where {T<:AbstractFloat}
    coverPadded= padarray(cover,Pad(:symmetric,[16,16],[16,16]))
    xi=Array{Array{Float64,2},1}(3)
    R=Array{Array{Float64,2},1}(3)
    for f = 1:3
        R[f] = sameconv2(coverPadded.parent, flts[f])
        xi[f]= sameconv2(1./(abs.(R[f])+sigma),abs.(rotflts[f]))[16:end-17,16:end-17]
    end
    costs=xi[1].+xi[2].+xi[3]
    return(costs)
end

function suniwardcosts(cover::Array{T},stego::Array{T},sigma::T=T(1)) where {T<:AbstractFloat}
    df=cover-stego
    coverPadded= padarray(cover,Pad(:symmetric,[16,16],[16,16]))
    xi=Array{Array{Float64,2},1}(3)
    invr = padarray(zeros(eltype(cover),size(cover)),Pad(:symmetric,[16,16],[16,16]))
    for f = 1:3
        invr.parent.=1./(abs.(sameconv2(coverPadded.parent, flts[f]))+sigma)
        xi[f] = sameconv2(df,rotflts[f]).*invr.parent[17:end-16,17:end-16]
    end
    return(ceiledsumabs(xi[1],1e8)+ceiledsumabs(xi[2],1e8)+ceiledsumabs(xi[3],1e8))
end

suniwardcosts(cover,stego,sigma = 1) = suniwardcosts(Float64.(cover),Float64.(stego),Float64(sigma))



end