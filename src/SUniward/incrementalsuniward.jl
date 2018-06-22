type IncrementalSUniward{B,T<:AbstractFloat} <: AbstractImageDistortion
    image::Array{UInt8}
    invrs::Array{B,1}
    waveletplanes::Array{B,1}
    distortion::T
    height::Int
    width::Int
end

Base.similar(i::IncrementalSUniward) = IncrementalSUniward(deepcopy(i.image),i.invrs,deepcopy(i.waveletplanes),i.distortion,i.height,i.width)
IncrementalSUniward(cover::Array{UInt8};sigma=1.0,T::Type=Float64) = IncrementalSUniward(cover,1,sigma=sigma,T=T)[1]

function IncrementalSUniward(cover::Array{UInt8},n::Int;sigma=1.0,T::Type=Float64)
    paddedcover = padarray(T.(cover),Pad(:symmetric,[16,16],[16,16]))
    invrs=Array{typeof(paddedcover),1}(3);
    for f in 1:length(invrs)
        invrs[f] = padarray(zeros(T,size(cover)),Pad(:symmetric,[16,16],[16,16]))
        invrs[f].parent.=1./(abs.(sameconv2(paddedcover.parent, flts[f]))+T(sigma))
    end
    mn = min(map(minimum,invrs)...)
    mx = max(map(maximum,invrs)...)
    [IncrementalSUniward(deepcopy(cover),
        invrs,
        [padarray(zeros(T,size(cover)),Pad(:symmetric,[16,16],[16,16])) for i in 1:3],
        0.0,size(cover,1),size(cover,2)) for i in 1:n]
end

function setindex!(image::IncrementalSUniward,v,i)
    setindex!(image,v,ind2sub((image.height,image.width),i)...)
end

function setindex!(image::IncrementalSUniward,v,i,j)
    (v>255 || v<0) && error("Value $v outside the allowed range [0,255]")
    idxs=(max(1,i-8):min(image.height,i+7),max(1,j-8):min(image.width,j+7))
    delta=eltype(image.waveletplanes[1])(v)-eltype(image.waveletplanes[1])(image.image[i,j])
    for f in 1:3
        A=image.waveletplanes[f]
        invr=image.invrs[f]
        flt=rotflts[f]
        for l in max(1,10-j):min(image.width+9-j,16)
            @simd for k in max(1,10-i):min(image.height+9-i,16)
                @inbounds image.distortion-=abs.(A[i-9+k,j-9+l])
                δ = delta*flt[k,l]*invr[i-9+k,j-9+l]
                @inbounds A[i-9+k,j-9+l]+=δ
                @inbounds image.distortion+=abs.(A[i-9+k,j-9+l])
            end
        end
    end
    image.image[i,j]=v
end

function tryvalue(image::IncrementalSUniward,v,i)
    tryvalue(image,v,ind2sub((image.height,image.width),i)...)
end

function tryvalue(image::IncrementalSUniward{B,T},v,i,j) where {B,T}
    if v>255 || v<0
        return(typemax(T))
    end
    distortion=image.distortion
    delta=v-image.image[i,j]
    for f in 1:3
        A=image.waveletplanes[f]
        invr=image.invrs[f]
        flt=rotflts[f]
        for l in max(1,10-j):min(image.width+9-j,16)
            @simd for k in max(1,10-i):min(image.height+9-i,16)
                @inbounds distortion-=abs.(A[i-9+k,j-9+l])
                @inbounds distortion+=abs.(A[i-9+k,j-9+l]+delta*flt[k,l]*invr[i-9+k,j-9+l])
            end
        end
    end
    return(distortion)
end

function trypmone(image::IncrementalSUniward,i)
    trypmone(image,ind2sub((image.height,image.width),i)...)
end

function trypmone(image::IncrementalSUniward,i,j)
    dp=0.0
    dm=0.0
    d=image.distortion
    for f in 1:3
        A=image.waveletplanes[f]
        invr=image.invrs[f]
        flt=rotflts[f]
        for l in max(1,10-j):min(image.width+9-j,16)
            @simd for k in max(1,10-i):min(image.height+9-i,16)
                @inbounds d-=abs.(A[i-9+k,j-9+l])
                @inbounds dp+=abs.(A[i-9+k,j-9+l]+flt[k,l]*invr[i-9+k,j-9+l])
                @inbounds dm+=abs.(A[i-9+k,j-9+l]-flt[k,l]*invr[i-9+k,j-9+l])
            end
        end
    end
    dp = (image[i,j] == 255) ? typemax(Float64) : dp
    dm = (image[i,j] == 0) ? typemax(Float64) : dm
    return(d+dp,d+dm)
end

function copy!(dest::IncrementalSUniward,src::IncrementalSUniward)
    if dest.width!=src.width || dest.height!=src.height
        error("source and destination image has to have the same height")
    end
    copy!(dest.image,src.image)
    for f in 1:3
        dest.invrs[f]=src.invrs[f]
        copy!(dest.waveletplanes[f].parent,src.waveletplanes[f].parent)
    end
    dest.distortion=src.distortion
    return(dest)
end