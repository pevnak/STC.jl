type MinChanges
    image::Matrix{UInt8}
    cover::Matrix{UInt8}
    distortion::Int
end

MinChanges(image::Array{UInt8}) = MinChanges(deepcopy(image),deepcopy(image),0)
getindex(image::MinChanges,i::Int64) = image.image[i]

function setindex!(a::MinChanges,v,i::Int64)
    a.distortion -= a.image[i] != a.cover[i]
    a.image[i] = v
    a.distortion += a.image[i] != a.cover[i]
end

function LSBcostfun(image::MinChanges,pidx)
    if image.image[pidx]&0x01==0x01
        return(image.distortion + 1,0)
    else
        return(0,image.distortion + 1)
    end    
end

function copy!(dest::MinChanges,src::MinChanges)
    copy!(dest.image,src.image)
    dest.distortion = src.distortion
    return(dest)
end


readlsb(image::MinChanges,embpath) = readlsb(image.image,embpath)

function readlsb(img::Matrix{UInt8},embpath)
    y=zeros(UInt8,length(embpath))
    for (i,ii) in enumerate(embpath)
        y[i] = img[ii] & 0x1
    end
    return(y)
end

function matchlsb!(img::MinChanges,mbit,pidx::Int)
    original=img[pidx]
    if img[pidx]&0x1!=mbit
        if original==255
            img[pidx]=original-1
            return
        end
        if original==0
            img[pidx]=original+1
            return
        end
        if rand()>0.5
            img[pidx]=original-1
        else
            img[pidx]=original+1
        end
    end
    @assert(mod(img[pidx],2) == mbit)
end