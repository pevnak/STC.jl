type SUniwardAdd <: AbstractImageDistortion
    image::Matrix{UInt8}
    rho::Matrix{Float64}
    distortion::Float64
    height::Int
    width::Int
end

SUniwardAdd(cover::Array{UInt8}) = SUniwardAdd(cover,suniwardcosts(Float64.(cover)))

function SUniwardAdd(cover::Array{UInt8},rho::Array{Float64})
  SUniwardAdd(deepcopy(cover),rho,0,size(cover)...)
end

setindex!(image::SUniwardAdd,v,i) = setindex!(image,v,ind2sub((image.height,image.width),i)...)

function setindex!(image::SUniwardAdd,v,i,j)
    image.distortion = tryvalue(image,v,i,j)
    image.image[i,j] = v
end

tryvalue(image::SUniwardAdd,v,i) = tryvalue(image,v,ind2sub((image.height,image.width),i)...)

function tryvalue(image::SUniwardAdd,v,i,j)
    if v == image.image[i,j]
        return(image.distortion)
    end
    abs(v-image.image[i,j])>1 && error("distortion can be calculated only for changes +-1")
    (v>255 || v<0) && return(typemax(Float64))

    image.distortion + image.rho[i,j]
end

trypmone(image::SUniwardAdd,i) = trypmone(image,ind2sub((image.height,image.width),i)...)

function trypmone(image::SUniwardAdd,i,j)
    v = image.image[i,j]
    return(tryvalue(image,v+1,i,j),tryvalue(image,v-1,i,j))
end

function copy!(dest::SUniwardAdd,src::SUniwardAdd)
    if dest.width!=src.width || dest.height!=src.height
        error("source and destination image has to have the same height")
    end
    copy!(dest.image,src.image) 
    dest.rho = src.rho
    dest.distortion = src.distortion
    dest.height = src.height
    dest.width = src.width
    return(dest)
end