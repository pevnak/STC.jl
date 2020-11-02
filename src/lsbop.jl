function LSBcostfun(stego::T,pidx) where {T<:AbstractImageDistortion}
    c0=stego.distortion
    c1=min(trypmone(stego,pidx)...)
    (stego[pidx]&0x1 == 0x0) ? (c0,c1) : (c1,c0)
end


function LSBcostfun2(stego::T,pidx) where {T<:AbstractImageDistortion}
    d = LSBcostfun(stego,pidx)
    [d[1],d[2]]
end

function LSBcostfun3(stego::T,pidx) where {T<:AbstractImageDistortion}
    d0 = stego.distortion
    dp,dm = trypmone(stego,pidx)
    m = mod(stego[pidx],3)
    if m == 0
        !((mod(stego[pidx],3)==0) && (mod(stego[pidx]+1,3)==1) && (mod(stego[pidx]-1,3)==2)) && println(stego[pidx])
        return[d0,dp,dm]
    elseif m == 1
        !((mod(stego[pidx]-1,3)==0) && (mod(stego[pidx],3)==1) && (mod(stego[pidx]+1,3)==2)) && println(stego[pidx])
        return[dm,d0,dp]
    else
        !((mod(stego[pidx]+1,3)==0) && (mod(stego[pidx]-1,3)==1) && (mod(stego[pidx],3)==2)) && println(stego[pidx])
        return[dp,dm,d0]
    end
end

LSBcostfun(stego::T,pidx,q) where {T<:AbstractImageDistortion} = (q==2) ? LSBcostfun2(stego,pidx) : LSBcostfun3(stego,pidx);


function matchlsb!(img::T,mbits,embpath::Array{Int,1}) where {T<:AbstractImageDistortion}
    for i in 1:length(mbits)
        matchlsb!(img,mbits[i],embpath[i])
    end
    return(img)
end

function matchlsb2!(img::T,mbit,pidx::Int) where {T<:AbstractImageDistortion}
    original=img[pidx]
    # println("matchlsb mbit = $mbit")
    if img[pidx]&0x1!=mbit
        if original==255
            img[pidx]=original-1
            return
        end
        if original==0
            img[pidx]=original+1
            return
        end
        (dp,dm)=trypmone(img,pidx)
        if abs(dp-dm)<1e-8
            img[pidx]+=rand([-1,+1])
        elseif dp>dm
            img[pidx]=original-1
        else
            img[pidx]=original+1
        end
    end
    @assert(mod(img[pidx],2) == mbit)
end

function matchlsb3!(img::T,mbit,pidx::Int) where {T<:AbstractImageDistortion}
    original=img[pidx]
    m = mod(original,3)
    if m != mbit 
        if mod(original+1,3) == mbit
            img[pidx]=original+1
        elseif mod(original-1,3) == mbit
            img[pidx]=original-1
        else 
            error("cannot embed")
        end
    end
    @assert(mod(img[pidx],3) == mbit)
end

matchlsb!(img::T,mbit,pidx::Int) where {T<:AbstractImageDistortion} = matchlsb2!(img,mbit,pidx)
matchlsb!(img::T,mbit,pidx::Int,q) where {T<:AbstractImageDistortion} = (q==2) ? matchlsb2!(img,mbit,pidx) : matchlsb3!(img,mbit,pidx);

readlsb(img::T,embpath) where {T<:AbstractImageDistortion} = readlsbq(img,embpath,2)

function readlsbq(img::T,embpath,q) where {T<:AbstractImageDistortion}
    y=zeros(UInt8,length(embpath))
    for (i,ii) in enumerate(embpath)
        y[i]=Int(mod(img.image[ii],q))
    end
    return(y)
end