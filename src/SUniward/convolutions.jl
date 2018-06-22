

function reflectedsetindex!(b::OffsetArrays.OffsetArray,k,l,v)
    b[k,l]=v    
    s=(size(b.parent,1) + 2*(b.offsets[1]),size(b.parent,2) + 2*(b.offsets[2]))
    #the four rules below resolve all the basic cases
    if k<=-b.offsets[1]
        b[1-k,l]=v
    end
    if l<=-b.offsets[2]
        b[k,1-l]=v
    end

    if k>s[1]+b.offsets[1]
       b[2*s[1]-k+1,l]=v 
    end
    
    if l>s[2]+b.offsets[2]
       b[k,2*s[2]-l+1]=v 
    end

    #the next four resolves remaining corner cases
    if k<=-b.offsets[1]&&l<=-b.offsets[2]
        b[1-k,1-l]=v
    end

    if k>s[1]+b.offsets[1] && l>s[2]+b.offsets[2]
       b[2*s[1]-k+1,2*s[2]-l+1]=v 
    end

    if k<=-b.offsets[1] && l>s[2]+b.offsets[2]
       b[1-k,2*s[2]-l+1]=v 
    end

    if k>s[1]+b.offsets[1] && l<=-b.offsets[2]
       b[2*s[1]-k+1,1-l]=v 
    end
end

"""
    function sameconv2(A,f)

    returns convolution of A with f that mimics matlab 'conv2(A,f,'same')
"""
function sameconv2(A,f)
    (ul,vl)=(cld(size(f,1),2),cld(size(f,2),2))
    (up,vp)=(fld(size(f,1),2),fld(size(f,2),2))
    paddedA=padarray(A,Fill(0,[ul,vl],[ul,vl]))
    r=zeros(eltype(A),size(A))
    f=rot180(f)
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            r[i,j] = sum(view(paddedA,i-ul+1:i+up,j-vl+1:j+vp) .* f)
        end
    end
    return(r)
end


"""
    function oneitemconv2!(A,f,i,j,v)

    returns updates the result of convolution A (A=B*f), after B is changed at coordinates (i,j) by v (B[i,j]+=v)
"""
function oneitemconv2!(A,f,i,j,v)
    (ul,vl)=(cld(size(f,1),2),cld(size(f,2),2))
    (up,vp)=(fld(size(f,1),2),fld(size(f,2),2))
    A[i-up:i+ul-1,j-vp:j+vl-1].+=v.*f
    return(A)
end

"""
    function oneitemconv2!(A,invr,f,i,j,v)

    returns updates the result of convolution A (A=(B*f)./invr), after B is changed at coordinates (i,j) by v (B[i,j]+=v)
"""
function oneitemconv2!(A,invr,f,i,j,v)
    (ul,vl)=(cld(size(f,1),2),cld(size(f,2),2))
    (up,vp)=(fld(size(f,1),2),fld(size(f,2),2))
    for l in 1:size(f,2)
        @simd for k in 1:size(f,1)
            @inbounds A[i-up+k-1,j-vp+l-1]+=v*f[k,l]*invr[i-up+k-1,j-vp+l-1]
        end
    end
    return(A)
end

function oneitemconv2!(A,invr,f,i,j,v,ul,vl,up,vp)
    for l in 1:size(f,2)
        @simd for k in 1:size(f,1)
            @inbounds A[i-up+k-1,j-vp+l-1]+=v*f[k,l]*invr[i-up+k-1,j-vp+l-1]
        end
    end
    return(A)
end

"""
function ceiledsumabs(a,c)

returns sumabs(min(a,c))
"""
function ceiledsumabs(a,c)
r=zero(eltype(a))
for v in a
    r+=min(c,abs(v))
end
return(r)
end
