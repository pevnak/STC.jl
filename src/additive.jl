"""
    costructtrellis(x,message,hhat,h,costfun)

hhat --- submatrix stored in vector, where each item is a integer representation of a column considered as a binary number
h --- height if the matrix
message --- message that should be embedded
x --- cover signal
"""    
function additivecoding(cover,message,hhat,h,costfun,embpath)
    x=readlsb(cover,embpath);
    img=mod.(deepcopy(cover),8)

    path=zeros(UInt8,2^h,length(x));
    numofblocks=length(message);

    wght=zeros(2^h);
    newwght=zeros(2^h);
    wght[2:end]=Inf64;
    indx=1  #this indexes the cover
    indm=1  #this indexes the message bits
    y=zeros(UInt8,length(x))
    for i in 1:numofblocks #iterate over the number of subblocks, each encodes one message bit
        hh=min(h,length(message)-indm+1)  #this is effective height of the submatrix, which is needed
        for j in 1:length(hhat) #iterate over columns of the H matrix
            # println("pidx $(embpath[indx])")
            pixelcosts=costfun(cover,embpath[indx])    #determine the cost of the pixel for being zero
            for k in 0:2^hh-1 # iterate over each state of the matrix
                # there are two possibilities to reach the state k. Either we code zero, 
                # in which case the price (weight) is equal to weight of the same state in previous iteration + possible change,
                # or we code one and the state is reached from other state determined by the column of the hhat matrix.
                instates = [k+1,(k ⊻ hhat[j])&(2^hh-1)+1]
                w0 = wght[instates[1]] + x[indx]*pixelcosts[1]
                w1 = wght[instates[2]] + (1-x[indx])*pixelcosts[2]
                ww = [w0,w1]
                # println("instates = $(instates-1) wght = $(wght[instates]) w = $ww")

                # choose the less costly path to the state
                path[k+1,indx] = (w1 < w0) ? 0x1 : 0x0 
                newwght[k+1] = min(w0, w1)
            end
            indx+=1
            (wght,newwght)=swap(wght,newwght)
        end

        # once we have traversed the trellis, we need to determine which path through trellises codes the message and
        # remove the others
        for j = 0:2^(hh-1)-1
            wght[j+1] = wght[2*j + message[indm]+1]
        end
        wght[2^(hh-1)+1:end]=Inf64;
        newwght[2^(hh-1)+1:end]=Inf64;
        indm+=1
    end
    # println(map(Int,y))

    yy=walkbacktrellis(view(path,:,1:indx-1),view(x,1:indx-1),view(message,1:indm-1),hhat,h);
    # println(map(Int,yy))
    stego=deepcopy(cover);
    matchlsb!(stego,yy,embpath)
    embedding_cost = wght[1]

    #check that the extracted message match the embedded message
    extm=extract(readlsb(stego,embpath),hhat,h)[1:length(message)];
    println(@sprintf("embedding costs by stc = %f",embedding_cost))
    # @assert(any(extm.!=message)==false)
    return(stego,embedding_cost)
end

additivecoding(cover,message,hhat,h,costfun) = additivecoding(cover,message,hhat,h,costfun,randperm(length(cover)))

function matchlsb!(img::Array{UInt8,2},mbits::Vector,embpath::Vector)
    for i in 1:length(mbits)
        matchlsb!(img,mbits[i],embpath[i])
    end
    return(img)
end

function matchlsb!(img::Array{UInt8,2},m::T,i::S) where {T<:Integer,S<:Integer}
    if mod(img[i],2) != m 
        if img[i] == 255 
            img[i] = 254
        elseif img[i] == 0
            img[i] = 1
        else
            img[i] += rand([-1,+1])
        end 
    end
end