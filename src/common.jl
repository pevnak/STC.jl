"""
    walkbacktrellis(x,message,h)

    Return the the stego signal with embedded message
"""    
function walkbacktrellis(path,x,message,hhat,h)
    numofblocks=length(message);
    indx=length(x)
    indm=length(message)
    state=0
    nstates=2^h-1;

    y=zeros(UInt8,length(x))  # the stego signal
    @inbounds for i = numofblocks:-1:1
        hh=(length(message)-indm+1<h)?hh=2^(length(message)-indm+1)-1:nstates;
        state = (2*state + message[indm])&nstates
        indm-=1
        @inbounds for j = length(hhat):-1:1
            y[indx] = path[state+1,indx]
            state = state ⊻ (y[indx]*(hhat[j]&hh)) 
            indx-=1
        end
    end
    return(y)
end


"""
    extract(y,hhat,h)

    Extract the message from stego object y. 

"""
function extract(y,hhat,h)
    numofblocks=div(length(y),length(hhat))
    state=0
    indy=1
    exm=zeros(Int,numofblocks)
    for b in 1:numofblocks
        for i=1:length(hhat)
            state=state ⊻ (y[indy]*hhat[i])
            indy+=1
        end
        exm[b]=mod(state,2)
        state>>=1
    end
    return(exm)
end

function swap(a,b)
    return(b,a)
end