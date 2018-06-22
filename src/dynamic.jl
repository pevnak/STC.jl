using Base.Threads 
"""
    instates(k,w,q,hh) 

    k --- index of the current state (0 -- hh-1)
    w --- current coding word (column of hhat matrix)
    q --- base of the coding (2,3,4,5)?
    hh --- height of w
"""
instates(k,w,q,hh) = [(mod(k - i*w,q^hh) + 1) for i in 0:q-1]

"""
    (y,wy) = findbeststate(inds,wght,w)

    finds the word with the lowest cost and the word to code
"""
function findbeststate(inds,wght,w)
    ww = [isinf(wght[inds[i]]) ? wght[inds[i]] : w[inds[i],i] for i in 1:length(inds)]
    # println("instates = $(inds-1) wght = $(wght[inds]) w = $ww  ")
    y = 1
    wy = Inf64
    for i in 1:length(inds)
        t = isinf(wght[inds[i]]) ? wght[inds[i]] : w[inds[i],i]
        if t<wy
            wy = t;
            y = i
        end
    end
    (y,wy)
end

function updatestate!(k,hhat,q,hh,w,wght,newwght,stegos,newstegos,pidx)
    inds = instates(k,hhat,q,hh)
    y,wy = findbeststate(inds,wght,w)
    newwght[k+1] = wy
    if !isinf(wy)
        copy!(newstegos[k+1],stegos[inds[y]])
        matchlsb!(newstegos[k+1],y-1,pidx,q)
    end
end

"""

    variablecoding(cover,stegos,message,hhat,h,embpath,q=2)

    embed the message `m` into cover object `cover` using sub-matrix `hhat` with height `h`


"""
function variablecoding(cover,stegos,message,hhat,h,embpath,q=2)
    newstegos=deepcopy(stegos)
    numofblocks=length(message);

    w=zeros(q^h,q);        #holds the costs of current "state"
    wght=zeros(q^h);        #holds the costs of current "state"
    newwght=zeros(q^h);     #contained for holding costs of currently generated states
    wght[2:end]=Inf64;      
    indx=1  #this indexes the cover
    indm=1  #this indexes the message bits
    meter = Progress(numofblocks,5)
    for i in 1:numofblocks #iterate over sub-blocks, each coding a single message
        hh=min(h,length(message)-indm+1)  #determine  effective height of the sub-matrix
        for j in 1:length(hhat) #iterate over columns of the sub-block

            # determine costs of reaching state k from different stats
            for k in 1:q^hh 
                w[k,:] = LSBcostfun(stegos[k],embpath[indx],q)
            end

            # iterate over each state of the matrix
            for k in 0:q^hh-1 
                updatestate!(k,hhat[j],q,hh,w,wght,newwght,stegos,newstegos,embpath[indx])
            end
            indx+=1
            (stegos,newstegos) = swap(stegos,newstegos)
            (wght,newwght) = swap(wght,newwght)
        end

        # After processing a sub-block, remove states that does not code the message
        for j = 0:q^(hh-1)-1
            midx=q*j + message[indm]+1
            wght[j+1] = wght[midx]
            aa = stegos[midx]; stegos[midx] = stegos[j+1]; stegos[j+1] = aa;
        end
        # println()
        wght[q^(hh-1)+1:end]=Inf64;
        indm+=1
        next!(meter)
    end
    embedding_cost = wght[1]

    # check that the extracted message match the embedded message
    if !isinf(embedding_cost)
        info(@sprintf("embedding costs by stc = %f",embedding_cost))
        extm=extractq(readlsbq(stegos[1],embpath,q),hhat,h,q)[1:length(message)];
        sum(abs.(extm.!=message)) !=0 && warn("message has not been correctly extracted");
    else
        warn("embedding has failed");
    end
    return(deepcopy(stegos[1]),stegos[1].distortion)
end

"""
    extractq(y,hhat,h)

    Extract the message from stego object y. 

"""
function extractq(y,hhat,h,q)
    numofblocks=div(length(y),length(hhat))
    state, indy = 0, 1
    exm=zeros(Int,numofblocks)
    for b in 1:numofblocks
        hh=min(h,numofblocks-b+1)
        for i=1:length(hhat)
            state = mod(state + y[indy]*hhat[i],q^hh)
            indy += 1
        end
        exm[b] = mod(state,q)
        state = div(state,q)
    end
    return(exm)
end

 variablecoding(cover,stegos,message,hhat,h) = variablecoding(cover,stegos,message,hhat,h,randperm(length(cover)))
